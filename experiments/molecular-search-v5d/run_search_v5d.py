#!/usr/bin/env python3
"""V5d-1: high-throughput evolutionary molecular search."""
from __future__ import annotations

import argparse
import copy
import gzip
import hashlib
import json
import math
import platform
import statistics
import subprocess
import sys
from collections import Counter, defaultdict
from concurrent.futures import ProcessPoolExecutor, as_completed
from importlib.metadata import version
from pathlib import Path
from typing import Any, Iterable

import yaml

HERE = Path(__file__).resolve().parent
EXPERIMENTS = HERE.parent
ENTRY = EXPERIMENTS / "entry-evaluator"
V5B = EXPERIMENTS / "property-predictor-v5b"
V5C = EXPERIMENTS / "molecular-search-v5c"
for path in (ENTRY, V5B, V5C, HERE):
    if str(path) not in sys.path:
        sys.path.insert(0, str(path))

from applicability import (  # noqa: E402
    calibrate_metric,
    canonical_smiles,
)
from evaluator import write_summary  # noqa: E402
from graph_mutator import (  # noqa: E402
    generate_mutations,
    structural_bucket,
)
from property_provider import build_property_provider  # noqa: E402
from run_search import (  # noqa: E402
    _coolant_config,
    _descriptor_check,
    _load_calibration_records,
    _sha256,
    _stable_candidate_id,
    _water_reference,
    _write_csv,
)
from run_search_v5c2 import (  # noqa: E402
    _calibration_atomic_numbers,
    _entry_candidate,
    _prescreen_candidate,
    _prescreen_sort_key,
    _water_property_reference,
)
from thermo_prescreen import (  # noqa: E402
    compare_to_reference,
    passes_property_gate,
    property_screen,
)


def _validate_config(config: dict[str, Any]) -> None:
    ht = config["high_throughput"]
    positive_integer_keys = (
        "generations",
        "children_per_parent",
        "initial_parent_limit",
        "beam_width",
        "generation_unique_target",
        "total_unique_target",
        "coarse_eos_budget_per_generation",
        "full_prescreen_budget",
        "entry_rankable_budget",
        "entry_exploratory_budget",
    )
    for key in positive_integer_keys:
        value = int(ht[key])
        if value <= 0:
            raise ValueError(f"high_throughput.{key} must be positive")

    if int(ht["generation_unique_target"]) > int(ht["total_unique_target"]):
        raise ValueError(
            "generation_unique_target cannot exceed total_unique_target"
        )
    diversity = float(ht["diversity_fraction"])
    if not 0.0 <= diversity <= 1.0:
        raise ValueError("diversity_fraction must be in [0, 1]")

    operators = list(config["mutation"]["operators"])
    if not operators:
        raise ValueError("mutation operator list must be non-empty")

    mutation_atoms = {
        int(value)
        for value in config["mutation"]["allowed_atomic_numbers"]
    }
    descriptor_atoms = {
        int(value)
        for value in config["candidate_constraints"][
            "allowed_atomic_numbers"
        ]
    }
    if not mutation_atoms.issubset(descriptor_atoms):
        raise ValueError(
            "mutation elements must be allowed by candidate constraints"
        )

    rankable = set(config["domain_expansion"]["rankable_families"])
    exploratory = set(
        config["domain_expansion"]["exploratory_families"]
    )
    if rankable & exploratory:
        raise ValueError(
            "rankable and exploratory family policies must be disjoint"
        )

    coarse_states = (
        len(config["coarse_prescreen"]["temperature_grid_k"])
        * len(config["coarse_prescreen"]["pressure_grid_pa"])
    )
    full_states = (
        len(config["prescreen"]["temperature_grid_k"])
        * len(config["prescreen"]["pressure_grid_pa"])
    )
    if coarse_states >= full_states:
        raise ValueError(
            "coarse prescreen must contain fewer states than full prescreen"
        )


def _stable_hash(value: str) -> str:
    return hashlib.sha256(value.encode("utf-8")).hexdigest()


def _family_policy(config: dict[str, Any]) -> set[str]:
    return {
        str(value)
        for key in ("rankable_families", "exploratory_families")
        for value in config["domain_expansion"][key]
    }


def _candidate_record(
    smiles: str,
    *,
    parent_smiles: str,
    operator: str,
    generation: int,
    descriptor: dict[str, Any],
    family: str,
    bucket: tuple[str, int, int, int],
) -> dict[str, Any]:
    _, heavy_atoms, hetero_atoms, rings = bucket
    return {
        "candidate_id": _stable_candidate_id(smiles),
        "smiles": smiles,
        "generation": int(generation),
        "parent_smiles": parent_smiles,
        "mutation_operator": operator,
        "descriptor": descriptor,
        "structural_bucket": {
            "family": family,
            "heavy_atoms": heavy_atoms,
            "hetero_atoms": hetero_atoms,
            "rings": rings,
        },
        "provenance": [
            {
                "family": family,
                "operation": operator,
                "parameters": {
                    "generation": int(generation),
                    "parent_smiles": parent_smiles,
                },
                "raw_smiles": smiles,
            }
        ],
    }


def _load_json_candidates(path: Path | None) -> set[str]:
    if path is None or not path.is_file():
        return set()
    payload = json.loads(path.read_text(encoding="utf-8"))
    result: set[str] = set()
    for row in payload.get("candidates", []):
        smiles = row.get("smiles")
        if smiles:
            try:
                result.add(canonical_smiles(str(smiles)))
            except Exception:
                pass
    return result


def _load_v5c2_seed_rows(
    path: Path | None,
    *,
    limit: int,
) -> list[dict[str, Any]]:
    if path is None or not path.is_file() or limit <= 0:
        return []
    payload = json.loads(path.read_text(encoding="utf-8"))
    rows = list(payload.get("passes", []))
    rows.sort(
        key=lambda row: (
            -float(
                row.get("property_screen", {}).get(
                    "property_priority_score"
                )
                or -1.0e30
            ),
            str(row.get("candidate_id", "")),
        )
    )
    output: list[dict[str, Any]] = []
    for row in rows[:limit]:
        smiles = row.get("smiles")
        if not smiles:
            continue
        try:
            output.append(
                {
                    "smiles": canonical_smiles(str(smiles)),
                    "source": "v5c2_prescreen",
                }
            )
        except Exception:
            continue
    return output


def _initial_parents(
    calibration_records: list[dict[str, Any]],
    config: dict[str, Any],
    v5c2_prescreen_path: Path | None,
) -> list[str]:
    sources = config["seed_sources"]
    rows: list[dict[str, str]] = []
    if bool(sources.get("include_v5b_calibration", True)):
        rows.extend(
            {
                "smiles": canonical_smiles(str(record["smiles"])),
                "source": "v5b_calibration",
            }
            for record in calibration_records
        )
    if bool(sources.get("include_v5c2_prescreen", True)):
        rows.extend(
            _load_v5c2_seed_rows(
                v5c2_prescreen_path,
                limit=int(
                    sources.get(
                        "v5c2_prescreen_seed_limit",
                        600,
                    )
                ),
            )
        )

    unique: dict[str, str] = {}
    for row in rows:
        unique.setdefault(row["smiles"], row["source"])

    limit = int(
        config["high_throughput"]["initial_parent_limit"]
    )
    ordered = sorted(
        unique,
        key=lambda smiles: (
            unique[smiles] != "v5c2_prescreen",
            _stable_hash(smiles),
        ),
    )
    return ordered[:limit]


def _structural_select(
    rows: list[dict[str, Any]],
    budget: int,
) -> list[dict[str, Any]]:
    """Deterministic round-robin selection across structural buckets."""
    if budget <= 0:
        return []
    if len(rows) <= budget:
        return sorted(rows, key=lambda row: str(row["candidate_id"]))

    groups: dict[tuple[Any, ...], list[dict[str, Any]]] = defaultdict(list)
    for row in rows:
        structural = row.get("structural_bucket") or {}
        family = str(row["provenance"][0]["family"])
        bucket = (
            family,
            int(structural.get("heavy_atoms", row["descriptor"]["heavy_atoms"])),
            int(structural.get("hetero_atoms", 0)),
            int(structural.get("rings", row["descriptor"]["rings"])),
        )
        groups[bucket].append(row)

    for bucket_rows in groups.values():
        bucket_rows.sort(
            key=lambda row: _stable_hash(str(row["candidate_id"]))
        )

    buckets = sorted(groups, key=lambda value: tuple(map(str, value)))
    indices = {bucket: 0 for bucket in buckets}
    selected: list[dict[str, Any]] = []
    while len(selected) < budget:
        progressed = False
        for bucket in buckets:
            index = indices[bucket]
            rows_in_bucket = groups[bucket]
            if index >= len(rows_in_bucket):
                continue
            selected.append(rows_in_bucket[index])
            indices[bucket] = index + 1
            progressed = True
            if len(selected) >= budget:
                break
        if not progressed:
            break
    return selected


def _score_value(row: dict[str, Any]) -> float:
    value = row.get("coarse_property_screen", {}).get(
        "property_priority_score"
    )
    if isinstance(value, (int, float)) and math.isfinite(float(value)):
        return float(value)
    return -math.inf


def _beam_select(
    rows: list[dict[str, Any]],
    *,
    width: int,
    diversity_fraction: float,
) -> list[dict[str, Any]]:
    if width <= 0 or not rows:
        return []
    if len(rows) <= width:
        return sorted(
            rows,
            key=lambda row: (
                -_score_value(row),
                str(row["candidate_id"]),
            ),
        )

    diversity_fraction = max(0.0, min(1.0, float(diversity_fraction)))
    diversity_budget = int(round(width * diversity_fraction))
    score_budget = width - diversity_budget

    ranked = sorted(
        rows,
        key=lambda row: (
            -_score_value(row),
            str(row["candidate_id"]),
        ),
    )
    selected = ranked[:score_budget]
    selected_ids = {
        str(row["candidate_id"])
        for row in selected
    }

    remaining = [
        row
        for row in rows
        if str(row["candidate_id"]) not in selected_ids
    ]
    diverse = _structural_select(remaining, diversity_budget)
    selected.extend(diverse)

    if len(selected) < width:
        selected_ids = {
            str(row["candidate_id"])
            for row in selected
        }
        for row in ranked:
            if str(row["candidate_id"]) in selected_ids:
                continue
            selected.append(row)
            selected_ids.add(str(row["candidate_id"]))
            if len(selected) >= width:
                break

    selected.sort(
        key=lambda row: (
            -_score_value(row),
            str(row["candidate_id"]),
        )
    )
    return selected[:width]


def _coarse_config(config: dict[str, Any]) -> dict[str, Any]:
    value = copy.deepcopy(config)
    value["prescreen"] = copy.deepcopy(config["coarse_prescreen"])
    return value


def _coarse_candidate(task: dict[str, Any]) -> dict[str, Any]:
    candidate = task["candidate"]
    config = task["config"]
    water_properties = task["water_properties"]

    smiles = str(candidate["smiles"])
    try:
        coolant = _coolant_config(
            config,
            smiles=smiles,
            name=str(candidate["candidate_id"]),
        )
        provider = build_property_provider(coolant)
    except Exception as exc:
        return {
            **candidate,
            "coarse_status": "rejected",
            "coarse_rejection_stage": "gc_not_representable",
            "coarse_rejection_reason": str(exc),
        }

    try:
        metrics = property_screen(provider, config)
    except Exception as exc:
        return {
            **candidate,
            "coarse_status": "rejected",
            "coarse_rejection_stage": "property_failure",
            "coarse_rejection_reason": str(exc),
        }

    if (
        bool(config["storage"].get("require_liquid", True))
        and metrics.get("storage_phase") != "liquid"
    ):
        return {
            **candidate,
            "coarse_status": "rejected",
            "coarse_rejection_stage": "storage_not_liquid",
            "coarse_rejection_reason": (
                f"phase={metrics.get('storage_phase')}"
            ),
            "coarse_property_screen": metrics,
        }

    targets = compare_to_reference(metrics, water_properties)
    metrics = {**metrics, **targets}
    passed, reason = passes_property_gate(metrics, config)
    if not passed:
        return {
            **candidate,
            "coarse_status": "rejected",
            "coarse_rejection_stage": "property_target",
            "coarse_rejection_reason": reason,
            "coarse_property_screen": metrics,
        }

    return {
        **candidate,
        "coarse_status": "pass",
        "coarse_property_screen": metrics,
    }


def _run_parallel(
    function,
    tasks: list[dict[str, Any]],
    *,
    workers: int,
) -> list[dict[str, Any]]:
    if workers <= 1:
        return [function(task) for task in tasks]
    rows: list[dict[str, Any]] = []
    with ProcessPoolExecutor(max_workers=workers) as pool:
        futures = [pool.submit(function, task) for task in tasks]
        for future in as_completed(futures):
            rows.append(future.result())
    return rows


def _write_jsonl_gz(
    path: Path,
    rows: Iterable[dict[str, Any]],
) -> None:
    with gzip.open(path, "wt", encoding="utf-8") as handle:
        for row in rows:
            handle.write(
                json.dumps(row, sort_keys=True, separators=(",", ":"))
            )
            handle.write("\n")


def _distribution_summary(values: list[float]) -> dict[str, Any]:
    finite = [
        float(value)
        for value in values
        if math.isfinite(float(value))
    ]
    if not finite:
        return {
            "count": 0,
            "min": None,
            "median": None,
            "max": None,
        }
    return {
        "count": len(finite),
        "min": min(finite),
        "median": statistics.median(finite),
        "max": max(finite),
    }


def _approximate_diversity(
    rows: list[dict[str, Any]],
    *,
    sample_size: int = 400,
) -> dict[str, Any]:
    from rdkit import Chem, DataStructs
    from rdkit.Chem import rdFingerprintGenerator

    if len(rows) < 2:
        return {
            "sample_count": len(rows),
            "mean_tanimoto_distance": None,
        }

    sample = sorted(
        rows,
        key=lambda row: _stable_hash(str(row["candidate_id"])),
    )[:sample_size]
    generator = rdFingerprintGenerator.GetMorganGenerator(
        radius=2,
        fpSize=2048,
    )
    fingerprints = []
    for row in sample:
        mol = Chem.MolFromSmiles(str(row["smiles"]))
        if mol is not None:
            fingerprints.append(generator.GetFingerprint(mol))

    distances: list[float] = []
    for index, fp in enumerate(fingerprints):
        if index == 0:
            continue
        similarities = DataStructs.BulkTanimotoSimilarity(
            fp,
            fingerprints[:index],
        )
        distances.extend(1.0 - float(value) for value in similarities)

    return {
        "sample_count": len(fingerprints),
        "pair_count": len(distances),
        "mean_tanimoto_distance": (
            statistics.mean(distances)
            if distances
            else None
        ),
        "median_tanimoto_distance": (
            statistics.median(distances)
            if distances
            else None
        ),
    }


def _approximate_reference_novelty(
    rows: list[dict[str, Any]],
    calibration_records: list[dict[str, Any]],
    *,
    sample_size: int = 400,
) -> dict[str, Any]:
    from rdkit import Chem, DataStructs
    from rdkit.Chem import rdFingerprintGenerator

    if not rows or not calibration_records:
        return {
            "sample_count": 0,
            "nearest_similarity": _distribution_summary([]),
        }

    generator = rdFingerprintGenerator.GetMorganGenerator(
        radius=2,
        fpSize=2048,
    )
    reference_fps = []
    for record in calibration_records:
        mol = Chem.MolFromSmiles(str(record["smiles"]))
        if mol is not None:
            reference_fps.append(generator.GetFingerprint(mol))

    sample = sorted(
        rows,
        key=lambda row: _stable_hash(str(row["candidate_id"])),
    )[:sample_size]
    nearest: list[float] = []
    for row in sample:
        mol = Chem.MolFromSmiles(str(row["smiles"]))
        if mol is None or not reference_fps:
            continue
        fp = generator.GetFingerprint(mol)
        similarities = DataStructs.BulkTanimotoSimilarity(
            fp,
            reference_fps,
        )
        nearest.append(max(float(value) for value in similarities))

    return {
        "sample_count": len(nearest),
        "nearest_similarity": _distribution_summary(nearest),
        "fraction_below_0_45": (
            sum(value < 0.45 for value in nearest) / len(nearest)
            if nearest
            else None
        ),
        "fraction_below_0_25": (
            sum(value < 0.25 for value in nearest) / len(nearest)
            if nearest
            else None
        ),
    }


def _generate_structural_generation(
    parents: list[str],
    *,
    generation: int,
    config: dict[str, Any],
    seen: set[str],
    accepted_so_far: int = 0,
) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    ht = config["high_throughput"]
    mutation = config["mutation"]
    operators = [str(value) for value in mutation["operators"]]
    allowed_atomic_numbers = [
        int(value)
        for value in mutation["allowed_atomic_numbers"]
    ]
    allowed_families = _family_policy(config)
    per_parent = int(ht["children_per_parent"])
    generation_target = int(ht["generation_unique_target"])
    total_target = int(ht["total_unique_target"])

    generated_attempt_count = 0
    sanitizer_success_count = 0
    duplicate_count = 0
    descriptor_rejection_count = 0
    family_rejection_count = 0
    operator_counts: Counter[str] = Counter()
    family_counts: Counter[str] = Counter()
    rows: list[dict[str, Any]] = []

    for parent in parents:
        mutations = generate_mutations(
            parent,
            generation=generation,
            count=per_parent,
            seed=str(config["search"]["random_seed"]),
            allowed_atomic_numbers=allowed_atomic_numbers,
            operators=operators,
        )
        generated_attempt_count += per_parent
        sanitizer_success_count += len(mutations)

        for mutation_row in mutations:
            smiles = mutation_row.smiles
            if smiles in seen:
                duplicate_count += 1
                continue

            descriptor, failure = _descriptor_check(
                smiles,
                config["candidate_constraints"],
            )
            if failure is not None or descriptor is None:
                descriptor_rejection_count += 1
                seen.add(smiles)
                continue

            bucket = structural_bucket(smiles)
            family = bucket[0]
            if family not in allowed_families:
                family_rejection_count += 1
                seen.add(smiles)
                continue

            seen.add(smiles)
            row = _candidate_record(
                smiles,
                parent_smiles=mutation_row.parent_smiles,
                operator=mutation_row.operator,
                generation=generation,
                descriptor=descriptor,
                family=family,
                bucket=bucket,
            )
            rows.append(row)
            operator_counts[mutation_row.operator] += 1
            family_counts[family] += 1

            if len(rows) >= generation_target:
                break
            if accepted_so_far + len(rows) >= total_target:
                break
        if (
            len(rows) >= generation_target
            or accepted_so_far + len(rows) >= total_target
        ):
            break

    return rows, {
        "generation": generation,
        "parent_count": len(parents),
        "raw_mutation_attempt_count": generated_attempt_count,
        "sanitizer_success_count": sanitizer_success_count,
        "new_structural_candidate_count": len(rows),
        "duplicate_or_seen_count": duplicate_count,
        "descriptor_rejection_count": descriptor_rejection_count,
        "family_rejection_count": family_rejection_count,
        "operator_counts": dict(sorted(operator_counts.items())),
        "family_counts": dict(sorted(family_counts.items())),
    }


def _prescreen_csv_rows(
    rows: list[dict[str, Any]],
    *,
    metrics_key: str,
) -> list[dict[str, Any]]:
    output: list[dict[str, Any]] = []
    for row in rows:
        metrics = row.get(metrics_key) or {}
        output.append(
            {
                "candidate_id": row["candidate_id"],
                "smiles": row["smiles"],
                "generation": row.get("generation"),
                "family": row["provenance"][0]["family"],
                "mutation_operator": row.get("mutation_operator"),
                "property_priority_score": metrics.get(
                    "property_priority_score"
                ),
                "delta_h_q25_ratio_vs_water": metrics.get(
                    "delta_h_q25_ratio_vs_water"
                ),
                "delta_h_median_ratio_vs_water": metrics.get(
                    "delta_h_median_ratio_vs_water"
                ),
                "positive_delta_h_fraction": metrics.get(
                    "positive_delta_h_fraction"
                ),
                "storage_density_kg_m3": metrics.get(
                    "storage_density_kg_m3"
                ),
                "boiling_temperature_k": metrics.get(
                    "boiling_temperature_k"
                ),
            }
        )
    return output


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--config",
        type=Path,
        default=HERE / "config.yaml",
    )
    parser.add_argument(
        "--calibration-summary",
        type=Path,
        default=V5B / "property-v5b2-results" / "summary.json",
    )
    parser.add_argument(
        "--robustness-summary",
        type=Path,
        default=V5B / "property-v5b21-results" / "summary.json",
    )
    parser.add_argument(
        "--v5c1-generated",
        type=Path,
        default=V5C / "search-v5c1-results" / "generated_candidates.json",
    )
    parser.add_argument(
        "--v5c2-generated",
        type=Path,
        default=V5C / "search-v5c2-results" / "generated_candidates.json",
    )
    parser.add_argument(
        "--v5c2-prescreen",
        type=Path,
        default=V5C / "search-v5c2-results" / "prescreen_results.json",
    )
    parser.add_argument(
        "--base-config",
        type=Path,
        default=ENTRY / "config.yaml",
    )
    parser.add_argument(
        "--physical-study",
        type=Path,
        default=ENTRY / "physical-sensitivity.yaml",
    )
    parser.add_argument("--workers", type=int, default=None)
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("search-v5d1-results"),
    )
    args = parser.parse_args()

    if args.output_dir.exists() and any(args.output_dir.iterdir()):
        parser.error("output directory must be empty; use a new directory")
    for required in (
        args.config,
        args.calibration_summary,
        args.robustness_summary,
        args.base_config,
        args.physical_study,
    ):
        if not required.is_file():
            parser.error(f"missing required input: {required}")

    config = yaml.safe_load(args.config.read_text(encoding="utf-8"))
    try:
        _validate_config(config)
    except (KeyError, TypeError, ValueError) as exc:
        parser.error(f"invalid V5d configuration: {exc}")
    calibration_summary = json.loads(
        args.calibration_summary.read_text(encoding="utf-8")
    )
    robustness_summary = json.loads(
        args.robustness_summary.read_text(encoding="utf-8")
    )
    base = yaml.safe_load(args.base_config.read_text(encoding="utf-8"))
    physical = yaml.safe_load(
        args.physical_study.read_text(encoding="utf-8")
    )

    if bool(config["search"].get("require_v5c_ready", True)):
        if not bool(robustness_summary.get("study_complete")):
            parser.error("robustness summary is incomplete")
        if not bool(robustness_summary.get("v5c_ready")):
            parser.error("robustness summary does not authorize V5d-1")

    workers = (
        int(args.workers)
        if args.workers is not None
        else int(config["search"].get("workers", 16))
    )
    if workers <= 0:
        parser.error("--workers must be positive")

    calibration_records = _load_calibration_records(calibration_summary)
    source_candidate_count = robustness_summary.get("source_candidate_count")
    if (
        isinstance(source_candidate_count, int)
        and source_candidate_count != len(calibration_records)
    ):
        parser.error(
            "robustness/calibration candidate-count mismatch: "
            f"{source_candidate_count} != {len(calibration_records)}"
        )

    entry_calibration = calibrate_metric(
        calibration_records,
        error_key="entry_error",
        coverage=float(config["uncertainty"]["coverage_target"]),
        k_neighbors=int(config["uncertainty"]["k_neighbors"]),
        similarity_floor=float(config["uncertainty"]["similarity_floor"]),
    )
    if entry_calibration.get("conformal_factor") is None:
        parser.error("could not construct entry uncertainty model")
    calibration_atomic_numbers = _calibration_atomic_numbers(
        calibration_records
    )

    args.output_dir.mkdir(parents=True, exist_ok=True)

    water_entry = _water_reference(base, physical, config)
    water_mass_kg = float(water_entry["coolant_used_kg"])
    water_full = _water_property_reference(config)
    coarse_config = _coarse_config(config)
    water_coarse = _water_property_reference(coarse_config)

    parents = _initial_parents(
        calibration_records,
        config,
        args.v5c2_prescreen,
    )
    if not parents:
        parser.error("no initial molecular parents available")

    seen: set[str] = set(parents)
    seen.update(
        canonical_smiles(str(record["smiles"]))
        for record in calibration_records
    )
    prior_v5c1 = _load_json_candidates(args.v5c1_generated)
    prior_v5c2 = _load_json_candidates(args.v5c2_generated)
    seen.update(prior_v5c1)
    seen.update(prior_v5c2)

    ht = config["high_throughput"]
    generation_count = int(ht["generations"])
    total_unique_target = int(ht["total_unique_target"])
    coarse_budget = int(ht["coarse_eos_budget_per_generation"])
    beam_width = int(ht["beam_width"])
    diversity_fraction = float(ht["diversity_fraction"])

    all_structural: list[dict[str, Any]] = []
    all_coarse_pass: list[dict[str, Any]] = []
    coarse_rejections: list[dict[str, Any]] = []
    generation_reports: list[dict[str, Any]] = []

    for generation in range(1, generation_count + 1):
        structural_rows, generation_report = _generate_structural_generation(
            parents,
            generation=generation,
            config=config,
            seen=seen,
            accepted_so_far=len(all_structural),
        )
        all_structural.extend(structural_rows)

        coarse_selected = _structural_select(
            structural_rows,
            min(coarse_budget, len(structural_rows)),
        )
        coarse_tasks = [
            {
                "candidate": row,
                "config": coarse_config,
                "water_properties": water_coarse,
            }
            for row in coarse_selected
        ]
        coarse_rows = _run_parallel(
            _coarse_candidate,
            coarse_tasks,
            workers=workers,
        )
        coarse_pass = [
            row
            for row in coarse_rows
            if row["coarse_status"] == "pass"
        ]
        coarse_rejected = [
            row
            for row in coarse_rows
            if row["coarse_status"] != "pass"
        ]
        all_coarse_pass.extend(coarse_pass)
        coarse_rejections.extend(coarse_rejected)

        next_beam = _beam_select(
            coarse_pass,
            width=beam_width,
            diversity_fraction=diversity_fraction,
        )
        parents = [str(row["smiles"]) for row in next_beam]

        generation_report.update(
            {
                "coarse_selected_count": len(coarse_selected),
                "coarse_pass_count": len(coarse_pass),
                "coarse_rejection_count": len(coarse_rejected),
                "next_parent_count": len(parents),
                "best_coarse_property_priority_score": (
                    max(
                        (_score_value(row) for row in coarse_pass),
                        default=None,
                    )
                ),
                "global_seen_count": len(seen),
            }
        )
        generation_reports.append(generation_report)
        write_summary(
            args.output_dir / "progress.json",
            {
                "study_version": "v5d-1",
                "completed_generation": generation,
                "accepted_structural_count": len(all_structural),
                "global_seen_count": len(seen),
                "current_parent_count": len(parents),
                "current_parents": parents,
                "generations": generation_reports,
            },
        )
        print(
            json.dumps(
                {
                    "generation": generation,
                    "new_structures": len(structural_rows),
                    "coarse_evaluated": len(coarse_selected),
                    "coarse_pass": len(coarse_pass),
                    "next_parents": len(parents),
                    "global_seen": len(seen),
                    "best_property_score": generation_report[
                        "best_coarse_property_priority_score"
                    ],
                }
            ),
            flush=True,
        )

        if not parents:
            break
        if len(all_structural) >= total_unique_target:
            break

    _write_jsonl_gz(
        args.output_dir / "structural_candidates.jsonl.gz",
        all_structural,
    )
    write_summary(
        args.output_dir / "generation_report.json",
        {
            "generations": generation_reports,
            "initial_parent_count": len(
                _initial_parents(
                    calibration_records,
                    config,
                    args.v5c2_prescreen,
                )
            ),
            "prior_v5c1_excluded_count": len(prior_v5c1),
            "prior_v5c2_excluded_count": len(prior_v5c2),
        },
    )

    coarse_pass_unique = {
        str(row["smiles"]): row
        for row in all_coarse_pass
    }
    coarse_pass_rows = list(coarse_pass_unique.values())
    full_budget = int(ht["full_prescreen_budget"])
    full_selected = _beam_select(
        coarse_pass_rows,
        width=min(full_budget, len(coarse_pass_rows)),
        diversity_fraction=diversity_fraction,
    )

    full_tasks = [
        {
            "candidate": row,
            "config": config,
            "calibration_records": calibration_records,
            "calibration_atomic_numbers": sorted(
                calibration_atomic_numbers
            ),
            "entry_calibration": entry_calibration,
            "water_properties": water_full,
        }
        for row in full_selected
    ]
    full_rows = _run_parallel(
        _prescreen_candidate,
        full_tasks,
        workers=workers,
    )
    full_pass = [
        row
        for row in full_rows
        if row["status"] == "prescreen_pass"
    ]
    full_rejected = [
        row
        for row in full_rows
        if row["status"] != "prescreen_pass"
    ]

    rankable = sorted(
        [
            row
            for row in full_pass
            if row["lane"] == "rankable"
        ],
        key=_prescreen_sort_key,
    )
    exploratory = sorted(
        [
            row
            for row in full_pass
            if row["lane"] == "exploratory_domain_expansion"
        ],
        key=_prescreen_sort_key,
    )

    rankable_selected = rankable[
        : int(ht["entry_rankable_budget"])
    ]
    exploratory_selected = exploratory[
        : int(ht["entry_exploratory_budget"])
    ]
    entry_tasks = [
        {
            "row": row,
            "config": config,
            "base": base,
            "physical": physical,
            "water_mass_kg": water_mass_kg,
        }
        for row in rankable_selected + exploratory_selected
    ]
    entry_rows = _run_parallel(
        _entry_candidate,
        entry_tasks,
        workers=workers,
    )

    ranked = [
        row
        for row in entry_rows
        if row["status"] == "ranked"
    ]
    ranked.sort(
        key=lambda row: (
            float(row["conservative_coolant_kg"]),
            float(row["predicted_coolant_kg"]),
            str(row["candidate_id"]),
        )
    )
    for index, row in enumerate(ranked, start=1):
        row["rank"] = index

    exploratory_evaluated = [
        row
        for row in entry_rows
        if row["status"] == "exploratory_evaluated"
    ]
    exploratory_evaluated.sort(
        key=lambda row: (
            float(row["predicted_coolant_kg"]),
            str(row["candidate_id"]),
        )
    )
    entry_rejected = [
        row
        for row in entry_rows
        if row["status"] == "entry_rejected"
    ]

    conservative_winners = [
        row
        for row in ranked
        if row.get("beats_water_conservative") is True
    ]
    exploratory_below_water = [
        row
        for row in exploratory_evaluated
        if row.get("predicted_below_water") is True
    ]

    coarse_rejection_counts = Counter(
        str(row.get("coarse_rejection_stage", "unknown"))
        for row in coarse_rejections
    )
    full_rejection_counts = Counter(
        str(row.get("rejection_stage", "unknown"))
        for row in full_rejected + entry_rejected
    )
    family_counts = Counter(
        str(row["provenance"][0]["family"])
        for row in all_structural
    )
    operator_counts = Counter(
        str(row["mutation_operator"])
        for row in all_structural
    )
    mw_values = [
        float(row["descriptor"]["molecular_weight_g_mol"])
        for row in all_structural
    ]

    _write_csv(
        args.output_dir / "coarse_pass.csv",
        _prescreen_csv_rows(
            coarse_pass_rows,
            metrics_key="coarse_property_screen",
        ),
    )
    _write_csv(
        args.output_dir / "full_prescreen.csv",
        _prescreen_csv_rows(
            full_pass,
            metrics_key="property_screen",
        ),
    )
    _write_csv(
        args.output_dir / "ranking.csv",
        [
            {
                "rank": row["rank"],
                "candidate_id": row["candidate_id"],
                "smiles": row["smiles"],
                "generation": row.get("generation"),
                "family": row["provenance"][0]["family"],
                "mutation_operator": row.get("mutation_operator"),
                "property_priority_score": row[
                    "property_screen"
                ]["property_priority_score"],
                "predicted_coolant_kg": row["predicted_coolant_kg"],
                "entry_relative_uncertainty": row[
                    "entry_relative_uncertainty"
                ],
                "conservative_coolant_kg": row[
                    "conservative_coolant_kg"
                ],
                "conservative_ratio_vs_water": row[
                    "conservative_ratio_vs_water"
                ],
                "beats_water_conservative": row[
                    "beats_water_conservative"
                ],
            }
            for row in ranked
        ],
    )
    _write_csv(
        args.output_dir / "exploratory.csv",
        [
            {
                "candidate_id": row["candidate_id"],
                "smiles": row["smiles"],
                "generation": row.get("generation"),
                "family": row["provenance"][0]["family"],
                "mutation_operator": row.get("mutation_operator"),
                "lane_reason": row.get("lane_reason"),
                "property_priority_score": row[
                    "property_screen"
                ]["property_priority_score"],
                "predicted_coolant_kg": row["predicted_coolant_kg"],
                "predicted_ratio_vs_water": row[
                    "predicted_ratio_vs_water"
                ],
                "predicted_below_water": row["predicted_below_water"],
            }
            for row in exploratory_evaluated
        ],
    )

    summary = {
        "study_version": "v5d-1",
        "study_complete": True,
        "water_reference_kg": water_mass_kg,
        "search_funnel": {
            "initial_parent_count": generation_reports[0][
                "parent_count"
            ]
            if generation_reports
            else len(parents),
            "prior_v5c1_excluded_count": len(prior_v5c1),
            "prior_v5c2_excluded_count": len(prior_v5c2),
            "structural_unique_count": len(all_structural),
            "coarse_evaluated_count": sum(
                row["coarse_selected_count"]
                for row in generation_reports
            ),
            "coarse_pass_unique_count": len(coarse_pass_rows),
            "full_prescreen_selected_count": len(full_selected),
            "full_prescreen_pass_count": len(full_pass),
            "rankable_full_prescreen_count": len(rankable),
            "exploratory_full_prescreen_count": len(exploratory),
            "entry_rankable_evaluated_count": len(
                [
                    row
                    for row in entry_rows
                    if row["status"] in {"ranked", "entry_rejected"}
                    and row.get("lane") == "rankable"
                ]
            ),
            "entry_exploratory_evaluated_count": len(
                [
                    row
                    for row in entry_rows
                    if row["status"]
                    in {"exploratory_evaluated", "entry_rejected"}
                    and row.get("lane")
                    == "exploratory_domain_expansion"
                ]
            ),
        },
        "generation_reports": generation_reports,
        "structural_space": {
            "family_counts": dict(sorted(family_counts.items())),
            "operator_counts": dict(sorted(operator_counts.items())),
            "molecular_weight_g_mol": _distribution_summary(mw_values),
            "approximate_fingerprint_diversity": _approximate_diversity(
                all_structural
            ),
            "approximate_novelty_vs_v5b": (
                _approximate_reference_novelty(
                    all_structural,
                    calibration_records,
                )
            ),
        },
        "rejections": {
            "coarse": dict(sorted(coarse_rejection_counts.items())),
            "full_or_entry": dict(sorted(full_rejection_counts.items())),
        },
        "screening": {
            "ranked_candidate_count": len(ranked),
            "conservative_water_winner_count": len(
                conservative_winners
            ),
            "exploratory_evaluated_count": len(
                exploratory_evaluated
            ),
            "exploratory_predicted_below_water_count": len(
                exploratory_below_water
            ),
        },
        "best_rankable_candidate": ranked[0] if ranked else None,
        "conservative_water_winners": conservative_winners,
        "best_exploratory_candidate": (
            exploratory_evaluated[0]
            if exploratory_evaluated
            else None
        ),
        "exploratory_predicted_below_water": exploratory_below_water,
        "interpretation": (
            "V5d-1 is a high-throughput evolutionary screening search. "
            "Only rankable candidates with calibrated-domain uncertainty may "
            "count as conservative Water winners. Exploratory predictions "
            "remain hypotheses for domain expansion."
        ),
    }
    write_summary(args.output_dir / "summary.json", summary)
    write_summary(args.output_dir / "water_reference.json", water_entry)
    write_summary(
        args.output_dir / "full_prescreen_results.json",
        {
            "passes": full_pass,
            "rejections": full_rejected,
        },
    )
    write_summary(
        args.output_dir / "entry_results.json",
        {
            "ranked": ranked,
            "exploratory": exploratory_evaluated,
            "rejections": entry_rejected,
        },
    )

    try:
        git_commit = subprocess.check_output(
            ["git", "rev-parse", "HEAD"],
            cwd=HERE,
            text=True,
        ).strip()
        git_dirty = bool(
            subprocess.check_output(
                ["git", "status", "--porcelain"],
                cwd=HERE,
                text=True,
            ).strip()
        )
    except (OSError, subprocess.CalledProcessError):
        git_commit = None
        git_dirty = None

    manifest = {
        "study_version": "v5d-1",
        "git_commit": git_commit,
        "git_dirty": git_dirty,
        "python": platform.python_version(),
        "packages": {
            package: version(package)
            for package in (
                "feos",
                "rdkit",
                "CoolProp",
                "pymsis",
                "numpy",
                "PyYAML",
            )
        },
        "inputs": {
            "config_sha256": _sha256(args.config),
            "calibration_summary_sha256": _sha256(
                args.calibration_summary
            ),
            "robustness_summary_sha256": _sha256(
                args.robustness_summary
            ),
            "v5c1_generated_sha256": (
                _sha256(args.v5c1_generated)
                if args.v5c1_generated.is_file()
                else None
            ),
            "v5c2_generated_sha256": (
                _sha256(args.v5c2_generated)
                if args.v5c2_generated.is_file()
                else None
            ),
            "v5c2_prescreen_sha256": (
                _sha256(args.v5c2_prescreen)
                if args.v5c2_prescreen.is_file()
                else None
            ),
        },
        "search_contract": {
            "graph_generation": (
                "deterministic hash-seeded RDKit graph mutations"
            ),
            "hierarchy": (
                "structure -> stratified coarse EOS -> full T/P prescreen "
                "-> calibrated-domain gate -> limited entry evaluation"
            ),
            "rankable_water_win": (
                "requires finite uncertainty-conservative coolant mass below "
                "the same-run Water reference"
            ),
        },
    }
    write_summary(args.output_dir / "manifest.json", manifest)

    print(
        json.dumps(
            {
                "study_complete": True,
                "water_reference_kg": water_mass_kg,
                "search_funnel": summary["search_funnel"],
                "structural_space": summary["structural_space"],
                "rejections": summary["rejections"],
                "screening": summary["screening"],
                "best_rankable_candidate": (
                    {
                        "candidate_id": ranked[0]["candidate_id"],
                        "smiles": ranked[0]["smiles"],
                        "family": ranked[0]["provenance"][0]["family"],
                        "generation": ranked[0].get("generation"),
                        "predicted_coolant_kg": ranked[0][
                            "predicted_coolant_kg"
                        ],
                        "entry_relative_uncertainty": ranked[0][
                            "entry_relative_uncertainty"
                        ],
                        "conservative_coolant_kg": ranked[0][
                            "conservative_coolant_kg"
                        ],
                        "conservative_ratio_vs_water": ranked[0][
                            "conservative_ratio_vs_water"
                        ],
                    }
                    if ranked
                    else None
                ),
                "best_exploratory_candidate": (
                    {
                        "candidate_id": exploratory_evaluated[0][
                            "candidate_id"
                        ],
                        "smiles": exploratory_evaluated[0]["smiles"],
                        "family": exploratory_evaluated[0][
                            "provenance"
                        ][0]["family"],
                        "generation": exploratory_evaluated[0].get(
                            "generation"
                        ),
                        "predicted_coolant_kg": exploratory_evaluated[0][
                            "predicted_coolant_kg"
                        ],
                        "predicted_ratio_vs_water": exploratory_evaluated[0][
                            "predicted_ratio_vs_water"
                        ],
                    }
                    if exploratory_evaluated
                    else None
                ),
            },
            indent=2,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
