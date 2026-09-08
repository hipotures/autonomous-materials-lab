#!/usr/bin/env python3
"""V5c-2: property-targeted search with controlled domain expansion."""
from __future__ import annotations

import argparse
import json
import math
import platform
import subprocess
import sys
from concurrent.futures import ProcessPoolExecutor, as_completed
from importlib.metadata import version
from pathlib import Path
from typing import Any

import yaml

HERE = Path(__file__).resolve().parent
ENTRY = HERE.parent / "entry-evaluator"
V5B = HERE.parent / "property-predictor-v5b"
for path in (ENTRY, V5B):
    if str(path) not in sys.path:
        sys.path.insert(0, str(path))

from applicability import (  # noqa: E402
    assess_domain,
    calibrate_metric,
    canonical_smiles,
    domain_to_dict,
    predict_relative_uncertainty,
)
from evaluator import evaluate, write_summary  # noqa: E402
from generator_v5c2 import generate_candidates  # noqa: E402
from property_provider import build_property_provider  # noqa: E402
from run_search import (  # noqa: E402
    SUCCESS_STATUSES,
    _base_entry_config,
    _coolant_config,
    _descriptor_check,
    _load_calibration_records,
    _sha256,
    _stable_candidate_id,
    _water_reference,
    _write_csv,
    conservative_score,
)
from thermo_prescreen import (  # noqa: E402
    compare_to_reference,
    passes_property_gate,
    property_screen,
)


def _atomic_numbers(smiles: str) -> set[int]:
    from rdkit import Chem

    molecule = Chem.MolFromSmiles(smiles)
    if molecule is None:
        raise ValueError(f"invalid SMILES: {smiles}")
    return {
        int(atom.GetAtomicNum())
        for atom in molecule.GetAtoms()
    }


def _calibration_atomic_numbers(
    records: list[dict[str, Any]],
) -> set[int]:
    result: set[int] = set()
    for record in records:
        result.update(_atomic_numbers(str(record["smiles"])))
    return result


def _candidate_pool(
    config: dict[str, Any],
    calibration_records: list[dict[str, Any]],
    prior_smiles: set[str] | None = None,
) -> tuple[list[dict[str, Any]], dict[str, int]]:
    generated = generate_candidates(config)
    known = {
        canonical_smiles(str(record["smiles"]))
        for record in calibration_records
    }

    prior_smiles = prior_smiles or set()
    by_smiles: dict[str, dict[str, Any]] = {}
    invalid_count = 0
    known_count = 0
    prior_count = 0
    duplicate_count = 0

    for raw in generated:
        try:
            smiles = canonical_smiles(raw.smiles)
        except Exception:
            invalid_count += 1
            continue
        if smiles in known:
            known_count += 1
            continue
        if smiles in prior_smiles:
            prior_count += 1
            continue

        provenance = {
            "family": raw.family,
            "operation": raw.operation,
            "parameters": raw.parameters,
            "raw_smiles": raw.smiles,
        }
        if smiles in by_smiles:
            duplicate_count += 1
            by_smiles[smiles]["provenance"].append(provenance)
            continue

        by_smiles[smiles] = {
            "candidate_id": _stable_candidate_id(smiles),
            "smiles": smiles,
            "provenance": [provenance],
        }

    candidates = [
        by_smiles[key]
        for key in sorted(by_smiles)
    ]
    maximum = int(
        config["search"].get(
            "maximum_candidates_after_dedup",
            800,
        )
    )
    if len(candidates) > maximum:
        raise ValueError(
            f"generated {len(candidates)} candidates after dedup, "
            f"exceeding configured maximum {maximum}"
        )

    return candidates, {
        "raw_generated_count": len(generated),
        "invalid_generated_count": invalid_count,
        "known_reference_removed_count": known_count,
        "prior_search_removed_count": prior_count,
        "duplicate_removed_count": duplicate_count,
        "novel_candidate_count": len(candidates),
    }


def _load_prior_smiles(path: Path | None) -> set[str]:
    if path is None or not path.is_file():
        return set()
    payload = json.loads(path.read_text(encoding="utf-8"))
    result: set[str] = set()
    for row in payload.get("candidates", []):
        smiles = row.get("smiles")
        if smiles:
            result.add(canonical_smiles(str(smiles)))
    return result


def _water_property_reference(
    config: dict[str, Any],
) -> dict[str, Any]:
    coolant = _coolant_config(
        config,
        smiles=None,
        name="water-property-reference",
        reference_coolprop_name="Water",
    )
    provider = build_property_provider(coolant)
    return property_screen(provider, config)


def _family_is_configured(
    candidate: dict[str, Any],
    config: dict[str, Any],
) -> bool:
    family = str(candidate["provenance"][0]["family"])
    configured = {
        str(value)
        for key in ("rankable_families", "exploratory_families")
        for value in config["domain_expansion"][key]
    }
    return family in configured


def _family_is_rankable(
    candidate: dict[str, Any],
    config: dict[str, Any],
) -> bool:
    family = str(candidate["provenance"][0]["family"])
    return family in {
        str(value)
        for value in config["domain_expansion"]["rankable_families"]
    }


def _prescreen_candidate(task: dict[str, Any]) -> dict[str, Any]:
    candidate = task["candidate"]
    config = task["config"]
    calibration_records = task["calibration_records"]
    calibration_atomic_numbers = set(
        int(value)
        for value in task["calibration_atomic_numbers"]
    )
    entry_calibration = task["entry_calibration"]
    water_properties = task["water_properties"]

    smiles = str(candidate["smiles"])
    descriptor, descriptor_failure = _descriptor_check(
        smiles,
        config["candidate_constraints"],
    )
    if descriptor_failure is not None:
        return {
            **candidate,
            "status": "rejected",
            "rejection_stage": "descriptor_constraint",
            "rejection_reason": descriptor_failure,
        }

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
            "status": "rejected",
            "rejection_stage": "gc_not_representable",
            "rejection_reason": str(exc),
            "descriptor": descriptor,
        }

    try:
        metrics = property_screen(provider, config)
    except Exception as exc:
        return {
            **candidate,
            "status": "rejected",
            "rejection_stage": "property_prescreen_failure",
            "rejection_reason": str(exc),
            "descriptor": descriptor,
        }

    if (
        bool(config["storage"].get("require_liquid", True))
        and metrics.get("storage_phase") != "liquid"
    ):
        return {
            **candidate,
            "status": "rejected",
            "rejection_stage": "storage_not_liquid",
            "rejection_reason": (
                f"predicted phase={metrics.get('storage_phase')}"
            ),
            "descriptor": descriptor,
            "property_screen": metrics,
        }

    if (
        bool(
            config["storage"].get(
                "require_saturation_at_storage_pressure",
                True,
            )
        )
        and not bool(metrics.get("saturation_supported"))
    ):
        return {
            **candidate,
            "status": "rejected",
            "rejection_stage": "saturation_unsupported",
            "rejection_reason": metrics.get(
                "saturation_failure_reason"
            ),
            "descriptor": descriptor,
            "property_screen": metrics,
        }

    targets = compare_to_reference(
        metrics,
        water_properties,
    )
    metrics = {
        **metrics,
        **targets,
    }
    property_pass, property_failure = passes_property_gate(
        metrics,
        config,
    )
    if not property_pass:
        return {
            **candidate,
            "status": "rejected",
            "rejection_stage": "property_target",
            "rejection_reason": property_failure,
            "descriptor": descriptor,
            "property_screen": metrics,
        }

    if not _family_is_configured(candidate, config):
        return {
            **candidate,
            "status": "rejected",
            "rejection_stage": "unconfigured_family",
            "rejection_reason": str(
                candidate["provenance"][0]["family"]
            ),
            "descriptor": descriptor,
            "property_screen": metrics,
        }

    domain = assess_domain(
        smiles,
        calibration_records,
        in_domain_similarity=float(
            config["uncertainty"]["in_domain_similarity"]
        ),
        edge_similarity=float(
            config["uncertainty"]["edge_similarity"]
        ),
        minimum_neighbors=int(
            config["uncertainty"]["minimum_neighbors"]
        ),
    )
    candidate_elements = set(
        int(value)
        for value in descriptor["atomic_numbers"]
    )
    element_supported = candidate_elements.issubset(
        calibration_atomic_numbers
    )
    family_supported = _family_is_rankable(
        candidate,
        config,
    )

    uncertainty = None
    if (
        domain.status == "in_domain"
        and element_supported
        and family_supported
    ):
        uncertainty = predict_relative_uncertainty(
            smiles,
            calibration_records,
            entry_calibration,
            error_key="entry_error",
            k_neighbors=int(
                config["uncertainty"]["k_neighbors"]
            ),
            similarity_floor=float(
                config["uncertainty"]["similarity_floor"]
            ),
        )

    multiplier = float(
        config["search"].get(
            "conservative_score_multiplier",
            1.0,
        )
    )
    finite_bound = (
        isinstance(uncertainty, (int, float))
        and math.isfinite(float(uncertainty))
        and float(uncertainty) >= 0.0
        and multiplier * float(uncertainty) < 1.0
    )

    if (
        domain.status == "in_domain"
        and element_supported
        and family_supported
        and finite_bound
    ):
        lane = "rankable"
        lane_reason = "validated_elements_and_in_domain"
    elif bool(
        config["search"].get(
            "allow_exploratory_domain_expansion",
            True,
        )
    ):
        lane = "exploratory_domain_expansion"
        reasons: list[str] = []
        if not element_supported:
            reasons.append("new_element_relative_to_calibration")
        if not family_supported:
            reasons.append("new_family_relative_to_calibration")
        if domain.status != "in_domain":
            reasons.append(f"domain_{domain.status}")
        if domain.status == "in_domain" and not finite_bound:
            reasons.append("uncertainty_not_finitely_rankable")
        lane_reason = ",".join(reasons) or "exploratory"
    else:
        return {
            **candidate,
            "status": "rejected",
            "rejection_stage": "applicability_domain",
            "rejection_reason": (
                f"domain={domain.status}; "
                f"element_supported={element_supported}"
            ),
            "descriptor": descriptor,
            "property_screen": metrics,
            "domain": domain_to_dict(domain),
        }

    return {
        **candidate,
        "status": "prescreen_pass",
        "lane": lane,
        "lane_reason": lane_reason,
        "descriptor": descriptor,
        "property_screen": metrics,
        "domain": domain_to_dict(domain),
        "calibration_element_supported": element_supported,
        "calibration_family_supported": family_supported,
        "entry_relative_uncertainty": (
            float(uncertainty)
            if isinstance(uncertainty, (int, float))
            else None
        ),
    }


def _entry_candidate(task: dict[str, Any]) -> dict[str, Any]:
    row = task["row"]
    config = task["config"]
    base = task["base"]
    physical = task["physical"]
    water_mass_kg = float(task["water_mass_kg"])

    smiles = str(row["smiles"])
    coolant = _coolant_config(
        config,
        smiles=smiles,
        name=str(row["candidate_id"]),
    )
    entry_config = _base_entry_config(
        base,
        physical,
        config,
    )
    entry_config["name"] = (
        f"v5c2-{row['candidate_id']}"
    )
    entry_config["coolant"] = coolant

    try:
        entry = evaluate(entry_config)
    except Exception as exc:
        return {
            **row,
            "status": "entry_rejected",
            "rejection_stage": "entry_exception",
            "rejection_reason": str(exc),
        }

    if entry.get("status") not in SUCCESS_STATUSES:
        return {
            **row,
            "status": "entry_rejected",
            "rejection_stage": "entry_failure",
            "rejection_reason": (
                f"{entry.get('status')}: "
                f"{entry.get('failure_reason')}"
            ),
            "entry": entry,
        }

    predicted = entry.get("coolant_used_kg")
    if (
        not isinstance(predicted, (int, float))
        or not math.isfinite(float(predicted))
        or float(predicted) <= 0.0
    ):
        return {
            **row,
            "status": "entry_rejected",
            "rejection_stage": "entry_invalid_score",
            "rejection_reason": "invalid coolant_used_kg",
            "entry": entry,
        }
    predicted = float(predicted)

    result = {
        **row,
        "status": (
            "ranked"
            if row["lane"] == "rankable"
            else "exploratory_evaluated"
        ),
        "predicted_coolant_kg": predicted,
        "water_reference_kg": water_mass_kg,
        "predicted_ratio_vs_water": predicted / water_mass_kg,
        "predicted_below_water": predicted < water_mass_kg,
        "entry": entry,
    }

    if row["lane"] == "rankable":
        uncertainty = float(row["entry_relative_uncertainty"])
        multiplier = float(
            config["search"].get(
                "conservative_score_multiplier",
                1.0,
            )
        )
        conservative = conservative_score(
            predicted,
            uncertainty,
            multiplier=multiplier,
        )
        result.update(
            {
                "conservative_coolant_kg": conservative,
                "conservative_ratio_vs_water": (
                    conservative / water_mass_kg
                ),
                "beats_water_conservative": (
                    conservative < water_mass_kg
                ),
            }
        )
    else:
        result.update(
            {
                "conservative_coolant_kg": None,
                "conservative_ratio_vs_water": None,
                "beats_water_conservative": None,
            }
        )

    return result


def _prescreen_sort_key(row: dict[str, Any]):
    metrics = row["property_screen"]
    priority = metrics.get("property_priority_score")
    latent = metrics.get("latent_proxy_ratio_vs_water")
    density = metrics.get("storage_density_kg_m3")
    return (
        -float(priority)
        if isinstance(priority, (int, float))
        else float("inf"),
        -float(latent)
        if isinstance(latent, (int, float))
        else float("inf"),
        -float(density)
        if isinstance(density, (int, float))
        else float("inf"),
        str(row["candidate_id"]),
    )


def _run_parallel(function, tasks: list[dict[str, Any]], workers: int):
    if workers <= 1:
        return [function(task) for task in tasks]
    rows: list[dict[str, Any]] = []
    with ProcessPoolExecutor(max_workers=workers) as pool:
        futures = [
            pool.submit(function, task)
            for task in tasks
        ]
        for future in as_completed(futures):
            rows.append(future.result())
    return rows


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--config",
        type=Path,
        default=HERE / "config-v5c2.yaml",
    )
    parser.add_argument(
        "--calibration-summary",
        type=Path,
        default=(
            V5B
            / "property-v5b2-results"
            / "summary.json"
        ),
    )
    parser.add_argument(
        "--robustness-summary",
        type=Path,
        default=(
            V5B
            / "property-v5b21-results"
            / "summary.json"
        ),
    )
    parser.add_argument(
        "--prior-generated-candidates",
        type=Path,
        default=HERE / "search-v5c1-results" / "generated_candidates.json",
        help=(
            "Optional V5c-1 generated-candidate artifact; if present, "
            "those structures are excluded from V5c-2."
        ),
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
        default=Path("search-v5c2-results"),
    )
    args = parser.parse_args()

    if args.output_dir.exists() and any(args.output_dir.iterdir()):
        parser.error(
            "output directory must be empty; use a new directory"
        )
    for path in (
        args.config,
        args.calibration_summary,
        args.robustness_summary,
        args.base_config,
        args.physical_study,
    ):
        if not path.is_file():
            parser.error(f"missing required input: {path}")

    config = yaml.safe_load(
        args.config.read_text(encoding="utf-8")
    )
    calibration_summary = json.loads(
        args.calibration_summary.read_text(encoding="utf-8")
    )
    robustness_summary = json.loads(
        args.robustness_summary.read_text(encoding="utf-8")
    )
    base = yaml.safe_load(
        args.base_config.read_text(encoding="utf-8")
    )
    physical = yaml.safe_load(
        args.physical_study.read_text(encoding="utf-8")
    )

    if bool(config["search"].get("require_v5c_ready", True)):
        if not bool(robustness_summary.get("study_complete")):
            parser.error("robustness summary is incomplete")
        if not bool(robustness_summary.get("v5c_ready")):
            parser.error(
                "robustness summary does not authorize V5c-2"
            )

    workers = (
        int(args.workers)
        if args.workers is not None
        else int(config["search"].get("workers", 16))
    )
    if workers <= 0:
        parser.error("--workers must be positive")

    calibration_records = _load_calibration_records(
        calibration_summary
    )
    source_candidate_count = robustness_summary.get(
        "source_candidate_count"
    )
    if (
        isinstance(source_candidate_count, int)
        and source_candidate_count != len(calibration_records)
    ):
        parser.error(
            "robustness/calibration candidate-count mismatch: "
            f"{source_candidate_count} != "
            f"{len(calibration_records)}"
        )

    calibration_elements = _calibration_atomic_numbers(
        calibration_records
    )
    entry_calibration = calibrate_metric(
        calibration_records,
        error_key="entry_error",
        coverage=float(
            config["uncertainty"]["coverage_target"]
        ),
        k_neighbors=int(
            config["uncertainty"]["k_neighbors"]
        ),
        similarity_floor=float(
            config["uncertainty"]["similarity_floor"]
        ),
    )
    if entry_calibration.get("conformal_factor") is None:
        parser.error(
            "could not construct production entry uncertainty model"
        )

    args.output_dir.mkdir(parents=True, exist_ok=True)

    water_entry = _water_reference(
        base,
        physical,
        config,
    )
    water_mass_kg = float(
        water_entry["coolant_used_kg"]
    )
    water_properties = _water_property_reference(config)

    prior_smiles = _load_prior_smiles(
        args.prior_generated_candidates
    )
    candidates, generation_stats = _candidate_pool(
        config,
        calibration_records,
        prior_smiles,
    )
    write_summary(
        args.output_dir / "generated_candidates.json",
        {
            "generation_stats": generation_stats,
            "candidates": candidates,
        },
    )
    write_summary(
        args.output_dir / "water_property_reference.json",
        water_properties,
    )

    prescreen_tasks = [
        {
            "candidate": candidate,
            "config": config,
            "calibration_records": calibration_records,
            "calibration_atomic_numbers": sorted(
                calibration_elements
            ),
            "entry_calibration": entry_calibration,
            "water_properties": water_properties,
        }
        for candidate in candidates
    ]
    prescreen_rows = _run_parallel(
        _prescreen_candidate,
        prescreen_tasks,
        workers,
    )

    prescreen_pass = [
        row
        for row in prescreen_rows
        if row["status"] == "prescreen_pass"
    ]
    rejected = [
        row
        for row in prescreen_rows
        if row["status"] != "prescreen_pass"
    ]
    write_summary(
        args.output_dir / "prescreen_results.json",
        {
            "passes": prescreen_pass,
            "rejections": rejected,
        },
    )

    rankable = sorted(
        [
            row
            for row in prescreen_pass
            if row["lane"] == "rankable"
        ],
        key=_prescreen_sort_key,
    )
    exploratory = sorted(
        [
            row
            for row in prescreen_pass
            if row["lane"]
            == "exploratory_domain_expansion"
        ],
        key=_prescreen_sort_key,
    )

    rankable_limit = int(
        config["search"]["full_entry_rankable_limit"]
    )
    exploratory_limit = int(
        config["search"]["full_entry_exploratory_limit"]
    )
    rankable_selected = rankable[:rankable_limit]
    exploratory_selected = exploratory[:exploratory_limit]

    entry_tasks = [
        {
            "row": row,
            "config": config,
            "base": base,
            "physical": physical,
            "water_mass_kg": water_mass_kg,
        }
        for row in (
            rankable_selected
            + exploratory_selected
        )
    ]
    entry_rows = _run_parallel(
        _entry_candidate,
        entry_tasks,
        workers,
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
    rejected.extend(entry_rejected)

    rejection_counts: dict[str, int] = {}
    for row in rejected:
        stage = str(
            row.get("rejection_stage", "unknown")
        )
        rejection_counts[stage] = (
            rejection_counts.get(stage, 0) + 1
        )

    conservative_winners = [
        row
        for row in ranked
        if row.get("beats_water_conservative") is True
    ]
    predicted_rankable_below_water = [
        row
        for row in ranked
        if row.get("predicted_below_water") is True
    ]
    exploratory_below_water = [
        row
        for row in exploratory_evaluated
        if row.get("predicted_below_water") is True
    ]

    exploratory_family_summary: dict[str, dict[str, Any]] = {}
    for row in exploratory:
        family = str(row["provenance"][0]["family"])
        metrics = row["property_screen"]
        item = exploratory_family_summary.setdefault(
            family,
            {
                "prescreen_pass_count": 0,
                "best_property_priority_score": None,
                "entry_evaluated_count": 0,
                "best_predicted_coolant_kg": None,
                "best_predicted_ratio_vs_water": None,
            },
        )
        item["prescreen_pass_count"] += 1
        priority = metrics.get("property_priority_score")
        if isinstance(priority, (int, float)):
            current = item["best_property_priority_score"]
            if current is None or float(priority) > float(current):
                item["best_property_priority_score"] = float(priority)

    for row in exploratory_evaluated:
        family = str(row["provenance"][0]["family"])
        item = exploratory_family_summary.setdefault(
            family,
            {
                "prescreen_pass_count": 0,
                "best_property_priority_score": None,
                "entry_evaluated_count": 0,
                "best_predicted_coolant_kg": None,
                "best_predicted_ratio_vs_water": None,
            },
        )
        item["entry_evaluated_count"] += 1
        predicted = float(row["predicted_coolant_kg"])
        current = item["best_predicted_coolant_kg"]
        if current is None or predicted < float(current):
            item["best_predicted_coolant_kg"] = predicted
            item["best_predicted_ratio_vs_water"] = float(
                row["predicted_ratio_vs_water"]
            )

    prescreen_csv: list[dict[str, Any]] = []
    for row in prescreen_pass:
        metrics = row["property_screen"]
        prescreen_csv.append(
            {
                "candidate_id": row["candidate_id"],
                "smiles": row["smiles"],
                "family": row["provenance"][0]["family"],
                "lane": row["lane"],
                "lane_reason": row["lane_reason"],
                "domain_status": row["domain"]["status"],
                "nearest_similarity": row["domain"][
                    "nearest_similarity"
                ],
                "calibration_element_supported": row[
                    "calibration_element_supported"
                ],
                "calibration_family_supported": row[
                    "calibration_family_supported"
                ],
                "property_priority_score": metrics[
                    "property_priority_score"
                ],
                "delta_h_q25_ratio_vs_water": metrics[
                    "delta_h_q25_ratio_vs_water"
                ],
                "delta_h_median_ratio_vs_water": metrics[
                    "delta_h_median_ratio_vs_water"
                ],
                "latent_proxy_ratio_vs_water": metrics[
                    "latent_proxy_ratio_vs_water"
                ],
                "positive_delta_h_fraction": metrics[
                    "positive_delta_h_fraction"
                ],
                "storage_density_kg_m3": metrics[
                    "storage_density_kg_m3"
                ],
                "boiling_temperature_k": metrics[
                    "boiling_temperature_k"
                ],
                "entry_relative_uncertainty": row.get(
                    "entry_relative_uncertainty"
                ),
            }
        )
    _write_csv(
        args.output_dir / "prescreen.csv",
        prescreen_csv,
    )

    ranking_csv = [
        {
            "rank": row["rank"],
            "candidate_id": row["candidate_id"],
            "smiles": row["smiles"],
            "family": row["provenance"][0]["family"],
            "property_priority_score": row[
                "property_screen"
            ]["property_priority_score"],
            "predicted_coolant_kg": row[
                "predicted_coolant_kg"
            ],
            "entry_relative_uncertainty": row[
                "entry_relative_uncertainty"
            ],
            "conservative_coolant_kg": row[
                "conservative_coolant_kg"
            ],
            "water_reference_kg": row[
                "water_reference_kg"
            ],
            "conservative_ratio_vs_water": row[
                "conservative_ratio_vs_water"
            ],
            "beats_water_conservative": row[
                "beats_water_conservative"
            ],
            "nearest_similarity": row["domain"][
                "nearest_similarity"
            ],
        }
        for row in ranked
    ]
    _write_csv(
        args.output_dir / "ranking.csv",
        ranking_csv,
    )

    exploratory_csv = [
        {
            "candidate_id": row["candidate_id"],
            "smiles": row["smiles"],
            "family": row["provenance"][0]["family"],
            "lane_reason": row["lane_reason"],
            "domain_status": row["domain"]["status"],
            "nearest_similarity": row["domain"][
                "nearest_similarity"
            ],
            "property_priority_score": row[
                "property_screen"
            ]["property_priority_score"],
            "predicted_coolant_kg": row[
                "predicted_coolant_kg"
            ],
            "predicted_ratio_vs_water": row[
                "predicted_ratio_vs_water"
            ],
            "predicted_below_water": row[
                "predicted_below_water"
            ],
        }
        for row in exploratory_evaluated
    ]
    _write_csv(
        args.output_dir / "exploratory.csv",
        exploratory_csv,
    )

    rejection_csv = [
        {
            "candidate_id": row["candidate_id"],
            "smiles": row["smiles"],
            "family": row["provenance"][0]["family"],
            "rejection_stage": row.get(
                "rejection_stage"
            ),
            "rejection_reason": row.get(
                "rejection_reason"
            ),
        }
        for row in rejected
    ]
    _write_csv(
        args.output_dir / "rejections.csv",
        rejection_csv,
    )

    summary = {
        "study_version": "v5c-2",
        "study_complete": True,
        "v5c_authorized_by_robustness": bool(
            robustness_summary.get("v5c_ready")
        ),
        "calibration_atomic_numbers": sorted(
            calibration_elements
        ),
        "water_reference": {
            "coolant_used_kg": water_mass_kg,
            "status": water_entry.get("status"),
            "property_screen": water_properties,
        },
        "generation": generation_stats,
        "prescreen": {
            "pass_count": len(prescreen_pass),
            "rankable_pass_count": len(rankable),
            "exploratory_pass_count": len(exploratory),
            "rankable_selected_for_entry": len(
                rankable_selected
            ),
            "exploratory_selected_for_entry": len(
                exploratory_selected
            ),
        },
        "screening": {
            "ranked_candidate_count": len(ranked),
            "exploratory_evaluated_count": len(
                exploratory_evaluated
            ),
            "rejected_candidate_count": len(rejected),
            "rejection_counts": rejection_counts,
            "predicted_rankable_below_water_count": len(
                predicted_rankable_below_water
            ),
            "conservative_water_winner_count": len(
                conservative_winners
            ),
            "exploratory_predicted_below_water_count": len(
                exploratory_below_water
            ),
            "exploratory_family_summary": exploratory_family_summary,
        },
        "ranking": ranked,
        "best_rankable_candidate": (
            ranked[0]
            if ranked
            else None
        ),
        "conservative_water_winners": conservative_winners,
        "exploratory_family_summary": exploratory_family_summary,
        "exploratory_evaluated": exploratory_evaluated,
        "best_exploratory_candidate": (
            exploratory_evaluated[0]
            if exploratory_evaluated
            else None
        ),
        "exploratory_predicted_below_water": (
            exploratory_below_water
        ),
        "interpretation": (
            "Only rankable candidates may claim a screening Water win. "
            "Exploratory candidates expose promising model-supported "
            "directions that require new reference compounds and uncertainty "
            "recalibration before ranking."
        ),
    }
    write_summary(
        args.output_dir / "summary.json",
        summary,
    )
    write_summary(
        args.output_dir / "water_reference.json",
        water_entry,
    )
    write_summary(
        args.output_dir / "production_uncertainty_model.json",
        {
            "calibration_records": calibration_records,
            "calibration_atomic_numbers": sorted(
                calibration_elements
            ),
            "entry_calibration": entry_calibration,
            "domain_settings": config["uncertainty"],
            "source_calibration_summary": str(
                args.calibration_summary
            ),
            "source_robustness_summary": str(
                args.robustness_summary
            ),
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
        "study_version": "v5c-2",
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
            "prior_generated_candidates": (
                str(args.prior_generated_candidates)
                if args.prior_generated_candidates.is_file()
                else None
            ),
            "prior_generated_candidates_sha256": (
                _sha256(args.prior_generated_candidates)
                if args.prior_generated_candidates.is_file()
                else None
            ),
            "calibration_summary_sha256": _sha256(
                args.calibration_summary
            ),
            "robustness_summary_sha256": _sha256(
                args.robustness_summary
            ),
        },
        "search_contract": {
            "rankable": (
                "candidate elements are covered by V5b calibration, "
                "structural domain is in_domain, and screening uncertainty "
                "has a finite conservative bound"
            ),
            "exploratory_domain_expansion": (
                "property/entry predictions are diagnostic only and cannot "
                "count as Water winners"
            ),
            "property_target": (
                "q25 usable enthalpy over the configured T/P grid is the "
                "primary pre-entry priority metric"
            ),
        },
    }
    write_summary(
        args.output_dir / "manifest.json",
        manifest,
    )

    print(json.dumps(
        {
            "study_complete": True,
            "water_reference_kg": water_mass_kg,
            "generation": generation_stats,
            "prescreen": summary["prescreen"],
            "rejection_counts": rejection_counts,
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
            "best_rankable_candidate": (
                {
                    "candidate_id": ranked[0]["candidate_id"],
                    "smiles": ranked[0]["smiles"],
                    "family": ranked[0]["provenance"][0][
                        "family"
                    ],
                    "property_priority_score": ranked[0][
                        "property_screen"
                    ]["property_priority_score"],
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
                    "smiles": exploratory_evaluated[0][
                        "smiles"
                    ],
                    "family": exploratory_evaluated[0][
                        "provenance"
                    ][0]["family"],
                    "lane_reason": exploratory_evaluated[0][
                        "lane_reason"
                    ],
                    "property_priority_score": (
                        exploratory_evaluated[0][
                            "property_screen"
                        ]["property_priority_score"]
                    ),
                    "predicted_coolant_kg": (
                        exploratory_evaluated[0][
                            "predicted_coolant_kg"
                        ]
                    ),
                    "predicted_ratio_vs_water": (
                        exploratory_evaluated[0][
                            "predicted_ratio_vs_water"
                        ]
                    ),
                }
                if exploratory_evaluated
                else None
            ),
        },
        indent=2,
    ))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
