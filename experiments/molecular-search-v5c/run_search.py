#!/usr/bin/env python3
"""V5c-1: constrained molecular generation and uncertainty-aware entry ranking."""
from __future__ import annotations

import argparse
import copy
import csv
import hashlib
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
from generator import generate_candidates  # noqa: E402
from property_provider import build_property_provider  # noqa: E402
from run_physical_sensitivity import apply_values  # noqa: E402

SUCCESS_STATUSES = {"terminal_velocity", "terminal_altitude"}


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _stable_candidate_id(smiles: str) -> str:
    digest = hashlib.sha256(smiles.encode("utf-8")).hexdigest()[:12]
    return f"v5c-{digest}"


def conservative_score(
    predicted_coolant_kg: float,
    relative_uncertainty: float,
    *,
    multiplier: float = 1.0,
) -> float:
    predicted = float(predicted_coolant_kg)
    uncertainty = float(relative_uncertainty)
    multiplier = float(multiplier)
    if (
        not math.isfinite(predicted)
        or predicted <= 0.0
        or not math.isfinite(uncertainty)
        or uncertainty < 0.0
        or not math.isfinite(multiplier)
        or multiplier < 0.0
    ):
        raise ValueError("invalid conservative-score input")
    effective_uncertainty = multiplier * uncertainty
    if effective_uncertainty >= 1.0:
        raise ValueError(
            "relative uncertainty >= 1 has no finite conservative upper bound"
        )
    # V5b error is |predicted-reference| / reference <= u.
    # Therefore predicted >= reference * (1-u), hence:
    # reference <= predicted / (1-u).
    return predicted / (1.0 - effective_uncertainty)


def _write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    if not rows:
        path.write_text("", encoding="utf-8")
        return
    fields = sorted({key for row in rows for key in row})
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def _load_calibration_records(summary: dict[str, Any]) -> list[dict[str, Any]]:
    records: list[dict[str, Any]] = []
    for row in summary.get("candidates", []):
        candidate_id = row.get("candidate_id")
        smiles = row.get("smiles")
        error = row.get("entry_error")
        if (
            candidate_id
            and smiles
            and isinstance(error, (int, float))
            and math.isfinite(float(error))
        ):
            records.append(
                {
                    "candidate_id": str(candidate_id),
                    "smiles": canonical_smiles(str(smiles)),
                    "entry_error": abs(float(error)),
                }
            )
    if len(records) < 3:
        raise ValueError(
            "calibration summary does not contain enough entry-error records"
        )
    return records


def _base_entry_config(
    base: dict[str, Any],
    physical: dict[str, Any],
    search: dict[str, Any],
) -> dict[str, Any]:
    config = apply_values(
        base,
        physical.get("fixed_overrides", {}) or {},
    )
    config = apply_values(
        config,
        physical.get("nominal_values", {}) or {},
    )
    config = copy.deepcopy(config)
    config["chemistry"] = {"mode": "disabled"}
    config.setdefault("heating", {})["backend"] = str(
        search["entry"].get(
            "heating_backend",
            "brandis_johnston_2014",
        )
    )
    config["vehicle"]["couple_coolant_mass_to_trajectory"] = False
    return config


def _coolant_config(
    search: dict[str, Any],
    *,
    smiles: str | None,
    name: str,
    reference_coolprop_name: str | None = None,
) -> dict[str, Any]:
    storage = search["storage"]
    cooling = search["cooling"]
    coolant: dict[str, Any] = {
        "storage_temperature_k": float(storage["temperature_k"]),
        "storage_pressure_pa": float(storage["pressure_pa"]),
        "max_exit_temperature_k": float(
            cooling["max_exit_temperature_k"]
        ),
        "wall_to_fluid_approach_k": float(
            cooling["wall_to_fluid_approach_k"]
        ),
        "cooled_area_m2": float(cooling["cooled_area_m2"]),
        "available_mass_kg": None,
        "injection_pressure_margin": 1.15,
        "porous_delta_p_pa": 100000.0,
        "max_injection_pressure_pa": 20000000.0,
    }
    if reference_coolprop_name is not None:
        coolant["coolprop_name"] = reference_coolprop_name
    else:
        if smiles is None:
            raise ValueError("structure-derived coolant requires SMILES")
        coolant["property_provider"] = {
            "type": "feos",
            "model": "gc_pcsaft_joback",
            "name": name,
            "smiles": smiles,
        }
    return coolant


def _water_reference(
    base: dict[str, Any],
    physical: dict[str, Any],
    search: dict[str, Any],
) -> dict[str, Any]:
    config = _base_entry_config(base, physical, search)
    config["name"] = "v5c-water-reference"
    config["coolant"] = _coolant_config(
        search,
        smiles=None,
        name="water-reference",
        reference_coolprop_name="Water",
    )
    result = evaluate(config)
    if result.get("status") not in SUCCESS_STATUSES:
        raise RuntimeError(
            "water reference did not reach a comparable terminal condition: "
            f"{result.get('status')}: {result.get('failure_reason')}"
        )
    mass = result.get("coolant_used_kg")
    if (
        not isinstance(mass, (int, float))
        or not math.isfinite(float(mass))
        or float(mass) <= 0.0
    ):
        raise RuntimeError("water reference returned invalid coolant mass")
    return result


def _descriptor_check(
    smiles: str,
    constraints: dict[str, Any],
) -> tuple[dict[str, Any] | None, str | None]:
    from rdkit import Chem
    from rdkit.Chem import Descriptors

    molecule = Chem.MolFromSmiles(smiles)
    if molecule is None:
        return None, "invalid_smiles"

    allowed = {
        int(value)
        for value in constraints["allowed_atomic_numbers"]
    }
    atomic_numbers = {
        int(atom.GetAtomicNum())
        for atom in molecule.GetAtoms()
    }
    if not atomic_numbers.issubset(allowed):
        return None, "disallowed_element"

    formal_charge = sum(
        int(atom.GetFormalCharge())
        for atom in molecule.GetAtoms()
    )
    if bool(constraints.get("require_neutral", True)) and formal_charge != 0:
        return None, "non_neutral"

    radical_electrons = sum(
        int(atom.GetNumRadicalElectrons())
        for atom in molecule.GetAtoms()
    )
    if radical_electrons != 0:
        return None, "radical"

    molecular_weight = float(Descriptors.MolWt(molecule))
    heavy_atoms = int(molecule.GetNumHeavyAtoms())
    rings = int(molecule.GetRingInfo().NumRings())

    if not (
        float(constraints["min_molecular_weight_g_mol"])
        <= molecular_weight
        <= float(constraints["max_molecular_weight_g_mol"])
    ):
        return None, "molecular_weight"
    if not (
        int(constraints["min_heavy_atoms"])
        <= heavy_atoms
        <= int(constraints["max_heavy_atoms"])
    ):
        return None, "heavy_atom_count"
    if rings > int(constraints["max_rings"]):
        return None, "ring_count"

    return {
        "molecular_weight_g_mol": molecular_weight,
        "heavy_atoms": heavy_atoms,
        "rings": rings,
        "formal_charge": formal_charge,
        "atomic_numbers": sorted(atomic_numbers),
    }, None


def _evaluate_candidate(task: dict[str, Any]) -> dict[str, Any]:
    candidate = task["candidate"]
    search = task["search"]
    base = task["base"]
    physical = task["physical"]
    calibration_records = task["calibration_records"]
    entry_calibration = task["entry_calibration"]
    water_mass_kg = float(task["water_mass_kg"])

    smiles = str(candidate["smiles"])
    candidate_id = str(candidate["candidate_id"])
    descriptor, descriptor_failure = _descriptor_check(
        smiles,
        search["candidate_constraints"],
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
            search,
            smiles=smiles,
            name=candidate_id,
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

    storage_t = float(search["storage"]["temperature_k"])
    storage_p = float(search["storage"]["pressure_pa"])
    try:
        storage_state = provider.state(storage_t, storage_p)
        storage_phase = provider.phase(storage_t, storage_p)
    except Exception as exc:
        return {
            **candidate,
            "status": "rejected",
            "rejection_stage": "property_model_failure",
            "rejection_reason": str(exc),
            "descriptor": descriptor,
        }

    if (
        bool(search["storage"].get("require_liquid", True))
        and storage_phase != "liquid"
    ):
        return {
            **candidate,
            "status": "rejected",
            "rejection_stage": "storage_not_liquid",
            "rejection_reason": f"predicted phase={storage_phase}",
            "descriptor": descriptor,
            "storage_density_kg_m3": storage_state.density_kg_m3,
        }

    if bool(
        search["storage"].get(
            "require_saturation_at_storage_pressure",
            True,
        )
    ):
        saturation = provider.saturation_at_pressure(storage_p)
        if not saturation.supported:
            return {
                **candidate,
                "status": "rejected",
                "rejection_stage": "saturation_unsupported",
                "rejection_reason": saturation.failure_reason,
                "descriptor": descriptor,
                "storage_density_kg_m3": storage_state.density_kg_m3,
            }
        saturation_temperature_k = saturation.bubble_temperature_k
    else:
        saturation_temperature_k = None

    domain = assess_domain(
        smiles,
        calibration_records,
        in_domain_similarity=float(
            search["uncertainty"]["in_domain_similarity"]
        ),
        edge_similarity=float(
            search["uncertainty"]["edge_similarity"]
        ),
        minimum_neighbors=int(
            search["uncertainty"]["minimum_neighbors"]
        ),
    )
    if (
        bool(search["search"].get("require_in_domain", True))
        and domain.status != "in_domain"
    ):
        return {
            **candidate,
            "status": "rejected",
            "rejection_stage": "applicability_domain",
            "rejection_reason": domain.status,
            "descriptor": descriptor,
            "domain": domain_to_dict(domain),
            "storage_density_kg_m3": storage_state.density_kg_m3,
            "saturation_temperature_k": saturation_temperature_k,
        }

    uncertainty = predict_relative_uncertainty(
        smiles,
        calibration_records,
        entry_calibration,
        error_key="entry_error",
        k_neighbors=int(search["uncertainty"]["k_neighbors"]),
        similarity_floor=float(
            search["uncertainty"]["similarity_floor"]
        ),
    )
    if (
        bool(
            search["search"].get(
                "require_screening_uncertainty",
                True,
            )
        )
        and uncertainty is None
    ):
        return {
            **candidate,
            "status": "rejected",
            "rejection_stage": "uncertainty_unavailable",
            "rejection_reason": "no entry-error uncertainty bound",
            "descriptor": descriptor,
            "domain": domain_to_dict(domain),
        }

    config = _base_entry_config(base, physical, search)
    config["name"] = f"v5c-{candidate_id}"
    config["coolant"] = coolant
    try:
        entry = evaluate(config)
    except Exception as exc:
        return {
            **candidate,
            "status": "rejected",
            "rejection_stage": "entry_exception",
            "rejection_reason": str(exc),
            "descriptor": descriptor,
            "domain": domain_to_dict(domain),
        }

    if entry.get("status") not in SUCCESS_STATUSES:
        return {
            **candidate,
            "status": "rejected",
            "rejection_stage": "entry_failure",
            "rejection_reason": (
                f"{entry.get('status')}: "
                f"{entry.get('failure_reason')}"
            ),
            "descriptor": descriptor,
            "domain": domain_to_dict(domain),
        }

    predicted_mass = entry.get("coolant_used_kg")
    if (
        not isinstance(predicted_mass, (int, float))
        or not math.isfinite(float(predicted_mass))
        or float(predicted_mass) <= 0.0
    ):
        return {
            **candidate,
            "status": "rejected",
            "rejection_stage": "entry_invalid_score",
            "rejection_reason": "invalid coolant_used_kg",
            "descriptor": descriptor,
            "domain": domain_to_dict(domain),
        }

    predicted_mass = float(predicted_mass)
    uncertainty_value = float(uncertainty or 0.0)
    multiplier = float(
        search["search"].get(
            "conservative_score_multiplier",
            1.0,
        )
    )
    try:
        conservative_mass = conservative_score(
            predicted_mass,
            uncertainty_value,
            multiplier=multiplier,
        )
    except ValueError as exc:
        return {
            **candidate,
            "status": "rejected",
            "rejection_stage": "uncertainty_too_large",
            "rejection_reason": str(exc),
            "descriptor": descriptor,
            "domain": domain_to_dict(domain),
            "entry_relative_uncertainty": uncertainty_value,
            "predicted_coolant_kg": predicted_mass,
        }
    storage_volume = (
        predicted_mass / storage_state.density_kg_m3
        if storage_state.density_kg_m3 > 0.0
        else None
    )
    return {
        **candidate,
        "status": "ranked",
        "descriptor": descriptor,
        "domain": domain_to_dict(domain),
        "entry_relative_uncertainty": uncertainty_value,
        "predicted_coolant_kg": predicted_mass,
        "conservative_coolant_kg": conservative_mass,
        "water_reference_kg": water_mass_kg,
        "predicted_ratio_vs_water": predicted_mass / water_mass_kg,
        "conservative_ratio_vs_water": conservative_mass / water_mass_kg,
        "predicted_advantage_fraction_vs_water": (
            water_mass_kg - predicted_mass
        ) / water_mass_kg,
        "conservative_advantage_fraction_vs_water": (
            water_mass_kg - conservative_mass
        ) / water_mass_kg,
        "beats_water_predicted": predicted_mass < water_mass_kg,
        "beats_water_conservative": conservative_mass < water_mass_kg,
        "storage_density_kg_m3": storage_state.density_kg_m3,
        "predicted_storage_volume_m3": storage_volume,
        "saturation_temperature_k": saturation_temperature_k,
        "entry_status": entry.get("status"),
        "entry": entry,
        "provider": provider.metadata(),
    }


def _candidate_pool(
    search: dict[str, Any],
    calibration_records: list[dict[str, Any]],
) -> tuple[list[dict[str, Any]], dict[str, int]]:
    generated = generate_candidates(search)
    known = {
        canonical_smiles(record["smiles"])
        for record in calibration_records
    }

    by_smiles: dict[str, dict[str, Any]] = {}
    invalid_count = 0
    known_count = 0
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
        search["search"].get(
            "maximum_candidates_after_dedup",
            500,
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
        "duplicate_removed_count": duplicate_count,
        "novel_candidate_count": len(candidates),
    }


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
        default=Path("search-v5c1-results"),
    )
    args = parser.parse_args()

    if args.output_dir.exists() and any(args.output_dir.iterdir()):
        parser.error("output directory must be empty; use a new directory")
    for path in (
        args.config,
        args.calibration_summary,
        args.robustness_summary,
        args.base_config,
        args.physical_study,
    ):
        if not path.is_file():
            parser.error(f"missing required input: {path}")

    search = yaml.safe_load(args.config.read_text(encoding="utf-8"))
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

    if bool(search["search"].get("require_v5c_ready", True)):
        if not bool(robustness_summary.get("study_complete")):
            parser.error(
                "robustness summary is incomplete"
            )
        if not bool(robustness_summary.get("v5c_ready")):
            parser.error(
                "robustness summary does not authorize V5c: "
                "v5c_ready != true"
            )

    workers = (
        int(args.workers)
        if args.workers is not None
        else int(search["search"].get("workers", 16))
    )
    if workers <= 0:
        parser.error("--workers must be positive")

    calibration_records = _load_calibration_records(calibration_summary)
    source_candidate_count = robustness_summary.get(
        "source_candidate_count"
    )
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
        coverage=float(
            search["uncertainty"]["coverage_target"]
        ),
        k_neighbors=int(search["uncertainty"]["k_neighbors"]),
        similarity_floor=float(
            search["uncertainty"]["similarity_floor"]
        ),
    )
    if entry_calibration.get("conformal_factor") is None:
        parser.error("could not construct production entry uncertainty model")

    args.output_dir.mkdir(parents=True, exist_ok=True)

    water = _water_reference(base, physical, search)
    water_mass = float(water["coolant_used_kg"])

    candidates, generation_stats = _candidate_pool(
        search,
        calibration_records,
    )
    write_summary(
        args.output_dir / "generated_candidates.json",
        {
            "generation_stats": generation_stats,
            "candidates": candidates,
        },
    )

    tasks = [
        {
            "candidate": candidate,
            "search": search,
            "base": base,
            "physical": physical,
            "calibration_records": calibration_records,
            "entry_calibration": entry_calibration,
            "water_mass_kg": water_mass,
        }
        for candidate in candidates
    ]

    rows: list[dict[str, Any]] = []
    if workers == 1:
        for task in tasks:
            rows.append(_evaluate_candidate(task))
    else:
        with ProcessPoolExecutor(max_workers=workers) as pool:
            futures = [
                pool.submit(_evaluate_candidate, task)
                for task in tasks
            ]
            for future in as_completed(futures):
                rows.append(future.result())

    ranked = [
        row
        for row in rows
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

    rejected = [
        row
        for row in rows
        if row["status"] != "ranked"
    ]
    rejected.sort(
        key=lambda row: (
            str(row.get("rejection_stage")),
            str(row["candidate_id"]),
        )
    )

    rejection_counts: dict[str, int] = {}
    for row in rejected:
        key = str(row.get("rejection_stage", "unknown"))
        rejection_counts[key] = rejection_counts.get(key, 0) + 1

    conservative_winners = [
        row
        for row in ranked
        if row["beats_water_conservative"]
    ]
    predicted_winners = [
        row
        for row in ranked
        if row["beats_water_predicted"]
    ]

    ranking_csv = []
    for row in ranked:
        ranking_csv.append(
            {
                "rank": row["rank"],
                "candidate_id": row["candidate_id"],
                "smiles": row["smiles"],
                "family": row["provenance"][0]["family"],
                "predicted_coolant_kg": row[
                    "predicted_coolant_kg"
                ],
                "entry_relative_uncertainty": row[
                    "entry_relative_uncertainty"
                ],
                "conservative_coolant_kg": row[
                    "conservative_coolant_kg"
                ],
                "water_reference_kg": row["water_reference_kg"],
                "predicted_ratio_vs_water": row[
                    "predicted_ratio_vs_water"
                ],
                "conservative_ratio_vs_water": row[
                    "conservative_ratio_vs_water"
                ],
                "conservative_advantage_fraction_vs_water": row[
                    "conservative_advantage_fraction_vs_water"
                ],
                "beats_water_conservative": row[
                    "beats_water_conservative"
                ],
                "storage_density_kg_m3": row[
                    "storage_density_kg_m3"
                ],
                "predicted_storage_volume_m3": row[
                    "predicted_storage_volume_m3"
                ],
                "nearest_similarity": row["domain"][
                    "nearest_similarity"
                ],
                "domain_neighbor_count": row["domain"][
                    "neighbor_count"
                ],
            }
        )

    rejection_csv = [
        {
            "candidate_id": row["candidate_id"],
            "smiles": row["smiles"],
            "family": row["provenance"][0]["family"],
            "rejection_stage": row.get("rejection_stage"),
            "rejection_reason": row.get("rejection_reason"),
        }
        for row in rejected
    ]
    _write_csv(args.output_dir / "ranking.csv", ranking_csv)
    _write_csv(args.output_dir / "rejections.csv", rejection_csv)

    summary = {
        "study_version": "v5c-1",
        "study_complete": True,
        "v5c_authorized_by_robustness": bool(
            robustness_summary.get("v5c_ready")
        ),
        "water_reference": {
            "coolant_used_kg": water_mass,
            "status": water.get("status"),
        },
        "generation": generation_stats,
        "calibration": {
            "reference_candidate_count": len(calibration_records),
            "coverage_target": float(
                search["uncertainty"]["coverage_target"]
            ),
            "entry_conformal_factor": entry_calibration.get(
                "conformal_factor"
            ),
        },
        "screening": {
            "ranked_candidate_count": len(ranked),
            "rejected_candidate_count": len(rejected),
            "rejection_counts": rejection_counts,
            "predicted_water_winner_count": len(predicted_winners),
            "conservative_water_winner_count": len(
                conservative_winners
            ),
        },
        "best_candidate": ranked[0] if ranked else None,
        "conservative_water_winners": conservative_winners,
        "ranking": ranked,
        "interpretation": (
            "V5c-1 is a constrained screening search. A conservative water "
            "win is a trigger for higher-fidelity validation, not evidence of "
            "flight-level superiority."
        ),
    }
    write_summary(args.output_dir / "summary.json", summary)
    write_summary(
        args.output_dir / "water_reference.json",
        water,
    )
    write_summary(
        args.output_dir / "production_uncertainty_model.json",
        {
            "calibration_records": calibration_records,
            "entry_calibration": entry_calibration,
            "domain_settings": search["uncertainty"],
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
        "study_version": "v5c-1",
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
        },
        "ranking_contract": (
            "conservative_coolant_kg = predicted_coolant_kg / "
            "(1 - conservative_score_multiplier * "
            "screening_entry_relative_uncertainty), for effective u < 1"
        ),
        "candidate_domain": (
            "novel relative to V5b calibration references; neutral C/O "
            "molecules from deterministic supported-family templates"
        ),
    }
    write_summary(args.output_dir / "manifest.json", manifest)

    print(json.dumps(
        {
            "study_complete": summary["study_complete"],
            "water_reference_kg": water_mass,
            "generation": generation_stats,
            "ranked_candidate_count": len(ranked),
            "rejection_counts": rejection_counts,
            "predicted_water_winner_count": len(predicted_winners),
            "conservative_water_winner_count": len(
                conservative_winners
            ),
            "best_candidate": (
                {
                    "candidate_id": ranked[0]["candidate_id"],
                    "smiles": ranked[0]["smiles"],
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
                    "beats_water_conservative": ranked[0][
                        "beats_water_conservative"
                    ],
                }
                if ranked
                else None
            ),
        },
        indent=2,
    ))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
