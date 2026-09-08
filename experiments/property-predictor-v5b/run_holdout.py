#!/usr/bin/env python3
"""V5b-1: SMILES -> GC-PC-SAFT/Joback -> holdout properties -> entry score."""
from __future__ import annotations

import argparse
import copy
import csv
import hashlib
import json
import math
import platform
import statistics
import subprocess
import sys
from importlib.metadata import version
from pathlib import Path
from typing import Any

import yaml

HERE = Path(__file__).resolve().parent
ENTRY = HERE.parent / "entry-evaluator"
if str(ENTRY) not in sys.path:
    sys.path.insert(0, str(ENTRY))

from evaluator import evaluate, write_summary  # noqa: E402
from property_provider import build_property_provider, state_to_dict  # noqa: E402
from run_physical_sensitivity import apply_values  # noqa: E402

SUCCESS_STATUSES = {"terminal_velocity", "terminal_altitude"}


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _rel_error(predicted: float | None, reference: float | None) -> float | None:
    if predicted is None or reference is None:
        return None
    if not math.isfinite(predicted) or not math.isfinite(reference):
        return None
    if abs(reference) < 1.0e-12:
        return None
    return (predicted - reference) / reference


def _candidate_coolant(
    candidate: dict[str, Any],
    study: dict[str, Any],
    *,
    reference: bool,
) -> dict[str, Any]:
    coolant: dict[str, Any] = {
        "storage_temperature_k": float(candidate["storage_temperature_k"]),
        "storage_pressure_pa": float(candidate["storage_pressure_pa"]),
        "max_exit_temperature_k": float(candidate["max_exit_temperature_k"]),
        "wall_to_fluid_approach_k": float(
            study["entry"].get("wall_to_fluid_approach_k", 25.0)
        ),
        "cooled_area_m2": float(study["entry"].get("cooled_area_m2", 12.0)),
        "available_mass_kg": None,
        "injection_pressure_margin": 1.15,
        "porous_delta_p_pa": 100000.0,
        "max_injection_pressure_pa": 20000000.0,
    }
    if reference:
        coolant["coolprop_name"] = str(candidate["reference_coolprop_name"])
    else:
        coolant["property_provider"] = {
            "type": "feos",
            "model": "gc_pcsaft_joback",
            "name": str(candidate["id"]),
            "smiles": str(candidate["smiles"]),
        }
    return coolant


def _entry_config(
    base: dict[str, Any],
    physical: dict[str, Any],
    study: dict[str, Any],
    candidate: dict[str, Any],
    *,
    reference: bool,
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
    mode = "reference" if reference else "prediction"
    config["name"] = f"v5b1-{candidate['id']}-{mode}"
    config["coolant"] = _candidate_coolant(candidate, study, reference=reference)
    config["chemistry"] = {"mode": "disabled"}
    config.setdefault("heating", {})["backend"] = str(
        study["entry"].get("heating_backend", "brandis_johnston_2014")
    )
    config["vehicle"]["couple_coolant_mass_to_trajectory"] = False
    return config


def _predict_candidate(
    candidate: dict[str, Any],
    study: dict[str, Any],
    base: dict[str, Any],
    physical: dict[str, Any],
) -> dict[str, Any]:
    coolant = _candidate_coolant(candidate, study, reference=False)
    provider = build_property_provider(coolant)
    storage_t = float(candidate["storage_temperature_k"])
    storage_p = float(candidate["storage_pressure_pa"])
    storage_state = provider.state(storage_t, storage_p)
    storage_h = storage_state.enthalpy_j_kg

    grid: list[dict[str, Any]] = []
    validation = study["validation"]
    for temperature_k in validation["temperatures_k"]:
        for pressure_pa in validation["pressures_pa"]:
            row: dict[str, Any] = {
                "temperature_k": float(temperature_k),
                "pressure_pa": float(pressure_pa),
                "supported": False,
                "failure_reason": None,
            }
            try:
                state = provider.state(float(temperature_k), float(pressure_pa))
                row.update(state_to_dict(state))
                row["delta_h_from_storage_j_kg"] = (
                    state.enthalpy_j_kg - storage_h
                )
                row["supported"] = True
            except Exception as exc:
                row["failure_reason"] = str(exc)
            grid.append(row)

    saturation: list[dict[str, Any]] = []
    for pressure_pa in validation.get("saturation_pressures_pa", []):
        sat = provider.saturation_at_pressure(float(pressure_pa))
        saturation.append(
            {
                "pressure_pa": float(pressure_pa),
                "temperature_k": sat.bubble_temperature_k,
                "supported": sat.supported,
                "failure_reason": sat.failure_reason,
            }
        )

    entry_summary = evaluate(
        _entry_config(
            base,
            physical,
            study,
            candidate,
            reference=False,
        )
    )
    return {
        "candidate_id": candidate["id"],
        "smiles": candidate["smiles"],
        "provider": provider.metadata(),
        "storage_state": state_to_dict(storage_state),
        "grid": grid,
        "saturation": saturation,
        "entry": entry_summary,
    }


def _reference_candidate(
    candidate: dict[str, Any],
    study: dict[str, Any],
    base: dict[str, Any],
    physical: dict[str, Any],
) -> dict[str, Any]:
    coolant = _candidate_coolant(candidate, study, reference=True)
    provider = build_property_provider(coolant)
    storage_t = float(candidate["storage_temperature_k"])
    storage_p = float(candidate["storage_pressure_pa"])
    storage_state = provider.state(storage_t, storage_p)
    storage_h = storage_state.enthalpy_j_kg

    grid: list[dict[str, Any]] = []
    validation = study["validation"]
    for temperature_k in validation["temperatures_k"]:
        for pressure_pa in validation["pressures_pa"]:
            row: dict[str, Any] = {
                "temperature_k": float(temperature_k),
                "pressure_pa": float(pressure_pa),
                "supported": False,
                "failure_reason": None,
            }
            try:
                state = provider.state(float(temperature_k), float(pressure_pa))
                row.update(state_to_dict(state))
                row["delta_h_from_storage_j_kg"] = (
                    state.enthalpy_j_kg - storage_h
                )
                row["supported"] = True
            except Exception as exc:
                row["failure_reason"] = str(exc)
            grid.append(row)

    saturation: list[dict[str, Any]] = []
    for pressure_pa in validation.get("saturation_pressures_pa", []):
        sat = provider.saturation_at_pressure(float(pressure_pa))
        saturation.append(
            {
                "pressure_pa": float(pressure_pa),
                "temperature_k": sat.bubble_temperature_k,
                "supported": sat.supported,
                "failure_reason": sat.failure_reason,
            }
        )

    entry_summary = evaluate(
        _entry_config(
            base,
            physical,
            study,
            candidate,
            reference=True,
        )
    )
    return {
        "candidate_id": candidate["id"],
        "reference_coolprop_name": candidate["reference_coolprop_name"],
        "provider": provider.metadata(),
        "storage_state": state_to_dict(storage_state),
        "grid": grid,
        "saturation": saturation,
        "entry": entry_summary,
    }


def _compare_candidate(
    prediction: dict[str, Any],
    reference: dict[str, Any],
) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    ref_grid = {
        (row["temperature_k"], row["pressure_pa"]): row
        for row in reference["grid"]
    }
    property_rows: list[dict[str, Any]] = []

    for pred in prediction["grid"]:
        key = (pred["temperature_k"], pred["pressure_pa"])
        ref = ref_grid[key]
        row = {
            "candidate_id": prediction["candidate_id"],
            "smiles": prediction["smiles"],
            "temperature_k": key[0],
            "pressure_pa": key[1],
            "prediction_supported": pred["supported"],
            "reference_supported": ref["supported"],
            "prediction_failure_reason": pred["failure_reason"],
            "reference_failure_reason": ref["failure_reason"],
            "predicted_phase": pred.get("phase"),
            "reference_phase": ref.get("phase"),
            "phase_match": (
                pred.get("phase") == ref.get("phase")
                if pred["supported"] and ref["supported"]
                else None
            ),
        }
        for field, output in (
            ("delta_h_from_storage_j_kg", "delta_h"),
            ("density_kg_m3", "density"),
            ("cp_j_kg_k", "cp"),
        ):
            p = pred.get(field)
            q = ref.get(field)
            row[f"predicted_{output}"] = p
            row[f"reference_{output}"] = q
            row[f"{output}_relative_error"] = _rel_error(p, q)
        property_rows.append(row)

    pred_sat = {
        row["pressure_pa"]: row
        for row in prediction["saturation"]
    }
    ref_sat = {
        row["pressure_pa"]: row
        for row in reference["saturation"]
    }
    saturation_rows = []
    for pressure_pa, pred in pred_sat.items():
        ref = ref_sat[pressure_pa]
        saturation_rows.append(
            {
                "pressure_pa": pressure_pa,
                "predicted_temperature_k": pred["temperature_k"],
                "reference_temperature_k": ref["temperature_k"],
                "relative_error": _rel_error(
                    pred["temperature_k"], ref["temperature_k"]
                ),
                "prediction_supported": pred["supported"],
                "reference_supported": ref["supported"],
            }
        )

    pred_entry = prediction["entry"]
    ref_entry = reference["entry"]
    comparable = (
        pred_entry.get("status") in SUCCESS_STATUSES
        and pred_entry.get("status") == ref_entry.get("status")
        and isinstance(pred_entry.get("coolant_used_kg"), (int, float))
        and isinstance(ref_entry.get("coolant_used_kg"), (int, float))
        and float(ref_entry["coolant_used_kg"]) > 0.0
    )
    entry_error = (
        _rel_error(
            float(pred_entry["coolant_used_kg"]),
            float(ref_entry["coolant_used_kg"]),
        )
        if comparable
        else None
    )
    comparison = {
        "candidate_id": prediction["candidate_id"],
        "smiles": prediction["smiles"],
        "prediction_status": pred_entry.get("status"),
        "reference_status": ref_entry.get("status"),
        "predicted_coolant_kg": (
            pred_entry.get("coolant_used_kg")
            if pred_entry.get("status") in SUCCESS_STATUSES
            else None
        ),
        "reference_coolant_kg": (
            ref_entry.get("coolant_used_kg")
            if ref_entry.get("status") in SUCCESS_STATUSES
            else None
        ),
        "entry_comparable": comparable,
        "entry_coolant_mass_relative_error": entry_error,
        "entry_coolant_mass_abs_relative_error": (
            abs(entry_error) if entry_error is not None else None
        ),
        "predicted_storage_density_kg_m3": prediction[
            "storage_state"
        ]["density_kg_m3"],
        "reference_storage_density_kg_m3": reference[
            "storage_state"
        ]["density_kg_m3"],
        "storage_density_relative_error": _rel_error(
            prediction["storage_state"]["density_kg_m3"],
            reference["storage_state"]["density_kg_m3"],
        ),
        "saturation": saturation_rows,
    }
    return property_rows, comparison


def _write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    fields = sorted({key for row in rows for key in row})
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def _summary(
    comparisons: list[dict[str, Any]],
    property_rows: list[dict[str, Any]],
) -> dict[str, Any]:
    def values(key: str) -> list[float]:
        return [
            abs(float(row[key]))
            for row in property_rows
            if isinstance(row.get(key), (int, float))
            and math.isfinite(float(row[key]))
        ]

    entry_errors = [
        float(row["entry_coolant_mass_abs_relative_error"])
        for row in comparisons
        if isinstance(row.get("entry_coolant_mass_abs_relative_error"), (int, float))
    ]
    delta_h = values("delta_h_relative_error")
    density = values("density_relative_error")
    cp = values("cp_relative_error")

    return {
        "study_complete": len(comparisons) > 0
        and all(row["entry_comparable"] for row in comparisons),
        "candidate_count": len(comparisons),
        "entry_comparable_count": sum(
            bool(row["entry_comparable"]) for row in comparisons
        ),
        "reference_revealed_after_prediction": True,
        "structure_only_prediction": True,
        "holdout_error": {
            "entry_coolant_mass_abs_relative_error": {
                "median": statistics.median(entry_errors) if entry_errors else None,
                "max": max(entry_errors) if entry_errors else None,
            },
            "delta_h_abs_relative_error": {
                "median": statistics.median(delta_h) if delta_h else None,
                "max": max(delta_h) if delta_h else None,
            },
            "density_abs_relative_error": {
                "median": statistics.median(density) if density else None,
                "max": max(density) if density else None,
            },
            "cp_abs_relative_error": {
                "median": statistics.median(cp) if cp else None,
                "max": max(cp) if cp else None,
            },
        },
        "candidates": comparisons,
        "uncertainty_status": (
            "empirical holdout error measured; calibrated per-candidate uncertainty "
            "is deferred to V5b-2"
        ),
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--study",
        type=Path,
        default=HERE / "benchmark.yaml",
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
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("property-v5b1-results"),
    )
    args = parser.parse_args()

    if args.output_dir.exists() and any(args.output_dir.iterdir()):
        parser.error("output directory must be empty; use a new directory")
    args.output_dir.mkdir(parents=True, exist_ok=True)

    study = yaml.safe_load(args.study.read_text(encoding="utf-8"))
    base = yaml.safe_load(args.base_config.read_text(encoding="utf-8"))
    physical = yaml.safe_load(args.physical_study.read_text(encoding="utf-8"))

    predictions: list[dict[str, Any]] = []
    prediction_failures: list[dict[str, Any]] = []
    for candidate in study["holdouts"]:
        try:
            result = _predict_candidate(candidate, study, base, physical)
            predictions.append(result)
            print(
                f"predicted {candidate['id']}: "
                f"{result['entry'].get('status')} "
                f"coolant={result['entry'].get('coolant_used_kg')}",
                flush=True,
            )
        except Exception as exc:
            prediction_failures.append(
                {
                    "candidate_id": candidate["id"],
                    "smiles": candidate["smiles"],
                    "failure_reason": str(exc),
                }
            )
            print(f"prediction failed {candidate['id']}: {exc}", flush=True)

    # Persist prediction artifacts before the hidden-reference phase.
    prediction_artifact = {
        "reference_properties_used_during_prediction": False,
        "predictions": predictions,
        "failures": prediction_failures,
    }
    write_summary(
        args.output_dir / "predictions_before_reference.json",
        prediction_artifact,
    )

    references: dict[str, dict[str, Any]] = {}
    reference_failures: list[dict[str, Any]] = []
    for candidate in study["holdouts"]:
        try:
            references[candidate["id"]] = _reference_candidate(
                candidate, study, base, physical
            )
        except Exception as exc:
            reference_failures.append(
                {
                    "candidate_id": candidate["id"],
                    "failure_reason": str(exc),
                }
            )

    property_rows: list[dict[str, Any]] = []
    comparisons: list[dict[str, Any]] = []
    for prediction in predictions:
        reference = references.get(prediction["candidate_id"])
        if reference is None:
            continue
        rows, comparison = _compare_candidate(prediction, reference)
        property_rows.extend(rows)
        comparisons.append(comparison)

    _write_csv(args.output_dir / "property_comparison.csv", property_rows)
    _write_csv(args.output_dir / "entry_comparison.csv", comparisons)
    report = _summary(comparisons, property_rows)
    report["prediction_failures"] = prediction_failures
    report["reference_failures"] = reference_failures
    write_summary(args.output_dir / "summary.json", report)

    parameter_dir = ENTRY / "parameters" / "v5b"
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
        "study_version": "v5b-1",
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
        "study_sha256": _sha256(args.study),
        "parameter_files": {
            path.name: _sha256(path)
            for path in sorted(parameter_dir.glob("*.json"))
        },
        "prediction_contract": (
            "SMILES plus pinned GC-PC-SAFT/Joback group parameters only; "
            "CoolProp reference evaluated after prediction artifact is written"
        ),
        "enthalpy_comparison": (
            "compare heat uptake delta_h relative to each provider's own "
            "storage reference; never subtract absolute enthalpies across backends"
        ),
        "important_limitation": (
            "V5b-1 measures empirical holdout error for two known pure liquids. "
            "It does not provide calibrated per-candidate uncertainty, mixture "
            "prediction, transport validation, decomposition chemistry, or "
            "system-level storage/feed penalties."
        ),
    }
    write_summary(args.output_dir / "manifest.json", manifest)
    print(json.dumps(report, indent=2))
    return 0 if report["study_complete"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
