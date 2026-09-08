#!/usr/bin/env python3
"""V5a: validate and score a known binary liquid-mixture composition sweep."""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import platform
import subprocess
from concurrent.futures import ProcessPoolExecutor, as_completed
from importlib.metadata import version
from pathlib import Path
from typing import Any

import yaml

from evaluator import evaluate, write_summary
from property_provider import (
    build_property_provider,
    coolprop_runtime_info,
    state_to_dict,
)
from run_physical_sensitivity import apply_values

HERE = Path(__file__).resolve().parent
SUCCESS_STATUSES = {"terminal_velocity", "terminal_altitude"}


def fraction_slug(value: float) -> str:
    return f"{value:.6f}".rstrip("0").rstrip(".").replace(".", "p")


def validate_fraction(value: float) -> float:
    value = float(value)
    if not math.isfinite(value) or not 0.0 <= value <= 1.0:
        raise ValueError("mixture fractions must be finite and in [0, 1]")
    return value


def provider_block(study: dict[str, Any], fraction_a: float) -> dict[str, Any]:
    return {
        **study.get("property_provider", {"type": "thermo", "model": "NRTL"}),
        "components": list(study["components"]),
        "composition_basis": str(study.get("composition_basis", "mole")),
        "fractions": [fraction_a, 1.0 - fraction_a],
    }


def coolant_block(study: dict[str, Any], fraction_a: float) -> dict[str, Any]:
    return {
        "storage_temperature_k": float(study["storage_temperature_k"]),
        "storage_pressure_pa": float(study["storage_pressure_pa"]),
        "max_exit_temperature_k": float(study["max_exit_temperature_k"]),
        "wall_to_fluid_approach_k": float(
            study.get("wall_to_fluid_approach_k", 25.0)
        ),
        "cooled_area_m2": float(study.get("cooled_area_m2", 12.0)),
        "available_mass_kg": None,
        "injection_pressure_margin": 1.15,
        "porous_delta_p_pa": 100000.0,
        "max_injection_pressure_pa": 20000000.0,
        "property_provider": provider_block(study, fraction_a),
    }


def build_entry_config(
    base: dict[str, Any],
    physical: dict[str, Any],
    study: dict[str, Any],
    fraction_a: float,
) -> dict[str, Any]:
    config = apply_values(
        base,
        physical.get("fixed_overrides", {}) or {},
    )
    config = apply_values(
        config,
        physical.get("nominal_values", {}) or {},
    )

    config["name"] = (
        f"v5a-{study['name']}-xa-{fraction_slug(fraction_a)}"
    )
    config["coolant"] = coolant_block(study, fraction_a)
    config["chemistry"] = {"mode": "disabled"}
    config.setdefault("heating", {})["backend"] = str(
        study.get("heating_backend", "brandis_johnston_2014")
    )
    config["vehicle"]["couple_coolant_mass_to_trajectory"] = False
    return config


def validate_candidate(
    study: dict[str, Any],
    fraction_a: float,
) -> tuple[dict[str, Any], list[dict[str, Any]]]:
    coolant = coolant_block(study, fraction_a)
    provider = build_property_provider(coolant)
    meta = provider.metadata()

    storage_t = float(study["storage_temperature_k"])
    storage_p = float(study["storage_pressure_pa"])

    result: dict[str, Any] = {
        "candidate_id": (
            f"{study['name']}-xa-{fraction_slug(fraction_a)}"
        ),
        "fraction_component_a": fraction_a,
        "fraction_component_b": 1.0 - fraction_a,
        "requested_components": list(study["components"]),
        "provider_identity": provider.identity,
        "provider_components": meta["components"],
        "provider_mole_fractions": meta["mole_fractions"],
        "composition_basis": meta["composition_basis"],
        "storage_supported": False,
        "storage_liquid": False,
        "storage_phase": None,
        "storage_density_kg_m3": None,
        "storage_enthalpy_j_kg": None,
        "storage_cp_j_kg_k": None,
        "storage_viscosity_pa_s": None,
        "storage_conductivity_w_m_k": None,
        "bubble_temperature_k": None,
        "dew_temperature_k": None,
        "saturation_supported": False,
        "state_points_requested": 0,
        "state_points_supported": 0,
        "transport_viscosity_supported": 0,
        "transport_conductivity_supported": 0,
        "failure_reason": None,
    }

    try:
        storage = provider.state(storage_t, storage_p)
        result.update(
            {
                "storage_supported": True,
                "storage_liquid": storage.phase
                in {"liquid", "supercritical_liquid"},
                "storage_phase": storage.phase,
                "storage_density_kg_m3": storage.density_kg_m3,
                "storage_enthalpy_j_kg": storage.enthalpy_j_kg,
                "storage_cp_j_kg_k": storage.cp_j_kg_k,
                "storage_viscosity_pa_s": storage.viscosity_pa_s,
                "storage_conductivity_w_m_k": (
                    storage.conductivity_w_m_k
                ),
            }
        )
    except Exception as exc:
        result["failure_reason"] = f"storage state: {exc}"

    saturation = provider.saturation_at_pressure(storage_p)
    result["bubble_temperature_k"] = saturation.bubble_temperature_k
    result["dew_temperature_k"] = saturation.dew_temperature_k
    result["saturation_supported"] = saturation.supported
    if (
        result["failure_reason"] is None
        and saturation.failure_reason is not None
    ):
        result["saturation_note"] = saturation.failure_reason

    validation = study.get("validation", {}) or {}
    temperatures = [
        float(value)
        for value in validation.get(
            "temperatures_k",
            [storage_t],
        )
    ]
    pressures = [
        float(value)
        for value in validation.get(
            "pressures_pa",
            [storage_p],
        )
    ]

    grid_rows: list[dict[str, Any]] = []
    for temperature_k in temperatures:
        for pressure_pa in pressures:
            result["state_points_requested"] += 1
            row = {
                "candidate_id": result["candidate_id"],
                "fraction_component_a": fraction_a,
                "temperature_k": temperature_k,
                "pressure_pa": pressure_pa,
                "supported": False,
                "failure_reason": None,
            }
            try:
                state = provider.state(
                    temperature_k,
                    pressure_pa,
                )
                row.update(state_to_dict(state))
                row["supported"] = True
                result["state_points_supported"] += 1
                if state.viscosity_pa_s is not None:
                    result["transport_viscosity_supported"] += 1
                if state.conductivity_w_m_k is not None:
                    result[
                        "transport_conductivity_supported"
                    ] += 1
            except Exception as exc:
                row["failure_reason"] = str(exc)
            grid_rows.append(row)

    return result, grid_rows


def continuity_report(
    validations: list[dict[str, Any]],
    components: list[str],
) -> dict[str, Any]:
    ordered = sorted(
        validations,
        key=lambda row: row["fraction_component_a"],
    )
    endpoint_b = ordered[0]
    endpoint_a = ordered[-1]

    endpoint_pass = (
        endpoint_b["fraction_component_a"] == 0.0
        and endpoint_b["provider_components"] == [components[1]]
        and endpoint_a["fraction_component_a"] == 1.0
        and endpoint_a["provider_components"] == [components[0]]
    )

    metrics: dict[str, Any] = {}
    for key in (
        "storage_density_kg_m3",
        "storage_cp_j_kg_k",
        "storage_viscosity_pa_s",
        "storage_conductivity_w_m_k",
        "bubble_temperature_k",
        "dew_temperature_k",
    ):
        def finite(value):
            return isinstance(value, (int, float)) and math.isfinite(value)

        changes = []
        missing = []
        for left, right in zip(ordered, ordered[1:]):
            x0, x1 = left["fraction_component_a"], right["fraction_component_a"]
            y0, y1 = left.get(key), right.get(key)
            if not (finite(y0) and finite(y1)):
                missing.append({"x0": x0, "x1": x1})
                continue
            changes.append({"x0": x0, "x1": x1,
                            "relative_change": abs(y1 - y0) / max(abs(y0), abs(y1), 1e-30)})
        metrics[key] = {
            "supported_compositions": sum(finite(row.get(key)) for row in ordered),
            "max_adjacent_relative_change": max((x["relative_change"] for x in changes), default=None),
            "adjacent_changes": changes,
            "unsupported_intervals": missing,
        }
    metrics["storage_enthalpy_j_kg"] = {
        "comparison": "not compared: absolute enthalpy references differ between thermo and CoolProp; compare delta_h"
    }

    return {
        "endpoint_pure_state_pass": endpoint_pass,
        "component_a": components[0],
        "component_b": components[1],
        "composition_count": len(ordered),
        "storage_liquid_count": sum(
            bool(row["storage_liquid"])
            for row in ordered
        ),
        "all_storage_states_liquid": all(
            bool(row["storage_liquid"])
            for row in ordered
        ),
        "all_thermo_grid_points_supported": all(
            row["state_points_requested"]
            == row["state_points_supported"]
            for row in ordered
        ),
        "metrics": metrics,
        "interpretation": (
            "Smoothness metrics are diagnostics only. Non-ideal mixtures can "
            "be physically non-monotonic; V5a does not hard-reject on shape."
        ),
    }


def preflight_failure_reason(validation: dict[str, Any]) -> str | None:
    reasons = []
    if not validation["storage_liquid"]:
        reasons.append("unsupported or non-liquid storage state")
    if validation["state_points_requested"] != validation["state_points_supported"]:
        reasons.append("unsupported required property grid points")
    return "; ".join(reasons) or None


def run_entry_task(task: dict[str, Any]) -> dict[str, Any]:
    try:
        summary = evaluate(task["config"])
    except Exception as exc:
        summary = {
            "status": "exception",
            "failure_reason": str(exc),
        }
    return {
        "candidate_id": task["candidate_id"],
        "fraction_component_a": task["fraction_component_a"],
        "fraction_component_b": 1.0 - task["fraction_component_a"],
        **summary,
    }


def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    fields = sorted({key for row in rows for key in row})
    with path.open(
        "w",
        newline="",
        encoding="utf-8",
    ) as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=fields,
        )
        writer.writeheader()
        writer.writerows(rows)


def valid_score(row: dict[str, Any]) -> bool:
    mass = row.get("coolant_used_kg")
    return (row.get("status") in SUCCESS_STATUSES
            and not row.get("failure_reason")
            and isinstance(mass, (int, float)) and math.isfinite(mass) and mass > 0)


def entry_summary(rows: list[dict[str, Any]], reference_fraction: float) -> dict[str, Any]:
    reference = next((row for row in rows if math.isclose(
        float(row["fraction_component_a"]), reference_fraction, rel_tol=0, abs_tol=1e-12)), None)
    reference_mass = float(reference["coolant_used_kg"]) if reference and valid_score(reference) else None
    comparable = []
    for row in rows:
        row["score_coolant_kg"] = float(row["coolant_used_kg"]) if valid_score(row) else None
        row["coolant_mass_ratio_vs_reference"] = None
        if (reference_mass is not None and valid_score(row)
                and row["status"] == reference["status"]):
            row["coolant_mass_ratio_vs_reference"] = row["score_coolant_kg"] / reference_mass
            comparable.append(row)
    best = min(comparable, key=lambda row: row["score_coolant_kg"], default=None)
    failed = [{"candidate_id": row["candidate_id"],
               "fraction_component_a": row["fraction_component_a"], "status": row.get("status"),
               "consumed_mass_diagnostic_kg": row.get("coolant_used_kg"),
               "failure_reason": row.get("failure_reason") or "missing valid completed-run score"}
              for row in rows if not valid_score(row)]
    return {
        "reference_fraction_component_a": reference_fraction,
        "reference_candidate_id": reference["candidate_id"] if reference else None,
        "reference_coolant_kg": reference_mass,
        "comparable_candidate_count": len(comparable),
        "beats_reference_count": sum(row["coolant_mass_ratio_vs_reference"] < 0.99 for row in comparable),
        "beat_threshold_ratio": 0.99,
        "failed_candidate_count": len(failed), "failed_candidates": failed,
        "best_candidate": None if best is None else {key: best[key] for key in (
            "candidate_id", "fraction_component_a", "score_coolant_kg", "coolant_mass_ratio_vs_reference")},
        "interpretation": "provisional model ranking; caloric validation pending; consumed mass on failure is not a score",
    }


def endpoint_model_comparison(study: dict[str, Any]) -> dict[str, Any]:
    """Compare heat uptake, never absolute enthalpy, across the backend boundary."""
    rows = []
    epsilon = 1e-6
    for pure_fraction, near_fraction in ((0.0, epsilon), (1.0, 1.0 - epsilon)):
        pure = build_property_provider(coolant_block(study, pure_fraction))
        near = build_property_provider(coolant_block(study, near_fraction))
        for pressure in study.get("validation", {}).get("pressures_pa", [101325.0]):
            row = {"pure_fraction_component_a": pure_fraction,
                   "near_fraction_component_a": near_fraction,
                   "outlet_temperature_k": float(study["max_exit_temperature_k"]),
                   "outlet_pressure_pa": float(pressure), "supported": False}
            try:
                deltas = [p.enthalpy_j_kg(row["outlet_temperature_k"], float(pressure))
                          - p.enthalpy_j_kg(float(study["storage_temperature_k"]), float(study["storage_pressure_pa"]))
                          for p in (pure, near)]
                if any(not math.isfinite(h) or h <= 0 for h in deltas):
                    raise ValueError("invalid endpoint heat uptake")
                row.update(supported=True, coolprop_delta_h_j_kg=deltas[0],
                           thermo_near_pure_delta_h_j_kg=deltas[1],
                           relative_difference=(deltas[1] - deltas[0]) / deltas[0])
            except Exception as exc:
                row["failure_reason"] = str(exc)
            rows.append(row)
    return {"rows": rows, "interpretation": "backend-limit diagnostic, not calorimetric validation; no offset correction applied"}


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--base-config",
        type=Path,
        default=HERE / "config.yaml",
    )
    parser.add_argument(
        "--physical-study",
        type=Path,
        default=HERE / "physical-sensitivity.yaml",
    )
    parser.add_argument(
        "--study",
        type=Path,
        default=HERE / "mixture-v5a.yaml",
    )
    parser.add_argument("--workers", type=int, default=16)
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("mixture-v5a-results"),
    )
    args = parser.parse_args()

    if args.workers <= 0:
        parser.error("--workers must be positive")

    if (
        args.output_dir.exists()
        and any(args.output_dir.iterdir())
    ):
        parser.error(
            "output directory must be empty; use a new directory"
        )
    args.output_dir.mkdir(parents=True, exist_ok=True)

    base = yaml.safe_load(
        args.base_config.read_text(encoding="utf-8")
    )
    physical = yaml.safe_load(
        args.physical_study.read_text(encoding="utf-8")
    )
    study = yaml.safe_load(
        args.study.read_text(encoding="utf-8")
    )

    components = [
        str(value)
        for value in study["components"]
    ]
    if len(components) != 2:
        parser.error("V5a requires exactly two components")

    fractions = sorted(
        {
            validate_fraction(value)
            for value in study["fraction_component_a"]
        }
    )
    if 0.0 not in fractions or 1.0 not in fractions:
        parser.error(
            "V5a sweep must include pure endpoints 0.0 and 1.0"
        )

    reference_fraction = validate_fraction(
        study.get(
            "reference_fraction_component_a",
            1.0,
        )
    )
    if reference_fraction not in fractions:
        parser.error(
            "reference_fraction_component_a must be in the sweep"
        )

    try:
        commit = subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=HERE, text=True).strip()
        dirty = bool(subprocess.check_output(["git", "status", "--porcelain"], cwd=HERE, text=True).strip())
    except (OSError, subprocess.CalledProcessError):
        commit, dirty = None, None
    manifest = {
        "study_version": "v5a-nrtl", "git_commit": commit, "git_dirty": dirty,
        "python": platform.python_version(),
        "packages": {p: version(p) for p in ("CoolProp", "thermo", "chemicals", "fluids", "scipy", "pymsis", "numpy", "PyYAML")},
        "coolprop_runtime": coolprop_runtime_info(),
        "input_configs": {"base": base, "physical": physical, "study": study},
        "source_sha256": {p.name: hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(HERE.glob("*.py"))},
        "config_sha256": {str(p.resolve()): hashlib.sha256(p.read_bytes()).hexdigest()
                          for p in (args.base_config, args.physical_study, args.study)},
        "providers": [build_property_provider(coolant_block(study, x)).metadata() for x in fractions],
        "physical_validation_complete": False,
        "important_limitation": "NRTL reproduces a documented VLE example; caloric fit validity is unestablished. Ideal gas vapor, additive liquid volumes, no mixture transport or chemistry. Ranking is provisional.",
    }
    write_summary(args.output_dir / "manifest.json", manifest)
    endpoint_comparison = endpoint_model_comparison(study)
    write_summary(args.output_dir / "endpoint_model_comparison.json", endpoint_comparison)

    validations: list[dict[str, Any]] = []
    grid_rows: list[dict[str, Any]] = []
    for fraction_a in fractions:
        validation, states = validate_candidate(
            study,
            fraction_a,
        )
        validations.append(validation)
        grid_rows.extend(states)

    continuity = continuity_report(
        validations,
        components,
    )
    write_csv(
        args.output_dir / "property_validation.csv",
        validations,
    )
    write_csv(
        args.output_dir / "property_state_grid.csv",
        grid_rows,
    )
    write_summary(
        args.output_dir / "continuity_report.json",
        continuity,
    )

    if not continuity["endpoint_pure_state_pass"]:
        write_summary(args.output_dir / "summary.json", {
            "study_complete": False, "numerical_sweep_complete": False,
            "physical_validation_complete": False, "property_validation": continuity,
            "failure_reason": "pure endpoint identity check failed", "entry_scoring": None,
        })
        return 2

    eligible_fractions = set()
    entry_rows: list[dict[str, Any]] = []
    for validation in validations:
        reason = preflight_failure_reason(validation)
        if reason is None:
            eligible_fractions.add(validation["fraction_component_a"])
        else:
            entry_rows.append({"candidate_id": validation["candidate_id"],
                               "fraction_component_a": validation["fraction_component_a"],
                               "fraction_component_b": validation["fraction_component_b"],
                               "status": "property_validation_failure", "coolant_used_kg": None,
                               "failure_reason": reason})
            print(f"{validation['candidate_id']}: property_validation_failure; score=n/a; {reason}", flush=True)

    tasks = [
        {
            "candidate_id": (
                f"{study['name']}-xa-{fraction_slug(fraction_a)}"
            ),
            "fraction_component_a": fraction_a,
            "config": build_entry_config(
                base,
                physical,
                study,
                fraction_a,
            ),
        }
        for fraction_a in fractions if fraction_a in eligible_fractions
    ]

    if tasks:
        with ProcessPoolExecutor(
            max_workers=min(args.workers, len(tasks))
        ) as pool:
            futures = [
                pool.submit(run_entry_task, task)
                for task in tasks
            ]
            for future in as_completed(futures):
                row = future.result()
                entry_rows.append(row)
                mass = row.get("coolant_used_kg")
                mass_text = (
                    f"{mass:.3f} kg"
                    if isinstance(mass, (int, float))
                    else "n/a"
                )
                print(
                    f"{row['candidate_id']}: "
                    f"{row.get('status')} {'score' if valid_score(row) else 'consumed_before_failure'}={mass_text}",
                    flush=True,
                )

    entry_rows.sort(
        key=lambda row: float(
            row["fraction_component_a"]
        )
    )
    score = entry_summary(
        entry_rows,
        reference_fraction,
    )
    write_csv(
        args.output_dir / "entry_results.csv",
        entry_rows,
    )

    report = {
        "study_complete": all(valid_score(row) for row in entry_rows)
            and score["comparable_candidate_count"] == len(fractions),
        "numerical_sweep_complete": all(valid_score(row) for row in entry_rows)
            and score["comparable_candidate_count"] == len(fractions),
        "physical_validation_complete": False,
        "caloric_validation": "pending experimental delta_h / excess-enthalpy validation",
        "endpoint_model_comparison": endpoint_comparison,
        "property_validation": continuity,
        "entry_scoring": score,
        "candidate_count": len(entry_rows),
        "thermo_grid_supported_fraction": (
            sum(
                row["state_points_supported"]
                for row in validations
            )
            / max(
                1,
                sum(
                    row["state_points_requested"]
                    for row in validations
                ),
            )
        ),
        "viscosity_coverage_fraction": (
            sum(
                row["transport_viscosity_supported"]
                for row in validations
            )
            / max(
                1,
                sum(
                    row["state_points_requested"]
                    for row in validations
                ),
            )
        ),
        "conductivity_coverage_fraction": (
            sum(
                row["transport_conductivity_supported"]
                for row in validations
            )
            / max(
                1,
                sum(
                    row["state_points_requested"]
                    for row in validations
                ),
            )
        ),
    }
    write_summary(
        args.output_dir / "summary.json",
        report,
    )
    print(json.dumps(report, indent=2))
    return 0 if report["study_complete"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
