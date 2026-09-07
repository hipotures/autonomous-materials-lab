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
    require_v5a_coolprop,
    saturation_to_dict,
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
        "type": "coolprop",
        "backend": "HEOS",
        "components": list(study["components"]),
        "composition_basis": str(
            study.get("composition_basis", "mole")
        ),
        "fractions": [
            fraction_a,
            1.0 - fraction_a,
        ],
        "stability_algorithm": int(
            study.get("stability_algorithm", 1)
        ),
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
        values = [
            (row["fraction_component_a"], row.get(key))
            for row in ordered
            if isinstance(row.get(key), (int, float))
            and math.isfinite(float(row[key]))
        ]
        changes = []
        for (x0, y0), (x1, y1) in zip(values, values[1:]):
            scale = max(
                abs(float(y0)),
                abs(float(y1)),
                1e-30,
            )
            changes.append(
                {
                    "x0": x0,
                    "x1": x1,
                    "relative_change": (
                        abs(float(y1) - float(y0)) / scale
                    ),
                }
            )
        metrics[key] = {
            "supported_compositions": len(values),
            "max_adjacent_relative_change": (
                max(
                    (
                        item["relative_change"]
                        for item in changes
                    ),
                    default=None,
                )
            ),
        }

    enthalpy_values = [
        (
            row["fraction_component_a"],
            row.get("storage_enthalpy_j_kg"),
        )
        for row in ordered
        if isinstance(
            row.get("storage_enthalpy_j_kg"),
            (int, float),
        )
    ]
    enthalpy_changes = [
        abs(float(b[1]) - float(a[1]))
        for a, b in zip(
            enthalpy_values,
            enthalpy_values[1:],
        )
    ]
    metrics["storage_enthalpy_j_kg"] = {
        "supported_compositions": len(enthalpy_values),
        "max_adjacent_absolute_change_j_kg": (
            max(enthalpy_changes)
            if enthalpy_changes
            else None
        ),
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


def entry_summary(
    rows: list[dict[str, Any]],
    reference_fraction: float,
) -> dict[str, Any]:
    reference = next(
        (
            row
            for row in rows
            if math.isclose(
                float(row["fraction_component_a"]),
                reference_fraction,
                rel_tol=0.0,
                abs_tol=1e-12,
            )
        ),
        None,
    )
    reference_mass = (
        float(reference["coolant_used_kg"])
        if reference is not None
        and reference.get("status") in SUCCESS_STATUSES
        and isinstance(reference.get("coolant_used_kg"), (int, float))
        and float(reference["coolant_used_kg"]) > 0.0
        else None
    )

    comparable: list[dict[str, Any]] = []
    for row in rows:
        mass = row.get("coolant_used_kg")
        ratio = None
        if (
            reference_mass is not None
            and row.get("status") == reference.get("status")
            and row.get("status") in SUCCESS_STATUSES
            and isinstance(mass, (int, float))
            and float(mass) > 0.0
        ):
            ratio = float(mass) / reference_mass
            comparable.append(row)
        row["coolant_mass_ratio_vs_reference"] = ratio

    best = (
        min(
            comparable,
            key=lambda row: float(row["coolant_used_kg"]),
        )
        if comparable
        else None
    )

    failed = [
        {
            "candidate_id": row["candidate_id"],
            "fraction_component_a": row["fraction_component_a"],
            "status": row.get("status"),
            "coolant_used_kg": row.get("coolant_used_kg"),
            "failure_reason": row.get("failure_reason"),
        }
        for row in rows
        if row.get("status") not in SUCCESS_STATUSES
    ]

    return {
        "reference_fraction_component_a": reference_fraction,
        "reference_candidate_id": (
            None if reference is None else reference["candidate_id"]
        ),
        "reference_coolant_kg": reference_mass,
        "comparable_candidate_count": len(comparable),
        "beats_reference_count": sum(
            1
            for row in rows
            if isinstance(
                row.get("coolant_mass_ratio_vs_reference"),
                (int, float),
            )
            and float(row["coolant_mass_ratio_vs_reference"]) < 0.99
        ),
        "failed_candidate_count": len(failed),
        "failed_candidates": failed,
        "best_candidate": (
            None
            if best is None
            else {
                "candidate_id": best["candidate_id"],
                "fraction_component_a": best[
                    "fraction_component_a"
                ],
                "coolant_used_kg": best["coolant_used_kg"],
                "coolant_mass_ratio_vs_reference": best[
                    "coolant_mass_ratio_vs_reference"
                ],
            }
        ),
    }


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
        "--stability-algorithm",
        type=int,
        choices=[0, 1],
        default=None,
        help=(
            "Override CoolProp mixture PT stability solver: "
            "1=Michelsen (v8 default), 0=legacy Gernert."
        ),
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("mixture-v5a-results"),
    )
    args = parser.parse_args()

    if args.workers <= 0:
        parser.error("--workers must be positive")

    try:
        require_v5a_coolprop()
    except RuntimeError as exc:
        parser.error(str(exc))
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

    if args.stability_algorithm is not None:
        study["stability_algorithm"] = args.stability_algorithm

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
        print(
            "V5a stopped: pure endpoint collapse validation failed",
            flush=True,
        )
        return 2
    if not continuity["all_storage_states_liquid"]:
        non_liquid = [
            {
                "candidate_id": row["candidate_id"],
                "fraction_component_a": row["fraction_component_a"],
                "storage_phase": row["storage_phase"],
                "failure_reason": row["failure_reason"],
            }
            for row in validations
            if not row["storage_liquid"]
        ]
        print(
            "V5a stopped: at least one composition is not liquid at "
            "the configured storage state",
            flush=True,
        )
        print(
            json.dumps(
                {"non_liquid_storage_candidates": non_liquid},
                indent=2,
            ),
            flush=True,
        )
        return 2

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
        for fraction_a in fractions
    ]

    entry_rows: list[dict[str, Any]] = []
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
                f"{row.get('status')} coolant={mass_text}",
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

    try:
        commit = subprocess.check_output(
            ["git", "rev-parse", "HEAD"],
            cwd=HERE,
            text=True,
        ).strip()
    except (
        OSError,
        subprocess.CalledProcessError,
    ):
        commit = None

    manifest = {
        "study_version": "v5a",
        "git_commit": commit,
        "python": platform.python_version(),
        "packages": {
            package: version(package)
            for package in (
                "CoolProp",
                "pymsis",
                "numpy",
                "PyYAML",
            )
        },
        "coolprop_runtime": coolprop_runtime_info(),
        "study_file": str(args.study.resolve()),
        "study_sha256": hashlib.sha256(
            args.study.read_bytes()
        ).hexdigest(),
        "components": components,
        "composition_basis": study.get(
            "composition_basis",
            "mole",
        ),
        "fractions_component_a": fractions,
        "reference_fraction_component_a": reference_fraction,
        "heating_backend": study.get(
            "heating_backend",
            "brandis_johnston_2014",
        ),
        "chemistry_mode": "disabled",
        "mixture_stability_algorithm": int(
            study.get("stability_algorithm", 1)
        ),
        "mixture_stability_algorithm_name": (
            "Michelsen"
            if int(study.get("stability_algorithm", 1)) == 1
            else "legacy_Gernert"
        ),
        "important_limitation": (
            "V5a entry scoring uses mixture enthalpy in the coolant energy "
            "balance. Density/cp/viscosity/conductivity/phase-envelope data "
            "are recorded for validation but are not yet coupled to porous "
            "flow, film cooling, tank mass or decomposition chemistry."
        ),
    }
    write_summary(
        args.output_dir / "manifest.json",
        manifest,
    )

    report = {
        "study_complete": all(
            row.get("status") in SUCCESS_STATUSES
            for row in entry_rows
        ),
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
