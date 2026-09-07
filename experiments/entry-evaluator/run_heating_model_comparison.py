#!/usr/bin/env python3
"""Compare legacy and Brandis-Johnston Earth heating on matched V2 scenarios."""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import platform
import subprocess
from collections import Counter, defaultdict
from concurrent.futures import ProcessPoolExecutor, as_completed
from importlib.metadata import version
from pathlib import Path
from typing import Any

import yaml

from evaluator import evaluate, write_summary
from run_batch import deep_merge, set_dotted
from run_physical_sensitivity import (
    apply_values,
    latin_hypercube,
    quantiles,
    safe_name,
    validate_ranges,
)

HERE = Path(__file__).resolve().parent
DEFAULT_BACKENDS = ["legacy", "brandis_johnston_2014"]
SUCCESS_STATUSES = {"terminal_velocity", "terminal_altitude"}


def finite_positive(value: Any) -> bool:
    return (
        isinstance(value, (int, float))
        and math.isfinite(float(value))
        and float(value) > 0.0
    )


def run_one(task: dict[str, Any]) -> dict[str, Any]:
    config = task["config"]
    try:
        summary = evaluate(config)
    except Exception as exc:
        summary = {
            "status": "exception",
            "failure_reason": str(exc),
            "heating_backend": task["backend"],
        }
    return {
        "run_id": task["run_id"],
        "scenario_id": task["scenario_id"],
        "case": task["case"],
        "fluid_expected": task["fluid_expected"],
        "backend_requested": task["backend"],
        **{
            f"input__{safe_name(key)}": value
            for key, value in task["values"].items()
        },
        **summary,
    }


def _ratio(a: Any, b: Any) -> float | None:
    if not finite_positive(a) or not finite_positive(b):
        return None
    return float(a) / float(b)


def analyze(
    rows: list[dict[str, Any]],
    scenario_ids: list[str],
    fluids: list[str],
    backends: list[str],
) -> tuple[dict[str, Any], list[dict[str, Any]], list[dict[str, Any]]]:
    index = {
        (
            str(row["scenario_id"]),
            str(row.get("fluid_expected") or row.get("fluid")),
            str(row.get("heating_backend") or row.get("backend_requested")),
        ): row
        for row in rows
    }

    expected = len(scenario_ids) * len(fluids) * len(backends)
    exceptions = sum(1 for row in rows if row.get("status") == "exception")
    missing = expected - len(rows)

    paired_rows: list[dict[str, Any]] = []
    mass_ratios_by_fluid: dict[str, list[float]] = defaultdict(list)
    heat_ratios_by_fluid: dict[str, list[float]] = defaultdict(list)
    peak_ratios_by_fluid: dict[str, list[float]] = defaultdict(list)
    same_status_by_fluid: Counter[str] = Counter()
    comparable_model_pairs_by_fluid: Counter[str] = Counter()

    if "legacy" in backends and "brandis_johnston_2014" in backends:
        for scenario_id in scenario_ids:
            for fluid in fluids:
                legacy = index.get((scenario_id, fluid, "legacy"))
                bj = index.get((scenario_id, fluid, "brandis_johnston_2014"))
                if legacy is None or bj is None:
                    continue
                same_status = legacy.get("status") == bj.get("status")
                if same_status:
                    same_status_by_fluid[fluid] += 1
                mass_ratio = heat_ratio = peak_ratio = None
                comparable = (
                    same_status
                    and legacy.get("status") in SUCCESS_STATUSES
                    and finite_positive(legacy.get("coolant_used_kg"))
                    and finite_positive(bj.get("coolant_used_kg"))
                )
                if comparable:
                    comparable_model_pairs_by_fluid[fluid] += 1
                    mass_ratio = _ratio(
                        bj.get("coolant_used_kg"),
                        legacy.get("coolant_used_kg"),
                    )
                    heat_ratio = _ratio(
                        bj.get("vehicle_incident_heat_mj"),
                        legacy.get("vehicle_incident_heat_mj"),
                    )
                    peak_ratio = _ratio(
                        bj.get("peak_heat_flux_w_m2"),
                        legacy.get("peak_heat_flux_w_m2"),
                    )
                    if mass_ratio is not None:
                        mass_ratios_by_fluid[fluid].append(mass_ratio)
                    if heat_ratio is not None:
                        heat_ratios_by_fluid[fluid].append(heat_ratio)
                    if peak_ratio is not None:
                        peak_ratios_by_fluid[fluid].append(peak_ratio)

                paired_rows.append({
                    "scenario_id": scenario_id,
                    "fluid": fluid,
                    "legacy_status": legacy.get("status"),
                    "bj_status": bj.get("status"),
                    "same_status": same_status,
                    "legacy_coolant_kg": legacy.get("coolant_used_kg"),
                    "bj_coolant_kg": bj.get("coolant_used_kg"),
                    "bj_over_legacy_coolant_mass": mass_ratio,
                    "legacy_vehicle_incident_heat_mj": legacy.get("vehicle_incident_heat_mj"),
                    "bj_vehicle_incident_heat_mj": bj.get("vehicle_incident_heat_mj"),
                    "bj_over_legacy_incident_heat": heat_ratio,
                    "legacy_peak_heat_flux_w_m2": legacy.get("peak_heat_flux_w_m2"),
                    "bj_peak_heat_flux_w_m2": bj.get("peak_heat_flux_w_m2"),
                    "bj_over_legacy_peak_heat_flux": peak_ratio,
                    "legacy_radiative_energy_valid_fraction": legacy.get(
                        "vehicle_radiative_energy_valid_fraction"
                    ),
                    "bj_radiative_energy_valid_fraction": bj.get(
                        "vehicle_radiative_energy_valid_fraction"
                    ),
                    "bj_convective_energy_valid_fraction": bj.get(
                        "vehicle_convective_energy_valid_fraction"
                    ),
                })

    ranking_rows: list[dict[str, Any]] = []
    h2_ratios_by_backend: dict[str, list[float]] = defaultdict(list)
    h2_outcomes: dict[tuple[str, str], str] = {}
    validity_by_backend: dict[str, list[float]] = defaultdict(list)
    conv_validity_by_backend: dict[str, list[float]] = defaultdict(list)

    for scenario_id in scenario_ids:
        for backend in backends:
            h2 = index.get((scenario_id, "Hydrogen", backend))
            water = index.get((scenario_id, "Water", backend))
            if h2 is None or water is None:
                continue
            comparable = (
                h2.get("status") == water.get("status")
                and h2.get("status") in SUCCESS_STATUSES
                and finite_positive(h2.get("coolant_used_kg"))
                and finite_positive(water.get("coolant_used_kg"))
            )
            ratio = (
                float(h2["coolant_used_kg"]) / float(water["coolant_used_kg"])
                if comparable
                else None
            )
            if ratio is None:
                outcome = "not_comparable"
            elif ratio < 0.99:
                outcome = "h2_better"
                h2_ratios_by_backend[backend].append(ratio)
            elif ratio > 1.01:
                outcome = "water_better"
                h2_ratios_by_backend[backend].append(ratio)
            else:
                outcome = "unresolved"
                h2_ratios_by_backend[backend].append(ratio)
            h2_outcomes[(scenario_id, backend)] = outcome

            water_rad_valid = water.get("vehicle_radiative_energy_valid_fraction")
            if isinstance(water_rad_valid, (int, float)) and math.isfinite(float(water_rad_valid)):
                validity_by_backend[backend].append(float(water_rad_valid))
            water_conv_valid = water.get("vehicle_convective_energy_valid_fraction")
            if isinstance(water_conv_valid, (int, float)) and math.isfinite(float(water_conv_valid)):
                conv_validity_by_backend[backend].append(float(water_conv_valid))

            ranking_rows.append({
                "scenario_id": scenario_id,
                "backend": backend,
                "status": h2.get("status"),
                "h2_coolant_kg": h2.get("coolant_used_kg"),
                "water_coolant_kg": water.get("coolant_used_kg"),
                "h2_water_mass_ratio": ratio,
                "h2_vs_water": outcome,
                "water_radiative_energy_valid_fraction": water_rad_valid,
                "water_convective_energy_valid_fraction": water_conv_valid,
            })

    common: list[str] = []
    preserved_h2 = reversed_to_water = unresolved_or_missing = 0
    shifts: list[float] = []
    if "legacy" in backends and "brandis_johnston_2014" in backends:
        legacy_ratios = {
            row["scenario_id"]: row["h2_water_mass_ratio"]
            for row in ranking_rows
            if row["backend"] == "legacy" and row["h2_water_mass_ratio"] is not None
        }
        bj_ratios = {
            row["scenario_id"]: row["h2_water_mass_ratio"]
            for row in ranking_rows
            if row["backend"] == "brandis_johnston_2014"
            and row["h2_water_mass_ratio"] is not None
        }
        for scenario_id in scenario_ids:
            if scenario_id not in legacy_ratios or scenario_id not in bj_ratios:
                continue
            common.append(scenario_id)
            shifts.append(bj_ratios[scenario_id] - legacy_ratios[scenario_id])
            lo = h2_outcomes[(scenario_id, "legacy")]
            bo = h2_outcomes[(scenario_id, "brandis_johnston_2014")]
            if lo == bo == "h2_better":
                preserved_h2 += 1
            elif lo == "h2_better" and bo == "water_better":
                reversed_to_water += 1
            else:
                unresolved_or_missing += 1

    report = {
        "study_version": "v3",
        "study_complete": exceptions == 0 and missing == 0,
        "physical_validation": "not_established",
        "interpretation": (
            "Cross-model engineering-correlation comparison on matched V2 scenarios. "
            "Agreement between correlations is not independent flight validation."
        ),
        "scenario_count": len(scenario_ids),
        "fluid_count": len(fluids),
        "backend_count": len(backends),
        "run_count": len(rows),
        "expected_run_count": expected,
        "exception_count": exceptions,
        "missing_run_count": missing,
        "model_pair_status": {
            fluid: {
                "same_status_scenarios": same_status_by_fluid[fluid],
                "comparable_mass_scenarios": comparable_model_pairs_by_fluid[fluid],
            }
            for fluid in fluids
        },
        "bj_over_legacy": {
            fluid: {
                "coolant_mass_ratio": quantiles(mass_ratios_by_fluid[fluid]),
                "incident_heat_ratio": quantiles(heat_ratios_by_fluid[fluid]),
                "peak_heat_flux_ratio": quantiles(peak_ratios_by_fluid[fluid]),
            }
            for fluid in fluids
        },
        "h2_water": {
            "by_backend": {
                backend: {
                    "comparable_scenarios": len(h2_ratios_by_backend[backend]),
                    "mass_ratio_quantiles": quantiles(h2_ratios_by_backend[backend]),
                    "h2_better_count": sum(
                        h2_outcomes.get((scenario_id, backend)) == "h2_better"
                        for scenario_id in scenario_ids
                    ),
                    "water_better_count": sum(
                        h2_outcomes.get((scenario_id, backend)) == "water_better"
                        for scenario_id in scenario_ids
                    ),
                    "unresolved_count": sum(
                        h2_outcomes.get((scenario_id, backend)) == "unresolved"
                        for scenario_id in scenario_ids
                    ),
                }
                for backend in backends
            },
            "common_comparable_scenarios": len(common),
            "h2_better_both_count": preserved_h2,
            "legacy_h2_to_bj_water_reversal_count": reversed_to_water,
            "other_or_unresolved_count": unresolved_or_missing,
            "ratio_shift_bj_minus_legacy": quantiles(shifts),
        },
        "validity": {
            backend: {
                "radiative_energy_valid_fraction": quantiles(validity_by_backend[backend]),
                "convective_energy_valid_fraction": quantiles(conv_validity_by_backend[backend]),
            }
            for backend in backends
        },
    }
    return report, paired_rows, ranking_rows


def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    fields = sorted({key for row in rows for key in row})
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--base-config", type=Path, default=HERE / "config.yaml")
    parser.add_argument("--batch", type=Path, default=HERE / "batch.yaml")
    parser.add_argument("--study", type=Path, default=HERE / "physical-sensitivity.yaml")
    parser.add_argument("--workers", type=int, default=16)
    parser.add_argument("--samples", type=int)
    parser.add_argument("--seed", type=int)
    parser.add_argument("--backends", nargs="+", default=DEFAULT_BACKENDS)
    parser.add_argument("--output-dir", type=Path, default=Path("heating-model-v3"))
    args = parser.parse_args()

    if args.workers <= 0:
        parser.error("--workers must be positive")
    if len(set(args.backends)) != len(args.backends):
        parser.error("--backends must be unique")
    if args.output_dir.exists() and any(args.output_dir.iterdir()):
        parser.error("output directory must be empty; use a new directory")
    args.output_dir.mkdir(parents=True, exist_ok=True)
    (args.output_dir / "runs").mkdir(exist_ok=True)

    base = yaml.safe_load(args.base_config.read_text(encoding="utf-8"))
    batch = yaml.safe_load(args.batch.read_text(encoding="utf-8"))
    study = yaml.safe_load(args.study.read_text(encoding="utf-8"))

    ranges = validate_ranges(study["ranges"])
    samples = int(args.samples if args.samples is not None else study.get("samples", 64))
    seed = int(args.seed if args.seed is not None else study.get("seed", 20260907))
    fixed_overrides = study.get("fixed_overrides", {}) or {}
    nominal_values = study.get("nominal_values", {}) or {}

    sampled = latin_hypercube(ranges, samples, seed)
    scenarios = [{"scenario_id": "nominal", "values": dict(nominal_values)}] + [
        {"scenario_id": f"lhs-{i:03d}", "values": values}
        for i, values in enumerate(sampled)
    ]

    cases = [deep_merge(base, case) for case in batch["cases"]]
    prepared_cases = []
    for case in cases:
        cfg = apply_values(case, fixed_overrides)
        if cfg["vehicle"].get("couple_coolant_mass_to_trajectory", False):
            parser.error("V3 requires fixed trajectory mass")
        if cfg["coolant"].get("available_mass_kg") is not None:
            parser.error("V3 requires uncapped coolant scoring")
        prepared_cases.append(cfg)

    fluids = [str(case["coolant"]["coolprop_name"]) for case in prepared_cases]
    if len(fluids) != len(set(fluids)):
        parser.error("fluid names in batch.yaml must be unique")

    tasks: list[dict[str, Any]] = []
    for scenario in scenarios:
        values = {**fixed_overrides, **scenario["values"]}
        for case in prepared_cases:
            for backend in args.backends:
                cfg = apply_values(case, scenario["values"])
                set_dotted(cfg, "heating.backend", backend)
                tasks.append({
                    "run_id": f"run-{len(tasks):05d}",
                    "scenario_id": scenario["scenario_id"],
                    "case": cfg["name"],
                    "fluid_expected": str(cfg["coolant"]["coolprop_name"]),
                    "backend": backend,
                    "values": {**values, "heating.backend": backend},
                    "config": cfg,
                })

    try:
        commit = subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=HERE, text=True).strip()
    except (OSError, subprocess.CalledProcessError):
        commit = None

    manifest = {
        "study_version": "v3",
        "git_commit": commit,
        "python": platform.python_version(),
        "packages": {p: version(p) for p in ("CoolProp", "pymsis", "numpy", "PyYAML")},
        "source_sha256": {
            p.name: hashlib.sha256(p.read_bytes()).hexdigest()
            for p in sorted(HERE.glob("*.py"))
        },
        "seed": seed,
        "samples": samples,
        "backends": args.backends,
        "ranges": ranges,
        "fixed_overrides": fixed_overrides,
        "nominal_values": nominal_values,
        "scenarios": scenarios,
        "fluids": fluids,
    }
    write_summary(args.output_dir / "manifest.json", manifest)

    print(
        f"scenarios={len(scenarios)} fluids={len(fluids)} "
        f"backends={len(args.backends)} runs={len(tasks)} "
        f"workers={min(args.workers, len(tasks))}",
        flush=True,
    )

    rows: list[dict[str, Any]] = []
    with ProcessPoolExecutor(max_workers=min(args.workers, len(tasks))) as pool:
        futures = [pool.submit(run_one, task) for task in tasks]
        for future in as_completed(futures):
            row = future.result()
            rows.append(row)
            write_summary(args.output_dir / "runs" / f"{row['run_id']}.json", row)
            mass = row.get("coolant_used_kg")
            mass_text = f"{mass:.3f} kg" if isinstance(mass, (int, float)) else "n/a"
            print(
                f"{row['run_id']} {row['scenario_id']} "
                f"{row.get('fluid', row['fluid_expected'])} "
                f"{row.get('heating_backend', row['backend_requested'])}: "
                f"{row.get('status')} coolant={mass_text}",
                flush=True,
            )

    rows.sort(key=lambda row: row["run_id"])
    scenario_ids = [str(s["scenario_id"]) for s in scenarios]
    report, paired_rows, ranking_rows = analyze(rows, scenario_ids, fluids, args.backends)

    write_csv(args.output_dir / "all_runs.csv", rows)
    write_csv(args.output_dir / "model_pairs.csv", paired_rows)
    write_csv(args.output_dir / "backend_ranking.csv", ranking_rows)
    write_summary(args.output_dir / "heating_model_report.json", report)

    print(json.dumps({
        "study_complete": report["study_complete"],
        "run_count": report["run_count"],
        "h2_water": report["h2_water"],
        "validity": report["validity"],
    }, indent=2))
    return 0 if report["study_complete"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
