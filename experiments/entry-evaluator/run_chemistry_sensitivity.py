#!/usr/bin/env python3
"""V4: H2 ignition-delay constrained coolant sensitivity study."""
from __future__ import annotations

import argparse
import csv
import hashlib
import itertools
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
SUCCESS_STATUSES = {"terminal_velocity", "terminal_altitude"}


def finite_positive(value: Any) -> bool:
    return (
        isinstance(value, (int, float))
        and math.isfinite(float(value))
        and float(value) > 0.0
    )


def chemistry_id(residence: float, phi: float, safety: float) -> str:
    return (
        f"res-{residence:.6g}_phi-{phi:g}_safety-{safety:g}"
        .replace(".", "p")
    )


def ignition_table_metadata(path: Path) -> dict[str, Any]:
    pressures: list[float] = []
    phis: list[float] = []
    temperatures: list[float] = []
    with path.open("r", newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle):
            pressures.append(float(row["pressure_bar"]))
            phis.append(float(row["phi"]))
            temperatures.append(float(row["temperature_initial_K"]))
    if not pressures:
        raise ValueError("ignition CSV contains no rows")
    return {
        "pressure_min_bar": min(pressures),
        "pressure_max_bar": max(pressures),
        "phis": sorted(set(phis)),
        "temperature_min_k": min(temperatures),
        "temperature_max_k": max(temperatures),
    }


def run_one(task: dict[str, Any]) -> dict[str, Any]:
    try:
        summary = evaluate(task["config"])
    except Exception as exc:
        summary = {
            "status": "exception",
            "failure_reason": str(exc),
        }
    return {
        "run_id": task["run_id"],
        "scenario_id": task["scenario_id"],
        "role": task["role"],
        "chemistry_id": task.get("chemistry_id"),
        **{
            f"input__{safe_name(key)}": value
            for key, value in task["values"].items()
        },
        **summary,
    }


def execute_tasks(
    tasks: list[dict[str, Any]],
    workers: int,
    output_dir: Path,
) -> list[dict[str, Any]]:
    if not tasks:
        return []
    rows: list[dict[str, Any]] = []
    with ProcessPoolExecutor(
        max_workers=min(workers, len(tasks))
    ) as pool:
        futures = [pool.submit(run_one, task) for task in tasks]
        for future in as_completed(futures):
            row = future.result()
            rows.append(row)
            write_summary(
                output_dir / "runs" / f"{row['run_id']}.json",
                row,
            )
            mass = row.get("coolant_used_kg")
            mass_text = (
                f"{mass:.3f} kg"
                if isinstance(mass, (int, float))
                else "n/a"
            )
            print(
                f"{row['run_id']} {row['scenario_id']} "
                f"{row['role']} {row.get('chemistry_id') or '-'}: "
                f"{row.get('status')} coolant={mass_text}",
                flush=True,
            )
    rows.sort(key=lambda row: row["run_id"])
    return rows


def coverage_preflight(
    water_rows: list[dict[str, Any]],
    ignition_meta: dict[str, Any],
    requested_phis: list[float],
) -> dict[str, Any]:
    successful = [
        row
        for row in water_rows
        if row.get("status") in SUCCESS_STATUSES
    ]
    mins = [
        float(row["minimum_active_cooling_surface_pressure_bar"])
        for row in successful
        if isinstance(
            row.get("minimum_active_cooling_surface_pressure_bar"),
            (int, float),
        )
    ]
    maxs = [
        float(row["maximum_active_cooling_surface_pressure_bar"])
        for row in successful
        if isinstance(
            row.get("maximum_active_cooling_surface_pressure_bar"),
            (int, float),
        )
    ]

    required_min = min(mins) if mins else None
    required_max = max(maxs) if maxs else None
    missing_phis = [
        phi
        for phi in sorted(set(requested_phis))
        if phi not in ignition_meta["phis"]
    ]

    pressure_ok = (
        (required_min is None or required_min >= ignition_meta["pressure_min_bar"])
        and
        (required_max is None or required_max <= ignition_meta["pressure_max_bar"])
    )
    pass_gate = pressure_ok and not missing_phis

    return {
        "pass": pass_gate,
        "successful_physical_scenarios": len(successful),
        "required_active_cooling_pressure_min_bar": required_min,
        "required_active_cooling_pressure_max_bar": required_max,
        "table_pressure_min_bar": ignition_meta["pressure_min_bar"],
        "table_pressure_max_bar": ignition_meta["pressure_max_bar"],
        "requested_phis": sorted(set(requested_phis)),
        "table_phis": ignition_meta["phis"],
        "missing_phis": missing_phis,
    }


def analyze(
    rows: list[dict[str, Any]],
    scenario_ids: list[str],
    chemistry_cases: list[dict[str, float]],
) -> tuple[dict[str, Any], list[dict[str, Any]]]:
    water = {
        str(row["scenario_id"]): row
        for row in rows
        if row["role"] == "water"
    }
    hydrogen = {
        (str(row["scenario_id"]), str(row["chemistry_id"])): row
        for row in rows
        if row["role"] == "hydrogen"
    }

    detail: list[dict[str, Any]] = []
    ratios_by_case: dict[str, list[float]] = defaultdict(list)
    margins_by_case: dict[str, list[float]] = defaultdict(list)
    exit_t_by_case: dict[str, list[float]] = defaultdict(list)
    wins: Counter[str] = Counter()
    losses: Counter[str] = Counter()
    unresolved: Counter[str] = Counter()
    noncomparable: Counter[str] = Counter()

    case_lookup = {
        chemistry_id(
            case["residence_time_s"],
            case["phi"],
            case["safety_factor"],
        ): case
        for case in chemistry_cases
    }

    for scenario_id in scenario_ids:
        water_row = water.get(scenario_id)
        for cid, case in case_lookup.items():
            h2_row = hydrogen.get((scenario_id, cid))
            comparable = (
                water_row is not None
                and h2_row is not None
                and water_row.get("status") == h2_row.get("status")
                and water_row.get("status") in SUCCESS_STATUSES
                and finite_positive(water_row.get("coolant_used_kg"))
                and finite_positive(h2_row.get("coolant_used_kg"))
            )

            ratio = None
            outcome = "not_comparable"
            if comparable:
                ratio = (
                    float(h2_row["coolant_used_kg"])
                    / float(water_row["coolant_used_kg"])
                )
                ratios_by_case[cid].append(ratio)

                margin = h2_row.get("minimum_ignition_margin")
                if isinstance(margin, (int, float)) and math.isfinite(float(margin)):
                    margins_by_case[cid].append(float(margin))

                exit_t = h2_row.get("maximum_coolant_exit_temperature_k")
                if isinstance(exit_t, (int, float)) and math.isfinite(float(exit_t)):
                    exit_t_by_case[cid].append(float(exit_t))

                if ratio < 0.99:
                    outcome = "h2_better"
                    wins[cid] += 1
                elif ratio > 1.01:
                    outcome = "water_better"
                    losses[cid] += 1
                else:
                    outcome = "unresolved"
                    unresolved[cid] += 1
            else:
                noncomparable[cid] += 1

            detail.append({
                "scenario_id": scenario_id,
                "chemistry_id": cid,
                **case,
                "water_status": None if water_row is None else water_row.get("status"),
                "h2_status": None if h2_row is None else h2_row.get("status"),
                "water_coolant_kg": None if water_row is None else water_row.get("coolant_used_kg"),
                "h2_coolant_kg": None if h2_row is None else h2_row.get("coolant_used_kg"),
                "h2_water_mass_ratio": ratio,
                "outcome": outcome,
                "minimum_ignition_margin": None if h2_row is None else h2_row.get("minimum_ignition_margin"),
                "maximum_h2_exit_temperature_k": None if h2_row is None else h2_row.get("maximum_coolant_exit_temperature_k"),
                "failure_reason": None if h2_row is None else h2_row.get("failure_reason"),
            })

    per_case: dict[str, Any] = {}
    for cid, case in case_lookup.items():
        ratios = ratios_by_case[cid]
        per_case[cid] = {
            **case,
            "required_ignition_delay_s": (
                case["residence_time_s"] * case["safety_factor"]
            ),
            "comparable_scenarios": len(ratios),
            "h2_better_count": wins[cid],
            "water_better_count": losses[cid],
            "unresolved_count": unresolved[cid],
            "not_comparable_count": noncomparable[cid],
            "h2_better_fraction": wins[cid] / len(ratios) if ratios else None,
            "h2_water_mass_ratio": quantiles(ratios),
            "minimum_ignition_margin": quantiles(margins_by_case[cid]),
            "maximum_h2_exit_temperature_k": quantiles(exit_t_by_case[cid]),
        }

    comparable_cases = [
        (cid, data)
        for cid, data in per_case.items()
        if data["comparable_scenarios"] > 0
    ]
    robust = [
        (cid, data)
        for cid, data in comparable_cases
        if (
            data["h2_better_count"] == data["comparable_scenarios"]
            and data["not_comparable_count"] == 0
        )
    ]
    reversals = [
        (cid, data)
        for cid, data in comparable_cases
        if data["water_better_count"] > 0
    ]
    worst = max(
        comparable_cases,
        key=lambda item: item[1]["h2_water_mass_ratio"]["max"],
        default=None,
    )

    expected = len(scenario_ids) * (len(chemistry_cases) + 1)
    exceptions = sum(1 for row in rows if row.get("status") == "exception")
    report = {
        "study_version": "v4",
        "study_complete": exceptions == 0 and len(rows) == expected,
        "physical_validation": "not_established",
        "interpretation": (
            "Homogeneous ignition-delay constrained H2 screening. "
            "The Cantera lookup does not model mixing, shock-layer radicals, "
            "wall catalysis or reacting boundary-layer transport."
        ),
        "physical_scenario_count": len(scenario_ids),
        "chemistry_case_count": len(chemistry_cases),
        "run_count": len(rows),
        "expected_run_count": expected,
        "exception_count": exceptions,
        "chemistry_cases": per_case,
        "robust_h2_better_case_count": len(robust),
        "water_better_somewhere_case_count": len(reversals),
        "worst_h2_water_case": (
            None
            if worst is None
            else {"chemistry_id": worst[0], **worst[1]}
        ),
        "cases_with_water_better_somewhere": [
            {"chemistry_id": cid, **data}
            for cid, data in reversals
        ],
    }
    return report, detail


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
    parser.add_argument("--physical-study", type=Path, default=HERE / "physical-sensitivity.yaml")
    parser.add_argument("--chemistry-study", type=Path, default=HERE / "chemistry-sensitivity.yaml")
    parser.add_argument(
        "--ignition-csv",
        type=Path,
        default=HERE.parent / "hydrogen-ignition-delay" / "ignition_delay_v4.csv",
    )
    parser.add_argument("--workers", type=int, default=16)
    parser.add_argument(
        "--physical-samples",
        type=int,
        default=64,
        help="64 reproduces the V2 LHS plus nominal; 0 runs nominal only.",
    )
    parser.add_argument("--seed", type=int)
    parser.add_argument("--output-dir", type=Path, default=Path("chemistry-v4"))
    args = parser.parse_args()

    if args.workers <= 0:
        parser.error("--workers must be positive")
    if args.physical_samples < 0:
        parser.error("--physical-samples must be nonnegative")
    if not args.ignition_csv.exists():
        parser.error(
            f"ignition table not found: {args.ignition_csv}; "
            "generate it with ../hydrogen-ignition-delay/generate_v4_table.py"
        )
    if args.output_dir.exists() and any(args.output_dir.iterdir()):
        parser.error("output directory must be empty; use a new directory")

    args.output_dir.mkdir(parents=True, exist_ok=True)
    (args.output_dir / "runs").mkdir(exist_ok=True)

    base = yaml.safe_load(args.base_config.read_text(encoding="utf-8"))
    batch = yaml.safe_load(args.batch.read_text(encoding="utf-8"))
    physical = yaml.safe_load(args.physical_study.read_text(encoding="utf-8"))
    chemistry = yaml.safe_load(args.chemistry_study.read_text(encoding="utf-8"))

    ranges = validate_ranges(physical["ranges"])
    seed = int(args.seed if args.seed is not None else physical.get("seed", 20260907))
    fixed_overrides = physical.get("fixed_overrides", {}) or {}
    nominal_values = physical.get("nominal_values", {}) or {}

    sampled = (
        latin_hypercube(ranges, args.physical_samples, seed)
        if args.physical_samples > 0
        else []
    )
    scenarios = [{"scenario_id": "nominal", "values": dict(nominal_values)}] + [
        {"scenario_id": f"lhs-{i:03d}", "values": values}
        for i, values in enumerate(sampled)
    ]

    chemistry_cases = [
        {
            "residence_time_s": float(residence),
            "phi": float(phi),
            "safety_factor": float(safety),
        }
        for residence, phi, safety in itertools.product(
            chemistry["residence_time_s"],
            chemistry["phi"],
            chemistry["safety_factor"],
        )
    ]
    backend = str(chemistry.get("heating_backend", "brandis_johnston_2014"))

    merged_cases = [deep_merge(base, case) for case in batch["cases"]]
    h2_base = next(
        case for case in merged_cases
        if case["coolant"]["coolprop_name"] == "Hydrogen"
    )
    water_base = next(
        case for case in merged_cases
        if case["coolant"]["coolprop_name"] == "Water"
    )

    for original in (h2_base, water_base):
        checked = apply_values(original, fixed_overrides)
        if checked["vehicle"].get("couple_coolant_mass_to_trajectory", False):
            parser.error("V4 requires fixed trajectory mass")
        if checked["coolant"].get("available_mass_kg") is not None:
            parser.error("V4 requires uncapped coolant scoring")

    ignition_meta = ignition_table_metadata(args.ignition_csv)

    water_tasks: list[dict[str, Any]] = []
    h2_tasks: list[dict[str, Any]] = []
    for scenario in scenarios:
        physical_values = {**fixed_overrides, **scenario["values"]}

        water_cfg = apply_values(water_base, physical_values)
        set_dotted(water_cfg, "heating.backend", backend)
        set_dotted(water_cfg, "chemistry.mode", "disabled")
        water_tasks.append({
            "run_id": f"water-{len(water_tasks):04d}",
            "scenario_id": scenario["scenario_id"],
            "role": "water",
            "values": {**physical_values, "heating.backend": backend},
            "config": water_cfg,
        })

        for case in chemistry_cases:
            h2_cfg = apply_values(h2_base, physical_values)
            set_dotted(h2_cfg, "heating.backend", backend)
            set_dotted(h2_cfg, "chemistry.mode", "ignition_csv")
            set_dotted(h2_cfg, "chemistry.ignition_csv", str(args.ignition_csv.resolve()))
            set_dotted(h2_cfg, "chemistry.phi", case["phi"])
            set_dotted(h2_cfg, "chemistry.residence_time_s", case["residence_time_s"])
            set_dotted(h2_cfg, "chemistry.ignition_safety_factor", case["safety_factor"])
            set_dotted(h2_cfg, "chemistry.allow_nearest_phi", False)
            cid = chemistry_id(
                case["residence_time_s"],
                case["phi"],
                case["safety_factor"],
            )
            h2_tasks.append({
                "run_id": f"h2-{len(h2_tasks):05d}",
                "scenario_id": scenario["scenario_id"],
                "role": "hydrogen",
                "chemistry_id": cid,
                "values": {
                    **physical_values,
                    "heating.backend": backend,
                    "chemistry.residence_time_s": case["residence_time_s"],
                    "chemistry.phi": case["phi"],
                    "chemistry.ignition_safety_factor": case["safety_factor"],
                },
                "config": h2_cfg,
            })

    try:
        commit = subprocess.check_output(
            ["git", "rev-parse", "HEAD"], cwd=HERE, text=True
        ).strip()
    except (OSError, subprocess.CalledProcessError):
        commit = None

    manifest = {
        "study_version": "v4",
        "git_commit": commit,
        "python": platform.python_version(),
        "packages": {
            package: version(package)
            for package in ("CoolProp", "pymsis", "numpy", "PyYAML")
        },
        "ignition_csv": str(args.ignition_csv.resolve()),
        "ignition_csv_sha256": hashlib.sha256(args.ignition_csv.read_bytes()).hexdigest(),
        "ignition_table": ignition_meta,
        "seed": seed,
        "physical_samples": args.physical_samples,
        "physical_scenarios": scenarios,
        "chemistry_cases": chemistry_cases,
        "heating_backend": backend,
    }
    write_summary(args.output_dir / "manifest.json", manifest)

    print(
        f"preflight: physical_scenarios={len(scenarios)} "
        f"water_runs={len(water_tasks)} workers={min(args.workers, len(water_tasks))}",
        flush=True,
    )
    water_rows = execute_tasks(water_tasks, args.workers, args.output_dir)

    preflight = coverage_preflight(
        water_rows,
        ignition_meta,
        [case["phi"] for case in chemistry_cases],
    )
    write_summary(args.output_dir / "coverage_preflight.json", preflight)
    print(json.dumps({"coverage_preflight": preflight}, indent=2), flush=True)

    if not preflight["pass"]:
        write_csv(args.output_dir / "all_runs.csv", water_rows)
        print(
            "V4 stopped before H2 sweep because the ignition table does not "
            "cover the physical pressure/phi envelope.",
            flush=True,
        )
        return 2

    print(
        f"chemistry sweep: h2_runs={len(h2_tasks)} "
        f"workers={min(args.workers, len(h2_tasks))}",
        flush=True,
    )
    h2_rows = execute_tasks(h2_tasks, args.workers, args.output_dir)
    rows = sorted(water_rows + h2_rows, key=lambda row: row["run_id"])

    scenario_ids = [str(scenario["scenario_id"]) for scenario in scenarios]
    report, detail = analyze(rows, scenario_ids, chemistry_cases)
    report["coverage_preflight"] = preflight

    write_csv(args.output_dir / "all_runs.csv", rows)
    write_csv(args.output_dir / "chemistry_comparison.csv", detail)
    write_summary(args.output_dir / "chemistry_report.json", report)

    print(json.dumps({
        "study_complete": report["study_complete"],
        "run_count": report["run_count"],
        "robust_h2_better_case_count": report["robust_h2_better_case_count"],
        "water_better_somewhere_case_count": report["water_better_somewhere_case_count"],
        "worst_h2_water_case": report["worst_h2_water_case"],
    }, indent=2))
    return 0 if report["study_complete"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
