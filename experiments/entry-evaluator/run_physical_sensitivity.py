#!/usr/bin/env python3
"""Deterministic physical-sensitivity study for the V1 entry evaluator."""
from __future__ import annotations

import argparse
import copy
import csv
import hashlib
import json
import math
import platform
import random
import subprocess
from collections import Counter, defaultdict
from concurrent.futures import ProcessPoolExecutor, as_completed
from importlib.metadata import version
from pathlib import Path
from typing import Any

import numpy as np
import yaml

from evaluator import evaluate, write_summary
from run_batch import deep_merge, set_dotted

HERE = Path(__file__).resolve().parent
SUCCESS_STATUSES = {"terminal_velocity", "terminal_altitude"}


def validate_ranges(ranges: dict[str, Any]) -> dict[str, tuple[float, float]]:
    parsed: dict[str, tuple[float, float]] = {}
    for key, value in ranges.items():
        if not isinstance(value, list) or len(value) != 2:
            raise ValueError(f"range {key!r} must be [min, max]")
        low, high = map(float, value)
        if not all(math.isfinite(x) for x in (low, high)) or not low < high:
            raise ValueError(f"range {key!r} must contain finite min < max")
        parsed[str(key)] = (low, high)
    if not parsed:
        raise ValueError("at least one sensitivity range is required")
    return parsed


def latin_hypercube(
    ranges: dict[str, tuple[float, float]],
    samples: int,
    seed: int,
) -> list[dict[str, float]]:
    if samples <= 0:
        raise ValueError("samples must be positive")
    rng = random.Random(seed)
    keys = sorted(ranges)
    columns: dict[str, list[float]] = {}
    for key in keys:
        low, high = ranges[key]
        unit = [(i + rng.random()) / samples for i in range(samples)]
        rng.shuffle(unit)
        columns[key] = [low + u * (high - low) for u in unit]
    return [
        {key: columns[key][i] for key in keys}
        for i in range(samples)
    ]


def apply_values(config: dict[str, Any], values: dict[str, Any]) -> dict[str, Any]:
    result = copy.deepcopy(config)
    for key, value in values.items():
        set_dotted(result, key, value)
    return result


def safe_name(path: str) -> str:
    return path.replace(".", "__")


def run_one(task: dict[str, Any]) -> dict[str, Any]:
    config = task["config"]
    try:
        summary = evaluate(config)
    except Exception as exc:
        summary = {
            "status": "exception",
            "failure_reason": str(exc),
        }
    return {
        "run_id": task["run_id"],
        "scenario_id": task["scenario_id"],
        "case": task["case"],
        **{
            f"input__{safe_name(key)}": value
            for key, value in task["values"].items()
        },
        **summary,
    }


def finite_positive(value: Any) -> bool:
    return (
        isinstance(value, (int, float))
        and math.isfinite(float(value))
        and float(value) > 0.0
    )


def quantiles(values: list[float]) -> dict[str, float] | None:
    if not values:
        return None
    a = np.asarray(values, dtype=float)
    return {
        "min": float(np.min(a)),
        "p05": float(np.percentile(a, 5)),
        "median": float(np.percentile(a, 50)),
        "p95": float(np.percentile(a, 95)),
        "max": float(np.max(a)),
    }


def spearman(x: list[float], y: list[float]) -> float | None:
    if len(x) < 3 or len(x) != len(y):
        return None
    xa = np.asarray(x, dtype=float)
    ya = np.asarray(y, dtype=float)
    if not np.all(np.isfinite(xa)) or not np.all(np.isfinite(ya)):
        return None
    if np.ptp(xa) == 0 or np.ptp(ya) == 0:
        return None

    def ranks(a: np.ndarray) -> np.ndarray:
        order = np.argsort(a, kind="mergesort")
        out = np.empty(len(a), dtype=float)
        i = 0
        while i < len(a):
            j = i + 1
            while j < len(a) and a[order[j]] == a[order[i]]:
                j += 1
            rank = 0.5 * (i + j - 1) + 1.0
            out[order[i:j]] = rank
            i = j
        return out

    xr = ranks(xa)
    yr = ranks(ya)
    corr = np.corrcoef(xr, yr)[0, 1]
    return float(corr) if math.isfinite(float(corr)) else None


def analyze(
    rows: list[dict[str, Any]],
    scenarios: list[dict[str, Any]],
    case_names: list[str],
    margin_fraction: float,
) -> tuple[dict[str, Any], list[dict[str, Any]]]:
    by_scenario: dict[str, list[dict[str, Any]]] = defaultdict(list)
    for row in rows:
        by_scenario[str(row["scenario_id"])].append(row)

    scenario_lookup = {s["scenario_id"]: s for s in scenarios}
    scenario_rows: list[dict[str, Any]] = []
    ratios: list[float] = []
    h2_masses: list[float] = []
    water_masses: list[float] = []
    comparable_inputs: dict[str, list[float]] = defaultdict(list)
    status_counts: Counter[str] = Counter()
    all_four_best: Counter[str] = Counter()
    pairwise_wins: Counter[str] = Counter()
    pairwise_comparable: Counter[str] = Counter()
    radiation_validity: list[float] = []
    radiation_fraction: list[float] = []

    h2_wins = h2_losses = h2_unresolved = comparable = 0
    complete_all_four = 0

    for scenario_id in sorted(by_scenario):
        group = by_scenario[scenario_id]
        by_fluid = {str(r.get("fluid")): r for r in group if r.get("fluid")}
        scenario = scenario_lookup[scenario_id]
        water = by_fluid.get("Water")
        h2 = by_fluid.get("Hydrogen")

        representative = water or (group[0] if group else {})
        status_counts[str(representative.get("status", "missing"))] += 1

        row_out: dict[str, Any] = {
            "scenario_id": scenario_id,
            **{
                f"input__{safe_name(key)}": value
                for key, value in scenario["values"].items()
            },
            "water_status": water.get("status") if water else None,
            "h2_status": h2.get("status") if h2 else None,
            "water_coolant_kg": water.get("coolant_used_kg") if water else None,
            "h2_coolant_kg": h2.get("coolant_used_kg") if h2 else None,
        }

        if water:
            rv = water.get("vehicle_radiative_energy_valid_fraction")
            rf = water.get("vehicle_radiative_energy_fraction")
            if isinstance(rv, (int, float)) and math.isfinite(float(rv)):
                radiation_validity.append(float(rv))
                row_out["radiative_energy_valid_fraction"] = float(rv)
            else:
                row_out["radiative_energy_valid_fraction"] = None
            if isinstance(rf, (int, float)) and math.isfinite(float(rf)):
                radiation_fraction.append(float(rf))
                row_out["radiative_energy_fraction"] = float(rf)
            else:
                row_out["radiative_energy_fraction"] = None

        pair_ok = (
            water is not None
            and h2 is not None
            and water.get("status") == h2.get("status")
            and water.get("status") in SUCCESS_STATUSES
            and finite_positive(water.get("coolant_used_kg"))
            and finite_positive(h2.get("coolant_used_kg"))
        )
        if pair_ok:
            comparable += 1
            ratio = float(h2["coolant_used_kg"]) / float(water["coolant_used_kg"])
            ratios.append(ratio)
            h2_masses.append(float(h2["coolant_used_kg"]))
            water_masses.append(float(water["coolant_used_kg"]))
            row_out["h2_water_mass_ratio"] = ratio
            if ratio < 1.0 - margin_fraction:
                outcome = "h2_better"
                h2_wins += 1
            elif ratio > 1.0 + margin_fraction:
                outcome = "water_better"
                h2_losses += 1
            else:
                outcome = "unresolved"
                h2_unresolved += 1
            row_out["h2_vs_water"] = outcome
            for key, value in scenario["values"].items():
                comparable_inputs[key].append(float(value))
        else:
            row_out["h2_water_mass_ratio"] = None
            row_out["h2_vs_water"] = "not_comparable"

        complete = [
            r for r in group
            if r.get("status") in SUCCESS_STATUSES
            and finite_positive(r.get("coolant_used_kg"))
        ]
        if (
            len(complete) == len(case_names)
            and len({r.get("status") for r in complete}) == 1
        ):
            complete_all_four += 1
            ordered = sorted(complete, key=lambda r: float(r["coolant_used_kg"]))
            row_out["all_fluid_order"] = " < ".join(str(r["fluid"]) for r in ordered)
            all_four_best[str(ordered[0]["fluid"])] += 1
            for i, a in enumerate(ordered):
                for b in ordered[i + 1:]:
                    key = f"{a['fluid']}<{b['fluid']}"
                    pairwise_wins[key] += 1
            fluids = sorted(str(r["fluid"]) for r in complete)
            for i, a in enumerate(fluids):
                for b in fluids[i + 1:]:
                    pairwise_comparable[f"{a}|{b}"] += 1
        else:
            row_out["all_fluid_order"] = None

        scenario_rows.append(row_out)

    ratio_array = np.asarray(ratios, dtype=float) if ratios else np.asarray([])
    sensitivity: dict[str, Any] = {}
    for key in sorted(comparable_inputs):
        x = comparable_inputs[key]
        sensitivity[key] = {
            "spearman_vs_h2_water_ratio": spearman(x, ratios),
            "spearman_vs_h2_mass_kg": spearman(x, h2_masses),
            "spearman_vs_water_mass_kg": spearman(x, water_masses),
        }

    exceptions = sum(1 for r in rows if r.get("status") == "exception")
    missing_runs = len(scenarios) * len(case_names) - len(rows)
    report = {
        "study_version": "v2",
        "study_complete": exceptions == 0 and missing_runs == 0,
        "physical_validation": "not_established",
        "interpretation": (
            "Deterministic sensitivity screening over declared parameter ranges; "
            "ranges are not probability distributions and results are not mission qualification."
        ),
        "scenario_count": len(scenarios),
        "fluid_case_count": len(case_names),
        "run_count": len(rows),
        "expected_run_count": len(scenarios) * len(case_names),
        "exception_count": exceptions,
        "missing_run_count": missing_runs,
        "scenario_status_counts": dict(sorted(status_counts.items())),
        "h2_water": {
            "comparable_scenarios": comparable,
            "win_margin_fraction": margin_fraction,
            "h2_better_count": h2_wins,
            "water_better_count": h2_losses,
            "unresolved_count": h2_unresolved,
            "h2_better_fraction_of_comparable": (
                h2_wins / comparable if comparable else None
            ),
            "mass_ratio_quantiles": quantiles(ratios),
        },
        "all_four": {
            "comparable_scenarios": complete_all_four,
            "best_fluid_counts": dict(sorted(all_four_best.items())),
            "pairwise_order_counts": dict(sorted(pairwise_wins.items())),
        },
        "radiation": {
            "vehicle_radiative_energy_valid_fraction_quantiles": quantiles(
                radiation_validity
            ),
            "vehicle_radiative_energy_fraction_quantiles": quantiles(
                radiation_fraction
            ),
            "fraction_with_at_least_50pct_radiative_energy_in_nominal_range": (
                sum(v >= 0.5 for v in radiation_validity) / len(radiation_validity)
                if radiation_validity
                else None
            ),
            "fraction_with_at_least_80pct_radiative_energy_in_nominal_range": (
                sum(v >= 0.8 for v in radiation_validity) / len(radiation_validity)
                if radiation_validity
                else None
            ),
        },
        "parameter_sensitivity": sensitivity,
    }
    if ratio_array.size:
        worst_index = int(np.argmax(ratio_array))
        best_index = int(np.argmin(ratio_array))
        comparable_rows = [r for r in scenario_rows if r["h2_water_mass_ratio"] is not None]
        report["h2_water"]["best_case"] = comparable_rows[best_index]
        report["h2_water"]["worst_case"] = comparable_rows[worst_index]
    return report, scenario_rows


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
    parser.add_argument(
        "--study",
        type=Path,
        default=HERE / "physical-sensitivity.yaml",
    )
    parser.add_argument("--workers", type=int, default=16)
    parser.add_argument("--samples", type=int)
    parser.add_argument("--seed", type=int)
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("physical-sensitivity-results"),
    )
    args = parser.parse_args()

    if args.workers <= 0:
        parser.error("--workers must be positive")
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
    margin = float(study.get("ranking_margin_fraction", 0.01))
    if not 0.0 <= margin < 1.0:
        parser.error("ranking_margin_fraction must be in [0, 1)")

    fixed_overrides = study.get("fixed_overrides", {}) or {}
    nominal_values = study.get("nominal_values", {}) or {}
    sampled = latin_hypercube(ranges, samples, seed)
    scenarios = [
        {"scenario_id": "nominal", "values": dict(nominal_values)}
    ] + [
        {"scenario_id": f"lhs-{i:03d}", "values": values}
        for i, values in enumerate(sampled)
    ]

    cases = [deep_merge(base, c) for c in batch["cases"]]
    names = [str(c["name"]) for c in cases]
    if len(names) < 2 or len(names) != len(set(names)):
        parser.error("batch must contain at least two uniquely named fluid cases")

    prepared_cases = []
    for case in cases:
        cfg = apply_values(case, fixed_overrides)
        if cfg["vehicle"].get("couple_coolant_mass_to_trajectory", False):
            parser.error("V2 requires fixed trajectory mass")
        if cfg["coolant"].get("available_mass_kg") is not None:
            parser.error("V2 requires uncapped coolant scoring")
        prepared_cases.append(cfg)

    tasks: list[dict[str, Any]] = []
    for scenario in scenarios:
        values = {**fixed_overrides, **scenario["values"]}
        for case in prepared_cases:
            cfg = apply_values(case, scenario["values"])
            task = {
                "run_id": f"run-{len(tasks):05d}",
                "scenario_id": scenario["scenario_id"],
                "case": cfg["name"],
                "values": values,
                "config": cfg,
            }
            tasks.append(task)

    try:
        commit = subprocess.check_output(
            ["git", "rev-parse", "HEAD"], cwd=HERE, text=True
        ).strip()
    except (OSError, subprocess.CalledProcessError):
        commit = None

    manifest = {
        "study_version": "v2",
        "git_commit": commit,
        "python": platform.python_version(),
        "packages": {
            p: version(p)
            for p in ("CoolProp", "pymsis", "numpy", "PyYAML")
        },
        "source_sha256": {
            p.name: hashlib.sha256(p.read_bytes()).hexdigest()
            for p in sorted(HERE.glob("*.py"))
        },
        "seed": seed,
        "samples": samples,
        "ranges": ranges,
        "fixed_overrides": fixed_overrides,
        "nominal_values": nominal_values,
        "ranking_margin_fraction": margin,
        "scenarios": scenarios,
        "fluid_cases": names,
    }
    write_summary(args.output_dir / "manifest.json", manifest)

    print(
        f"scenarios={len(scenarios)} fluids={len(prepared_cases)} "
        f"runs={len(tasks)} workers={min(args.workers, len(tasks))}",
        flush=True,
    )

    rows: list[dict[str, Any]] = []
    with ProcessPoolExecutor(max_workers=min(args.workers, len(tasks))) as pool:
        futures = [pool.submit(run_one, task) for task in tasks]
        for future in as_completed(futures):
            row = future.result()
            rows.append(row)
            write_summary(
                args.output_dir / "runs" / f"{row['run_id']}.json",
                row,
            )
            mass = row.get("coolant_used_kg")
            mass_text = f"{mass:.3f} kg" if isinstance(mass, (int, float)) else "n/a"
            print(
                f"{row['run_id']} {row['scenario_id']} "
                f"{row.get('fluid', row['case'])}: "
                f"{row.get('status')} coolant={mass_text}",
                flush=True,
            )

    rows.sort(key=lambda r: r["run_id"])
    report, scenario_rows = analyze(rows, scenarios, names, margin)

    write_csv(args.output_dir / "all_runs.csv", rows)
    write_csv(args.output_dir / "scenario_summary.csv", scenario_rows)
    write_summary(args.output_dir / "physical_sensitivity_report.json", report)
    write_summary(
        args.output_dir / "scenario_summary.json",
        scenario_rows,
    )

    print(
        json.dumps(
            {
                "study_complete": report["study_complete"],
                "comparable_scenarios": report["h2_water"]["comparable_scenarios"],
                "h2_better_fraction": report["h2_water"][
                    "h2_better_fraction_of_comparable"
                ],
                "mass_ratio_quantiles": report["h2_water"]["mass_ratio_quantiles"],
                "radiative_validity": report["radiation"][
                    "vehicle_radiative_energy_valid_fraction_quantiles"
                ],
            },
            indent=2,
        )
    )
    return 0 if report["study_complete"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
