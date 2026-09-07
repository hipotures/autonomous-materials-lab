#!/usr/bin/env python3
"""Time-step and angular-grid verification with an explicit ranking stability gate."""
from __future__ import annotations

import argparse
import copy
import csv
import hashlib
import itertools
import json
import math
import platform
import subprocess
from concurrent.futures import ProcessPoolExecutor, as_completed
from importlib.metadata import version
from pathlib import Path

import yaml

from evaluator import evaluate, write_summary
from run_batch import deep_merge

HERE = Path(__file__).resolve().parent
# Relative tolerances plus absolute floors, chosen before running the study.
METRICS = {
    "coolant_used_kg": (0.01, 0.01),
    "peak_heat_flux_w_m2": (0.01, 1.0),
    "minimum_altitude_km": (0.001, 0.01),
    "peak_dynamic_pressure_pa": (0.01, 1.0),
    "peak_wall_temperature_k": (0.0, 0.5),
    "peak_coolant_flow_kg_s": (0.01, 0.01),
    "vehicle_incident_heat_mj": (0.01, 0.001),
}


def compare(coarse: dict, fine: dict) -> dict:
    """Invalid or differently terminated runs cannot pass numerical verification."""
    valid = coarse.get("status") == fine.get("status") == "terminal_velocity"
    details = {}
    for key, (rtol, atol) in METRICS.items():
        a, b = coarse.get(key), fine.get(key)
        finite = all(isinstance(x, (float, int)) and math.isfinite(x) for x in (a, b))
        if finite:
            a, b = float(a), float(b)
        delta = abs(a - b) if finite else None
        tolerance = atol + rtol * abs(b) if finite else None
        passed = finite and delta <= tolerance
        details[key] = {"absolute_difference": delta, "tolerance": tolerance,
                        "relative_difference": delta / abs(b) if finite and b != 0 else None,
                        "pass": passed}
        valid = valid and passed
    return {"pass": valid, "same_successful_termination": coarse.get("status") == fine.get("status") == "terminal_velocity",
            "metrics": details}


def run_one(task: dict) -> dict:
    config = task["config"]
    try:
        summary = evaluate(config)
    except Exception as exc:
        summary = {"status": "exception", "failure_reason": str(exc)}
    return {"id": task["id"], "profile": task["profile"], "case": task["case"],
            "dt_s": config["numerics"]["dt_s"], "zones": task["zones"],
            "kind": task["kind"], **summary}


def pair_order(a: dict, b: dict) -> str:
    if any(r.get("status") != "terminal_velocity" or not math.isfinite(r.get("coolant_used_kg", math.nan))
           or r["coolant_used_kg"] <= 0 for r in (a, b)):
        return "unavailable"
    # Treat small gaps as unresolved; do not invent a strict ordering for a tie.
    delta = a["coolant_used_kg"] - b["coolant_used_kg"]
    tolerance = 0.01 + 0.01 * max(a["coolant_used_kg"], b["coolant_used_kg"])
    if abs(delta) <= tolerance:
        return "unresolved"
    return "a_less" if delta < 0 else "b_less"


def analyze(rows: list[dict], case_names: list[str], profile_names: list[str], dts: list[float]) -> dict:
    index = {(r["case"], r["profile"], r["dt_s"], r["kind"]): r for r in rows}
    comparisons = []
    for case, profile in itertools.product(case_names, profile_names):
        for coarse, fine in zip(dts, dts[1:]):
            result = compare(index[case, profile, coarse, "time"], index[case, profile, fine, "time"])
            comparisons.append({"case": case, "profile": profile, "kind": "time",
                                "coarse_dt_s": coarse, "fine_dt_s": fine, **result})
        comparisons.append({"case": case, "profile": profile, "kind": "angular_grid",
                            **compare(index[case, profile, dts[-1], "time"],
                                      index[case, profile, dts[-1], "refined_grid"])})
    pairs = []
    for a, b in itertools.combinations(case_names, 2):
        orders = []
        for profile in profile_names:
            for dt in dts:
                orders.append(pair_order(index[a, profile, dt, "time"], index[b, profile, dt, "time"]))
            orders.append(pair_order(index[a, profile, dts[-1], "refined_grid"],
                                     index[b, profile, dts[-1], "refined_grid"]))
        stable = len(set(orders)) == 1 and orders[0] in {"a_less", "b_less"}
        pairs.append({"a": a, "b": b, "orders_observed": sorted(set(orders)), "stable": stable})
    numerical_pass = all(c["pass"] for c in comparisons)
    ranking_stable = bool(pairs) and all(p["stable"] for p in pairs)
    return {"numerical_pass": numerical_pass, "ranking_stable": ranking_stable,
            "gate_pass": numerical_pass and ranking_stable,
            "comparisons": comparisons, "ranking_pairs": pairs,
            "interpretation": "Numerical consistency under tested assumptions only; not physical validation or mission qualification."}


def write_csv(path: Path, rows: list[dict]) -> None:
    fields = sorted({key for row in rows for key in row})
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--base-config", type=Path, default=HERE / "config.yaml")
    parser.add_argument("--batch", type=Path, default=HERE / "batch.yaml")
    parser.add_argument("--study", type=Path, default=HERE / "convergence.yaml")
    parser.add_argument("--workers", type=int, default=4)
    parser.add_argument("--output-dir", type=Path, default=Path("convergence-results"))
    args = parser.parse_args()
    base = yaml.safe_load(args.base_config.read_text())
    batch = yaml.safe_load(args.batch.read_text())
    study = yaml.safe_load(args.study.read_text())
    dts = [float(x) for x in study["dt_s"]]
    if len(dts) < 2 or any(not math.isfinite(t) or t <= 0 for t in dts) or dts != sorted(set(dts), reverse=True):
        parser.error("study dt_s must contain at least two distinct finite positive steps in decreasing order")
    if args.workers <= 0:
        parser.error("workers must be positive")
    profiles = study["surface_profiles"]
    if not profiles:
        parser.error("at least one surface profile is required")
    # Deliberately ignore batch.matrix: this study fixes one entry for every fluid.
    cases = [deep_merge(base, c) for c in batch["cases"]]
    names = [str(c["name"]) for c in cases]
    if len(names) < 2 or len(names) != len(set(names)):
        parser.error("at least two uniquely named cases are required")
    fixed = []
    for c in cases:
        c.update({"entry": deep_merge(c["entry"], study["entry"])})
        if c["vehicle"].get("couple_coolant_mass_to_trajectory", False) or c["coolant"].get("available_mass_kg") is not None:
            parser.error("ranking study requires fixed trajectory mass and uncapped coolant")
        shared = {k: v for k, v in c.items() if k not in {"name", "coolant", "chemistry", "numerics", "surface"}}
        # Delivery limits, area and numerical horizon must also match.
        shared["delivery"] = {k: v for k, v in c["coolant"].items() if k not in {
            "coolprop_name", "storage_temperature_k", "storage_pressure_pa", "max_exit_temperature_k"}}
        shared["numerics"] = {k: v for k, v in c["numerics"].items() if k != "dt_s"}
        fixed.append(shared)
    if any(c != fixed[0] for c in fixed[1:]):
        parser.error("all fluid cases must share vehicle, entry, wall, delivery and numerical settings")
    if args.output_dir.exists() and any(args.output_dir.iterdir()):
        parser.error("output directory must be empty; use a new directory to preserve previous runs")
    args.output_dir.mkdir(parents=True, exist_ok=True)
    tasks = []
    for c, (profile, surface) in itertools.product(cases, profiles.items()):
        for dt, kind in [(t, "time") for t in dts] + [(dts[-1], "refined_grid")]:
            cfg = copy.deepcopy(c)
            cfg["surface"] = copy.deepcopy(surface)
            zones = int(surface.get("zones", 1))
            if kind == "refined_grid" and surface["model"] != "uniform":
                zones *= 2
                cfg["surface"]["zones"] = zones
            cfg["numerics"]["dt_s"] = dt
            task = {"id": f"run-{len(tasks):04d}", "case": cfg["name"], "profile": profile,
                    "kind": kind, "zones": zones, "config": cfg}
            tasks.append(task)
    try:
        commit = subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=HERE, text=True).strip()
    except (OSError, subprocess.CalledProcessError):
        commit = None
    manifest = {"git_commit": commit, "python": platform.python_version(),
                "packages": {p: version(p) for p in ("CoolProp", "pymsis", "numpy", "PyYAML")},
                "source_sha256": {p.name: hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(HERE.glob("*.py"))},
                "study": study, "metric_tolerances": METRICS, "tasks": tasks}
    write_summary(args.output_dir / "manifest.json", manifest)
    rows = []
    print(f"cases={len(tasks)} workers={min(args.workers, len(tasks))}; batch.matrix is ignored", flush=True)
    with ProcessPoolExecutor(max_workers=min(args.workers, len(tasks))) as pool:
        futures = [pool.submit(run_one, t) for t in tasks]
        for future in as_completed(futures):
            row = future.result()
            rows.append(row)
            write_summary(args.output_dir / f"{row['id']}.json", row)
            print(f"{row['id']} {row['case']} {row['profile']} dt={row['dt_s']} {row['kind']}: "
                  f"{row['status']} mass={row.get('coolant_used_kg')}", flush=True)
    rows.sort(key=lambda r: r["id"])
    report = analyze(rows, names, list(profiles), dts)
    write_summary(args.output_dir / "convergence_report.json", report)
    write_csv(args.output_dir / "all_runs.csv", rows)
    ranking = []
    for row in rows:
        if row["dt_s"] != dts[-1] or row["kind"] != "refined_grid":
            continue
        water = next((w for w in rows if w.get("fluid") == "Water" and w["profile"] == row["profile"]
                      and w["dt_s"] == dts[-1] and w["kind"] == "refined_grid"), None)
        comparable = (water is not None and row["status"] == water["status"] == "terminal_velocity"
                      and water.get("coolant_used_kg", 0) > 0)
        ranking.append({**row, **{f"entry_{k}": v for k, v in study["entry"].items()}, "mass_ratio_vs_water": row["coolant_used_kg"] / water["coolant_used_kg"] if comparable else None,
                        "numerical_study_pass": report["numerical_pass"], "ranking_stable": report["ranking_stable"],
                        "physical_validation": "not_established"})
    write_csv(args.output_dir / "ranking.csv", ranking)
    print(json.dumps({k: report[k] for k in ("numerical_pass", "ranking_stable", "gate_pass")}, indent=2))
    return 0 if report["gate_pass"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
