#!/usr/bin/env python3
"""Generate the dense H2/air ignition-delay table used by entry-evaluator V4."""
from __future__ import annotations

import argparse
import json
import platform
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

import cantera as ct

from ignition_delay import evaluate_case, write_csv

HERE = Path(__file__).resolve().parent


def default_temperatures() -> list[float]:
    return [float(value) for value in range(600, 1101, 25)]


def default_pressures() -> list[float]:
    return [
        0.001, 0.003, 0.01, 0.03, 0.05, 0.1, 0.2,
        0.3, 0.5, 0.75, 1.0, 1.5, 2.0, 3.0,
        5.0, 8.0, 12.0, 15.0, 20.0, 30.0,
    ]


def run_case(task: tuple):
    return evaluate_case(*task)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--workers", type=int, default=16)
    parser.add_argument(
        "--output",
        type=Path,
        default=HERE / "ignition_delay_v4.csv",
    )
    parser.add_argument(
        "--manifest",
        type=Path,
        default=HERE / "ignition_delay_v4.manifest.json",
    )
    parser.add_argument("--mechanism", default="h2o2.yaml")
    parser.add_argument("--max-time-s", type=float, default=1.0)
    parser.add_argument("--ignition-rise-k", type=float, default=200.0)
    parser.add_argument("--advance-limit-k", type=float, default=5.0)
    parser.add_argument(
        "--temperatures-k",
        nargs="+",
        type=float,
        default=default_temperatures(),
    )
    parser.add_argument(
        "--pressures-bar",
        nargs="+",
        type=float,
        default=default_pressures(),
    )
    parser.add_argument(
        "--phi",
        nargs="+",
        type=float,
        default=[2.0, 4.0, 8.0, 16.0],
    )
    args = parser.parse_args()

    if args.workers <= 0:
        parser.error("--workers must be positive")
    if args.max_time_s <= 0.0:
        parser.error("--max-time-s must be positive")

    for name, values in (
        ("temperatures", args.temperatures_k),
        ("pressures", args.pressures_bar),
        ("phi", args.phi),
    ):
        if not values or any(value <= 0.0 for value in values):
            parser.error(f"{name} must contain positive values")

    tasks = [
        (
            args.mechanism,
            temperature_k,
            pressure_bar,
            phi,
            "constant-pressure",
            args.max_time_s,
            args.ignition_rise_k,
            args.advance_limit_k,
        )
        for temperature_k in args.temperatures_k
        for pressure_bar in args.pressures_bar
        for phi in args.phi
    ]

    print(
        f"Cantera {ct.__version__} | cases={len(tasks)} | "
        f"workers={min(args.workers, len(tasks))}",
        flush=True,
    )

    results = []
    with ProcessPoolExecutor(
        max_workers=min(args.workers, len(tasks))
    ) as pool:
        futures = [
            pool.submit(run_case, task)
            for task in tasks
        ]
        for index, future in enumerate(
            as_completed(futures),
            start=1,
        ):
            results.append(future.result())
            if index % 50 == 0 or index == len(tasks):
                print(
                    f"completed {index}/{len(tasks)}",
                    flush=True,
                )

    results.sort(
        key=lambda result: (
            result.temperature_initial_k,
            result.pressure_bar,
            result.phi,
        )
    )
    write_csv(args.output, results)

    manifest = {
        "study": "V4 H2 homogeneous ignition-delay lookup",
        "python": platform.python_version(),
        "cantera": ct.__version__,
        "mechanism": args.mechanism,
        "reactor": "constant-pressure",
        "max_time_s": args.max_time_s,
        "ignition_rise_k": args.ignition_rise_k,
        "advance_limit_k": args.advance_limit_k,
        "temperatures_k": args.temperatures_k,
        "pressures_bar": args.pressures_bar,
        "phi": args.phi,
        "cases": len(tasks),
        "interpretation": (
            "Zero-dimensional homogeneous H2/air ignition screening. "
            "No mixing, shock-layer radicals, wall chemistry or CFD."
        ),
    }
    args.manifest.write_text(
        json.dumps(
            manifest,
            indent=2,
            sort_keys=True,
        )
        + "\n",
        encoding="utf-8",
    )

    print(f"CSV: {args.output}")
    print(f"manifest: {args.manifest}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
