#!/usr/bin/env python3
from __future__ import annotations

import argparse
import copy
import csv
import itertools
import json
import re
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path
from typing import Any

import yaml

from evaluator import evaluate, write_summary


def deep_merge(
    base: dict[str, Any],
    override: dict[str, Any],
) -> dict[str, Any]:
    result = copy.deepcopy(base)
    for key, value in override.items():
        if (
            isinstance(value, dict)
            and isinstance(result.get(key), dict)
        ):
            result[key] = deep_merge(result[key], value)
        else:
            result[key] = copy.deepcopy(value)
    return result


def set_dotted(
    config: dict[str, Any],
    path: str,
    value: Any,
) -> None:
    keys = path.split(".")
    node = config
    for key in keys[:-1]:
        if key not in node or not isinstance(node[key], dict):
            node[key] = {}
        node = node[key]
    node[keys[-1]] = value


def slug(value: Any) -> str:
    text = str(value).replace("-", "m").replace(".", "p")
    return re.sub(r"[^A-Za-z0-9_]+", "_", text).strip("_")


def expand_cases(
    base: dict[str, Any],
    batch: dict[str, Any],
) -> list[dict[str, Any]]:
    fluid_cases = batch.get("cases", [{}])
    matrix = batch.get("matrix", {}) or {}
    matrix_keys = list(matrix.keys())
    matrix_values = [matrix[key] for key in matrix_keys]
    combinations = (
        list(itertools.product(*matrix_values))
        if matrix_keys
        else [()]
    )

    expanded: list[dict[str, Any]] = []
    for case_override in fluid_cases:
        case = deep_merge(base, case_override)
        base_name = str(case.get("name", "case"))
        for values in combinations:
            item = copy.deepcopy(case)
            suffixes: list[str] = []
            for key, value in zip(matrix_keys, values):
                set_dotted(item, key, value)
                suffixes.append(
                    f"{key.split('.')[-1]}_{slug(value)}"
                )
            if suffixes:
                item["name"] = (
                    base_name
                    + "__"
                    + "__".join(suffixes)
                )
            expanded.append(item)
    return expanded


def worker(
    config: dict[str, Any],
    output_dir: str,
    write_history: bool,
) -> dict[str, Any]:
    name = str(
        config.get("name", "case")
    ).replace("/", "_")
    out = Path(output_dir)
    history = (
        out / f"{name}.history.csv"
        if write_history
        else None
    )
    summary = evaluate(config, history)
    write_summary(
        out / f"{name}.summary.json",
        summary,
    )
    return summary


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description=(
            "Run independent entry cases in parallel processes."
        )
    )
    p.add_argument(
        "--base-config",
        type=Path,
        default=Path(__file__).with_name("config.yaml"),
    )
    p.add_argument(
        "--batch",
        type=Path,
        default=Path(__file__).with_name("batch.yaml"),
    )
    p.add_argument(
        "--workers",
        type=int,
        default=16,
    )
    p.add_argument(
        "--output-dir",
        type=Path,
        default=Path("batch-results"),
    )
    p.add_argument(
        "--history",
        action="store_true",
        help="Write per-case trajectory history CSV files.",
    )
    return p.parse_args()


def main() -> None:
    args = parse_args()
    if args.workers <= 0:
        raise SystemExit("--workers must be positive")

    base = yaml.safe_load(
        args.base_config.read_text(encoding="utf-8")
    )
    batch = yaml.safe_load(
        args.batch.read_text(encoding="utf-8")
    )
    cases = expand_cases(base, batch)
    args.output_dir.mkdir(
        parents=True,
        exist_ok=True,
    )

    print(
        f"cases={len(cases)} "
        f"workers={min(args.workers, len(cases))}"
    )
    summaries: list[dict[str, Any]] = []
    with ProcessPoolExecutor(
        max_workers=min(args.workers, len(cases))
    ) as pool:
        futures = {
            pool.submit(
                worker,
                case,
                str(args.output_dir),
                args.history,
            ): case.get("name", "case")
            for case in cases
        }
        for future in as_completed(futures):
            name = futures[future]
            try:
                summary = future.result()
                summaries.append(summary)
                print(
                    f"{name}: {summary['status']} "
                    f"coolant={summary['coolant_used_kg']:.3f} kg"
                )
            except Exception as exc:
                summaries.append(
                    {
                        "case_name": name,
                        "status": "exception",
                        "failure_reason": str(exc),
                    }
                )
                print(f"{name}: EXCEPTION {exc}")

    summaries.sort(
        key=lambda row: str(row.get("case_name", ""))
    )
    keys = sorted(
        {
            key
            for row in summaries
            for key in row.keys()
        }
    )
    with (
        args.output_dir / "batch_summary.csv"
    ).open(
        "w",
        newline="",
        encoding="utf-8",
    ) as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=keys,
        )
        writer.writeheader()
        writer.writerows(summaries)

    (
        args.output_dir / "batch_summary.json"
    ).write_text(
        json.dumps(
            summaries,
            indent=2,
            sort_keys=True,
        )
        + "\n",
        encoding="utf-8",
    )
    print(
        "Batch summary: "
        f"{args.output_dir / 'batch_summary.csv'}"
    )


if __name__ == "__main__":
    main()
