#!/usr/bin/env python3
from __future__ import annotations

import argparse
import copy
import csv
import json
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
            result[key] = deep_merge(
                result[key],
                value,
            )
        else:
            result[key] = copy.deepcopy(value)
    return result


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
    base = yaml.safe_load(
        args.base_config.read_text(encoding="utf-8")
    )
    batch = yaml.safe_load(
        args.batch.read_text(encoding="utf-8")
    )
    cases = [
        deep_merge(base, item)
        for item in batch["cases"]
    ]
    args.output_dir.mkdir(
        parents=True,
        exist_ok=True,
    )

    summaries: list[dict[str, Any]] = []
    with ProcessPoolExecutor(
        max_workers=args.workers
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
