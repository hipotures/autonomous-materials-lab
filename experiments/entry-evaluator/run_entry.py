#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
from pathlib import Path

import yaml

from evaluator import evaluate, write_summary


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description=(
            "Run one low-fidelity Earth entry working-fluid evaluation."
        )
    )
    p.add_argument(
        "--config",
        type=Path,
        default=Path(__file__).with_name("config.yaml"),
    )
    p.add_argument(
        "--history",
        type=Path,
        default=Path("entry_history.csv"),
    )
    p.add_argument(
        "--summary",
        type=Path,
        default=Path("entry_summary.json"),
    )
    return p.parse_args()


def main() -> None:
    args = parse_args()
    config = yaml.safe_load(
        args.config.read_text(encoding="utf-8")
    )
    summary = evaluate(config, args.history)
    write_summary(args.summary, summary)
    print(json.dumps(summary, indent=2, sort_keys=True))
    print(f"\nHistory: {args.history}")
    print(f"Summary: {args.summary}")


if __name__ == "__main__":
    main()
