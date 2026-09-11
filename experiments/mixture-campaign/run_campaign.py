#!/usr/bin/env python3
"""Run the broad mission campaign, or the preserved property-only protocol."""
from __future__ import annotations
import argparse
import json
import sys

# Preserve the old import-level API for existing integrations and regressions.
from run_property_campaign import execute, references, latest_snapshot


def main():
    parser = argparse.ArgumentParser(add_help=False)
    parser.add_argument('--mode', choices=('mission', 'mission-v4', 'properties'), default='mission')
    options, remaining = parser.parse_known_args()
    if options.mode == 'mission':
        from runtime_tuning import tune_mission_argv
        tuned, info = tune_mission_argv(remaining)
        if info['injected']:
            printable = {**info, 'available_memory_gib': (
                round(info['available_memory_bytes'] / 1024**3, 2)
                if info['available_memory_bytes'] is not None else None)}
            print('[runtime-tuning] ' + json.dumps(printable, sort_keys=True), flush=True)
        from mission_batch_campaign import main as mission_main
        return mission_main(tuned)
    if options.mode == 'mission-v4':
        from mission_campaign import main as legacy_main
        return legacy_main(remaining)
    from run_property_campaign import main as property_main
    original = sys.argv
    try:
        sys.argv = [original[0], *remaining]
        return property_main()
    finally:
        sys.argv = original


if __name__ == '__main__':
    raise SystemExit(main())
