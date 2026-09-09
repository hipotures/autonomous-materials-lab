#!/usr/bin/env python3
"""Reanalyze a completed V5m-2 grid; publish verified gzip archives and small reports."""
from __future__ import annotations

import argparse
from contextlib import contextmanager
import json
from pathlib import Path
import platform
import sys
import tempfile
import zlib

from v5m21_review import (SOURCE_FILES, VERSION, calculate_review, compress_file, sha256,
                          validate_source, write_json, write_table)

HERE = Path(__file__).resolve().parent


@contextmanager
def source_lock(source: Path):
    """Cooperate with the original runner without creating or changing source files."""
    path = source / ".run.lock"
    if path.is_symlink():
        raise ValueError("source lock must not be a symlink")
    if not path.exists():
        yield
        return
    import fcntl
    with path.open("rb") as handle:
        try:
            fcntl.flock(handle, fcntl.LOCK_SH | fcntl.LOCK_NB)
        except BlockingIOError as exc:
            raise ValueError("source run is still active; review a completed run") from exc
        try:
            yield
        finally:
            fcntl.flock(handle, fcntl.LOCK_UN)


def execute(source: Path, output: Path, *, gain_margin: float = .001,
            max_report_bytes: int = 48 * 1024) -> dict:
    if source.is_symlink() or output.is_symlink():
        raise ValueError("input/output directories must not be symlinks")
    source, output = source.resolve(), output.resolve()
    if not source.is_dir():
        raise ValueError("input directory not found")
    if output == source or source in output.parents or output in source.parents:
        raise ValueError("review output must be separate from the source tree")
    if output.exists():
        raise ValueError("review directory already exists; use a new output directory")
    if not 4096 <= max_report_bytes <= 128 * 1024:
        raise ValueError("report byte limit must be between 4096 and 131072")
    output.parent.mkdir(parents=True, exist_ok=True)
    with source_lock(source):
        print("[verify] checking original artifact hashes and grid completeness", flush=True)
        source_summary, config, rows, source_hashes = validate_source(source)
        print(f"[review] {len(rows)} existing states; no model calls", flush=True)
        report = calculate_review(rows, config, gain_margin)
        # An interrupted review never publishes a misleading complete directory.
        with tempfile.TemporaryDirectory(prefix=".v5m21-building-", dir=output.parent) as staging_name:
            staging = Path(staging_name)
            inventory = []
            for label in ("global-formulations", "paired-conditions", "regime-summary", "model-regimes", "disagreements",
                          "pareto-water", "pareto-water-flow-proxy"):
                inventory.extend(write_table(staging, label, report[label],
                    text=label not in {"global-formulations", "paired-conditions"}, max_bytes=max_report_bytes))
            # Keep every failed state in small, per-pair text files as well as in
            # the full archive. Each shard remains below the configured byte cap.
            identities = sorted({row["pair_id"] for row in rows})
            for i, pair in enumerate(identities, 1):
                detail = [r for r in report["paired-conditions"] if r["pair_id"] == pair]
                compact = [{k: r.get(k) for k in (
                    "pair_id", "name", "cas_number", "additive_mass_fraction", "outlet_temperature_k",
                    "pressure_pa", "comparison_eligible", "failure_counts", "minimum_paired_delta_h_ratio",
                    "maximum_model_spread_delta_h_ratio", "maximum_storage_bubble_pressure_ratio",
                    "water_dominated_thermodynamic", "pareto_with_water",
                    "required_net_heat_reduction_for_mass_parity")} for r in detail]
                inventory.extend(write_table(staging, f"pair-{i:03d}", compact,
                                              text=True, max_bytes=max_report_bytes))
            print("[archive] compressing allowlisted CSV/JSON; originals remain untouched", flush=True)
            archive = []
            for name in SOURCE_FILES:
                if name in source_hashes:
                    item = compress_file(source / name, staging / ("archive-" + name + ".gz"), source_hashes[name])
                    archive.append({**item, "source_file": name})
            # Verify again after processing to detect concurrent noncooperative edits.
            for name, digest in source_hashes.items():
                if sha256(source / name) != digest:
                    raise ValueError("source changed during review: " + name)
            summary = {**report["summary"], "source_summary_reproduced": {
                "predicted_state_count": source_summary["predicted_state_count"],
                "complete_two_model_formulation_count": source_summary["complete_two_model_formulation_count"]},
                "source_files_archived": len(archive),
                "archived_original_bytes": sum(r["original_bytes"] for r in archive),
                "archive_compressed_bytes": sum(r["compressed_bytes"] for r in archive),
                "utf8_report_shard_count": sum(r["kind"] == "utf8_report" for r in inventory),
                "maximum_utf8_report_bytes": max_report_bytes, "original_outputs_modified": False,
                "gzip_auto_decode_by_connector_assumed": False}
            if summary["complete_two_model_formulation_count"] != source_summary["complete_two_model_formulation_count"]:
                raise ValueError("reconstructed comparison count differs from source summary")
            write_json(staging / "summary.json", summary)
            manifest = {"study_version": VERSION, "source_artifact_sha256": source_hashes,
                "review_code_sha256": {p.name: sha256(p) for p in
                    (Path(__file__), HERE / "v5m21_review.py", HERE / "v5m2_core.py")},
                "python": platform.python_version(), "zlib": zlib.ZLIB_RUNTIME_VERSION,
                "gzip": {"mtime": 0, "filename": "", "compresslevel": 6, "lossless": True,
                         "byte_reproducibility_scope": "same implementation and zlib version"},
                "gain_margin_fraction": gain_margin, "margin_is_uncertainty": False,
                "physical_model_recomputed": False, "source_path_recorded": source.name}
            write_json(staging / "manifest.json", manifest)
            # Small summary/manifest stay readable, and also receive compressed copies.
            for name in ("summary.json", "manifest.json"):
                inventory.append({**compress_file(staging/name, staging/(name+".gz")), "kind": "review_metadata_gzip"})
            write_json(staging / "publication-index.json", {"archives": archive, "reports": inventory,
                "entry_point": "summary.json", "all_pairs_in_text_reports": True,
                "connector_note": "Gzip archives may require external decompression. Use bounded UTF-8 report shards for text-only tools."})
            text = (
                "# V5m-2.1: frozen-grid review\n\n"
                "Start with summary.json and publication-index.json. No thermodynamic model was rerun.\n\n"
                "report-regime-summary-* contains each pair's best sampled composition per T/P.\n"
                "report-model-regimes-* keeps models separate; report-disagreements-* shows paired\n"
                "phase fractions at the largest discrepancy. report-pair-* retains every paired\n"
                "composition/condition, including incomplete comparisons. report-pareto-water-*\n"
                "adds pure water explicitly to the existing descriptive Pareto objectives.\n\n"
                "All report tables and all allowlisted original CSV/JSON files have verified gzip\n"
                "copies. Original artifacts and their hashes were not changed. The compressed\n"
                "archives preserve original bytes; these text shards are derived reports.\n\n"
                "A 0.1% default reporting margin is only a sensitivity threshold, not calibrated\n"
                "uncertainty. Per-regime gains remain retrospective model hypotheses; selecting\n"
                "a favorable regime after seeing results is not new validation. No wetting,\n"
                "pore-flow, handling, phase or caloric accuracy has been validated by this step.\n"
            )
            (staging / "README.md").write_text(text, encoding="utf-8")
            write_json(staging / "artifact-hashes.json", {p.name: sha256(p) for p in sorted(staging.iterdir()) if p.is_file()})
            if output.exists():
                raise ValueError("output appeared during review; refusing to replace it")
            staging.rename(output)
    print(json.dumps(summary, indent=2, sort_keys=True), flush=True)
    return summary


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-dir", type=Path, default=Path("predictive-v5m2-results"))
    parser.add_argument("--output-dir", type=Path, default=Path("review-v5m21-results"))
    parser.add_argument("--gain-margin", type=float, default=.001,
                        help="Reporting sensitivity threshold, not uncertainty (default 0.001 = 0.1%%).")
    parser.add_argument("--max-report-bytes", type=int, default=48*1024)
    args = parser.parse_args()
    try:
        execute(args.input_dir, args.output_dir, gain_margin=args.gain_margin, max_report_bytes=args.max_report_bytes)
    except (ValueError, OSError, KeyError, TypeError) as exc:
        print(f"review failed: {exc}", file=sys.stderr)
        return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
