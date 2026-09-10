#!/usr/bin/env python3
"""Gzip every CSV report, including shards, without rerunning any experiment.

Small Markdown previews are explicitly incomplete; complete tables and report
shards are lossless gzip files. JSON summaries and publication indexes stay
readable. Repacking changes packaging only, never scientific JSON or task cache.
"""
from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import io
import json
from pathlib import Path
import re
from typing import Any

from campaign_store import atomic_bytes, encoded, exclusive_lock, read_json, write_json

POLICY = "gzip-csv-with-bounded-preview-v1"
PREVIEW_ROWS = 8
PREVIEW_BYTES = 12288


def sha256(raw: bytes) -> str:
    return hashlib.sha256(raw).hexdigest()


def artifact(name: str, raw: bytes) -> dict:
    return {"file": name, "bytes": len(raw), "sha256": sha256(raw)}


def gzip_bytes(name: str, raw: bytes) -> tuple[bytes, dict]:
    data = gzip.compress(raw, compresslevel=9, mtime=0)
    if gzip.decompress(data) != raw:
        raise ValueError("gzip round-trip failure: " + name)
    return data, {**artifact(name, data), "uncompressed_bytes": len(raw),
                  "uncompressed_sha256": sha256(raw)}


def csv_line(values) -> bytes:
    stream = io.StringIO(newline="")
    csv.writer(stream, lineterminator="\n").writerow(values)
    return stream.getvalue().encode("utf-8")


def preview(name: str, raw: bytes) -> bytes:
    """Bounded first-row preview, never a replacement for the full table."""
    reader = csv.DictReader(io.StringIO(raw.decode("utf-8"), newline=""))
    count, examples = 0, []
    for row in reader:
        count += 1
        if len(examples) < PREVIEW_ROWS:
            text = json.dumps(row, ensure_ascii=False, sort_keys=True)
            # Indented text prevents embedded markup from becoming instructions.
            text = text.replace("\r", " ").replace("\n", " ")
            if len(text.encode("utf-8")) > 1200:
                text = text.encode("utf-8")[:1150].decode("utf-8", errors="ignore") + " ... [preview row truncated]"
            examples.append("    " + text + "\n")
    head = (f"# {name}: limited preview\n\n"
            f"Full table: `table-{name}.csv.gz`\n\n"
            f"Total records: {count}. Shown: first {len(examples)} in input order. "
            "This is not a selected ranking. Preview rows may be truncated; "
            "the gzip table retains every original byte.\n\n")
    result = (head + "\n".join(examples)).encode("utf-8")
    if len(result) > PREVIEW_BYTES:
        raise ValueError("preview exceeds its byte budget")
    return result


def compressed_table(output: Path, name: str, rows: list[dict], limit: int) -> list[dict]:
    """Same CSV serialization and shard boundaries as the original table writer."""
    if not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9_.-]{0,180}", name):
        raise ValueError("unsafe table name")
    if isinstance(limit, bool) or not isinstance(limit, int) or limit < 1:
        raise ValueError("uncompressed shard limit must be a positive integer")
    fields = sorted({k for row in rows for k in row}) or ["status"]
    header = csv_line(fields)
    if len(header) > limit:
        raise ValueError("CSV header exceeds uncompressed shard byte limit")
    full, block = bytearray(header), bytearray(header)
    count, shard, records = 0, 1, []

    def flush():
        nonlocal block, count, shard
        name_out = f"report-{name}-{shard:04d}.csv.gz"
        data, record = gzip_bytes(name_out, bytes(block))
        atomic_bytes(output / name_out, data)
        records.append({**record, "rows": count})
        block, count, shard = bytearray(header), 0, shard + 1

    for row in rows:
        line = csv_line([json.dumps(row.get(k), sort_keys=True, allow_nan=False, separators=(",", ":"))
                         if isinstance(row.get(k), (dict, list)) else row.get(k) for k in fields])
        if len(header) + len(line) > limit:
            raise ValueError("single CSV row exceeds report byte limit: " + name)
        if len(block) + len(line) > limit:
            flush()
        block.extend(line)
        full.extend(line)
        count += 1
    if count or not rows:
        flush()
    data, record = gzip_bytes("table-" + name + ".csv.gz", bytes(full))
    atomic_bytes(output / record["file"], data)
    records.append(record)
    text = preview(name, bytes(full))
    name_out = "preview-" + name + ".md"
    atomic_bytes(output / name_out, text)
    records.append({**artifact(name_out, text), "preview_only": True})
    return records


def safe_file(root: Path, name: str) -> Path:
    if (not isinstance(name, str) or name in {"", ".", ".."} or
            Path(name).name != name or "/" in name or "\\" in name):
        raise ValueError("unsafe publication path: " + repr(name))
    path = root / name
    if path.is_symlink():
        raise ValueError("symlink in publication: " + name)
    return path


def verified(root: Path, record: dict) -> bytes:
    raw = safe_file(root, record["file"]).read_bytes()
    if len(raw) != record["bytes"] or sha256(raw) != record["sha256"]:
        raise ValueError("published file integrity mismatch: " + record["file"])
    return raw


def cleanup_originals(root: Path, index: dict, *, dry_run: bool) -> int:
    """Resume safely after interruption between committing the index and cleanup."""
    entries = {r["file"]: r for r in index["files"]}
    pending = []
    for old in index.get("csv_repack", {}).get("original_csv_records", []):
        path = safe_file(root, old["file"])
        if not path.exists():
            continue
        raw = verified(root, old)
        current = entries.get(old["file"] + ".gz")
        if current is None or gzip.decompress(verified(root, current)) != raw:
            raise ValueError("cannot remove CSV without verified gzip: " + old["file"])
        pending.append(path)
    if not dry_run:
        for path in pending:
            path.unlink()
    return len(pending)


def repack_snapshot(root: Path, *, dry_run: bool = False) -> dict:
    if root.is_symlink() or not root.is_dir():
        raise ValueError("snapshot must be an existing non-symlink directory")
    index_path = safe_file(root, "publication-index.json")
    source_bytes = index_path.read_bytes()
    index = json.loads(source_bytes)
    files = index["files"]
    names = [r["file"] for r in files]
    if len(names) != len(set(names)):
        raise ValueError("duplicate publication-index entry")
    for name in names:
        safe_file(root, name)
    plain = [r for r in files if r["file"].endswith(".csv")]
    allowed_old = {r["file"] for r in plain} | {
        r["file"] for r in index.get("csv_repack", {}).get("original_csv_records", [])}
    unknown = {p.name for p in root.glob("*.csv")} - allowed_old
    if unknown:
        raise ValueError("unindexed CSV files; refusing removal: " + ",".join(sorted(unknown)))
    if not plain:
        removed = cleanup_originals(root, index, dry_run=dry_run)
        return {"snapshot": root.name, "converted_csv_count": 0,
                "removed_plain_csv_count": removed, "already_compressed": True, "dry_run": dry_run}

    # Preflight every conversion and destination before writing anything.
    writes, replacements = {}, {}
    for old in plain:
        raw = verified(root, old)
        name = old["file"] + ".gz"
        data, record = gzip_bytes(name, raw)
        replacements[old["file"]] = {**old, **record}
        writes[name] = data
    previews = []
    for record in files:
        match = re.fullmatch(r"table-([A-Za-z0-9][A-Za-z0-9_.-]{0,180})\.csv\.gz", record["file"])
        if not match:
            continue
        raw = gzip.decompress(verified(root, record))
        if (record.get("uncompressed_sha256", sha256(raw)) != sha256(raw) or
                record.get("uncompressed_bytes", len(raw)) != len(raw)):
            raise ValueError("uncompressed table integrity mismatch: " + record["file"])
        name = "preview-" + match[1] + ".md"
        data = preview(match[1], raw)
        writes[name] = data
        previews.append({**artifact(name, data), "preview_only": True})
    backup_name = "publication-index.before-csv-repack.json.gz"
    backup, backup_record = gzip_bytes(backup_name, source_bytes)
    writes[backup_name] = backup
    additional = previews + [backup_record]
    if (set(replacements[r["file"]]["file"] for r in plain) | {r["file"] for r in additional}) & set(names):
        raise ValueError("destination already indexed; refusing conflicting migration")
    for name, data in writes.items():
        path = safe_file(root, name)
        if path.exists() and path.read_bytes() != data:
            raise ValueError("destination exists with different contents: " + name)

    updated = {**index, "files": [replacements.get(r["file"], r) for r in files] + additional,
               "csv_policy": POLICY, "gzip_connector_decode_assumed": False,
               "interpretation": "All CSV data are gzip; Markdown previews are bounded and incomplete.",
               "csv_repack": {"schema": "csv-packaging-migration-v1",
                              "source_publication_index_sha256": sha256(source_bytes),
                              "original_csv_records": plain,
                              "scientific_data_changed": False}}
    if not dry_run:
        for name, data in writes.items():
            atomic_bytes(root / name, data)
        # Commit new paths first. Any leftover original is safe to clean on retry.
        write_json(index_path, updated)
        cleanup_originals(root, updated, dry_run=False)
    return {"snapshot": root.name, "converted_csv_count": len(plain),
            "removed_plain_csv_count": len(plain), "already_compressed": False,
            "csv_input_bytes": sum(r["bytes"] for r in plain),
            "gzip_output_bytes": sum(r["bytes"] for r in replacements.values()),
            "scientific_data_changed": False, "dry_run": dry_run}


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--results-dir", type=Path,
                        default=Path(__file__).resolve().parent / "mission-v5m4-results")
    parser.add_argument("--all-runs", action="store_true", help="Repack all published run-* snapshots, not just latest.")
    parser.add_argument("--check", action="store_true", help="Verify and report; do not modify snapshots.")
    args = parser.parse_args(argv)
    try:
        root = args.results_dir.expanduser().absolute()
        if root.is_symlink() or not root.is_dir():
            raise ValueError("results directory is missing or is a symlink")
        with exclusive_lock(root):
            if args.all_runs:
                runs = sorted(p for p in root.glob("run-*") if p.is_dir())
            else:
                pointer = read_json(safe_file(root, "latest.json"))
                name = pointer["run_id"]
                if not isinstance(name, str) or not name.startswith("run-"):
                    raise ValueError("invalid latest snapshot name")
                runs = [safe_file(root, name)]
            if not runs:
                raise ValueError("no published snapshots")
            for run in runs:
                print(json.dumps(repack_snapshot(run, dry_run=args.check), sort_keys=True), flush=True)
        return 0
    except (OSError, ValueError, KeyError, TypeError, RuntimeError) as exc:
        parser.error(str(exc))


if __name__ == "__main__":
    raise SystemExit(main())
