"""CSV compression, index migration and recovery tests; no numerical stack."""
from __future__ import annotations

import csv
import gzip
import hashlib
import io
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
from unittest.mock import patch

HERE = Path(__file__).resolve().parents[2] / "mixture-campaign"
sys.path.insert(0, str(HERE))
import repack_csv_reports as mod
from campaign_store import encoded, read_json, write_json


class CsvPublicationTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name)

    def fixture(self):
        run = self.root / "run-fixture"
        run.mkdir()
        raw = b'id,text\r\n1,"two\nlines"\r\n2,"a,b"\r\n'
        report = "report-example-0001.csv"
        (run / report).write_bytes(raw)
        data, entry = mod.gzip_bytes("table-example.csv.gz", raw)
        (run / entry["file"]).write_bytes(data)
        write_json(run / "summary.json", {"scientific_result_sha256": "unchanged"})
        write_json(run / "publication-index.json", {"files": [
            {**mod.artifact(report, raw), "rows": 2}, entry]})
        write_json(self.root / "latest.json", {"run_id": run.name, "scientific_result_sha256": "unchanged"})
        return run, report, raw

    def test_new_writer_has_no_plain_csv(self):
        index = mod.compressed_table(self.root, "sample", [{"n": i, "s": "a" * 12} for i in range(30)], 100)
        self.assertFalse(list(self.root.glob("*.csv")))
        shards = [r for r in index if r["file"].startswith("report-")]
        self.assertGreater(len(shards), 1)
        self.assertEqual(sum(r["rows"] for r in shards), 30)
        for r in shards:
            raw = gzip.decompress((self.root / r["file"]).read_bytes())
            self.assertLessEqual(len(raw), 100)
            self.assertEqual(mod.sha256(raw), r["uncompressed_sha256"])

    def test_unicode_commas_quotes_newlines_and_nested_values(self):
        rows = [{"id": 1, "s": 'caf\u00e9, "quoted"\nnext', "j": [1, 2]}, {"id": 2, "s": "ok", "j": None}]
        mod.compressed_table(self.root, "unicode", rows, 300)
        with gzip.open(self.root / "table-unicode.csv.gz", "rt", encoding="utf-8", newline="") as h:
            loaded = list(csv.DictReader(h))
        self.assertEqual(loaded[0]["s"], rows[0]["s"])
        self.assertEqual(loaded[0]["j"], "[1,2]")

    def test_deterministic_bytes(self):
        a = self.root / "a"; b = self.root / "b"
        first = mod.compressed_table(a, "same", [{"x": 1}], 100)
        second = mod.compressed_table(b, "same", [{"x": 1}], 100)
        self.assertEqual(first, second)
        self.assertTrue(all((a/r["file"]).read_bytes() == (b/r["file"]).read_bytes() for r in first))

    def test_empty_table(self):
        index = mod.compressed_table(self.root, "empty", [], 100)
        self.assertEqual(gzip.decompress((self.root / "table-empty.csv.gz").read_bytes()), b"status\n")
        self.assertEqual(next(x["rows"] for x in index if x["file"].startswith("report-")), 0)

    def test_oversized_row_rejected(self):
        with self.assertRaises(ValueError):
            mod.compressed_table(self.root, "bad", [{"x": "a"*300}], 100)

    def test_unsafe_name_rejected(self):
        with self.assertRaises(ValueError):
            mod.compressed_table(self.root, "../escape", [], 100)

    def test_preview_is_bounded_and_explicitly_incomplete(self):
        raw = b"text\n" + (b"x"*10000+b"\n")*10
        data = mod.preview("bounded", raw)
        self.assertLessEqual(len(data), mod.PREVIEW_BYTES)
        self.assertIn(b"first 8", data)
        self.assertIn(b"preview row truncated", data)

    def test_migration_roundtrip_and_scientific_data_unchanged(self):
        run, name, raw = self.fixture()
        summary = (run / "summary.json").read_bytes()
        pointer = (self.root / "latest.json").read_bytes()
        old_index = (run / "publication-index.json").read_bytes()
        result = mod.repack_snapshot(run)
        self.assertEqual(result["converted_csv_count"], 1)
        self.assertFalse((run / name).exists())
        self.assertEqual(gzip.decompress((run/(name+".gz")).read_bytes()), raw)
        self.assertEqual((run/"summary.json").read_bytes(), summary)
        self.assertEqual((self.root/"latest.json").read_bytes(), pointer)
        self.assertEqual(gzip.decompress((run/"publication-index.before-csv-repack.json.gz").read_bytes()), old_index)
        index = read_json(run/"publication-index.json")
        self.assertEqual(index["csv_policy"], mod.POLICY)
        self.assertFalse(any(r["file"].endswith(".csv") for r in index["files"]))
        for record in index["files"]:
            self.assertEqual(mod.sha256((run/record["file"]).read_bytes()), record["sha256"])

    def test_dry_run_changes_no_snapshot_bytes(self):
        run, _, _ = self.fixture()
        before = {p.name:p.read_bytes() for p in run.iterdir()}
        mod.repack_snapshot(run, dry_run=True)
        self.assertEqual(before, {p.name:p.read_bytes() for p in run.iterdir()})

    def test_idempotent(self):
        run, _, _ = self.fixture()
        mod.repack_snapshot(run)
        before = {p.name:p.read_bytes() for p in run.iterdir()}
        result = mod.repack_snapshot(run)
        self.assertTrue(result["already_compressed"])
        self.assertEqual(before, {p.name:p.read_bytes() for p in run.iterdir()})

    def test_corrupt_report_rejected_before_writes(self):
        run, name, _ = self.fixture()
        (run/name).write_bytes(b"changed")
        before = {p.name:p.read_bytes() for p in run.iterdir()}
        with self.assertRaises(ValueError): mod.repack_snapshot(run)
        self.assertEqual(before, {p.name:p.read_bytes() for p in run.iterdir()})

    def test_corrupt_full_table_rejected(self):
        run, _, _ = self.fixture()
        (run/"table-example.csv.gz").write_bytes(b"changed")
        with self.assertRaises(ValueError): mod.repack_snapshot(run)
        self.assertFalse((run/"report-example-0001.csv.gz").exists())

    def test_different_existing_destination_preserved(self):
        run, name, _ = self.fixture()
        (run/(name+".gz")).write_bytes(b"unrelated")
        with self.assertRaises(ValueError): mod.repack_snapshot(run)
        self.assertTrue((run/name).exists())
        self.assertEqual((run/(name+".gz")).read_bytes(), b"unrelated")

    def test_identical_partial_destination_can_resume(self):
        run, name, raw = self.fixture()
        (run/(name+".gz")).write_bytes(mod.gzip_bytes(name+".gz", raw)[0])
        mod.repack_snapshot(run)
        self.assertFalse((run/name).exists())

    def test_interrupted_cleanup_can_resume(self):
        run, name, raw = self.fixture()
        with patch.object(mod, "cleanup_originals", side_effect=RuntimeError("interruption")):
            with self.assertRaises(RuntimeError): mod.repack_snapshot(run)
        self.assertTrue((run/name).exists())
        result = mod.repack_snapshot(run)
        self.assertEqual(result["removed_plain_csv_count"], 1)
        self.assertEqual(gzip.decompress((run/(name+".gz")).read_bytes()), raw)
        self.assertFalse((run/name).exists())

    def test_changed_original_after_interruption_is_not_deleted(self):
        run, name, _ = self.fixture()
        with patch.object(mod, "cleanup_originals", side_effect=RuntimeError("interruption")):
            with self.assertRaises(RuntimeError): mod.repack_snapshot(run)
        (run/name).write_bytes(b"new data")
        with self.assertRaises(ValueError): mod.repack_snapshot(run)
        self.assertEqual((run/name).read_bytes(), b"new data")

    def test_symlink_is_rejected(self):
        run, name, raw = self.fixture()
        (run/name).unlink()
        target = self.root/"target.csv"; target.write_bytes(raw)
        (run/name).symlink_to(target)
        with self.assertRaises(ValueError): mod.repack_snapshot(run)
        self.assertEqual(target.read_bytes(), raw)

    def test_unknown_csv_is_not_deleted(self):
        run, _, _ = self.fixture()
        (run/"user-data.csv").write_text("secret\n")
        with self.assertRaises(ValueError): mod.repack_snapshot(run)
        self.assertTrue((run/"user-data.csv").exists())

    def test_duplicate_index_rejected(self):
        run, _, _ = self.fixture()
        index = read_json(run/"publication-index.json")
        index["files"].append(index["files"][0])
        write_json(run/"publication-index.json", index)
        with self.assertRaises(ValueError): mod.repack_snapshot(run)

    def test_path_escape_rejected(self):
        run, _, _ = self.fixture()
        index = read_json(run/"publication-index.json")
        index["files"][0]["file"] = "../user.csv"
        write_json(run/"publication-index.json", index)
        with self.assertRaises(ValueError): mod.repack_snapshot(run)

    def test_cli_uses_latest_only(self):
        run, name, _ = self.fixture()
        old = self.root/"run-old"; old.mkdir(); (old/"other.csv").write_text("untouched")
        completed = subprocess.run([sys.executable, str(HERE/"repack_csv_reports.py"),
                                   "--results-dir", str(self.root)], capture_output=True, text=True)
        self.assertEqual(completed.returncode, 0, completed.stderr)
        self.assertEqual(json.loads(completed.stdout)["converted_csv_count"], 1)
        self.assertFalse((run/name).exists())
        self.assertEqual((old/"other.csv").read_text(), "untouched")

    def test_gitignore_accepts_gzip_preview_and_tracks_existing_deletions(self):
        subprocess.run(["git", "init", "-q", str(self.root)], check=True)
        (self.root/".gitignore").write_bytes((HERE/".gitignore").read_bytes())
        for prefix in ["mission-v5m4-results/run-test/", "campaign-v5m3-results/run-test/"]:
            allowed = ["report-x-0001.csv.gz", "table-x.csv.gz", "preview-x.md", "publication-index.before-csv-repack.json.gz"]
            blocked = ["report-x-0001.csv", "temp.csv", "workers/x.csv.gz"]
            for name in allowed + blocked:
                proc = subprocess.run(["git", "check-ignore", "--no-index", prefix+name], cwd=self.root, capture_output=True)
                self.assertEqual(proc.returncode, 0 if name in blocked else 1, prefix+name)
        run, name, _ = self.fixture()
        # Track a legacy snapshot before migrating; git add -A must stage deletion.
        subprocess.run(["git", "add", "--", run.name, "latest.json"], cwd=self.root, check=True)
        subprocess.run(["git", "-c", "user.name=Test", "-c", "user.email=test@example.invalid", "commit", "-qm", "fixture"], cwd=self.root, check=True)
        mod.repack_snapshot(run)
        subprocess.run(["git", "add", "-A", "--", run.name], cwd=self.root, check=True)
        status = subprocess.check_output(["git", "diff", "--cached", "--name-status"], cwd=self.root, text=True)
        self.assertIn("D\t"+run.name+"/"+name, status)
        self.assertIn("A\t"+run.name+"/"+name+".gz", status)


if __name__ == "__main__":
    unittest.main()
