"""Offline review/publication regression tests. All scientific fixtures are synthetic."""
from __future__ import annotations

import copy
import csv
import gzip
import io
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
from contextlib import redirect_stdout

HERE = Path(__file__).resolve().parents[2] / "formulation-screen-v5m"
sys.path.insert(0, str(HERE))
import v5m21_review as review
import review_predictive_grid as runner
from v5m2_core import MODELS, compare_models, digest


def config():
    return {"version": "v5m-2", "models": list(MODELS),
            "grid": {"outlet_temperatures_k": [350., 450.], "pressures_pa": [101325.]}}


def rows():
    return [{"pair_id": "synthetic-pair", "cas_number": "synthetic-cas", "name": "Synthetic pair",
             "smiles": "CCO", "model": model, "additive_mass_fraction": w,
             "inlet_temperature_k": 293.15, "outlet_temperature_k": t, "pressure_pa": 101325.,
             "status": "ok", "model_comparison_eligible": True,
             "same_model_water_delta_h_ratio": 1.1 if t == 350 else .9,
             "storage_bubble_pressure_ratio": 1.05, "viscosity_ratio_proxy": 1.1,
             "outlet": {"vapor_mole_fraction": 0. if t == 350 else 1., "liquid_phase_count": 1 if t == 350 else 0}}
            for model in MODELS for w in (.01, .03) for t in (350., 450.)]


def write_fixture(path):
    path.mkdir()
    rs = rows(); cfg = config()
    fields = sorted({k for r in rs for k in r})
    (path / "prediction-grid.csv").write_text(review.csv_text(rs, fields), encoding="utf-8", newline="")
    catalog = {"candidates": [{k: rs[0][k] for k in ("pair_id", "cas_number", "name", "smiles")}]}
    contract = {"config": cfg, "packages": {"fixture": "synthetic"}}
    objects = {
        "catalog.json": catalog,
        "grid-plan.json": {"models": list(MODELS), "composition_basis": "mass", "pairs": 1,
                           "mass_fractions": [.01, .03], "conditions_per_composition_model": 2},
        "refinement-plan.json": {"tasks": []},
        "manifest.json": {"version": "v5m-2", "contract": contract, "signature": digest(contract)},
        "summary.json": {"study_version": "v5m-2", "study_complete": True, "predicted_state_count": 8,
                         "state_status_counts": {"ok": 8}, "complete_two_model_formulation_count": 2},
        "predictions.json": {"test_fixture_only": True, "rows": rs},
    }
    for name, obj in objects.items():
        review.write_json(path / name, obj)
    rehash(path)
    return rs


def rehash(path):
    review.write_json(path / "artifact-hashes.json", {p.name: review.sha256(p) for p in path.iterdir()
                                                     if p.is_file() and p.name != "artifact-hashes.json"})


class CompressionTests(unittest.TestCase):
    def test_lossless_and_unchanged_input(self):
        with tempfile.TemporaryDirectory() as tmp:
            d = Path(tmp); p = d / "input.csv"; data = b'a,b\r\n1,"line\nsecond"\r\n'
            p.write_bytes(data)
            report = review.compress_file(p, d / "out.gz", review.sha256(p))
            self.assertEqual(gzip.decompress((d/"out.gz").read_bytes()), data)
            self.assertEqual(p.read_bytes(), data)
            self.assertTrue(report["round_trip_verified"])
    def test_repeatable_gzip_no_timestamp_or_name(self):
        with tempfile.TemporaryDirectory() as tmp:
            d = Path(tmp); p = d/"first.txt"; p.write_text("same"*100)
            q = d/"second.txt"; q.write_bytes(p.read_bytes())
            review.compress_file(p,d/"a.gz");review.compress_file(q,d/"b.gz")
            self.assertEqual((d/"a.gz").read_bytes(),(d/"b.gz").read_bytes())
    def test_hash_mismatch_no_archive(self):
        with tempfile.TemporaryDirectory() as tmp:
            d=Path(tmp);p=d/"a";p.write_text("changed")
            with self.assertRaisesRegex(ValueError,"changed"):
                review.compress_file(p,d/"a.gz","0"*64)
            self.assertFalse((d/"a.gz").exists())
            self.assertFalse((d/"a.gz.tmp").exists())
    def test_existing_destination_and_temp_preserved(self):
        with tempfile.TemporaryDirectory() as tmp:
            d=Path(tmp);p=d/"a";p.write_text("input")
            for name in ("a.gz", "b.gz.tmp"):
                (d/name).write_text("keep")
            for name in ("a.gz", "b.gz"):
                with self.assertRaises(ValueError):review.compress_file(p,d/name)
            self.assertEqual((d/"a.gz").read_text(),"keep")
            self.assertEqual((d/"b.gz.tmp").read_text(),"keep")
    def test_symlink_source_rejected(self):
        with tempfile.TemporaryDirectory() as tmp:
            d=Path(tmp);(d/"a").write_text("private");(d/"link").symlink_to(d/"a")
            with self.assertRaises(ValueError):review.compress_file(d/"link",d/"out.gz")
    def test_shards_preserve_all_rows_and_byte_limit(self):
        with tempfile.TemporaryDirectory() as tmp:
            d=Path(tmp);rs=[{"id":i,"text":"\u03bc"*100+",\nline"} for i in range(100)]
            items=review.write_table(d,"test",rs,text=True,max_bytes=1024)
            shards=[r for r in items if r["kind"]=="utf8_report"]
            self.assertGreater(len(shards),1)
            recovered=[]
            for item in shards:
                p=d/item["file"];self.assertLessEqual(p.stat().st_size,1024)
                with p.open(newline="") as f:recovered+=list(csv.DictReader(f))
            self.assertEqual([int(r["id"]) for r in recovered], list(range(100)))
            self.assertEqual(recovered[0]["text"],rs[0]["text"])
    def test_empty_table_has_header(self):
        with tempfile.TemporaryDirectory() as tmp:
            d=Path(tmp);review.write_table(d,"empty",[],text=True)
            self.assertEqual((d/"report-empty-0001.csv").read_text(),"status\n")
    def test_large_row_rejected(self):
        with tempfile.TemporaryDirectory() as tmp:
            with self.assertRaisesRegex(ValueError,"byte limit"):
                review.write_table(Path(tmp),"test",[{"x":"A"*5000}],text=True,max_bytes=1024)


class ReviewTests(unittest.TestCase):
    def test_add_water_dominates_global_loser(self):
        r=review.calculate_review(rows(),config())
        self.assertEqual(r["summary"]["global_pareto_count_including_water"],1)
        self.assertEqual(r["summary"]["water_dominated_global_thermodynamic_count"],2)
        self.assertEqual(r["pareto-water"][0]["name"],"water")
    def test_regime_gain_is_not_global_gain(self):
        s=review.calculate_review(rows(),config())["summary"]
        self.assertEqual(s["global_paired_gain_any_positive_count"],0)
        self.assertEqual(s["condition_results"][0]["paired_gain_above_margin_count"],2)
        self.assertEqual(s["condition_results"][1]["paired_gain_above_margin_count"],0)
    def test_one_model_failure_not_dropped(self):
        rs=rows();rs[0].update(status="inlet_not_single_liquid",model_comparison_eligible=False)
        result=review.calculate_review(rs,config())
        self.assertEqual(result["summary"]["complete_two_model_formulation_count"],1)
        self.assertEqual(result["summary"]["condition_results"][0]["eligible_formulation_count"],1)
    def test_pressure_tradeoff_not_called_gain(self):
        rs=rows()
        for row in rs:row.update(same_model_water_delta_h_ratio=.9,storage_bubble_pressure_ratio=.9)
        result=review.calculate_review(rs,config())
        self.assertEqual(result["summary"]["global_pareto_additive_count"],2)
        self.assertEqual(result["summary"]["global_paired_gain_any_positive_count"],0)
    def test_break_even_not_achieved(self):
        r=review.calculate_review(rows(),config())["global-formulations"][0]
        self.assertAlmostEqual(r["required_net_heat_reduction_for_mass_parity"],.1)
        self.assertAlmostEqual(r["mass_ratio_equal_net_heat"],1/.9)
        self.assertFalse(r["rankable"])
    def test_reporting_margin_does_not_hide_raw_positive(self):
        rs=rows()
        for r in rs:r["same_model_water_delta_h_ratio"]=1.00001
        s=review.calculate_review(rs,config())["summary"]
        self.assertEqual(s["global_paired_gain_any_positive_count"],2)
        self.assertEqual(s["global_paired_gain_above_margin_count"],0)
    def test_missing_viscosity_excluded_only_from_flow_front(self):
        rs=rows()
        for r in rs:r.update(storage_bubble_pressure_ratio=.9,viscosity_ratio_proxy=None)
        result=review.calculate_review(rs,config())
        self.assertGreater(len(result["pareto-water"]),1)
        self.assertEqual(len(result["pareto-water-flow-proxy"]),1)
    def test_disagreement_uses_same_composition(self):
        rs=rows();rs[0]["same_model_water_delta_h_ratio"]=1.3
        result=review.calculate_review(rs,config())
        r=result["disagreements"][0]
        self.assertEqual(r["additive_mass_fraction"],.01)
        self.assertAlmostEqual(r["spread_delta_h_ratio"],.2)
        self.assertAlmostEqual(r["chemsep_nrtl_delta_h_ratio"],1.3)
        self.assertAlmostEqual(r["unifac_dortmund_delta_h_ratio"],1.1)
    def test_nan_and_boolean_margin_rejected(self):
        for margin in (float("nan"),float("inf"),True,-1,1):
            with self.assertRaises(ValueError):review.calculate_review(rows(),config(),margin)
    def test_json_duplicate_or_nonfinite_rejected(self):
        for text in ('{"a":1,"a":2}','{"a":NaN}'):
            with self.assertRaises(ValueError):review.json_loads(text)


class EndToEndTests(unittest.TestCase):
    def setUp(self):
        self.tmp=tempfile.TemporaryDirectory();self.addCleanup(self.tmp.cleanup)
        self.root=Path(self.tmp.name);self.source=self.root/"source";write_fixture(self.source)
    def test_plan_and_hashes_validate(self):
        _,_,rs,hashes=review.validate_source(self.source)
        self.assertEqual(len(rs),8);self.assertIn("predictions.json",hashes)
    def test_changed_file_rejected(self):
        p=self.source/"prediction-grid.csv";p.write_text(p.read_text()+"garbage")
        with self.assertRaisesRegex(ValueError,"hash mismatch"):review.validate_source(self.source)
    def test_missing_state_not_cherry_picked(self):
        rs=rows()[:-1]
        (self.source/"prediction-grid.csv").write_text(review.csv_text(rs,sorted(rs[0])))
        rehash(self.source)
        with self.assertRaisesRegex(ValueError,"grid plan mismatch"):review.validate_source(self.source)
    def test_duplicate_state_rejected(self):
        rs=rows()+[rows()[0]]
        (self.source/"prediction-grid.csv").write_text(review.csv_text(rs,sorted(rs[0])))
        rehash(self.source)
        with self.assertRaisesRegex(ValueError,"duplicate prediction"):review.validate_source(self.source)
    def test_incorrect_summary_count_rejected(self):
        s=review.load_json(self.source/"summary.json");s["predicted_state_count"]=10
        review.write_json(self.source/"summary.json",s);rehash(self.source)
        with self.assertRaisesRegex(ValueError,"summary/grid"):review.validate_source(self.source)
    def test_smoke_run_rejected(self):
        s=review.load_json(self.source/"summary.json");s["study_complete"]=False
        review.write_json(self.source/"summary.json",s);rehash(self.source)
        with self.assertRaisesRegex(ValueError,"non-smoke"):review.validate_source(self.source)
    def test_end_to_end_no_source_changes_no_secret_exports(self):
        (self.source/".env").write_text("secret")
        (self.source/"workers").mkdir();(self.source/"workers"/"x.json").write_text("secret")
        old={p.name:p.read_bytes() for p in self.source.iterdir() if p.is_file()}
        output=self.root/"review"
        with redirect_stdout(io.StringIO()):s=runner.execute(self.source,output,max_report_bytes=4096)
        self.assertTrue(s["study_complete"]);self.assertEqual(s["new_thermodynamic_states_computed"],0)
        self.assertEqual(old,{p.name:p.read_bytes() for p in self.source.iterdir() if p.is_file()})
        self.assertFalse(any("env" in p.name or "worker" in p.name for p in output.iterdir()))
        self.assertEqual(gzip.decompress((output/"archive-prediction-grid.csv.gz").read_bytes()),old["prediction-grid.csv"])
        idx=review.load_json(output/"publication-index.json")
        for item in idx["reports"]:
            if item["kind"]=="utf8_report":self.assertLessEqual((output/item["file"]).stat().st_size,4096)
        hashes=review.load_json(output/"artifact-hashes.json")
        self.assertTrue(all(review.sha256(output/n)==h for n,h in hashes.items()))
    def test_output_reuse_rejected(self):
        output=self.root/"review";output.mkdir()
        with self.assertRaisesRegex(ValueError,"already exists"):runner.execute(self.source,output)
    def test_no_output_nested_in_source(self):
        with self.assertRaisesRegex(ValueError,"separate"):runner.execute(self.source,self.source/"review")
    def test_failure_does_not_publish_partial_directory(self):
        output=self.root/"review"
        with redirect_stdout(io.StringIO()),self.assertRaises(ValueError):
            runner.execute(self.source,output,gain_margin=float("nan"))
        self.assertFalse(output.exists())
    def test_source_lock_refuses_active_writer(self):
        import fcntl
        p=self.source/".run.lock"
        with p.open("w") as f:
            fcntl.flock(f,fcntl.LOCK_EX|fcntl.LOCK_NB)
            with self.assertRaisesRegex(ValueError,"still active"):
                with runner.source_lock(self.source):pass
    def test_cli_help(self):
        result = subprocess.run([sys.executable, str(HERE/"review_predictive_grid.py"), "--help"],
                                capture_output=True, text=True)
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertIn("0.1%", result.stdout)
    def test_synthetic_cli_full_run(self):
        result=subprocess.run([sys.executable,str(HERE/"review_predictive_grid.py"),"--input-dir",str(self.source),
                               "--output-dir",str(self.root/"out")],capture_output=True,text=True)
        self.assertEqual(result.returncode,0,result.stderr)
        self.assertIn('"new_thermodynamic_states_computed": 0',result.stdout)


class IgnorePolicyTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory(); self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name)
        subprocess.run(["git", "init", "-q", str(self.root)], check=True)
        (self.root / ".gitignore").write_bytes((HERE / ".gitignore").read_bytes())
    def ignored(self, path):
        p = self.root / path; p.parent.mkdir(parents=True, exist_ok=True); p.write_text("fixture")
        return subprocess.run(["git", "check-ignore", "-q", "--", path], cwd=self.root).returncode == 0
    def test_summary_is_visible(self):
        self.assertFalse(self.ignored("review-v5m21-results/summary.json"))
    def test_exact_archive_is_visible(self):
        self.assertFalse(self.ignored("review-v5m21-results/archive-predictions.json.gz"))
    def test_report_shards_are_visible(self):
        self.assertFalse(self.ignored("review-v5m21-results/report-pair-001-0001.csv"))
    def test_unknown_archive_is_ignored(self):
        self.assertTrue(self.ignored("review-v5m21-results/archive-credentials.json.gz"))
    def test_nested_report_is_ignored(self):
        self.assertTrue(self.ignored("review-v5m21-results/workers/report-secret.csv"))
    def test_temporary_builds_are_ignored(self):
        self.assertTrue(self.ignored(".v5m21-building-anything/private.json"))
    def test_old_publication_unchanged(self):
        self.assertFalse(self.ignored("predictive-v5m2-results/predictions.json"))
        self.assertTrue(self.ignored("predictive-v5m2-results/workers/data.json"))


if __name__ == "__main__":
    unittest.main()
