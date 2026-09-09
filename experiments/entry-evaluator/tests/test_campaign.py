"""Synthetic architecture tests; never counted as physical or thermodynamic evidence."""
from __future__ import annotations
from copy import deepcopy
from contextlib import redirect_stdout
import gzip
import importlib.util
import io
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
from unittest.mock import patch

import yaml
HERE = Path(__file__).resolve().parents[2] / "mixture-campaign"
sys.path.insert(0, str(HERE))
import campaign_backend as b
import campaign_config as cc
import campaign_design as d
import campaign_publication as pub
import campaign_studies as s
import campaign_workflow as w
import run_campaign as cli
from campaign_store import Store, read_json, write_json, digest, encoded, exclusive_lock

CONFIG = yaml.safe_load((HERE / "campaign.yaml").read_text())


def compact():
    c = deepcopy(CONFIG)
    c["coarse"].update(log_points=2, linear_points=2, temperatures_k=[350., 450.], pressures_pa=[101325., 200000.])
    c["adaptive"].update(max_rounds=2, max_new_points_per_pair=4, points_per_round=2)
    c["execution"].update(workers=1, batch_points=9)
    return c


def pair(n="a", bias=0.):
    return {"pair_id": "synthetic-"+n, "cas_number": n, "name": "Synthetic fixture "+n,
            "smiles": "CCO", "inchi_key": "SYNTHETIC-TEST-ONLY",
            "models": {m: {"bias": bias} for m in b.MODELS}}


class FakeDriver:
    def __init__(self, *args):
        self.environment = {"synthetic_fixture": True}
        self.code = {"synthetic_test_backend": "1"}
        self.calls = []
        self.pairs = [pair()]
        self.fail = False
    def catalog(self, config, historical):
        self.calls.append("catalog")
        return {"candidates": deepcopy(self.pairs), "rejections": []}
    def freeze(self, p, model, config):
        self.calls.append("freeze")
        return {"pair": b.model_pair(p, model), "model": model, "config": b.numerical_config(config),
                "selected_correlations": {"synthetic": "fixed-envelope"}}
    def states(self, frozen, points):
        self.calls.append(("states", len(points)))
        if self.fail: raise RuntimeError("synthetic numerical failure")
        p, model = frozen["pair"], frozen["model"]
        rows = []
        for x, t, pressure in points:
            vapor = .1 if t < 380 and pressure < 150000 else 0.
            ratio = 1-.5*x + (3*x if vapor else 0.) + p["models"][model]["bias"]
            if model == b.MODELS[1]: ratio -= x*.1
            rows.append({"pair_id": p["pair_id"], "cas_number": p["cas_number"], "name": p["name"], "smiles": p["smiles"],
                "model": model, "additive_mass_fraction": x, "outlet_temperature_k": t, "pressure_pa": pressure,
                "inlet_temperature_k": frozen["config"]["grid"]["inlet_temperature_k"], "status": "ok", "model_comparison_eligible": True,
                "same_model_water_delta_h_ratio": ratio, "delta_h_j_kg": ratio*100000.,
                "storage_bubble_pressure_ratio": 1+x, "viscosity_ratio_proxy": 1+x,
                "outlet": {"vapor_mole_fraction": vapor, "liquid_phase_count": 1},
                "evidence_kind": "synthetic_test_only", "endpoint_relative_error": 0.})
        return rows


class StoreTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory(); self.root = Path(self.temp.name)
    def tearDown(self): self.temp.cleanup()
    def store(self, **kw): return Store(self.root/"cache", self.root/"journal.jsonl", **kw)
    def test_cache_hit_does_not_execute(self):
        a = self.store().run("a", {}, {}, lambda: 1)
        b1 = self.store().run("a", {}, {}, lambda: self.fail("executed on hit"))
        self.assertEqual(a.result_hash, b1.result_hash)
    def test_code_change_invalidates(self):
        a = self.store().run("a", {}, {"source": "1"}, lambda: 1)
        b1 = self.store().run("a", {}, {"source": "2"}, lambda: 2)
        self.assertNotEqual(a.key, b1.key)
    def test_dependency_change_invalidates(self):
        st = self.store(); a = st.run("parent", 1, {}, lambda: 1)
        first = st.run("child", {}, {}, lambda: 2, dependencies={"a": a})
        a2 = st.run("parent", 2, {}, lambda: 2)
        second = st.run("child", {}, {}, lambda: 3, dependencies={"a": a2})
        self.assertNotEqual(first.key, second.key)
    def test_failure_retry_retains_immutable_old_outcome(self):
        def fail(): raise ValueError("test")
        old = self.store().run("task", {}, {}, fail)
        self.assertEqual(old.outcome["status"], "failed")
        self.assertEqual(self.store().run("task", {}, {}, lambda: 3).outcome["status"], "failed")
        new = self.store(retry_failures=True).run("task", {}, {}, lambda: 3)
        self.assertEqual(new.key, old.key); self.assertNotEqual(new.result_hash, old.result_hash)
        self.assertTrue((self.root/"cache"/"objects"/(old.result_hash+".json.gz")).exists())
    def test_retry_changes_downstream_even_when_parent_key_same(self):
        st = self.store(); spec = st.spec("parent", {}, {})
        a = st.put(spec, error="failure")
        child1 = st.run("child", {}, {}, lambda: None, dependencies={"a": a})
        b1 = self.store(retry_failures=True).run("parent", {}, {}, lambda: 3)
        child2 = self.store().run("child", {}, {}, lambda: 3, dependencies={"a": b1})
        self.assertNotEqual(child1.key, child2.key)
    def test_canonical_tuple_list_roundtrip(self):
        a = self.store().run("a", {"points": [(1, 2)]}, {}, lambda: {"values": (1, 2)})
        b1 = self.store().run("a", {"points": [[1, 2]]}, {}, lambda: None)
        self.assertEqual(a.data, b1.data)
    def test_input_mutation_cannot_change_saved_spec(self):
        inp = {"x": 1}
        def fn(): inp["x"] = 2; return 1
        a = self.store().run("a", inp, {}, fn)
        self.assertEqual(a.spec["inputs"]["x"], 1)
    def test_nan_is_failed_outcome(self):
        a = self.store().run("bad", {}, {}, lambda: float("nan"))
        self.assertEqual(a.outcome["status"], "failed")
    def test_corrupt_blob_rejected(self):
        a = self.store().run("a", {}, {}, lambda: 1)
        p = self.root/"cache"/"objects"/(a.result_hash+".json.gz")
        p.write_bytes(gzip.compress(encoded({"wrong": 1}), mtime=0))
        with self.assertRaises(ValueError): self.store().run("a", {}, {}, lambda: 2)
    def test_keyboard_interrupt_not_cached_as_success(self):
        def interrupted(): raise KeyboardInterrupt()
        with self.assertRaises(KeyboardInterrupt): self.store().run("a", {}, {}, interrupted)
        self.assertEqual(self.store().run("a", {}, {}, lambda: 1).data, 1)
    def test_lock_conflict(self):
        with exclusive_lock(self.root/"lock"):
            with self.assertRaises(RuntimeError):
                with exclusive_lock(self.root/"lock"): pass


class DesignTests(unittest.TestCase):
    def setUp(self):
        self.config = compact(); self.driver = FakeDriver()
        self.freeze = self.driver.freeze(pair(), b.MODELS[0], self.config)
    def nodes(self, pts):
        rows = []
        for model in b.MODELS:
            rows += self.driver.states(self.driver.freeze(pair(), model, self.config), pts)
        return d.paired(rows, b.MODELS)
    def test_default_grid_matches_prior_144_coordinates(self):
        self.assertEqual(len(cc.basic_points(CONFIG)), 144)
    def test_config_boolean_rejected(self):
        c=compact(); c["execution"]["workers"]=True
        with self.assertRaises(ValueError): cc.validate(c)
    def test_scope_endpoints_required(self):
        c=compact(); c["coarse"]["temperatures_k"]=[350., 400.]
        with self.assertRaises(ValueError): cc.validate(c)
    def test_invalid_pressure_scope_rejected(self):
        c=compact(); c["scope"]["pressure_pa"][1]=1e7
        with self.assertRaises(ValueError): cc.validate(c)
    def test_duplicate_states_rejected(self):
        rows=self.driver.states(self.freeze, [[.1,350.,101325.]])
        with self.assertRaises(ValueError): d.paired(rows+rows, b.MODELS)
    def test_single_model_not_eligible(self):
        rows=self.driver.states(self.freeze, [[.1,350.,101325.]])
        self.assertFalse(d.paired(rows, b.MODELS)[0]["eligible"])
    def test_phase_change_triggers_temperature_refinement(self):
        nodes=self.nodes([[.1,350.,101325.],[.1,450.,101325.]])
        edges=d.triggered_edges(nodes,self.config["adaptive"])
        self.assertTrue(any(e["axis"]=="temperature_k" for e in edges))
    def test_pressure_midpoint_is_geometric(self):
        edges=d.triggered_edges(self.nodes([[.1,350.,101325.],[.1,350.,200000.]]),self.config["adaptive"])
        self.assertAlmostEqual(edges[0]["point"][2], (101325.*200000.)**.5, places=5)
    def test_budget_cap_and_no_duplicates(self):
        nodes=self.nodes([list(q) for q in cc.basic_points(self.config)])
        plan=d.refinement_plan(nodes,self.config["adaptive"],1)
        self.assertEqual(len(plan["points"]),1)
        self.assertNotIn(tuple(plan["points"][0]),{tuple(n["point"]) for n in nodes})
    def test_zero_budget_keeps_unresolved_count(self):
        nodes=self.nodes([[.1,350.,101325.],[.1,450.,101325.]])
        p=d.refinement_plan(nodes,self.config["adaptive"],0)
        self.assertFalse(p["points"]); self.assertGreater(p["unresolved_edge_count"],0)
    def test_all_new_points_inside_scope(self):
        nodes=self.nodes([list(q) for q in cc.basic_points(self.config)])
        for q in d.refinement_plan(nodes,self.config["adaptive"],100)["points"]:
            for axis,x in zip(d.AXES,q): self.assertTrue(self.config["scope"][axis][0]<=x<=self.config["scope"][axis][1])
    def test_local_gain_not_global_gain(self):
        nodes=self.nodes([[.1,350.,101325.],[.1,450.,101325.]])
        r=d.comparisons(nodes,self.config["adaptive"],self.config["scope"])
        self.assertEqual(r["paired_gain_point_count"],1)
        self.assertFalse(r["continuum_convergence_claimed"])
    def test_margins_are_not_uncertainty(self):
        r=d.comparisons(self.nodes([[.1,350.,101325.]]),self.config["adaptive"],self.config["scope"])
        self.assertFalse(r["empirically_validated"])
    def test_model_freeze_independent_of_adaptive_grid(self):
        a=b.numerical_config(self.config)
        self.config["adaptive"]["max_new_points_per_pair"]=500
        self.config["coarse"]["log_points"]=50
        self.assertEqual(a,b.numerical_config(self.config))


class SuiteTests(unittest.TestCase):
    def test_cycle_rejected(self):
        reg={"a":s.Study("a","1",("b",),(),lambda d,c:{}),"b":s.Study("b","1",("a",),(),lambda d,c:{})}
        with self.assertRaises(ValueError): s.ordered_studies(["a"],reg)
    def test_unknown_study_rejected(self):
        with self.assertRaises(ValueError): s.ordered_studies(["unknown"],s.registry())
    def test_default_order_dependency_valid(self):
        order=s.ordered_studies(CONFIG["studies"],s.registry())
        self.assertEqual(len(order),4)
    def test_missing_references_are_not_validated(self):
        out=s.reference_audit({"references":{"records":[]},"reference_predictions":{"rows":[]}}, {"models":b.MODELS})
        self.assertEqual(out["status"],"missing_external_references")
        self.assertFalse(out["physical_tps_validation"])
    def test_reference_inlet_mismatch_not_compared(self):
        ref={"id":"r","source":"test","additive_mass_fraction":.1,"outlet_temperature_k":350.,"pressure_pa":101325.,"inlet_temperature_k":280.,"value":1.}
        driver=FakeDriver(); rows=driver.states(driver.freeze(pair(),b.MODELS[0],compact()),[[.1,350.,101325.]])
        out=s.reference_audit({"references":{"records":[ref]},"reference_predictions":{"rows":rows}},{"models":b.MODELS})
        self.assertTrue(all(x["predicted_value"] is None for x in out["comparisons"]))


class PublicationTests(unittest.TestCase):
    def test_gzip_roundtrip(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp); raw=b'a,b\r\n1,2\r\n'
            pub.compressed(root,"x.gz",raw)
            self.assertEqual(gzip.decompress((root/"x.gz").read_bytes()),raw)
    def test_utf8_shards_preserve_every_row(self):
        with tempfile.TemporaryDirectory() as tmp:
            rows=[{"id":i,"text":"é"*12} for i in range(20)]
            index=pub.table(Path(tmp),"test",rows,120)
            parts=[p for p in index if p["file"].endswith(".csv")]
            self.assertEqual(sum(p["rows"] for p in parts),20)
            self.assertTrue(all(p["bytes"]<=120 for p in parts))
    def test_oversized_row_is_not_truncated(self):
        with tempfile.TemporaryDirectory() as tmp:
            with self.assertRaises(ValueError): pub.table(Path(tmp),"test",[{"text":"x"*200}],50)
    def test_diff_covers_removed_pairs(self):
        changes=pub.scientific_diff({"pairs":[{"pair_id":"a"}]},[{"pair_id":"b"}])
        self.assertEqual({x["change"] for x in changes},{"added","removed_or_out_of_scope"})
    def test_reference_model_labels_rejected(self):
        with tempfile.TemporaryDirectory() as tmp:
            p=Path(tmp)/"r.json"; write_json(p,{"schema":"mixture-campaign-measurements-v1","observations":[{"id":"x","source":"x","evidence_kind":"model"}]})
            with self.assertRaises(ValueError): cli.references(p)


class CampaignTests(unittest.TestCase):
    def setUp(self):
        self.tmp=tempfile.TemporaryDirectory(); self.root=Path(self.tmp.name)
        self.c=compact(); self.driver=FakeDriver()
    def tearDown(self): self.tmp.cleanup()
    def run_it(self, **kw):
        with redirect_stdout(io.StringIO()):
            return cli.execute(self.c,self.root/"cache",self.root/"results",driver_factory=lambda *a:self.driver,**kw)
    def test_end_to_end_snapshot_and_missing_evidence(self):
        summary,path=self.run_it()
        self.assertTrue(summary["suite_execution_complete"])
        self.assertGreater(summary["state_task_count"],0)
        self.assertFalse(summary["physical_validation_complete"])
        self.assertEqual(summary["pairs"][0]["reference_status"],"missing_external_references")
        self.assertTrue((path/"publication-index.json").is_file())
    def test_rerun_reuses_all_states(self):
        first,p1=self.run_it(); second,p2=self.run_it()
        self.assertNotEqual(p1,p2)
        self.assertEqual(second["executed_state_count"],0)
        self.assertEqual(first["scientific_result_sha256"],second["scientific_result_sha256"])
    def test_report_setting_does_not_rerun_solver(self):
        self.run_it(); self.c["publication"]["shard_bytes"]=24000
        summary,_=self.run_it()
        self.assertEqual(summary["executed_state_count"],0)
    def test_cold_rebuild_equals_incremental(self):
        first,_=self.run_it()
        with redirect_stdout(io.StringIO()):
            second,_=cli.execute(self.c,self.root/"clean",self.root/"fresh",driver_factory=lambda *a:self.driver)
        self.assertEqual(first["scientific_result_sha256"],second["scientific_result_sha256"])
    def test_parameter_fix_recomputes_and_preserves_old_snapshot(self):
        first,p1=self.run_it(); original=(p1/"scientific-results.json.gz").read_bytes()
        self.driver.pairs=[pair(bias=.2)]; self.c["catalog"]["test_revision"]=2
        second,_=self.run_it()
        self.assertGreater(second["executed_state_count"],0)
        self.assertNotEqual(first["scientific_result_sha256"],second["scientific_result_sha256"])
        self.assertEqual((p1/"scientific-results.json.gz").read_bytes(),original)
    def test_new_pair_reuses_old_pair_states(self):
        first,_=self.run_it(); self.driver.pairs.append(pair("b")); self.c["additional_cas"]=["b"]
        second,_=self.run_it()
        self.assertEqual(second["processed_pair_count"],2)
        self.assertEqual(second["cached_state_count"],first["state_task_count"])
    def test_forced_recompute_matches_cached_science(self):
        first,_=self.run_it(); second,_=self.run_it(recompute=True)
        self.assertGreater(second["executed_state_count"],0)
        self.assertEqual(first["scientific_result_sha256"],second["scientific_result_sha256"])
        self.assertEqual(second["same_key_changed_outcome_count"],0)
    def test_smoke_does_not_replace_full_latest(self):
        _,path=self.run_it(); self.run_it(limit_pairs=1)
        self.assertEqual(read_json(self.root/"results"/"latest.json")["run_id"],path.name)
    def test_backend_failure_stays_explicit(self):
        self.driver.fail=True
        summary,_=self.run_it()
        self.assertGreater(summary["numerical_failure_count"],0)
        self.assertEqual(summary["pairs"][0]["qualification"],"no_eligible_paired_states")
    def test_retry_failure_runs_dependents(self):
        self.driver.fail=True; self.run_it(); self.driver.fail=False
        summary,_=self.run_it(retry_failures=True)
        self.assertEqual(summary["numerical_failure_count"],0)
        self.assertGreater(summary["executed_state_count"],0)
    def test_new_study_backfills_old_pairs_without_solver(self):
        self.driver.pairs=[pair(),pair("b")]; self.run_it()
        text='from campaign_studies import Study\ndef extra(deps, settings):\n    return {"rows":len(deps["observations"]["rows"])}\ndef register_studies(reg):\n    reg["extra"] = Study("extra", "1", ("observations",), (), extra)\n'
        path=self.root/"extension_test.py"; path.write_text(text)
        sys.path.insert(0,str(self.root))
        try:
            self.c["plugins"]=["extension_test"]; self.c["studies"].append("extra")
            summary,_=self.run_it()
        finally:
            sys.path.remove(str(self.root)); sys.modules.pop("extension_test",None)
        self.assertEqual(summary["executed_state_count"],0)
        self.assertTrue(all(p["required_studies"]["extra"]["status"]=="complete" for p in summary["pairs"]))
    def test_failed_study_cannot_leave_complete_qualification(self):
        def fail(deps,settings): raise RuntimeError("synthetic study bug")
        reg=s.registry(); reg["broken"]=s.Study("broken","1",("observations",),(),fail)
        with patch.object(s,"registry",return_value=reg):
            self.c["studies"].append("broken"); summary,_=self.run_it()
        self.assertFalse(summary["suite_execution_complete"])
        self.assertEqual(summary["pairs"][0]["qualification"],"study_failed")


REAL_STACK = all(importlib.util.find_spec(m) is not None for m in ("thermo","chemicals","CoolProp"))
@unittest.skipUnless(REAL_STACK, "real thermo/chemicals/CoolProp stack unavailable")
class RealBackendTests(unittest.TestCase):
    def test_ethanol_frozen_context_and_state(self):
        c=compact(); cfg={"catalog":c["catalog"],"grid":{"inlet_temperature_k":293.15}}
        cat=b.worker({"operation":"catalog","config":cfg,"additional_cas":[]})
        p=next(p for p in cat["candidates"] if p["cas_number"]=="64-17-5")
        frozen=b.worker({"operation":"freeze","pair":b.model_pair(p,b.MODELS[0]),"model":b.MODELS[0],"config":b.numerical_config(c)})
        result=b.worker({"operation":"states","frozen":frozen,"points":[[.01,350.,101325.]]})
        self.assertEqual(result["rows"][0]["status"],"ok")
    def test_acetone_transition_both_models(self):
        c=compact(); cat=b.worker({"operation":"catalog","config":{"catalog":c["catalog"],"grid":{"inlet_temperature_k":293.15}},"additional_cas":[]})
        p=next(p for p in cat["candidates"] if p["cas_number"]=="67-64-1")
        for m in b.MODELS:
            frozen=b.worker({"operation":"freeze","pair":b.model_pair(p,m),"model":m,"config":b.numerical_config(c)})
            rows=b.worker({"operation":"states","frozen":frozen,"points":[[.3,350.,101325.],[.3,350.,200000.]]})["rows"]
            self.assertTrue(all(r["status"]=="ok" for r in rows))
            self.assertGreater(rows[0]["same_model_water_delta_h_ratio"],1)
            self.assertLess(rows[1]["same_model_water_delta_h_ratio"],1)




class AdditionalRegressionTests(unittest.TestCase):
    def test_runtime_change_invalidates_numeric_spec(self):
        with tempfile.TemporaryDirectory() as tmp:
            store=Store(Path(tmp)/"cache",Path(tmp)/"j")
            a=store.run("a",{}, {},lambda:1,environment={"version":"1"})
            b1=store.run("a",{}, {},lambda:1,environment={"version":"2"})
            self.assertNotEqual(a.key,b1.key)
    def test_historical_ids_rechecked(self):
        class HistoricalDriver(FakeDriver):
            def catalog(self,c,h):
                self.historical=deepcopy(h)
                return super().catalog(c,h)
        driver=HistoricalDriver()
        with tempfile.TemporaryDirectory() as tmp, redirect_stdout(io.StringIO()):
            root=Path(tmp)
            for _ in range(2): cli.execute(compact(),root/"cache",root/"out",driver_factory=lambda *a:driver)
            self.assertEqual(driver.historical,["a"])
    def test_empty_catalog_is_reported_not_vacuous_success(self):
        driver=FakeDriver(); driver.pairs=[]
        with tempfile.TemporaryDirectory() as tmp, redirect_stdout(io.StringIO()):
            root=Path(tmp); summary,_=cli.execute(compact(),root/"c",root/"r",driver_factory=lambda *a:driver)
            self.assertFalse(summary["suite_execution_complete"])
            self.assertEqual(summary["run_status"],"no_eligible_pairs")
    def test_subprocess_timeout_is_a_failure(self):
        driver=b.Driver.__new__(b.Driver)
        with tempfile.TemporaryDirectory() as tmp:
            driver.work=Path(tmp); driver.timeout=.01
            with patch.object(b.subprocess,"run",side_effect=subprocess.TimeoutExpired("test",.01)):
                with self.assertRaisesRegex(RuntimeError,"worker_timeout"): driver.call({"operation":"test"})
    def test_worker_wrong_pair_is_rejected(self):
        driver=b.Driver.__new__(b.Driver)
        frozen={"pair":{"pair_id":"a"},"model":"m"}
        driver.call=lambda req:{"rows":[{"additive_mass_fraction":.1,"outlet_temperature_k":350.,"pressure_pa":1e5,"pair_id":"WRONG","model":"m"}]}
        with self.assertRaises(ValueError):driver.states(frozen,[[.1,350.,1e5]])
    def test_new_reference_conditions_trigger_only_missing_states(self):
        driver=FakeDriver(); c=compact()
        ref={"id":"test","source":"synthetic-unit-test-only","value":100000.,"inlet_temperature_k":293.15,
             "outlet_temperature_k":366.,"pressure_pa":120000.,"additive_mass_fraction":.077,"additive_cas":"a"}
        with tempfile.TemporaryDirectory() as tmp, redirect_stdout(io.StringIO()):
            root=Path(tmp); cli.execute(c,root/"c",root/"r",driver_factory=lambda *a:driver)
            summary,path=cli.execute(c,root/"c",root/"r",driver_factory=lambda *a:driver,reference_records=[ref])
            self.assertEqual(summary["executed_state_count"],2)
            dossier=read_json(path/"pair-0001.json.gz")
            self.assertEqual(dossier["studies"]["reference_audit"]["data"]["status"],"compared_not_certified")
    def test_planner_failure_marks_suite_incomplete(self):
        with tempfile.TemporaryDirectory() as tmp, redirect_stdout(io.StringIO()):
            root=Path(tmp); driver=FakeDriver()
            with patch.object(d,"refinement_plan",side_effect=ValueError("synthetic planner bug")):
                summary,_=cli.execute(compact(),root/"c",root/"r",driver_factory=lambda *a:driver)
            self.assertFalse(summary["suite_execution_complete"])
            self.assertEqual(summary["pairs"][0]["adaptive_stop_reason"],"planner_failed")
    def test_publication_ignore_rules_with_git(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp); subprocess.run(["git","init","-q",str(root)],check=True)
            (root/".gitignore").write_bytes((HERE/".gitignore").read_bytes())
            allowed=["campaign-v5m3-results/latest.json","campaign-v5m3-results/run-test/summary.json",
                     "campaign-v5m3-results/run-test/report-pair-0001.csv","campaign-v5m3-results/run-test/pair-0001.json.gz"]
            blocked=[".campaign-cache/objects/x.json.gz","campaign-v5m3-results/.building-run-test/x",
                     "campaign-v5m3-results/run-test/workers/x.json","campaign-v5m3-results/run-test/secret.env",
                     "campaign-v5m3-results/run-test/summary.json/nested"]
            for name in allowed+blocked:
                proc=subprocess.run(["git","check-ignore","--no-index",name],cwd=root,capture_output=True)
                self.assertEqual(proc.returncode,0 if name in blocked else 1,name)

if __name__ == "__main__": unittest.main()
