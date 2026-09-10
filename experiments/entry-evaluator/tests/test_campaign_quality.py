"""Quality-suite regression tests. Synthetic fixtures are not physical evidence."""
from __future__ import annotations

from concurrent.futures import ThreadPoolExecutor
from contextlib import redirect_stdout
from copy import deepcopy
import hashlib
import io
import math
from pathlib import Path
import subprocess
import tempfile
import unittest
from unittest.mock import patch

from test_campaign import FakeDriver, compact, pair, REAL_STACK, HERE
import campaign_backend as backend
import campaign_config as config_module
import campaign_design as design
import campaign_quality as quality
import campaign_studies as studies
import campaign_water as water
import run_campaign as cli
from campaign_store import Store, digest, read_json, write_json


class WaterTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.root = Path(self.tmp.name)
        self.config, self.driver = compact(), FakeDriver()
        config_module.validate(self.config)
        self.store = Store(self.root/"cache", self.root/"events")
        self.service = water.SharedWaterAudit(self.store, self.driver, self.config)
        self.frozen = {m: self.store.run("fixture.freeze", {"model": m}, {},
            lambda m=m: self.driver.freeze(pair(), m, self.config)) for m in backend.MODELS}
    def tearDown(self):
        self.tmp.cleanup()
    def raw(self):
        return self.driver.states(self.frozen[backend.MODELS[0]].data, [[.1, 350., 101325.]])
    def bundle(self):
        return self.service.bundle(self.frozen, [[350., 101325.]])
    def test_context_is_shared_across_activity_models_and_additives(self):
        c1 = water.context(self.frozen[backend.MODELS[0]].data)
        c2 = water.context(self.driver.freeze(pair("different"), backend.MODELS[1], self.config))
        self.assertEqual(c1, c2)
        self.assertNotIn("pair_id", c1)
    def test_different_selected_correlations_do_not_share_context(self):
        frozen = deepcopy(self.frozen[backend.MODELS[0]].data)
        c1 = water.context(frozen)
        frozen["selected_correlations"]["pure_methods"]["VaporPressures"][0]["method"] = "different"
        self.assertNotEqual(c1, water.context(frozen))
    def test_identical_control_is_computed_once_under_concurrency(self):
        with ThreadPoolExecutor(max_workers=4) as pool:
            list(pool.map(lambda _: self.bundle(), range(8)))
        self.assertEqual(sum(call[1] for call in self.driver.calls if isinstance(call, tuple) and call[0] == "water_states"), 1)
    def test_control_tasks_have_no_mixture_freeze_dependency(self):
        self.bundle()
        a = next(t for t in self.store.graph() if t["spec"]["kind"] == "reference.water_state")
        self.assertEqual(a["spec"]["dependencies"], {})
        self.assertIn("pure_methods", a["spec"]["inputs"]["context"])
    def test_preflight_contains_scope_endpoints(self):
        points = water.control_grid(self.config)
        self.assertIn((350.,101325.), points)
        self.assertIn((450.,200000.), points)
    def test_good_control_does_not_change_raw_state(self):
        raw = self.raw(); before = deepcopy(raw)
        rows = water.annotate(raw, self.bundle().data, self.config)
        self.assertEqual(raw, before)
        self.assertTrue(rows[0]["model_comparison_eligible"])
        self.assertEqual(rows[0]["same_model_water_delta_h_ratio"], raw[0]["same_model_water_delta_h_ratio"])
    def test_reference_tolerance_blocks_without_erasing_solver_result(self):
        bundle = deepcopy(self.bundle().data)
        bundle["controls"][0]["outcome"]["data"].update(endpoint_relative_error=.06)
        row = water.annotate(self.raw(), bundle, self.config)[0]
        self.assertEqual(row["status"], "ok")
        self.assertFalse(row["model_comparison_eligible"])
        self.assertIn("water_reference_tolerance_exceeded", row["comparison_blockers"])
    def test_mismatched_embedded_denominator_blocks_reuse(self):
        raw = self.raw(); raw[0]["same_model_water_delta_h_j_kg"] *= 1.01
        row = water.annotate(raw, self.bundle().data, self.config)[0]
        self.assertIn("water_control_context_mismatch", row["comparison_blockers"])
        self.assertFalse(row["model_comparison_eligible"])
    def test_missing_control_fails_closed(self):
        row = water.annotate(self.raw(), {"controls": [], "model_contexts": {}, "context_errors": {}}, self.config)[0]
        self.assertFalse(row["model_comparison_eligible"])
        self.assertIn("water_control_unavailable", row["comparison_blockers"])
    def test_water_phase_mismatch_blocks_comparison(self):
        bundle = deepcopy(self.bundle().data)
        bundle["controls"][0]["outcome"]["data"]["status"] = "water_phase_mismatch"
        row = water.annotate(self.raw(), bundle, self.config)[0]
        self.assertFalse(row["model_comparison_eligible"])
    def test_raw_ineligible_state_cannot_be_promoted(self):
        raw = self.raw(); raw[0]["model_comparison_eligible"] = False
        self.assertFalse(water.annotate(raw, self.bundle().data, self.config)[0]["model_comparison_eligible"])
    def test_duplicate_worker_coordinates_are_failed_tasks(self):
        self.driver.water_states = lambda ctx, points: [{"coordinate": points[0], "status": "ok"}]*2
        self.service.ensure(water.context(self.frozen[backend.MODELS[0]].data), [[350.,101325.]])
        self.assertTrue(all(a.outcome["status"] == "failed" for a in self.service.controls.values()))
    def test_corrupt_control_cache_is_fatal_not_a_measurement_gap(self):
        self.bundle()
        artifact = next(iter(self.service.controls.values()))
        (self.root/"cache"/"objects"/(artifact.result_hash+".json.gz")).write_bytes(b"corrupt")
        fresh = water.SharedWaterAudit(Store(self.root/"cache", self.root/"fresh"), self.driver, self.config)
        with self.assertRaises((ValueError, OSError)):
            fresh.bundle(self.frozen, [[350.,101325.]])
    def test_worker_timeout_is_not_a_success(self):
        driver = water.AuditDriver.__new__(water.AuditDriver)
        driver.work, driver.timeout = self.root/"workers", .01
        with patch.object(water.subprocess, "run", side_effect=subprocess.TimeoutExpired("test", .01)):
            with self.assertRaisesRegex(RuntimeError, "water_worker_timeout"):
                driver.water_states({}, [[350., 101325.]])
    def test_worker_failure_can_be_retried_without_rewriting_old_object(self):
        self.driver.water_states = lambda *_: (_ for _ in ()).throw(ValueError("fixture failure"))
        old = self.bundle()
        a = next(iter(self.service.controls.values()))
        self.assertEqual(a.outcome["status"], "failed")
        driver = FakeDriver()
        store = Store(self.root/"cache", self.root/"retry", retry_failures=True)
        service = water.SharedWaterAudit(store, driver, self.config)
        new = service.bundle(self.frozen, [[350., 101325.]])
        self.assertNotEqual(new.result_hash, old.result_hash)
        self.assertTrue((self.root/"cache"/"objects"/(a.result_hash+".json.gz")).is_file())


class DecisionTests(unittest.TestCase):
    def setUp(self):
        self.config = compact(); config_module.validate(self.config)
        self.policy = self.config["adaptive"]
    def node(self, x, low=1.02, high=1.12, *, blocked=False, vapor=False):
        return {"point": [x,350.,101325.], "eligible": not blocked,
            "minimum_ratio": None if blocked else low, "maximum_ratio": None if blocked else high,
            "spread": None if blocked else high-low,
            "phase_signatures": {m: ["ok",1,vapor] for m in backend.MODELS},
            "models": {m: {"comparison_blockers": ["water_reference_tolerance_exceeded"] if blocked else []} for m in backend.MODELS}}
    def test_small_persistent_disagreement_goes_to_evidence_not_refinement(self):
        report = design.edge_analysis([self.node(.1), self.node(.105)], self.policy)
        self.assertEqual(len(report["evidence_limited_edges"]), 1)
        self.assertFalse(report["refinable_edges"])
    def test_large_disagreement_interval_is_still_refined(self):
        self.assertTrue(design.triggered_edges([self.node(.1),self.node(.2)],self.policy))
    def test_changing_disagreement_is_not_declared_persistent(self):
        report = design.edge_analysis([self.node(.1),self.node(.105,high=1.3)],self.policy)
        self.assertFalse(report["evidence_limited_edges"])
        self.assertTrue(report["refinable_edges"])
    def test_gain_boundary_is_preserved_inside_evidence_limited_edge(self):
        report = design.edge_analysis([self.node(.1,low=.999,high=1.099),self.node(.105,low=1.005,high=1.105)],self.policy)
        self.assertTrue(report["evidence_limited_edges"])
        self.assertIn("paired_gain_boundary", report["refinable_edges"][0]["reasons"])
    def test_phase_change_survives_water_reference_block(self):
        report = design.edge_analysis([self.node(.1,blocked=True),self.node(.105,blocked=True,vapor=True)],self.policy)
        self.assertTrue(report["refinable_edges"])
        self.assertTrue(report["reference_blocked_edges"])
        self.assertFalse(report["evidence_limited_edges"])
    def test_water_block_without_phase_change_does_not_generate_gain_grid(self):
        nodes = [self.node(.1,blocked=True),self.node(.105,blocked=True)]
        plan = design.refinement_plan(nodes,self.policy,100)
        self.assertFalse(plan["points"])
        self.assertEqual(design.stop_assessment(plan,100,False)["primary_reason"],"reference_blocked")
    def test_resolved_phase_boundary_has_distinct_stop_status(self):
        nodes=[self.node(.1,high=1.021),self.node(.1001,high=1.021,vapor=True)]
        plan=design.refinement_plan(nodes,self.policy,10)
        self.assertEqual(design.stop_assessment(plan,10,False)["primary_reason"],"sampled_boundaries_resolved")
    def test_simultaneous_budget_reference_and_evidence_causes_are_retained(self):
        plan={"unresolved_edge_count":1,"reference_blocked_point_count":2,"evidence_limited_edge_count":3}
        stop=design.stop_assessment(plan,0,True)
        self.assertEqual(stop["causes"],["point_budget","reference_blocked","model_disagreement_requires_evidence"])
        self.assertFalse(stop["continuum_convergence_claimed"])
    def test_round_limit_is_not_point_budget(self):
        stop=design.stop_assessment({"unresolved_edge_count":1},20,True)
        self.assertEqual(stop["primary_reason"],"max_rounds")
    def test_invalid_quality_configuration_is_rejected(self):
        for key, value in (("temperature_points", True),("pressure_points",1),("consistency_relative_tolerance",float("nan"))):
            c=compact(); c["water_audit"][key]=value
            with self.assertRaises(ValueError): config_module.validate(c)


class QualityCampaignTests(unittest.TestCase):
    def setUp(self):
        self.tmp=tempfile.TemporaryDirectory(); self.root=Path(self.tmp.name)
        self.config,self.driver=compact(),FakeDriver()
    def tearDown(self): self.tmp.cleanup()
    def execute(self, **kwargs):
        with redirect_stdout(io.StringIO()):
            return cli.execute(self.config,self.root/"cache",self.root/"out",driver_factory=lambda *a:self.driver,**kwargs)
    def test_two_pairs_share_water_states(self):
        self.driver.pairs=[pair("a"),pair("b")]; self.config["execution"]["workers"]=2
        summary,path=self.execute()
        records=read_json(path/"water-controls.json.gz")
        self.assertEqual(len(records["contexts"]),1)
        self.assertEqual(summary["water_control_state_count"],len(records["controls"]))
        self.assertTrue(all(p["required_studies"]["evidence_needs"]["status"]=="complete" for p in summary["pairs"]))
    def test_quality_policy_change_reuses_mixture_states(self):
        self.execute(); self.config["water_audit"]["temperature_points"]=5
        summary,_=self.execute()
        self.assertEqual(summary["executed_state_count"],0)
    def test_evidence_budget_change_reuses_all_numerical_states(self):
        self.execute(); self.config["evidence"]["max_points_per_pair"]=2
        summary,path=self.execute()
        events=read_json(path/"execution-events.json.gz")
        self.assertEqual(summary["executed_state_count"],0)
        self.assertFalse(any(e["kind"]=="reference.water_state" and e["action"]=="executed" for e in events))
    def test_old_snapshots_are_immutable_after_quality_change(self):
        _,p1=self.execute(); raw=(p1/"scientific-results.json.gz").read_bytes()
        self.config["evidence"]["max_points_per_pair"]=2
        _,p2=self.execute()
        self.assertEqual((p1/"scientific-results.json.gz").read_bytes(),raw)
        self.assertNotEqual(p1,p2)
    def test_no_reference_is_fabricated(self):
        summary,path=self.execute()
        needs=read_json(path/"evidence-needs.json.gz")
        self.assertGreater(len(needs["requests"]),0)
        self.assertTrue(all(r["reference_value"] is None and r["sources"]==[] for r in needs["requests"]))
        self.assertFalse(summary["rankable_promotions"])
    def test_exact_declared_measurement_is_available_but_not_certified(self):
        _,path=self.execute()
        request=next(r for r in read_json(path/"evidence-needs.json.gz")["requests"] if r["kind"]=="mixture_caloric_reference")
        record={**request["conditions"],"id":"fixture-measurement", "property":"delta_h_j_kg", "unit":"J/kg",
                "evidence_kind":"measured","source":"synthetic-test-fixture-only", "value":123456.}
        _,path=self.execute(reference_records=[record])
        request=next(r for r in read_json(path/"evidence-needs.json.gz")["requests"] if r["id"]==request["id"])
        self.assertEqual(request["status"],"source_linked_measurement_available_not_certified")
        self.assertFalse(request["satisfied"])
        self.assertEqual(request["sources"][0]["value"],123456.)
    def test_reference_at_different_inlet_does_not_close_need(self):
        _,path=self.execute()
        request=next(r for r in read_json(path/"evidence-needs.json.gz")["requests"] if r["kind"]=="mixture_caloric_reference")
        record={**request["conditions"],"id":"fixture", "property":"delta_h_j_kg", "unit":"J/kg",
            "evidence_kind":"measured","source":"fixture", "value":123456.,"inlet_temperature_k":280.}
        _,path=self.execute(reference_records=[record])
        r=next(r for r in read_json(path/"evidence-needs.json.gz")["requests"] if r["id"]==request["id"])
        self.assertEqual(r["status"],"pending_external_reference")
    def test_water_failure_is_reported_and_can_be_retried(self):
        self.driver.water_states=lambda *_: (_ for _ in ()).throw(ValueError("fixture control failure"))
        summary,_=self.execute()
        self.assertGreater(summary["water_control_task_failure_count"],0)
        self.assertFalse(summary["pairs"][0]["rankable"])
        self.driver=FakeDriver()
        summary,_=self.execute(retry_failures=True)
        self.assertEqual(summary["water_control_task_failure_count"],0)
    def test_quality_artifacts_have_verified_publication_hashes(self):
        _,path=self.execute()
        index=read_json(path/"publication-index.json")
        files={r["file"]:r for r in index["files"]}
        for name in ("water-controls.json.gz","evidence-needs.json.gz"):
            self.assertEqual(hashlib.sha256((path/name).read_bytes()).hexdigest(),files[name]["sha256"])
        self.assertTrue(all(r["uncompressed_bytes"]<=self.config["publication"]["shard_bytes"] for r in files.values() if r["file"].startswith("report-")))
    def test_unmodified_backend_signature_matches_original_archive(self):
        self.assertEqual(backend.code_signature()["campaign_backend.py"],"a367d05a57b4e4494be1fc315c1425687507d1a50d036932823318a7c44b6ec4")
    def test_control_failure_does_not_mutate_raw_cache(self):
        self.driver.water_states=lambda *_: (_ for _ in ()).throw(ValueError("fixture"))
        _,path=self.execute()
        graph=read_json(path/"task-graph.json.gz")
        a=next(t for t in graph if t["spec"]["kind"]=="model.state")
        raw=read_json(self.root/"cache"/"objects"/(a["result_hash"]+".json.gz"))["outcome"]["data"]
        self.assertTrue(raw["model_comparison_eligible"])
        self.assertNotIn("comparison_blockers",raw)
    def test_new_artifacts_allowed_by_gitignore(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp); subprocess.run(["git","init","-q",str(root)],check=True)
            (root/".gitignore").write_bytes((HERE/".gitignore").read_bytes())
            for name in ("water-controls.json.gz","evidence-needs.json.gz","report-evidence-needs-0001.csv.gz"):
                p=subprocess.run(["git","check-ignore","--no-index","campaign-v5m3-results/run-test/"+name],cwd=root,capture_output=True)
                self.assertEqual(p.returncode,1)


@unittest.skipUnless(REAL_STACK,"real thermo/chemicals/CoolProp stack unavailable")
class RealWaterControlTests(unittest.TestCase):
    def test_shared_pure_limit_matches_both_binary_models(self):
        c=compact()
        catalog=backend.worker({"operation":"catalog","config":{"catalog":c["catalog"],"grid":{"inlet_temperature_k":293.15}},"additional_cas":[]})
        for cas in ("67-64-1","64-17-5"):
            p=next(p for p in catalog["candidates"] if p["cas_number"]==cas)
            for model in backend.MODELS:
                frozen=backend.worker({"operation":"freeze","pair":backend.model_pair(p,model),"model":model,"config":backend.numerical_config(c)})
                points=[[.1,350.,101325.],[.1,362.5,101325.],[.1,375.,200000.],[.1,400.,101325.]]
                raw=backend.worker({"operation":"states","frozen":frozen,"points":points})["rows"]
                controls=water.calculate({"context":water.context(frozen),"coordinates":[q[1:] for q in points]})["rows"]
                for row,check in zip(raw,controls):
                    self.assertEqual(check["status"],"ok")
                    self.assertEqual(row["status"],"ok")
                    self.assertTrue(math.isclose(check["model_delta_h_j_kg"],row["same_model_water_delta_h_j_kg"],rel_tol=1e-7))
                    self.assertAlmostEqual(check["endpoint_relative_error"],row["endpoint_relative_error"],places=7)
    def test_pinned_stack_water_worker_subprocess(self):
        backend.runtime()
        with tempfile.TemporaryDirectory() as tmp:
            c=compact(); d=water.AuditDriver(Path(tmp),60.)
            p=next(p for p in d.catalog(c,[])["candidates"] if p["cas_number"]=="67-64-1")
            frozen=d.freeze(p,backend.MODELS[0],c)
            self.assertEqual(d.water_states(water.context(frozen),[[350.,101325.]])[0]["status"],"ok")


if __name__ == "__main__":
    unittest.main()
