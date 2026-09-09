"""V5m-2 tests. All non-integration fixtures are explicitly synthetic."""
from __future__ import annotations

import copy
import importlib.util
import io
import json
from pathlib import Path
import subprocess
import sys
import tempfile
from types import SimpleNamespace
import unittest
from unittest.mock import Mock, patch
from contextlib import redirect_stdout

import yaml
HERE = Path(__file__).resolve().parents[2] / "formulation-screen-v5m"
sys.path.insert(0, str(HERE))
import v5m2_core as c
import v5m2_backend as b
import run_predictive_grid as r
import validate_predictive_grid as v
CONFIG = yaml.safe_load((HERE / "config-v5m2.yaml").read_text())


def pair_fixture():
    return {"pair_id": "synthetic-pair", "name": "Synthetic model-test pair", "cas_number": "64-17-5",
            "smiles": "CCO", "inchi_key": "LFQSCWFLJHTTHZ-UHFFFAOYSA-N", "models": dict.fromkeys(c.MODELS, {})}


def row_fixture(model="chemsep_nrtl", w=.01, **kwargs):
    return {"pair_id": "synthetic-pair", "name": "Synthetic model-test pair", "cas_number": "64-17-5",
            "smiles": "CCO", "model": model, "additive_mass_fraction": w,
            "outlet_temperature_k": 400., "pressure_pa": 101325., "status": "ok",
            "model_comparison_eligible": True, "same_model_water_delta_h_ratio": .99,
            "storage_bubble_pressure_ratio": 1.1, "viscosity_ratio_proxy": 1.2, **kwargs}


def state_fixture(vf=0., liquids=1, h=0., rho=900.):
    return {"enthalpy_j_kg": h, "density_kg_m3": rho, "vapor_mole_fraction": vf,
            "liquid_phase_count": liquids, "vapor_additive_mass_fraction": .03 if vf > 0 else None}


class CoreTests(unittest.TestCase):
    def test_default_configuration_valid(self): c.validate_config(CONFIG)
    def test_boolean_numeric_rejected(self):
        with self.assertRaises(ValueError): c.number(True, "x")
    def test_nonfinite_rejected(self):
        for x in [float("nan"), float("inf"), -float("inf")]:
            with self.assertRaises(ValueError): c.number(x, "x")
    def test_grid_count_and_endpoints(self):
        xs = c.mass_grid(CONFIG)
        self.assertEqual(len(xs), 24)
        self.assertEqual(xs[0], .0001); self.assertEqual(xs[-1], .30)
        self.assertEqual(len(xs), len(set(xs)))
    def test_grid_mass_fraction_not_percentage(self):
        self.assertIn(.01, c.mass_grid(CONFIG))
        self.assertAlmostEqual(.01*100., 1.)
    def test_grid_is_deterministic(self): self.assertEqual(c.mass_grid(CONFIG), c.mass_grid(copy.deepcopy(CONFIG)))
    def test_invalid_mass_domain(self):
        cfg = copy.deepcopy(CONFIG); cfg["grid"]["maximum_mass_fraction"] = 1.
        with self.assertRaises(ValueError): c.validate_config(cfg)
    def test_fraction_range(self):
        for w in [-.1, 1.1]:
            with self.assertRaises(ValueError): c.mass_to_mole(w, [18., 46.])
    def test_mass_mole_roundtrip(self):
        for w in [0., .001, .03, .5, 1.]:
            self.assertAlmostEqual(c.mole_to_mass(c.mass_to_mole(w, [18.01528, 46.06844]), [18.01528, 46.06844])[1], w)
    def test_invalid_mole_sum(self):
        with self.assertRaises(ValueError): c.mole_to_mass([.4, .4], [18., 46.])
    def test_invalid_molecular_weight(self):
        with self.assertRaises(ValueError): c.mass_to_mole(.1, [18., -46.])
    def test_unknown_model_rejected(self):
        cfg = copy.deepcopy(CONFIG); cfg["models"] = ["ideal"]
        with self.assertRaises(ValueError): c.validate_config(cfg)
    def test_pressure_scope_rejected(self):
        cfg = copy.deepcopy(CONFIG); cfg["grid"]["pressures_pa"] = [1e7]
        with self.assertRaises(ValueError): c.validate_config(cfg)
    def test_negative_timeout_rejected(self):
        cfg = copy.deepcopy(CONFIG); cfg["pair_timeout_s"] = -1
        with self.assertRaises(ValueError): c.validate_config(cfg)
    def test_nan_config_rejected(self):
        cfg = copy.deepcopy(CONFIG); cfg["grid"]["outlet_temperatures_k"] = [float("nan")]
        with self.assertRaises(ValueError): c.validate_config(cfg)
    def test_digest_rejects_nonfinite(self):
        with self.assertRaises(ValueError): c.digest({"x": float("nan")})


class CatalogTests(unittest.TestCase):
    def test_automatic_water_partners(self):
        table = {"64-17-5 7732-18-5": {}, "7732-18-5 64-17-5": {}, "67-56-1 7732-18-5": {}, "1-1-1 2-2-2": {}}
        self.assertEqual(b.water_partners(table), ["64-17-5", "67-56-1"])
    def db(self):
        db = SimpleNamespace(metadata={"ChemSep NRTL": {"source": "synthetic fixture"}})
        db.has_ip_specific = lambda *args: True
        db.get_ip_specific = lambda name, order, field: .3 if field == "alphaij" else (100. if order[0] == c.WATER_CAS else 200.)
        return db
    def test_directional_nrtl(self):
        p = b.nrtl_parameters(self.db(), "64-17-5")
        self.assertEqual(p["component_order"], [c.WATER_CAS, "64-17-5"])
        self.assertEqual(p["tau_bs_k"], [[0., 100.], [200., 0.]])
    def test_missing_nrtl_not_zero_filled(self):
        db = self.db(); db.has_ip_specific = lambda *args: False
        with self.assertRaisesRegex(ValueError, "missing directed"): b.nrtl_parameters(db, "64-17-5")
    def test_nonfinite_nrtl_blocked(self):
        db = self.db(); db.get_ip_specific = lambda *args: float("nan")
        with self.assertRaises(ValueError): b.nrtl_parameters(db, "64-17-5")
    def groups(self):
        sg = {1: SimpleNamespace(main_group_id=1, R=1., Q=1.), 2: SimpleNamespace(main_group_id=2, R=1., Q=1.)}
        return [{1: 1}, {2: 1}], sg, {1: {2: (100., .1, .0)}, 2: {1: (50., .0, .0)}}
    def test_complete_directed_unifac(self):
        result = b.complete_group_interactions(*self.groups())
        self.assertEqual(result["missing_interactions"], [])
        self.assertEqual(len(result["interaction_sha256"]), 64)
    def test_missing_reverse_unifac_rejected(self):
        groups, sg, ip = self.groups(); ip[2] = {}
        with self.assertRaisesRegex(ValueError, "missing directed"): b.complete_group_interactions(groups, sg, ip)
    def test_unknown_subgroup_rejected(self):
        groups, sg, ip = self.groups(); groups[0] = {99: 1}
        with self.assertRaises(ValueError): b.complete_group_interactions(groups, sg, ip)
    def test_no_group_assignment_rejected(self):
        _, sg, ip = self.groups()
        with self.assertRaises(ValueError): b.complete_group_interactions([{}, {1: 1}], sg, ip)
    def test_missing_metadata_false_sentinel(self):
        for item in [False, None]:
            with self.assertRaisesRegex(ValueError, "missing_local_identity"):
                b.molecular_metadata(item, "64-17-5", CONFIG["catalog"])
    def test_metadata_consistency(self):
        item = SimpleNamespace(CASs="64-17-5", InChI_key="LFQSCWFLJHTTHZ-UHFFFAOYSA-N", smiles="CCO",
                               MW=46.069, common_name="ethanol", iupac_name="ethanol")
        result = b.molecular_metadata(item, "64-17-5", CONFIG["catalog"])
        self.assertFalse(result["commercial_availability_verified"])
    def test_metadata_structure_mismatch(self):
        item = SimpleNamespace(CASs="64-17-5", InChI_key="WRONG", smiles="CCO", MW=46.069)
        with self.assertRaisesRegex(ValueError, "inchikey_mismatch"): b.molecular_metadata(item, "64-17-5", CONFIG["catalog"])
    def test_electrolyte_outside_model(self):
        item = SimpleNamespace(CASs="x", smiles="[Na+]")
        with self.assertRaises(ValueError): b.molecular_metadata(item, "x", CONFIG["catalog"])


class BoundsAndFlashTests(unittest.TestCase):
    def obj(self):
        return SimpleNamespace(method="BOUNDED", T_limits={"BOUNDED": (280., 450.)}, tabular_data={},
                               test_method_validity=lambda t, m: 280. <= t <= 450.)
    def test_bounds_enforced(self):
        b.in_bounds(self.obj(), 300.)
        with self.assertRaises(ValueError): b.in_bounds(self.obj(), 500.)
    def test_unknown_bounds_rejected(self):
        obj = self.obj(); obj.T_limits = {}
        with self.assertRaises(ValueError): b.in_bounds(obj, 300.)
    def test_tabular_bounds_override(self):
        obj = self.obj(); obj.tabular_data = {"BOUNDED": ([300., 310.], [1., 2.])}
        with self.assertRaises(ValueError): b.in_bounds(obj, 320.)
    def test_method_choice_disables_extrapolation(self):
        obj = self.obj(); obj.all_methods = ["BOUNDED"]; obj.ranked_methods = ["BOUNDED"]
        meta = b.pin_correlation(obj, [298.15, 400.])
        self.assertIsNone(obj.extrapolation); self.assertFalse(obj.tabular_extrapolation_permitted)
        self.assertEqual(meta["method"], "BOUNDED")
    def flash(self, two_phase=False):
        liq = SimpleNamespace(zs=[.9, .1])
        if not two_phase:
            return SimpleNamespace(phases=[liq], betas=[1.], H=lambda: -1000., V=lambda: .00002,
                                   Cp=lambda: 80., liquids=[liq], gas=None, VF=0.)
        gas = SimpleNamespace(zs=[.1, .9])
        return SimpleNamespace(phases=[gas, liq], betas=[.5, .5], H=lambda: 1000., V=lambda: .01,
                               Cp=lambda: (_ for _ in ()).throw(AssertionError("two-phase Cp called")),
                               liquids=[liq], gas=gas, VF=.5)
    def test_single_phase_unit_conversion(self):
        row = c.summarize_flash(self.flash(), [.9, .1], [18., 46.], 1e-7)
        self.assertAlmostEqual(row["enthalpy_j_kg"], -1000./.0208)
        self.assertAlmostEqual(row["density_kg_m3"], .0208/.00002)
    def test_two_phase_cp_not_invented(self):
        row = c.summarize_flash(self.flash(True), [.5, .5], [18., 46.], 1e-7)
        self.assertIsNone(row["cp_j_kg_k"])
    def test_phase_mass_composition_converted(self):
        row = c.summarize_flash(self.flash(True), [.5, .5], [18., 46.], 1e-7)
        self.assertAlmostEqual(row["vapor_additive_mass_fraction"], .9*46/(.1*18+.9*46))
    def test_species_balance_failure(self):
        with self.assertRaisesRegex(ValueError, "species balance"):
            c.summarize_flash(self.flash(), [.8, .2], [18., 46.], 1e-7)
    def test_phase_fraction_failure(self):
        result = self.flash(); result.betas = [.5]
        with self.assertRaises(ValueError): c.summarize_flash(result, [.9, .1], [18., 46.], 1e-7)
    def test_invalid_phase_composition(self):
        result = self.flash(); result.phases[0].zs = [.9, -.1]
        with self.assertRaises(ValueError): c.summarize_flash(result, [.9, .1], [18., 46.], 1e-7)


class EvaluationTests(unittest.TestCase):
    def backend(self):
        obj = b.BinaryModel.__new__(b.BinaryModel)
        obj.pair, obj.model, obj.config = pair_fixture(), "chemsep_nrtl", CONFIG
        obj.state = lambda w, t, p: state_fixture(h=(t-293.15)*(100. if w == 0 else 95.), vf=0. if t < 350 else 1., liquids=1 if t < 350 else 0)
        obj.bubble = lambda w, t: {"bubble_pressure_pa": 1000*(1+w)}
        obj.viscosity_ratio = lambda w, t: 1.+w
        return obj
    def test_delta_h_and_mass_parity(self):
        obj = self.backend()
        with patch.object(b, "water_reference", return_value={"status": "ok", "delta_h_j_kg": (400-293.15)*100}):
            row = obj.evaluate(.01, 400., 101325.)
        self.assertEqual(row["status"], "ok")
        self.assertAlmostEqual(row["same_model_water_delta_h_ratio"], .95)
        self.assertAlmostEqual(row["mass_ratio_same_net_heat"], 1/.95)
        self.assertAlmostEqual(row["additional_net_heat_reduction_for_mass_parity"], .05)
        self.assertFalse(row["surface_enhancement_predicted"])
    def test_endpoint_discrepancy_blocks_comparison(self):
        with patch.object(b, "water_reference", return_value={"status": "ok", "delta_h_j_kg": 100.}):
            row = self.backend().evaluate(.01, 400., 101325.)
        self.assertFalse(row["model_comparison_eligible"])
    def test_inlet_two_liquids_rejected(self):
        obj = self.backend(); obj.state = lambda *args: state_fixture(liquids=2)
        self.assertEqual(obj.evaluate(.01, 400., 101325.)["status"], "inlet_not_single_liquid")
    def test_inlet_vapor_rejected(self):
        obj = self.backend(); obj.state = lambda *args: state_fixture(vf=.5)
        self.assertEqual(obj.evaluate(.01, 400., 101325.)["status"], "inlet_not_single_liquid")
    def test_negative_heat_window_not_scored(self):
        obj = self.backend(); obj.state = lambda *args: state_fixture(h=0.)
        self.assertEqual(obj.evaluate(.01, 400., 101325.)["status"], "model_state_rejected")
    def test_water_failure_no_fallback(self):
        with patch.object(b, "water_reference", return_value={"status": "water_reference_failed"}):
            row = self.backend().evaluate(.01, 400., 101325.)
        self.assertEqual(row["status"], "water_reference_unavailable")


class SelectionTests(unittest.TestCase):
    def aggregate(self, rows=None, n=1):
        return c.compare_models(rows or [row_fixture(m) for m in c.MODELS], n, CONFIG)[0]
    def test_both_models_required(self):
        self.assertFalse(self.aggregate([row_fixture()])["comparison_eligible"])
    def test_failures_not_dropped(self):
        rows = [row_fixture(), row_fixture(c.MODELS[1], status="worker_timeout")]
        self.assertFalse(self.aggregate(rows)["comparison_eligible"])
    def test_missing_condition_blocks_complete_grid(self): self.assertFalse(self.aggregate(n=2)["comparison_eligible"])
    def test_model_spread_not_uncertainty(self):
        row = self.aggregate([row_fixture(), row_fixture(c.MODELS[1], same_model_water_delta_h_ratio=.98)])
        self.assertAlmostEqual(row["maximum_model_spread_delta_h_ratio"], .01)
        self.assertIsNone(row["calibrated_uncertainty"])
        self.assertFalse(row["experimental_winner"])
    def test_duplicates_not_silently_overwritten(self):
        with self.assertRaises(ValueError): self.aggregate([row_fixture(), row_fixture()])
    def test_missing_viscosity_kept_missing(self):
        row = self.aggregate([row_fixture(viscosity_ratio_proxy=None), row_fixture(c.MODELS[1])])
        self.assertIsNone(row["maximum_viscosity_ratio_proxy"])
        self.assertEqual(c.pareto([row], c.FLOW_OBJECTIVES), [])
    def test_pareto_tradeoff_and_dominance(self):
        rows = [{"comparison_eligible": True, "h": 1., "p": 1.},
                {"comparison_eligible": True, "h": .9, "p": 1.1},
                {"comparison_eligible": True, "h": 1.1, "p": 1.2}]
        front = c.pareto(rows, {"h": "max", "p": "min"})
        self.assertEqual(front, [rows[0], rows[2]])
    def test_refinement_deterministic_new_midpoints(self):
        xs = c.mass_grid(CONFIG)
        rows = [{**self.aggregate(), "additive_mass_fraction": xs[3]}]
        fine = c.refinement_grid(xs, rows, CONFIG)
        self.assertEqual(fine, c.refinement_grid(xs, rows, CONFIG))
        self.assertFalse(set(fine) & set(xs))
        self.assertTrue(all(xs[0] < w < xs[-1] for w in fine))
    def test_shortlist_avoids_twenty_aliases_same_pair(self):
        row = self.aggregate(); rows = [{**row, "additive_mass_fraction": w} for w in [.01, .02, .03]]
        selected = c.shortlist(rows, CONFIG)
        self.assertEqual(len(selected), 2)
        self.assertEqual({r["acquisition_role"] for r in selected}, {"thermodynamic_tradeoff", "model_disagreement"})


class RunnerTests(unittest.TestCase):
    def task(self):
        cfg = copy.deepcopy(CONFIG); cfg["grid"]["outlet_temperatures_k"] = [400.]; cfg["grid"]["pressures_pa"] = [101325.]
        return {"pair": pair_fixture(), "config": cfg, "mass_fractions": [.01, .02], "stage": "coarse"}
    def test_failure_rows_cover_entire_requested_grid(self):
        task = self.task()
        rows = r.failure_rows(task["pair"], task["config"], task["mass_fractions"], "worker_timeout", "test", "coarse")
        self.assertEqual(len(rows), 4)
        self.assertTrue(all(not x["model_comparison_eligible"] for x in rows))
    def test_checkpoint_roundtrip(self):
        with tempfile.TemporaryDirectory() as folder:
            path = Path(folder)/"x.json"
            self.assertIsNone(r.checkpoint(path, "sig"))
            r.checkpoint(path, "sig", {"a": 1})
            self.assertEqual(r.checkpoint(path, "sig"), {"a": 1})
    def test_checkpoint_signature_mismatch(self):
        with tempfile.TemporaryDirectory() as folder:
            path = Path(folder)/"x.json"; r.checkpoint(path, "sig", {"a": 1})
            with self.assertRaises(ValueError): r.checkpoint(path, "other")
    def test_checkpoint_corruption_detected(self):
        with tempfile.TemporaryDirectory() as folder:
            path = Path(folder)/"x.json"; r.checkpoint(path, "sig", {"a": 1})
            payload = r.load_json(path); payload["data"]["a"] = 2; r.write_json(path, payload)
            with self.assertRaises(ValueError): r.checkpoint(path, "sig")
    def test_timeout_returned_as_records(self):
        with tempfile.TemporaryDirectory() as folder, patch.object(r.subprocess, "run", side_effect=subprocess.TimeoutExpired("test", 1)):
            out = r.run_worker(self.task(), Path(folder), 1)
        self.assertEqual(len(out["rows"]), 4)
        self.assertTrue(all(row["status"] == "worker_timeout" for row in out["rows"]))
    def test_worker_nonzero_return_no_success(self):
        with tempfile.TemporaryDirectory() as folder, patch.object(r.subprocess, "run", return_value=SimpleNamespace(returncode=1)):
            out = r.run_worker(self.task(), Path(folder), 1)
        self.assertTrue(all(row["status"] == "worker_failed" for row in out["rows"]))
    def test_missing_model_does_not_skip_rows(self):
        task = self.task(); task["pair"]["models"] = {}
        out = r.pair_calculation(task, factory=lambda *args: self.fail("unexpected model call"))
        self.assertEqual(len(out["rows"]), 4)
        self.assertTrue(all(row["status"] == "model_unavailable" for row in out["rows"]))
    def test_stage_resume_does_not_recompute(self):
        calls = []
        def worker(task, folder, timeout):
            calls.append(task)
            return {"rows": r.failure_rows(task["pair"], task["config"], task["mass_fractions"], "worker_timeout", "test", task["stage"]), "model_metadata": {}}
        with tempfile.TemporaryDirectory() as folder, redirect_stdout(io.StringIO()):
            first = r.execute_stage([self.task()], Path(folder), "sig", CONFIG, worker=worker)
            second = r.execute_stage([self.task()], Path(folder), "sig", CONFIG, worker=worker)
            self.assertEqual(first, second); self.assertEqual(len(calls), 1)
            r.execute_stage([self.task()], Path(folder), "sig", CONFIG, worker=worker, retry_failures=True)
            self.assertEqual(len(calls), 2)
    def test_json_rejects_nan(self):
        with tempfile.TemporaryDirectory() as folder:
            path = Path(folder)/"x.json"; path.write_text('{"a": NaN}')
            with self.assertRaises(ValueError): r.load_json(path)
    def test_csv_empty_still_has_header(self):
        with tempfile.TemporaryDirectory() as folder:
            path = Path(folder)/"x.csv"; r.write_csv(path, [])
            self.assertEqual(path.read_text().strip(), "status")


class ValidationTests(unittest.TestCase):
    def ref(self):
        return {"id": "synthetic-test-only", "source": "synthetic fixture for unit tests",
                "evidence_kind": "measured", "water_cas": c.WATER_CAS, "additive_cas": "64-17-5",
                "composition_basis": "liquid_mole", "liquid_additive_mole_fraction": .1,
                "phase": "homogeneous_liquid", "property": "excess_enthalpy_j_mol", "unit": "J/mol",
                "temperature_k": 300., "pressure_pa": 101325., "value": -100.}
    def load(self, rows): return v.validate_records({"schema": "binary-mixture-measurements-v1", "observations": rows})
    def test_signed_excess_enthalpy_accepted(self): self.assertEqual(self.load([self.ref()])[0]["value"], -100.)
    def test_model_reference_forbidden(self):
        row = self.ref(); row["evidence_kind"] = "model_prediction"
        with self.assertRaises(ValueError): self.load([row])
    def test_unlabeled_phase_rejected(self):
        row = self.ref(); row["phase"] = "unknown"
        with self.assertRaises(ValueError): self.load([row])
    def test_no_source_rejected(self):
        row = self.ref(); row["source"] = ""
        with self.assertRaises(ValueError): self.load([row])
    def test_mass_not_mole_cannot_silently_convert(self):
        row = self.ref(); row["composition_basis"] = "mass"
        with self.assertRaises(ValueError): self.load([row])
    def test_zero_HE_relative_error_null(self):
        row = self.ref(); row["value"] = 0.
        output, summary = v.compare_records(self.load([row]), lambda *args: 10.)
        self.assertIsNone(output[0]["relative_error"])
        self.assertFalse(summary["blind_holdout_claimed"])
    def test_failure_counted(self):
        def prediction(*args): raise ValueError("failed state")
        output, summary = v.compare_records(self.load([self.ref()]), prediction)
        self.assertEqual(sum(row["failed"] for row in summary["pair_property_results"]), 2)
    def test_duplicate_reference_id(self):
        with self.assertRaises(ValueError): self.load([self.ref(), self.ref()])


@unittest.skipUnless(all(importlib.util.find_spec(p) for p in ("thermo", "chemicals", "CoolProp")),
                     "real thermo/chemicals/CoolProp packages unavailable")
class IntegrationTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        from importlib.metadata import version
        for name, wanted in (("thermo", "0.6.1"), ("chemicals", "1.5.2"), ("CoolProp", "8.0.0")):
            if version(name) != wanted: raise unittest.SkipTest("integration requires pinned project libraries")
        cls.catalog = b.discover_catalog(CONFIG)
    def test_real_automatic_catalog_contains_ethanol(self):
        pairs = self.catalog["candidates"]
        self.assertGreater(len(pairs), 5)
        ethanol = next(p for p in pairs if p["cas_number"] == "64-17-5")
        self.assertEqual(set(ethanol["models"]), set(c.MODELS))
    def test_real_both_models_and_units(self):
        pair = next(p for p in self.catalog["candidates"] if p["cas_number"] == "64-17-5")
        for model in c.MODELS:
            obj = b.BinaryModel(pair, model, CONFIG)
            row = obj.evaluate(.01, 400., 101325.)
            self.assertEqual(row["status"], "ok", row)
            self.assertTrue(row["model_comparison_eligible"], row)
            self.assertGreater(row["delta_h_j_kg"], 1e6)
            self.assertLess(row["delta_h_j_kg"], 5e6)
            self.assertLess(row["inlet"]["component_balance_error"], 1e-7)
    def test_real_pure_endpoint_activity_limit(self):
        pair = next(p for p in self.catalog["candidates"] if p["cas_number"] == "64-17-5")
        for model in c.MODELS:
            obj = b.BinaryModel(pair, model, CONFIG)
            self.assertAlmostEqual(obj.bubble(0., 293.15)["water_activity"], 1., places=8)
    def test_real_worker_roundtrip(self):
        pair = next(p for p in self.catalog["candidates"] if p["cas_number"] == "64-17-5")
        cfg = copy.deepcopy(CONFIG); cfg["grid"]["outlet_temperatures_k"] = [400.]; cfg["grid"]["pressures_pa"] = [101325.]
        task = {"pair": pair, "config": cfg, "mass_fractions": [.01], "stage": "coarse"}
        with tempfile.TemporaryDirectory() as folder:
            output = r.run_worker(task, Path(folder), 90)
        self.assertEqual(len(output["rows"]), 2)
        self.assertTrue(all(row["status"] == "ok" for row in output["rows"]), output)


if __name__ == "__main__": unittest.main()
