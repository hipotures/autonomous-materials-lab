"""V5b-3 regression tests; reference/FeOS integrations run when installed."""
from __future__ import annotations

import copy
import importlib.util
import io
import json
from pathlib import Path
import subprocess
import sys
import tempfile
from contextlib import redirect_stdout
from types import SimpleNamespace
import unittest
from unittest.mock import patch

import yaml
V5B = Path(__file__).resolve().parents[2] / "property-predictor-v5b"
sys.path.insert(0, str(V5B))
import v5b3_audit as a
import v5b3_metrics as m
import v5b3_predict as p
import run_v5b3 as runner

CONFIG = yaml.safe_load((V5B / "benchmark-v5b3.yaml").read_text())
KEY = "LFQSCWFLJHTTHZ-UHFFFAOYSA-N"


def fixture():
    def constant(prop, value, method="IUPAC"):
        return {"value": value, "unit": a.UNITS[prop], "method": method, "source_family": "spoofed"}
    def series(prop, method, values):
        return {"unit": a.UNITS[prop], "method": method, "source_family": "spoofed",
                "points": [{"temperature_k": t, "value": v} for t, v in values]}
    return {"candidate_id": "ethanol-test", "smiles": "CCO", "canonical_smiles": "CCO", "inchi_key": KEY,
            "selected_cas_number": "64-17-5", "family": "oxygenated", "reference_quality": "high",
            "packet": {"cas_number": "64-17-5", "molecular_weight_g_mol": 46.069,
                       "constants": {a.TC: constant(a.TC, 514.0), a.PC: constant(a.PC, 6.14e6)},
                       "series": {a.PS: series(a.PS, "ANTOINE_POLING", [(300.0, 10000.0), (350.0, 80000.0)]),
                                  a.CP: series(a.CP, "POLING_CONST", [(293.15, 2400.0), (350.0, 2400.0)]),
                                  a.HV: series(a.HV, "CRC_HVAP_TB", [(300.0, 9e5), (350.0, 8e5)])}}}


class FakeBackend:
    def identity(self, cas):
        return {"cas": cas, "inchi_key": KEY, "smiles": "CCO", "formula": "C2H6O", "mw": 46.069}
    def constant(self, cas, prop, method):
        return {a.TC: 514.0, a.PC: 6.14e6, a.TB: 351.44}[prop]
    def alternatives(self, cas, prop):
        return []
    def native_point(self, cas, prop, method, packet):
        return (351.44, 38600.0) if prop == a.HV else (298.15, 110.5)
    def series_value(self, cas, prop, method, t, packet):
        if not 280.0 <= t <= 400.0:
            raise ValueError("outside_reference_bounds")
        return {300.0: 10000.0, 350.0: 80000.0}[t], [280.0, 400.0]


class AuditTests(unittest.TestCase):
    def setUp(self):
        self.target, self.backend = fixture(), FakeBackend()
    def audit(self):
        return a.audit_target(self.target, CONFIG["audit"], self.backend)
    def test_poling_aliases_are_one_source_group(self):
        self.assertEqual({a.source_group(x) for x in ("ANTOINE_POLING", "ANTOINE_EXTENDED_POLING", "WAGNER_POLING", "POLING_CONST")}, {"POLING"})
    def test_allowlist_rejects_unknown_predictors(self):
        for method in ("JOBACK", "RACKETT", "UNKNOWN", "DIRECT_NIST_V5E2"):
            self.assertFalse(a.allowed(a.TC, method))
        self.assertFalse(a.allowed(a.CP, "ZABRANSKY_NEW_GUESSED_METHOD"))
    def test_valid_cas_and_checksum(self):
        self.assertTrue(a.valid_cas("64-17-5"))
        self.assertFalse(a.valid_cas("64-17-6"))
    def test_exact_local_identity(self):
        result = self.audit()
        self.assertTrue(result["identity_verified"])
        self.assertEqual(result["molecule_group"], KEY[:14])
    def test_inchikey_mismatch_blocks_all_properties(self):
        self.target["inchi_key"] = "AAAAAAAAAAAAAA-UHFFFAOYSA-N"
        self.assertEqual(self.audit()["observations"], [])
    def test_conflicting_cas_metadata_blocks_all_properties(self):
        self.backend.identity = lambda cas: {"cas": cas, "inchi_key": "OTHER"}
        result = self.audit()
        self.assertFalse(result["identity_verified"])
        self.assertIn("cas_identity_mismatch", result["identity_reason"])
    def test_mw_mismatch_blocks_mass_conversion(self):
        self.target["packet"]["molecular_weight_g_mol"] = 90.0
        self.assertFalse(self.audit()["identity_verified"])
    def test_crc_watson_grid_replaced_by_native_point(self):
        rows = [x for x in self.audit()["observations"] if x["property"] == a.HV]
        self.assertEqual(len(rows), 1)
        self.assertEqual(rows[0]["temperature_k"], 351.44)
        self.assertAlmostEqual(rows[0]["reference_value"], 38600.0 * 1000 / 46.069)
        self.assertEqual(rows[0]["evidence_kind"], "native_reference_point")
    def test_cp_constant_is_single_native_diagnostic(self):
        rows = [x for x in self.audit()["observations"] if x["property"] == a.CP]
        self.assertEqual(len(rows), 1)
        self.assertEqual(rows[0]["temperature_k"], 298.15)
        self.assertFalse(rows[0]["benchmark_eligible"])
    def test_unknown_pressure_not_claimed_as_matched(self):
        self.assertFalse(a.physical_basis(a.RHO, "DIPPR_PERRY_8E")[1])
        self.assertFalse(a.physical_basis(a.CP, "COOLPROP")[1])
        self.assertTrue(a.physical_basis(a.RHO, "VDI_TABULAR")[1])
    def test_no_requirement_for_complete_packet(self):
        self.target["packet"]["constants"] = {}
        self.target["packet"]["series"] = {a.PS: self.target["packet"]["series"][a.PS]}
        self.assertEqual(len(self.audit()["observations"]), 2)
    def test_reference_outside_declared_range_rejected(self):
        self.target["packet"]["series"][a.PS]["points"] = [{"temperature_k": 450.0, "value": 1e5}]
        result = self.audit()
        self.assertFalse(any(x["property"] == a.PS for x in result["observations"]))
    def test_reproduction_mismatch_rejected(self):
        self.target["packet"]["constants"][a.TC]["value"] = 900.0
        self.assertFalse(any(x["property"] == a.TC for x in self.audit()["observations"]))
    def test_unit_mismatch_rejected(self):
        self.target["packet"]["constants"][a.PC]["unit"] = "bar"
        self.assertFalse(any(x["property"] == a.PC for x in self.audit()["observations"]))
    def test_constant_source_conflict_quarantined(self):
        self.backend.alternatives = lambda cas, prop: [{"method": "CRC", "source_group": "CRC", "value": 100.0}]
        self.assertFalse(any(x["property"] in a.CONSTANTS for x in self.audit()["observations"]))
    def test_unknown_measurement_uncertainty_stays_null(self):
        obs = self.audit()["observations"][0]
        self.assertIsNone(obs["provenance"]["measurement_uncertainty"])
        self.assertNotEqual(obs["provenance"]["source_group"], "spoofed")
    def test_native_hvap_does_not_call_watson(self):
        backend = a.LocalReferenceBackend.__new__(a.LocalReferenceBackend)
        backend.object = lambda *args: SimpleNamespace(all_methods={"CRC_HVAP_TB", "CRC_HVAP_298"},
            CRC_HVAP_TB_Tb=351.44, CRC_HVAP_TB_Hvap=38600.0, CRC_HVAP_298=42000.0)
        self.assertEqual(backend.native_point("64-17-5", a.HV, "CRC_HVAP_TB", {}), (351.44, 38600.0))
        self.assertEqual(backend.native_point("64-17-5", a.HV, "CRC_HVAP_298", {}), (298.15, 42000.0))
    def test_unknown_bounds_fail_closed(self):
        backend = a.LocalReferenceBackend.__new__(a.LocalReferenceBackend)
        backend.object = lambda *args: SimpleNamespace(all_methods={"VDI_PPDS"}, T_limits={}, tabular_data={})
        with self.assertRaisesRegex(ValueError, "reference_bounds_unknown"):
            backend.series_value("64-17-5", a.PS, "VDI_PPDS", 300.0, {})


class PredictionTests(unittest.TestCase):
    def request(self):
        audited = a.audit_target(fixture(), CONFIG["audit"], FakeBackend())
        return runner.prediction_request(audited)
    def test_worker_never_receives_references_or_cas(self):
        request = self.request()
        p.validate_request(request)
        text = json.dumps(request)
        self.assertNotIn("reference_value", text)
        self.assertNotIn("64-17-5", text)
        self.assertNotIn("method", text)
    def test_reference_leakage_field_rejected(self):
        request = self.request()
        request["requests"][0]["reference_value"] = 100.0
        with self.assertRaises(ValueError):
            p.validate_request(request)
    def test_failed_property_does_not_destroy_other_predictions(self):
        class Adapter:
            def __init__(self, smiles): pass
            def metadata(self): return {}
            def predict(self, request):
                if request["property"] == a.PC:
                    raise ValueError("numeric failure")
                return 100.0
        result = p.run_prediction(self.request(), Adapter)
        self.assertEqual(sum(x["status"] == "prediction_failed" for x in result["results"]), 1)
        self.assertGreater(sum(x["status"] == "ok" for x in result["results"]), 0)
    def test_provider_failure_recorded_for_every_request(self):
        def factory(smiles): raise ValueError("out of GC domain")
        request = self.request()
        result = p.run_prediction(request, factory)
        self.assertEqual(len(result["results"]), len(request["requests"]))
        self.assertTrue(all(r["status"] == "prediction_failed" for r in result["results"]))
    def test_subprocess_timeout_is_explicit(self):
        with tempfile.TemporaryDirectory() as folder:
            with patch.object(runner.subprocess, "run", side_effect=subprocess.TimeoutExpired("worker", 1)):
                result = runner.predict_subprocess(self.request(), Path(folder), 1)
        self.assertTrue(all(r["status"] == "worker_timeout" for r in result["results"]))
    def test_nonfinite_prediction_rejected(self):
        with self.assertRaises(ValueError): p.finite_positive(float("nan"))
    def test_feos_adapter_converts_mass_properties_and_hvap(self):
        adapter = p.FeosBenchmarkAdapter.__new__(p.FeosBenchmarkAdapter)
        adapter.si = SimpleNamespace(KELVIN=1., PASCAL=1., KILOGRAM=1., METER=1., JOULE=1.)
        liquid = SimpleNamespace(pressure=lambda: 1e5, specific_enthalpy=lambda: 1e5,
                                 mass_density=lambda: 789., specific_isobaric_heat_capacity=lambda: 2400.)
        vapor = SimpleNamespace(specific_enthalpy=lambda: 9e5)
        adapter.equilibrium = lambda t: SimpleNamespace(liquid=liquid, vapor=vapor)
        for prop, expected in [(a.HV, 8e5), (a.PS, 1e5), (a.RHO, 789.), (a.CP, 2400.)]:
            self.assertEqual(adapter.predict({"property": prop, "temperature_k": 300.}), expected)


class MetricsTests(unittest.TestCase):
    def rows(self):
        observations = a.audit_target(fixture(), CONFIG["audit"], FakeBackend())["observations"]
        predictions = [{"observation_id": r["observation_id"], "predicted_value": 1.1*r["reference_value"], "status": "ok"} for r in observations]
        return m.comparison_rows(observations, predictions, CONFIG["statistics"])
    def test_temperatures_never_split_across_groups(self):
        self.assertEqual(len({r["split"] for r in self.rows()}), 1)
    def test_molecule_not_gridpoint_is_statistical_unit(self):
        rows = [r for r in self.rows() if r["property"] == a.PS]
        self.assertEqual(len(m.molecule_metrics(rows)), 1)
        self.assertAlmostEqual(m.molecule_metrics(rows)[0]["median_absolute_relative_error"], .1)
    def test_insufficient_data_never_creates_finite_envelope(self):
        report = m.summarize(self.rows(), CONFIG["statistics"])
        self.assertTrue(all(r["maximum_error_envelope"] is None for r in report["properties"].values()))
        self.assertEqual(report["rankable_promotions"], 0)
    def test_diagnostic_points_do_not_calibrate(self):
        report = m.summarize(self.rows(), CONFIG["statistics"])["properties"][a.CP]
        self.assertEqual(report["strict_molecule_group_count"], 0)
        self.assertEqual(report["diagnostic_only_group_count"], 1)
    def test_missing_predictions_are_counted_as_failures(self):
        rows = self.rows()
        for row in rows:
            row["absolute_relative_error"] = None
            row["status"] = "worker_timeout"
        self.assertTrue(all(r["failed_point_count"] == r["point_count"] for r in m.molecule_metrics(rows)))
    def test_failure_on_evaluation_is_not_dropped_from_coverage(self):
        row = self.rows()[0]
        rows = []
        for i in range(12):
            rows.append({**row, "molecule_group": str(i), "split": "calibration" if i < 9 else "evaluation",
                         "absolute_relative_error": .1 if i != 11 else None})
        report = m.summarize(rows, CONFIG["statistics"])["properties"][row["property"]]
        self.assertAlmostEqual(report["evaluation_coverage_including_failures"], 2/3)
    def test_duplicate_predictions_raise(self):
        obs = a.audit_target(fixture(), CONFIG["audit"], FakeBackend())["observations"]
        one = {"observation_id": obs[0]["observation_id"], "predicted_value": 1.}
        with self.assertRaises(ValueError): m.comparison_rows(obs, [one, one], CONFIG["statistics"])


class RunnerTests(unittest.TestCase):
    def test_default_config_valid(self): runner.validate_config(CONFIG)
    def test_nonfinite_timeout_rejected(self):
        config = copy.deepcopy(CONFIG); config["worker_timeout_s"] = float("nan")
        with self.assertRaises(ValueError): runner.validate_config(config)
    def test_duplicate_target_ids_rejected(self):
        with self.assertRaises(ValueError): runner.load_targets({"version": "v5e-2.1", "targets": [fixture(), fixture()]})
    def test_checkpoint_corruption_rejected(self):
        with tempfile.TemporaryDirectory() as folder:
            path = Path(folder) / "check.json"
            runner.checkpoint(path, "signature", {"a": 1})
            saved = runner.load_json(path); saved["data"]["a"] = 2
            runner.write_json(path, saved)
            with self.assertRaises(ValueError): runner.checkpoint(path, "signature")
    def test_checkpoint_other_configuration_rejected(self):
        with tempfile.TemporaryDirectory() as folder:
            path = Path(folder) / "check.json"
            runner.checkpoint(path, "one", {"a": 1})
            with self.assertRaises(ValueError): runner.checkpoint(path, "two")
    def test_end_to_end_artifacts_and_resume(self):
        calls = []
        def predictor(request, folder, timeout):
            self.assertTrue((folder.parents[1] / "audited-reference.json").is_file())
            p.validate_request(request); calls.append(1)
            return {"results": [{"observation_id": r["observation_id"], "predicted_value": 100., "status": "ok"} for r in request["requests"]]}
        with tempfile.TemporaryDirectory() as folder, redirect_stdout(io.StringIO()):
            out = Path(folder)
            first = runner.execute([fixture()], CONFIG, out, "sig", FakeBackend(), predictor=predictor)
            second = runner.execute([fixture()], CONFIG, out, "sig", FakeBackend(), predictor=predictor)
            self.assertEqual(first, second)
            self.assertEqual(len(calls), 1)
            self.assertTrue(first["study_complete"])
            self.assertEqual(first["rankable_promotions"], 0)
            self.assertTrue((out / "property-comparison.csv").is_file())
    def test_audit_only_does_not_predict(self):
        with tempfile.TemporaryDirectory() as folder, redirect_stdout(io.StringIO()):
            result = runner.execute([fixture()], CONFIG, Path(folder), "sig", FakeBackend(), audit_only=True,
                                    predictor=lambda *args: self.fail("must not predict"))
            self.assertFalse(result["study_complete"])
            self.assertEqual(result["run_status"], "audit_only")
    def test_aliases_do_not_duplicate_reference_observations(self):
        alias = copy.deepcopy(fixture()); alias["candidate_id"] = "another-id"
        with tempfile.TemporaryDirectory() as folder, redirect_stdout(io.StringIO()):
            result = runner.execute([fixture(), alias], CONFIG, Path(folder), "sig", FakeBackend(), audit_only=True)
            self.assertEqual(result["targets_with_strict_reference_count"], 1)


class InstalledReferenceTests(unittest.TestCase):
    @unittest.skipUnless(importlib.util.find_spec("chemicals") and importlib.util.find_spec("thermo"), "chemicals/thermo not installed")
    def test_actual_ethanol_native_reference_points(self):
        backend = a.LocalReferenceBackend()
        target = fixture()
        identity = a.verify_identity(target, backend, CONFIG["audit"]["identity_mw_relative_tolerance"])
        self.assertEqual(identity["inchi_key"], KEY)
        t, hvap = backend.native_point("64-17-5", a.HV, "CRC_HVAP_TB", target["packet"])
        self.assertTrue(340 < t < 370 and hvap > 20000)
        t, cp = backend.native_point("64-17-5", a.CP, "POLING_CONST", target["packet"])
        self.assertEqual(t, 298.15)
        self.assertGreater(cp, 50)
    @unittest.skipUnless(importlib.util.find_spec("feos"), "FeOS not installed")
    def test_actual_feos_structure_only_ethanol(self):
        adapter = p.FeosBenchmarkAdapter("CCO")
        for prop in (a.TC, a.PC, a.PS, a.HV, a.RHO, a.CP):
            self.assertGreater(adapter.predict({"property": prop, "temperature_k": 300.0}), 0)


if __name__ == "__main__":
    unittest.main()
