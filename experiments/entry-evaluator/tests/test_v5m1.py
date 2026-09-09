"""V5m-1 regression tests; fixtures below are synthetic, not physical evidence."""
from __future__ import annotations

import copy
import importlib.util
import io
import json
from pathlib import Path
import re
import subprocess
import sys
import tempfile
from contextlib import redirect_stdout
from types import SimpleNamespace
import unittest
from unittest.mock import patch

import yaml
HERE = Path(__file__).resolve().parents[2] / "formulation-screen-v5m"
sys.path.insert(0, str(HERE))
import formulations as f
import mixture_models as models
import run_formulations as runner
import analyze_bench as bench
CONFIG = yaml.safe_load((HERE / "config.yaml").read_text())


class FakeCP:
    """Synthetic backend with deliberately different enthalpy reference offsets."""
    def __init__(self): self.calls = []
    def PhaseSI(self, *args): return "liquid" if args[1] < 370 else "gas"
    def PropsSI(self, output, _, temperature, __, pressure, fluid):
        self.calls.append((output, temperature, pressure, fluid))
        w = float(re.search(r"\[([^]]+)\]", fluid)[1]) if "[" in fluid else 0.0
        cp = 4200 * (1 - 0.5 * w)
        return {"D": 1000 * (1 + w), "C": cp, "V": 0.001 * (1 + 5*w),
                "L": 0.6 * (1 - w), "H": -2e6 + 1e6*w + cp * (temperature - 293.15)}[output]


class FakeNRTL:
    def __init__(self, **kwargs):
        assert kwargs["composition_basis"] == "mass"
        self.kwargs = kwargs
    def enthalpy_j_kg(self, t, pressure): return -7e5 + 3800*(t - 293.15)
    def phase(self, t, pressure): return "gas" if t > 375 else "liquid"
    def saturation_at_pressure(self, p):
        return SimpleNamespace(bubble_temperature_k=370., dew_temperature_k=375., supported=True)


def compact_config():
    config = copy.deepcopy(CONFIG)
    config["additives"] = [copy.deepcopy(CONFIG["additives"][0])]
    config["additives"][0]["mass_fractions"] = [0.01]
    return config


def bench_fixture():
    config = compact_config()
    plan = {"version": "v5m-1", "runs": f.make_bench_plan(f.build_formulations(config), config)}
    protocol = f.protocol_template()
    protocol.update({k: k + "-test" for k in bench.MATCH_KEYS})
    protocol.update(wall_temperature_limit_k=400., max_feed_pressure_pa=200000.,
        target_incident_energy_j=10000., target_duration_s=100., max_hydraulic_resistance_ratio=1.2,
        max_residue_mg=1.)
    measured = []
    for r in plan["runs"]:
        measured.append({**{k: r[k] for k in ("run_id", "formulation_id", "block_id", "role", "replicate")},
            **{k: protocol[k] for k in bench.MATCH_KEYS}, **{k: "true" for k in bench.PASS_FIELDS},
            "evidence_kind": "measured", "coupon_id": r["run_id"], "preparation_record_id": "prep-test",
            "source_data_id": "synthetic-unit-test-fixture", "safety_review_id": "review-test",
            "mass_used_g": 8. if r["role"] == "candidate" else 10., "mass_uncertainty_g": .1,
            "incident_energy_j": 10000., "energy_uncertainty_j": 10., "duration_s": 100.,
            "wall_peak_k": 380., "wall_uncertainty_k": 1., "feed_pressure_max_pa": 150000.,
            "pressure_uncertainty_pa": 1000., "hydraulic_resistance_ratio": 1.05, "residue_mg": .1})
    return plan, measured, protocol


class FormulationTests(unittest.TestCase):
    def test_default_catalog_count(self):
        self.assertEqual(len(f.build_formulations(CONFIG)), 26)
    def test_one_water_reference_only(self):
        self.assertEqual(sum(r["formulation_id"] == "water" for r in f.build_formulations(CONFIG)), 1)
    def test_recipes_close_mass_balance(self):
        for r in f.build_formulations(CONFIG):
            self.assertAlmostEqual(r["additive_mass_g"] + r["water_mass_g"], 100.)
    def test_one_percent_is_one_gram_per_hundred_grams(self):
        r = f.build_formulations(compact_config())[1]
        self.assertAlmostEqual(r["additive_mass_g"], 1.)
    def test_nonfinite_fraction_rejected(self):
        c = compact_config(); c["additives"][0]["mass_fractions"] = [float("nan")]
        with self.assertRaises(ValueError): f.build_formulations(c)
    def test_boolean_fraction_rejected(self):
        c = compact_config(); c["additives"][0]["mass_fractions"] = [True]
        with self.assertRaises(ValueError): f.build_formulations(c)
    def test_backend_cas_mismatch_rejected(self):
        c = compact_config(); c["additives"][0]["cas"] = "67-56-1"
        with self.assertRaises(ValueError): f.build_formulations(c)
    def test_no_implicit_suspension_enable(self):
        c = copy.deepcopy(CONFIG); c["additives"][-1].update(enabled=True, mass_fractions=[.001])
        with self.assertRaises(ValueError): f.build_formulations(c)
    def test_retention_is_mass_balance_not_solubility(self):
        r = f.retained_additive_stress(f.build_formulations(compact_config())[1], .9)
        self.assertAlmostEqual(r["retained_additive_mass_fraction"], 1/10.9)
        self.assertFalse(r["precipitation_predicted"])
    def test_retention_rejects_complete_water_loss(self):
        with self.assertRaises(ValueError): f.retained_additive_stress(f.build_formulations(CONFIG)[1], 1.)
    def test_break_even_not_invented_enhancement(self):
        r = f.heat_only_comparison(80., 100.)
        self.assertAlmostEqual(r["mass_ratio_equal_net_heat"], 1.25)
        self.assertAlmostEqual(r["additional_net_heat_reduction_needed_for_mass_parity"], .2)
        self.assertFalse(r["system_winner"])
    def test_equal_mass_and_volume_flux_are_distinct(self):
        row = {"formulation_id": "x", "temperature_k": 300., "pressure_pa": 1e5,
               "density_kg_m3": 1200., "viscosity_pa_s": .002}
        r = f.cold_ratios(row, {"density_kg_m3": 1000., "viscosity_pa_s": .001})
        self.assertAlmostEqual(r["pressure_drop_ratio_equal_volume_flux"], 2.)
        self.assertAlmostEqual(r["pressure_drop_ratio_equal_mass_flux"], 2/1.2)
    def test_bracketed_deterministic_plan(self):
        formulations = f.build_formulations(CONFIG)
        plan = f.make_bench_plan(formulations, CONFIG)
        self.assertEqual(plan, f.make_bench_plan(list(reversed(formulations)), CONFIG))
        for block in {r["block_id"] for r in plan}:
            rows = [r for r in plan if r["block_id"] == block]
            self.assertEqual(rows[0]["role"], "water_before")
            self.assertEqual(rows[-1]["role"], "water_after")
    def test_no_physical_execution_authorized(self):
        self.assertTrue(all(not r["physical_execution_authorized"] for r in f.build_formulations(CONFIG)))


class ModelsTests(unittest.TestCase):
    def setUp(self):
        self.cp = FakeCP(); self.m = models.MixtureModels(cp=self.cp, nrtl_factory=FakeNRTL)
        self.fs = f.build_formulations(CONFIG)
    def get(self, additive): return next(r for r in self.fs if r["additive_id"] == additive)
    def test_water_uses_true_pure_backend(self):
        r = self.m.cold(self.fs[0], 293.15, 101325.)
        self.assertEqual(r["backend"], "HEOS::Water")
    def test_explicit_incomp_mass_fraction_not_percent(self):
        r = self.m.cold(self.get("propylene_glycol"), 300., 101325.)
        self.assertEqual(r["backend"], "INCOMP::MPG[0.001]")
    def test_zero_fit_endpoint_is_explicit_not_a_recipe(self):
        probe = {**self.get("glycerol"), "additive_mass_fraction": 0.0}
        self.assertEqual(self.m.cold(probe, 300., 101325.)["status"], "outside_model_range")
        r = self.m.cold(probe, 300., 101325., fit_endpoint_probe=True)
        self.assertEqual(r["backend"], "INCOMP::MGL[0]")
        self.assertTrue(r["fit_endpoint_probe"])
    def test_glycerol_high_temperature_not_extrapolated(self):
        r = self.m.cold(self.get("glycerol"), 350., 101325.)
        self.assertEqual(r["status"], "outside_model_range"); self.assertFalse(self.cp.calls)
    def test_salt_concentration_range(self):
        r = self.m.cold({**self.get("sodium_chloride"), "additive_mass_fraction": .3}, 300., 101325.)
        self.assertEqual(r["status"], "outside_model_range")
    def test_unknown_wetting_is_null(self):
        self.assertIsNone(self.m.cold(self.get("ethanol"), 300., 101325.)["surface_tension_n_m"])
    def test_sds_is_measurement_only(self):
        self.assertEqual(self.m.cold(self.get("sds"), 300., 101325.)["status"], "measurement_required")
    def test_negative_absolute_enthalpy_is_valid(self):
        r = self.m.cold(self.get("glycerol"), 300., 101325.)
        self.assertLess(r["enthalpy_j_kg"], 0.); self.assertEqual(r["status"], "ok")
    def test_no_boiling_fallback_for_glycol(self):
        r = self.m.hot(self.get("propylene_glycol"), 293.15, 500., 101325.)
        self.assertEqual(r["status"], "unsupported_phase_change_model")
        self.assertIsNone(r["delta_h_j_kg"])
    def test_nrtl_delta_h_uses_same_backend_both_ends(self):
        r = self.m.hot(self.get("ethanol"), 293.15, 500., 101325.)
        self.assertAlmostEqual(r["delta_h_j_kg"], 3800*(500.-293.15))
        self.assertFalse(self.cp.calls)
    def test_backend_failure_not_zero_or_pure_water(self):
        self.cp.PropsSI = lambda *args: (_ for _ in ()).throw(ValueError("bad fit"))
        r = self.m.cold(self.get("glycerol"), 300., 101325.)
        self.assertIsNone(r["density_kg_m3"]); self.assertEqual(r["status"], "property_failed")


class RunnerTests(unittest.TestCase):
    def test_checkpoint_corruption_rejected(self):
        with tempfile.TemporaryDirectory() as d:
            p = Path(d)/"checkpoint.json"; runner.checkpoint(p, "a", {"x": 1})
            r = f.read_json(p); r["data"]["x"] = 2; f.write_json(p, r)
            with self.assertRaises(ValueError): runner.checkpoint(p, "a")
    def test_signature_change_rejected(self):
        with tempfile.TemporaryDirectory() as d:
            p = Path(d)/"checkpoint.json"; runner.checkpoint(p, "a", {"x": 1})
            with self.assertRaises(ValueError): runner.checkpoint(p, "b")
    def test_full_synthetic_screen_outputs_and_resume(self):
        with tempfile.TemporaryDirectory() as d, redirect_stdout(io.StringIO()):
            output = Path(d); m = models.MixtureModels(cp=FakeCP(), nrtl_factory=FakeNRTL)
            r = runner.execute(compact_config(), output, "test", m)
            self.assertEqual(r["unexpected_model_failure_count"], 0)
            self.assertEqual(r["experimental_winners"], 0)
            (output/"measurements-template.csv").write_text("do not replace")
            with patch.object(m, "cold", side_effect=AssertionError("recomputed")):
                runner.execute(compact_config(), output, "test", m)
            self.assertEqual((output/"measurements-template.csv").read_text(), "do not replace")
            self.assertTrue((output/"enthalpy-window.csv").is_file())
            self.assertTrue((output/"model-endpoint-comparison.csv").is_file())
    def test_plan_only_command_and_resume(self):
        with tempfile.TemporaryDirectory() as d:
            cmd = [sys.executable, str(HERE/"run_formulations.py"), "--plan-only", "--output-dir", d]
            self.assertEqual(subprocess.run(cmd, capture_output=True).returncode, 0)
            self.assertEqual(subprocess.run(cmd+["--resume"], capture_output=True).returncode, 0)
            self.assertNotEqual(subprocess.run(cmd, capture_output=True).returncode, 0)


class BenchTests(unittest.TestCase):
    def setUp(self): self.plan, self.data, self.p = bench_fixture()
    def report(self): return bench.analyze(self.plan, self.data, self.p)
    def candidate(self): return next(r for r in self.data if r["role"] == "candidate")
    def test_three_matched_replicates_can_be_promising(self):
        self.assertEqual(self.report()["bench_promising_count"], 1)
    def test_model_data_never_becomes_physical_evidence(self):
        self.candidate()["evidence_kind"] = "model"
        self.assertEqual(self.report()["bench_promising_count"], 0)
    def test_synthetic_data_label_rejected(self):
        for r in self.data: r["evidence_kind"] = "synthetic"
        self.assertEqual(self.report()["valid_candidate_comparison_count"], 0)
    def test_missing_water_control_blocks_comparison(self):
        self.data = [r for r in self.data if r["role"] != "water_after"]
        self.assertEqual(self.report()["bench_promising_count"], 0)
    def test_no_optional_stopping_or_dropped_failed_replicate(self):
        self.candidate()["completed"] = "false"
        self.assertEqual(self.report()["bench_promising_count"], 0)
    def test_uncertainty_can_erase_improvement(self):
        for r in self.data:
            r["mass_uncertainty_g"] = 2.0
        self.assertEqual(self.report()["bench_promising_count"], 0)
    def test_peak_temperature_bound_checked(self):
        self.candidate().update(wall_peak_k=399., wall_uncertainty_k=2.)
        self.assertEqual(self.report()["bench_promising_count"], 0)
    def test_pressure_bound_checked(self):
        self.candidate().update(feed_pressure_max_pa=199500., pressure_uncertainty_pa=1000.)
        self.assertEqual(self.report()["bench_promising_count"], 0)
    def test_low_exposure_cannot_fake_mass_saving(self):
        self.candidate()["incident_energy_j"] = 8000.
        self.assertEqual(self.report()["bench_promising_count"], 0)
    def test_residue_is_failure_not_bonus(self):
        self.candidate()["residue_mg"] = 5.
        self.assertEqual(self.report()["bench_promising_count"], 0)
    def test_clogging_rejects_apparent_thermal_gain(self):
        self.candidate()["hydraulic_resistance_ratio"] = 2.
        self.assertEqual(self.report()["bench_promising_count"], 0)
    def test_reused_coupon_confounds_experiment(self):
        self.data[1]["coupon_id"] = self.data[0]["coupon_id"]
        self.assertEqual(self.report()["bench_promising_count"], 0)
    def test_dirty_fluid_path_rejected(self):
        self.candidate()["fluid_path_clean"] = "false"
        self.assertEqual(self.report()["bench_promising_count"], 0)
    def test_unknown_hardware_blocks_analysis(self):
        self.p["geometry_id"] = None
        with self.assertRaises(ValueError): self.report()
    def test_nan_measurement_rejected(self):
        self.candidate()["mass_used_g"] = "nan"
        self.assertEqual(self.report()["bench_promising_count"], 0)
    def test_duplicate_run_id_rejected(self):
        self.data.append(self.data[0])
        with self.assertRaises(ValueError): self.report()
    def test_control_drift_blocks_winner(self):
        next(r for r in self.data if r["role"] == "water_after")["mass_used_g"] = 15.
        self.assertEqual(self.report()["bench_promising_count"], 0)
    def test_partial_submitted_data_does_not_complete_experiment(self):
        self.data = self.data[:3]
        self.assertEqual(self.report()["bench_promising_count"], 0)
    def test_zero_positive_property_rejected(self):
        self.candidate()["incident_energy_j"] = 0.
        self.assertEqual(self.report()["bench_promising_count"], 0)
    def test_null_provenance_id_is_not_a_filled_field(self):
        self.candidate()["source_data_id"] = None
        self.assertEqual(self.report()["bench_promising_count"], 0)
    def test_template_is_not_physical_evidence(self):
        self.data = [{"run_id": r["run_id"]} for r in self.plan["runs"]]
        self.assertEqual(self.report()["bench_promising_count"], 0)


class InstalledIntegrationTests(unittest.TestCase):
    def test_real_coolprop_cold_models(self):
        if importlib.util.find_spec("CoolProp") is None:
            self.skipTest("CoolProp not installed")
        m = models.MixtureModels()
        for row in f.build_formulations(CONFIG):
            if row["additive_id"] == "sds": continue
            for t in CONFIG["cold_temperatures_k"]:
                state = m.cold(row, t, CONFIG["pressure_pa"])
                self.assertEqual(state["status"], "ok", state)
    def test_existing_nrtl_same_backend_window(self):
        if any(importlib.util.find_spec(p) is None for p in ("CoolProp", "thermo", "chemicals")):
            self.skipTest("CoolProp/thermo/chemicals not installed")
        m = models.MixtureModels()
        row = f.build_formulations(compact_config())[1]
        result = m.hot(row, 293.15, 400., 101325.)
        self.assertEqual(result["status"], "ok", result)
        self.assertGreater(result["delta_h_j_kg"], 0.)


if __name__ == "__main__":
    unittest.main()
