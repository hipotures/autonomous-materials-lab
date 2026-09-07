from __future__ import annotations

import copy
import math
from pathlib import Path
import sys
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import CoolProp.CoolProp as CP
import yaml
import json
import numpy as np
from chemistry import ChemistryLimiter
from coolant import CoolantModel, CoolantStep
from evaluator import evaluate
from heating import HeatingState
from run_convergence import METRICS, analyze, compare, pair_order
from surface import Forebody, build_zones
from wall import WallState, step_wall


class ConstantEnthalpyCoolant:
    def evaluate(self, wall_temperature_k, pressure, q, *, area_m2):
        flux = q / 1000.0
        return CoolantStep(300, 1000, flux, flux * area_m2, pressure,
                           "disabled", None, True, None)


class SurfaceTests(unittest.TestCase):
    def test_area_and_analytic_integral(self):
        for theta in (30, 60, 90):
            mu = math.cos(math.radians(theta))
            for exponent in (0, 0.5, 1, 2):
                expected = (1 - mu ** (exponent + 1)) / ((exponent + 1) * (1 - mu))
                for count in (1, 8, 64):
                    zones = build_zones({"model": "cosine_cap", "theta_max_deg": theta,
                                         "zones": count, "convective_exponent": exponent}, 12)
                    self.assertAlmostEqual(sum(z.area_m2 for z in zones), 12)
                    self.assertAlmostEqual(sum(z.area_m2 * z.convective_factor for z in zones), 12 * expected)

    def test_invalid_geometry(self):
        for cfg in ({"model": "unknown"}, {"model": "cosine_cap", "zones": 0},
                    {"model": "cosine_cap", "zones": 1.5}, {"model": "cosine_cap", "zones": True},
                    {"model": "cosine_cap", "theta_max_deg": 91},
                    {"model": "cosine_cap", "convective_exponent": float("nan")}):
            with self.assertRaises(ValueError):
                build_zones(cfg, 12)
        for area in (0, -1, math.nan, math.inf):
            with self.assertRaises(ValueError):
                build_zones({}, area)

    def wall_config(self):
        return {"initial_temperature_k": 300, "temperature_setpoint_k": 400,
                "areal_heat_capacity_j_m2_k": 100, "emissivity": 0,
                "backface_loss_w_m2": 0}

    def test_uniform_matches_single_wall_and_probe_does_not_add_mass(self):
        cfg = self.wall_config()
        surface = Forebody({"model": "uniform"}, 12, cfg)
        wall = WallState(300)
        for q in (5000, 20000, 30000, 0):
            old = step_wall(wall, q, 300, 0.5, cfg)
            result = surface.evaluate(HeatingState(q, 0, q, 1, False), 300, 1e5, 0.5, ConstantEnthalpyCoolant())
            self.assertAlmostEqual(result.mass_flow_kg_s, old.coolant_required_w_m2 * 12 / 1000)
            self.assertAlmostEqual(result.stagnation_wall.next_temperature_k, old.next_temperature_k)
            self.assertLess(abs(result.energy_residual_w), 1e-8)
            surface.commit(result)
            wall.temperature_k = old.next_temperature_k

    def test_ring_storage_not_a_global_area_multiplier(self):
        surface = Forebody({"model": "cosine_cap", "zones": 4}, 12, self.wall_config())
        # q(mu)=20000*mu, C*dT/dt=10000. Integral max(q-10000,0) dmu = 2500.
        result = surface.evaluate(HeatingState(20000, 0, 20000, 1, False), 300, 1e5, 1, ConstantEnthalpyCoolant())
        self.assertAlmostEqual(result.incident_power_w, 120000)
        self.assertAlmostEqual(result.coolant_power_w, 30000)
        self.assertAlmostEqual(result.mass_flow_kg_s, 30)
        self.assertAlmostEqual(result.peak_temperature_k, 400)
        self.assertAlmostEqual(result.energy_residual_w, 0, places=6)
        self.assertLess(result.zone_walls[0].next_temperature_k, 400)
        # Area-averaged heating applied to a single wall would request zero cooling.
        average = step_wall(WallState(300), 10000, 300, 1, self.wall_config())
        self.assertEqual(average.coolant_required_w_m2, 0)

    def test_convective_and_radiative_profiles_are_independent(self):
        s = Forebody({"model": "cosine_cap", "zones": 8,
                      "convective_exponent": 1, "radiative_exponent": 0}, 12, self.wall_config())
        r = s.evaluate(HeatingState(20000, 10000, 30000, 1, True), 300, 1e5, 1, ConstantEnthalpyCoolant())
        self.assertAlmostEqual(r.incident_power_w, 12 * (20000 / 2 + 10000))


class CoolantTests(unittest.TestCase):
    def test_exact_cache_is_order_independent_and_matches_propssi(self):
        cfg = {"coolprop_name": "Water", "storage_temperature_k": 293.15,
               "storage_pressure_pa": 101325, "max_exit_temperature_k": 1000, "cooled_area_m2": 12}
        a = CoolantModel(cfg, ChemistryLimiter({"mode": "disabled"}))
        b = CoolantModel(cfg, ChemistryLimiter({"mode": "disabled"}))
        points = [(749.96, 123400), (749.99, 123401)]  # Old rounded cache aliased these.
        forward = [a._h_out(*p) for p in points]
        reverse = [b._h_out(*p) for p in reversed(points)][::-1]
        self.assertEqual(forward, reverse)
        self.assertNotEqual(forward[0], forward[1])
        for (temperature, pressure), value in zip(points, forward):
            self.assertAlmostEqual(value, CP.PropsSI("H", "T", temperature, "P", pressure, "Water"), places=6)


class ConvergenceTests(unittest.TestCase):
    def row(self, **kwargs):
        row = {key: 100.0 for key in METRICS}
        row.update(status="terminal_velocity")
        row.update(kwargs)
        return row

    def test_compare_rejects_changed_status_nan_and_large_errors(self):
        self.assertTrue(compare(self.row(), self.row())["pass"])
        for changes in ({"status": "failed"}, {"status": "atmospheric_exit"},
                        {"coolant_used_kg": math.nan}, {"coolant_used_kg": 120},
                        {"peak_wall_temperature_k": 101}):
            self.assertFalse(compare(self.row(**changes), self.row())["pass"])
        self.assertFalse(compare(self.row(status="atmospheric_exit", coolant_used_kg=0),
                                 self.row(status="atmospheric_exit", coolant_used_kg=0))["pass"])

    def test_comparison_with_numpy_observations_is_json_serializable(self):
        result = compare(self.row(minimum_altitude_km=np.float64(100)), self.row())
        json.dumps(result, allow_nan=False)

    def test_ranking_zero_failure_and_ties_are_not_stable_wins(self):
        self.assertEqual(pair_order(self.row(coolant_used_kg=0), self.row()), "unavailable")
        self.assertEqual(pair_order(self.row(), self.row(coolant_used_kg=100.1)), "unresolved")
        self.assertEqual(pair_order(self.row(coolant_used_kg=50), self.row()), "a_less")

    def test_ranking_reversal_is_reported(self):
        rows = []
        for profile in ("broad", "narrow"):
            for case in ("a", "b"):
                mass = 50 if (profile == "broad") == (case == "a") else 100
                for dt, kind in ((0.1, "time"), (0.05, "time"), (0.05, "refined_grid")):
                    rows.append(self.row(case=case, profile=profile, dt_s=dt, kind=kind, coolant_used_kg=mass))
        report = analyze(rows, ["a", "b"], ["broad", "narrow"], [0.1, 0.05])
        self.assertTrue(report["numerical_pass"])
        self.assertFalse(report["ranking_stable"])
        self.assertFalse(report["gate_pass"])


class EvaluatorTests(unittest.TestCase):
    def test_skip_is_not_a_zero_mass_cooling_success(self):
        cfg = yaml.safe_load((Path(__file__).resolve().parents[1] / "config.yaml").read_text())
        cfg["entry"]["flight_path_angle_deg"] = -8
        cfg["numerics"]["dt_s"] = 0.1
        r = evaluate(cfg)
        self.assertEqual(r["status"], "atmospheric_exit")
        self.assertEqual(r["coolant_used_kg"], 0)
        self.assertGreater(r["rarefied_heating_disabled_time_s"], 0)
        self.assertFalse(compare(r, r)["pass"])


if __name__ == "__main__":
    unittest.main()
