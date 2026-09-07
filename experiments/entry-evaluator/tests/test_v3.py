from __future__ import annotations

import math
from pathlib import Path
import sys
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from atmosphere import AtmosphereState
from heating import (
    brandis_johnston_convective,
    brandis_johnston_radiative,
    total_heating,
)
from run_heating_model_comparison import analyze


def atmosphere(rho=5.0e-4, continuum=1.0):
    return AtmosphereState(
        altitude_m=50000.0,
        density_kg_m3=rho,
        temperature_k=250.0,
        pressure_pa=80.0,
        number_density_m3=1.0,
        mean_free_path_m=1.0e-5,
        knudsen=1.0e-5,
        continuum_factor=continuum,
        n2_m3=0.0,
        o2_m3=0.0,
        o_m3=0.0,
    )


class BrandisJohnstonTests(unittest.TestCase):
    def test_low_speed_convective_reference_value_and_units(self):
        q, valid = brandis_johnston_convective(atmosphere(), 9000.0, 0.5)
        self.assertTrue(valid)
        self.assertAlmostEqual(q, 4.903472502602395e6, delta=1e-6 * q)

    def test_high_speed_convective_reference_value_and_units(self):
        q, valid = brandis_johnston_convective(atmosphere(), 11000.0, 0.5)
        self.assertTrue(valid)
        self.assertAlmostEqual(q, 8.252634746140069e6, delta=1e-6 * q)

    def test_radiative_reference_value_and_velocity_domain(self):
        q, valid = brandis_johnston_radiative(atmosphere(), 11000.0, 0.5)
        self.assertTrue(valid)
        self.assertAlmostEqual(q, 3.6953436124299695e6, delta=1e-6 * q)
        q_low, valid_low = brandis_johnston_radiative(atmosphere(), 9000.0, 0.5)
        self.assertEqual(q_low, 0.0)
        self.assertFalse(valid_low)

    def test_validity_is_separate_from_algebraic_evaluation(self):
        q, valid = brandis_johnston_convective(atmosphere(rho=8.0e-3), 11000.0, 1.0)
        self.assertGreater(q, 0.0)
        self.assertFalse(valid)

    def test_backend_switch_is_explicit(self):
        bj = total_heating(
            atmosphere(), 11000.0, 0.5, {"backend": "brandis_johnston_2014"}
        )
        legacy = total_heating(atmosphere(), 11000.0, 0.5, {"backend": "legacy"})
        self.assertEqual(bj.backend, "brandis_johnston_2014")
        self.assertEqual(legacy.backend, "legacy")
        self.assertNotAlmostEqual(bj.total_external_w_m2, legacy.total_external_w_m2)
        with self.assertRaises(ValueError):
            total_heating(atmosphere(), 11000.0, 0.5, {"backend": "unknown"})


class ComparisonAnalysisTests(unittest.TestCase):
    def row(self, scenario, fluid, backend, mass, status="terminal_velocity", **extra):
        return {
            "scenario_id": scenario,
            "fluid_expected": fluid,
            "fluid": fluid,
            "heating_backend": backend,
            "status": status,
            "coolant_used_kg": mass,
            **extra,
        }

    def test_analysis_preserves_h2_win_and_reports_model_shift(self):
        rows = [
            self.row("s0", "Hydrogen", "legacy", 40.0),
            self.row("s0", "Water", "legacy", 100.0),
            self.row(
                "s0", "Hydrogen", "brandis_johnston_2014", 45.0,
                vehicle_radiative_energy_valid_fraction=0.8,
            ),
            self.row(
                "s0", "Water", "brandis_johnston_2014", 90.0,
                vehicle_radiative_energy_valid_fraction=0.8,
            ),
        ]
        report, paired, rankings = analyze(
            rows,
            scenario_ids=["s0"],
            fluids=["Hydrogen", "Water"],
            backends=["legacy", "brandis_johnston_2014"],
        )
        self.assertTrue(report["study_complete"])
        self.assertEqual(report["h2_water"]["common_comparable_scenarios"], 1)
        self.assertEqual(report["h2_water"]["h2_better_both_count"], 1)
        self.assertAlmostEqual(
            report["h2_water"]["ratio_shift_bj_minus_legacy"]["median"], 0.1
        )
        self.assertEqual(len(paired), 2)
        self.assertEqual(len(rankings), 2)


if __name__ == "__main__":
    unittest.main()
