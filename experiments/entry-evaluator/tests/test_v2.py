from __future__ import annotations

import copy
import math
from pathlib import Path
import sys
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from run_physical_sensitivity import (
    analyze,
    apply_values,
    latin_hypercube,
    validate_ranges,
)


class SamplingTests(unittest.TestCase):
    def test_latin_hypercube_is_deterministic_and_stratified(self):
        ranges = {"x": (0.0, 1.0), "y": (10.0, 20.0)}
        a = latin_hypercube(ranges, 8, 1234)
        b = latin_hypercube(ranges, 8, 1234)
        self.assertEqual(a, b)
        bins = sorted(int(row["x"] * 8) for row in a)
        self.assertEqual(bins, list(range(8)))
        self.assertTrue(all(10.0 <= row["y"] <= 20.0 for row in a))

    def test_invalid_ranges_rejected(self):
        for value in ([1], [2, 1], [1, math.inf], "bad"):
            with self.assertRaises(ValueError):
                validate_ranges({"x": value})

    def test_apply_values_does_not_mutate_base(self):
        base = {"a": {"b": 1}, "x": 2}
        original = copy.deepcopy(base)
        changed = apply_values(base, {"a.b": 7, "a.c": 9})
        self.assertEqual(base, original)
        self.assertEqual(changed["a"], {"b": 7, "c": 9})


class AnalysisTests(unittest.TestCase):
    def row(self, scenario, fluid, status, mass, **extra):
        return {
            "scenario_id": scenario,
            "fluid": fluid,
            "status": status,
            "coolant_used_kg": mass,
            **extra,
        }

    def test_analyze_counts_wins_losses_and_skips(self):
        scenarios = [
            {"scenario_id": "a", "values": {"x": 0.2}},
            {"scenario_id": "b", "values": {"x": 0.5}},
            {"scenario_id": "c", "values": {"x": 0.8}},
        ]
        rows = [
            self.row("a", "Hydrogen", "terminal_velocity", 50),
            self.row("a", "Water", "terminal_velocity", 100),
            self.row("b", "Hydrogen", "terminal_velocity", 120),
            self.row("b", "Water", "terminal_velocity", 100),
            self.row("c", "Hydrogen", "atmospheric_exit", 0),
            self.row("c", "Water", "atmospheric_exit", 0),
        ]
        report, per_scenario = analyze(
            rows, scenarios, ["h2", "water"], margin_fraction=0.01
        )
        self.assertEqual(report["h2_water"]["comparable_scenarios"], 2)
        self.assertEqual(report["h2_water"]["h2_better_count"], 1)
        self.assertEqual(report["h2_water"]["water_better_count"], 1)
        self.assertAlmostEqual(
            report["h2_water"]["h2_better_fraction_of_comparable"], 0.5
        )
        self.assertEqual(
            next(r for r in per_scenario if r["scenario_id"] == "c")[
                "h2_vs_water"
            ],
            "not_comparable",
        )

    def test_radiative_energy_metrics_are_summarized(self):
        scenarios = [
            {"scenario_id": "a", "values": {"x": 0.2}},
            {"scenario_id": "b", "values": {"x": 0.8}},
        ]
        rows = []
        for scenario, valid, fraction in (("a", 0.25, 0.7), ("b", 0.75, 0.9)):
            rows.extend(
                [
                    self.row(
                        scenario,
                        "Hydrogen",
                        "terminal_velocity",
                        50,
                    ),
                    self.row(
                        scenario,
                        "Water",
                        "terminal_velocity",
                        100,
                        vehicle_radiative_energy_valid_fraction=valid,
                        vehicle_radiative_energy_fraction=fraction,
                    ),
                ]
            )
        report, _ = analyze(rows, scenarios, ["h2", "water"], 0.01)
        q = report["radiation"][
            "vehicle_radiative_energy_valid_fraction_quantiles"
        ]
        self.assertAlmostEqual(q["median"], 0.5)
        self.assertAlmostEqual(
            report["radiation"][
                "fraction_with_at_least_50pct_radiative_energy_in_nominal_range"
            ],
            0.5,
        )
        self.assertEqual(report["physical_validation"], "not_established")


if __name__ == "__main__":
    unittest.main()
