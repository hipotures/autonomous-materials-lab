from __future__ import annotations

import csv
from pathlib import Path
import sys
import tempfile
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from chemistry import ChemistryLimiter
from run_chemistry_sensitivity import analyze, coverage_preflight


def write_table(path: Path) -> None:
    fields = [
        "temperature_initial_K",
        "pressure_bar",
        "phi",
        "max_time_s",
        "tau_dTdt_s",
        "status",
    ]
    rows = []
    for pressure, scale in (
        (1.0, 1.0),
        (3.0, 0.5),
    ):
        for temperature, tau in (
            (700.0, 1.0),
            (800.0, 0.1),
            (900.0, 0.01),
        ):
            rows.append(
                {
                    "temperature_initial_K": temperature,
                    "pressure_bar": pressure,
                    "phi": 4.0,
                    "max_time_s": 10.0,
                    "tau_dTdt_s": tau * scale,
                    "status": "ignited",
                }
            )

    with path.open(
        "w",
        newline="",
        encoding="utf-8",
    ) as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=fields,
        )
        writer.writeheader()
        writer.writerows(rows)


class ChemistryLimiterTests(unittest.TestCase):
    def make_limiter(
        self,
        path: Path,
        **overrides,
    ) -> ChemistryLimiter:
        config = {
            "mode": "ignition_csv",
            "ignition_csv": str(path),
            "phi": 4.0,
            "residence_time_s": 0.01,
            "ignition_safety_factor": 10.0,
        }
        config.update(overrides)
        return ChemistryLimiter(config)

    def test_temperature_interpolation_finds_threshold(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "table.csv"
            write_table(path)
            limiter = self.make_limiter(path)
            result = limiter.limit(1e5, 900.0)

            self.assertTrue(result.feasible)
            self.assertAlmostEqual(
                result.max_exit_temperature_k,
                800.0,
                places=8,
            )
            self.assertAlmostEqual(
                result.estimated_ignition_delay_s,
                0.1,
                places=8,
            )

    def test_log_pressure_interpolation(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "table.csv"
            write_table(path)
            limiter = self.make_limiter(path)

            pressure_pa = (3.0 ** 0.5) * 1e5
            tau = limiter.estimate_delay(
                pressure_pa,
                800.0,
            )

            self.assertIsNotNone(tau)
            self.assertAlmostEqual(
                tau,
                (0.1 * 0.05) ** 0.5,
                places=10,
            )

    def test_pressure_outside_table_is_not_clamped(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "table.csv"
            write_table(path)
            limiter = self.make_limiter(path)

            result = limiter.limit(
                5e5,
                900.0,
            )

            self.assertFalse(result.feasible)
            self.assertIn(
                "outside ignition table range",
                result.failure_reason,
            )

    def test_missing_phi_is_rejected(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "table.csv"
            write_table(path)
            limiter = self.make_limiter(
                path,
                phi=8.0,
            )

            result = limiter.limit(
                1e5,
                900.0,
            )

            self.assertFalse(result.feasible)
            self.assertIn(
                "not present",
                result.failure_reason,
            )


class V4AnalysisTests(unittest.TestCase):
    def test_analysis_detects_water_better_case(self):
        chemistry_cases = [
            {
                "residence_time_s": 0.001,
                "phi": 4.0,
                "safety_factor": 10.0,
            }
        ]
        rows = [
            {
                "scenario_id": "nominal",
                "role": "water",
                "status": "terminal_velocity",
                "coolant_used_kg": 100.0,
            },
            {
                "scenario_id": "nominal",
                "role": "hydrogen",
                "chemistry_id": (
                    "res-0p001_phi-4_safety-10"
                ),
                "status": "terminal_velocity",
                "coolant_used_kg": 120.0,
                "minimum_ignition_margin": 1.2,
                "maximum_coolant_exit_temperature_k": 760.0,
            },
        ]

        report, detail = analyze(
            rows,
            ["nominal"],
            chemistry_cases,
        )

        self.assertTrue(report["study_complete"])
        self.assertEqual(
            report["water_better_somewhere_case_count"],
            1,
        )
        self.assertEqual(
            detail[0]["outcome"],
            "water_better",
        )


class CoveragePreflightTests(unittest.TestCase):
    def test_preflight_rejects_pressure_gap_before_h2_sweep(self):
        rows = [
            {
                "status": "terminal_velocity",
                "minimum_active_cooling_surface_pressure_bar": 0.01,
                "maximum_active_cooling_surface_pressure_bar": 12.15,
            }
        ]
        meta = {
            "pressure_min_bar": 0.001,
            "pressure_max_bar": 12.0,
            "phis": [2.0, 4.0, 8.0, 16.0],
        }
        result = coverage_preflight(rows, meta, [2.0, 4.0, 8.0, 16.0])
        self.assertFalse(result["pass"])
        self.assertAlmostEqual(
            result["required_active_cooling_pressure_max_bar"],
            12.15,
        )

    def test_preflight_accepts_extended_pressure_table(self):
        rows = [
            {
                "status": "terminal_velocity",
                "minimum_active_cooling_surface_pressure_bar": 0.01,
                "maximum_active_cooling_surface_pressure_bar": 12.15,
            }
        ]
        meta = {
            "pressure_min_bar": 0.001,
            "pressure_max_bar": 30.0,
            "phis": [2.0, 4.0, 8.0, 16.0],
        }
        result = coverage_preflight(rows, meta, [2.0, 4.0, 8.0, 16.0])
        self.assertTrue(result["pass"])


if __name__ == "__main__":
    unittest.main()
