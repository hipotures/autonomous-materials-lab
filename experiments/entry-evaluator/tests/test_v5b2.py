from __future__ import annotations

import importlib.util
from pathlib import Path
import sys
import unittest

import yaml

ENTRY = Path(__file__).resolve().parents[1]
V5B = ENTRY.parent / "property-predictor-v5b"
sys.path.insert(0, str(V5B))

HAS_RDKIT = importlib.util.find_spec("rdkit") is not None

from failure_taxonomy import SUCCESS, expected_prediction_outcome  # noqa: E402
from applicability import (  # noqa: E402
    assess_domain,
    calibrate_metric,
    conformal_quantile,
    deterministic_split,
    local_error_scale,
    predict_relative_uncertainty,
)


class V5b2CalibrationMathTests(unittest.TestCase):
    def test_deterministic_split_is_stable_and_disjoint(self):
        ids = [f"c{i}" for i in range(20)]
        a_cal, a_eval = deterministic_split(
            ids,
            calibration_fraction=0.7,
            seed="fixed",
        )
        b_cal, b_eval = deterministic_split(
            reversed(ids),
            calibration_fraction=0.7,
            seed="fixed",
        )
        self.assertEqual(a_cal, b_cal)
        self.assertEqual(a_eval, b_eval)
        self.assertFalse(set(a_cal) & set(a_eval))
        self.assertEqual(set(a_cal) | set(a_eval), set(ids))
        self.assertEqual(len(a_cal), 14)
        self.assertEqual(len(a_eval), 6)

    def test_conformal_quantile_uses_finite_sample_order_statistic(self):
        self.assertEqual(
            conformal_quantile(
                [0.1, 0.2, 0.3, 0.4],
                coverage=0.8,
            ),
            0.4,
        )
        self.assertIsNone(
            conformal_quantile([], coverage=0.9)
        )

    def test_benchmark_has_broad_holdouts_and_domain_probes(self):
        study = yaml.safe_load(
            (V5B / "benchmark-v5b2.yaml").read_text(
                encoding="utf-8"
            )
        )
        holdouts = study["holdouts"]
        self.assertGreaterEqual(len(holdouts), 20)
        supported = [
            row
            for row in holdouts
            if expected_prediction_outcome(row) == SUCCESS
        ]
        probes = [
            row
            for row in holdouts
            if expected_prediction_outcome(row) != SUCCESS
        ]
        self.assertGreaterEqual(len(supported), 15)
        self.assertGreaterEqual(len(probes), 3)
        self.assertEqual(
            len({row["id"] for row in holdouts}),
            len(holdouts),
        )


@unittest.skipUnless(HAS_RDKIT, "install requirements-v5b.txt")
class V5b2ApplicabilityTests(unittest.TestCase):
    def setUp(self):
        self.records = [
            {
                "candidate_id": "pentane",
                "smiles": "CCCCC",
                "entry_error": 0.04,
            },
            {
                "candidate_id": "hexane",
                "smiles": "CCCCCC",
                "entry_error": 0.03,
            },
            {
                "candidate_id": "heptane",
                "smiles": "CCCCCCC",
                "entry_error": 0.05,
            },
            {
                "candidate_id": "octane",
                "smiles": "CCCCCCCC",
                "entry_error": 0.06,
            },
        ]

    def test_homolog_is_in_domain_with_neighbors(self):
        domain = assess_domain(
            "CCCCCC",
            self.records,
            in_domain_similarity=0.45,
            edge_similarity=0.25,
            minimum_neighbors=2,
            exclude_id="hexane",
        )
        self.assertEqual(domain.status, "in_domain")
        self.assertTrue(domain.uncertainty_valid)
        self.assertGreaterEqual(domain.neighbor_count, 2)

    def test_unrelated_structure_has_no_supported_uncertainty(self):
        domain = assess_domain(
            "N",
            self.records,
            in_domain_similarity=0.45,
            edge_similarity=0.25,
            minimum_neighbors=2,
        )
        self.assertIn(domain.status, {"edge", "out_of_domain"})
        self.assertFalse(domain.uncertainty_valid)

    def test_metric_calibration_produces_positive_bound(self):
        calibration = calibrate_metric(
            self.records,
            error_key="entry_error",
            coverage=0.8,
            k_neighbors=3,
            similarity_floor=0.05,
        )
        self.assertIsNotNone(calibration["conformal_factor"])
        self.assertGreater(calibration["conformal_factor"], 0.0)
        value = predict_relative_uncertainty(
            "CCCCCC",
            self.records,
            calibration,
            error_key="entry_error",
            k_neighbors=3,
            similarity_floor=0.05,
            exclude_id="hexane",
        )
        self.assertIsNotNone(value)
        self.assertGreater(value, 0.0)

    def test_local_scale_uses_structural_neighbors(self):
        value = local_error_scale(
            "CCCCCC",
            self.records,
            error_key="entry_error",
            k_neighbors=2,
            similarity_floor=0.05,
            exclude_id="hexane",
        )
        self.assertIsNotNone(value)
        self.assertGreater(value, 0.0)
        self.assertLess(value, 0.1)


if __name__ == "__main__":
    unittest.main()
