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

from failure_taxonomy import (  # noqa: E402
    ENTRY_EVALUATOR_FAILURE,
    PROPERTY_MODEL_FAILURE,
    STORAGE_STATE_INFEASIBLE,
    STRUCTURE_OUT_OF_DOMAIN,
    SUCCESS,
    classify_prediction,
    expected_prediction_outcome,
)
from run_robustness import (  # noqa: E402
    _aggregate_repeated,
    _evaluate_split,
)


class V5b21FailureTaxonomyTests(unittest.TestCase):
    def test_structure_out_of_domain_is_distinct_from_storage_failure(self):
        self.assertEqual(
            classify_prediction(
                result=None,
                failure_reason=(
                    "SMILES is outside the pinned V5b GC-PC-SAFT/Joback "
                    "domain: N: Molecule cannot be built from groups!"
                ),
            ),
            STRUCTURE_OUT_OF_DOMAIN,
        )
        self.assertEqual(
            classify_prediction(
                result=None,
                failure_reason=(
                    "configured coolant storage state is not liquid: "
                    "feos:gc-pcsaft+joback:CO phase=gas"
                ),
            ),
            STORAGE_STATE_INFEASIBLE,
        )

    def test_property_and_entry_failures_have_separate_classes(self):
        self.assertEqual(
            classify_prediction(
                result=None,
                failure_reason=(
                    "FeOS TP state failed for candidate"
                ),
            ),
            PROPERTY_MODEL_FAILURE,
        )
        self.assertEqual(
            classify_prediction(
                result={
                    "entry": {
                        "status": "failed",
                        "failure_reason": "wall exceeded limit",
                    }
                },
                failure_reason=None,
            ),
            ENTRY_EVALUATOR_FAILURE,
        )
        self.assertEqual(
            classify_prediction(
                result={"entry": {"status": "terminal_velocity"}},
                failure_reason=None,
            ),
            SUCCESS,
        )

    def test_benchmark_probe_expectations_match_semantics(self):
        study = yaml.safe_load(
            (V5B / "benchmark-v5b2.yaml").read_text(
                encoding="utf-8"
            )
        )
        by_id = {
            row["id"]: row
            for row in study["holdouts"]
        }
        self.assertEqual(
            expected_prediction_outcome(
                by_id["methanol-domain-probe"]
            ),
            STORAGE_STATE_INFEASIBLE,
        )
        for candidate_id in (
            "water-domain-probe",
            "ammonia-domain-probe",
            "cyclopropane-domain-probe",
        ):
            self.assertEqual(
                expected_prediction_outcome(by_id[candidate_id]),
                STRUCTURE_OUT_OF_DOMAIN,
            )

    def test_supported_holdouts_have_family_labels(self):
        study = yaml.safe_load(
            (V5B / "benchmark-v5b2.yaml").read_text(
                encoding="utf-8"
            )
        )
        supported = [
            row
            for row in study["holdouts"]
            if expected_prediction_outcome(row) == SUCCESS
        ]
        self.assertGreaterEqual(len(supported), 20)
        self.assertTrue(
            all(row.get("challenge_family") for row in supported)
        )
        self.assertGreaterEqual(
            len(
                {
                    row["challenge_family"]
                    for row in supported
                }
            ),
            5,
        )


@unittest.skipUnless(HAS_RDKIT, "install requirements-v5b.txt")
class V5b21RepeatedSplitTests(unittest.TestCase):
    def setUp(self):
        self.records = [
            {
                "candidate_id": "pentane",
                "smiles": "CCCCC",
                "entry_error": 0.04,
                "storage_density_error": 0.03,
                "delta_h_error": 0.05,
                "cp_error": 0.04,
                "saturation_error": 0.02,
            },
            {
                "candidate_id": "hexane",
                "smiles": "CCCCCC",
                "entry_error": 0.03,
                "storage_density_error": 0.04,
                "delta_h_error": 0.04,
                "cp_error": 0.05,
                "saturation_error": 0.03,
            },
            {
                "candidate_id": "heptane",
                "smiles": "CCCCCCC",
                "entry_error": 0.05,
                "storage_density_error": 0.02,
                "delta_h_error": 0.06,
                "cp_error": 0.03,
                "saturation_error": 0.04,
            },
            {
                "candidate_id": "octane",
                "smiles": "CCCCCCCC",
                "entry_error": 0.06,
                "storage_density_error": 0.05,
                "delta_h_error": 0.07,
                "cp_error": 0.06,
                "saturation_error": 0.05,
            },
            {
                "candidate_id": "nonane",
                "smiles": "CCCCCCCCC",
                "entry_error": 0.04,
                "storage_density_error": 0.03,
                "delta_h_error": 0.05,
                "cp_error": 0.04,
                "saturation_error": 0.03,
            },
            {
                "candidate_id": "decane",
                "smiles": "CCCCCCCCCC",
                "entry_error": 0.05,
                "storage_density_error": 0.04,
                "delta_h_error": 0.06,
                "cp_error": 0.05,
                "saturation_error": 0.04,
            },
        ]
        self.cfg = {
            "calibration_fraction": 0.67,
            "coverage_target": 0.8,
            "k_neighbors": 3,
            "similarity_floor": 0.05,
            "in_domain_similarity": 0.45,
            "edge_similarity": 0.25,
            "minimum_neighbors": 1,
            "minimum_in_domain_evaluation_candidates": 1,
        }

    def test_repeated_split_is_deterministic(self):
        a = _evaluate_split(
            self.records,
            self.cfg,
            seed="repeat-001",
        )
        b = _evaluate_split(
            list(reversed(self.records)),
            self.cfg,
            seed="repeat-001",
        )
        self.assertEqual(
            a["calibration_candidate_ids"],
            b["calibration_candidate_ids"],
        )
        self.assertEqual(
            a["evaluation_candidate_ids"],
            b["evaluation_candidate_ids"],
        )
        self.assertEqual(a["coverage"], b["coverage"])

    def test_aggregate_reports_pooled_and_split_coverage(self):
        splits = [
            _evaluate_split(
                self.records,
                self.cfg,
                seed=f"repeat-{index:03d}",
            )
            for index in range(4)
        ]
        report = _aggregate_repeated(splits)
        self.assertEqual(report["split_count"], 4)
        self.assertGreater(report["valid_split_count"], 0)
        entry = report["metrics"]["entry_error"]
        self.assertIn("pooled_coverage", entry)
        self.assertIn("split_coverage", entry)
        self.assertGreater(
            entry["pooled_eligible_count"],
            0,
        )


class V5b21LanguageTests(unittest.TestCase):
    def test_uncertainty_code_does_not_claim_certification(self):
        paths = [
            V5B / "run_calibration.py",
            V5B / "run_robustness.py",
            V5B / "README.md",
        ]
        text = "\n".join(
            path.read_text(encoding="utf-8").lower()
            for path in paths
        )
        self.assertNotIn("certified", text)


if __name__ == "__main__":
    unittest.main()
