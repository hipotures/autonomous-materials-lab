from __future__ import annotations

import copy
from pathlib import Path
import sys
import unittest

import yaml

ENTRY = Path(__file__).resolve().parents[1]
ROOT = ENTRY.parent.parent
V5E = ROOT / "experiments" / "active-learning-v5e"
V5B = ROOT / "experiments" / "property-predictor-v5b"
for path in (V5E, V5B):
    sys.path.insert(0, str(path))

from acquisition import (  # noqa: E402
    hypothetical_domain_expansion,
    nearest_similarity,
    prepare_features,
    select_targets,
)
from run_acquisition import (  # noqa: E402
    _candidate_pool,
    _identity_metadata,
    _validate_config,
)


class V5eAcquisitionTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.config = yaml.safe_load(
            (V5E / "config.yaml").read_text(encoding="utf-8")
        )

    def _candidate(
        self,
        candidate_id: str,
        smiles: str,
        family: str,
        property_score: float,
    ):
        return {
            "candidate_id": candidate_id,
            "smiles": smiles,
            "lane": "exploratory_domain_expansion",
            "lane_reason": "domain_out_of_domain",
            "domain": {"status": "out_of_domain"},
            "property_screen": {
                "property_priority_score": property_score,
            },
            "provenance": [
                {
                    "family": family,
                    "operation": "test",
                    "parameters": {},
                    "raw_smiles": smiles,
                }
            ],
        }

    def test_nearest_similarity_recognizes_identical_reference(self):
        self.assertAlmostEqual(
            nearest_similarity("CCCO", ["CCO", "CCCO"]),
            1.0,
        )

    def test_prepare_features_rewards_property_and_novelty(self):
        candidates = [
            self._candidate("a", "CCCN", "amines", 0.9),
            self._candidate("b", "CCCCN", "amines", 0.5),
            self._candidate("c", "CC(=O)OC", "esters", 0.7),
        ]
        entry = {
            "CCCN": {
                "predicted_ratio_vs_water": 1.5,
            }
        }
        features = prepare_features(
            candidates,
            ["CCO", "CCCCC"],
            entry,
            self.config,
        )
        self.assertEqual(len(features), 3)
        by_id = {
            feature.candidate_id: feature
            for feature in features
        }
        self.assertGreater(
            by_id["a"].property_percentile,
            by_id["b"].property_percentile,
        )
        self.assertIsNotNone(
            by_id["a"].entry_percentile
        )
        self.assertGreaterEqual(
            by_id["a"].calibration_novelty,
            0.0,
        )

    def test_selection_is_deterministic_and_family_capped(self):
        config = copy.deepcopy(self.config)
        config["acquisition"].update(
            {
                "target_count": 6,
                "minimum_per_family": 1,
                "maximum_family_fraction": 0.5,
            }
        )
        candidates = [
            self._candidate(
                f"a{i}",
                "C" * (i + 2) + "N",
                "amines",
                0.9 - i * 0.03,
            )
            for i in range(5)
        ]
        candidates += [
            self._candidate(
                f"e{i}",
                "C" * (i + 2) + "C(=O)OC",
                "esters",
                0.8 - i * 0.03,
            )
            for i in range(5)
        ]
        candidates += [
            self._candidate(
                f"k{i}",
                "C" * (i + 2) + "C#C",
                "alkynes",
                0.7 - i * 0.03,
            )
            for i in range(5)
        ]

        features = prepare_features(
            candidates,
            ["CCO", "CCCCC"],
            {},
            config,
        )
        first = select_targets(features, config)
        second = select_targets(
            list(reversed(features)),
            config,
        )
        self.assertEqual(
            [row["candidate_id"] for row in first],
            [row["candidate_id"] for row in second],
        )
        self.assertEqual(len(first), 6)

        counts = {}
        for row in first:
            counts[row["family"]] = counts.get(row["family"], 0) + 1
        self.assertLessEqual(max(counts.values()), 3)
        self.assertTrue(
            {"amines", "esters", "alkynes"}.issubset(counts)
        )

    def test_domain_whatif_never_reduces_structural_coverage(self):
        candidates = [
            self._candidate("a", "CCCN", "amines", 0.9),
            self._candidate("b", "CCCCN", "amines", 0.8),
            self._candidate("c", "CC(=O)OC", "esters", 0.7),
            self._candidate("d", "CCC#C", "alkynes", 0.6),
        ]
        features = prepare_features(
            candidates,
            ["CCO", "CCCCC"],
            {},
            self.config,
        )
        config = copy.deepcopy(self.config)
        config["acquisition"]["target_count"] = 2
        config["acquisition"]["minimum_per_family"] = 0
        selected = select_targets(features, config)
        report = hypothetical_domain_expansion(
            features,
            ["CCO", "CCCCC"],
            selected,
            config,
        )
        self.assertGreaterEqual(
            report["after_if_all_selected_are_validated"][
                "fraction_ge_edge"
            ],
            report["before"]["fraction_ge_edge"],
        )
        self.assertGreaterEqual(
            report["after_if_all_selected_are_validated"][
                "fraction_ge_in_domain_similarity"
            ],
            report["before"][
                "fraction_ge_in_domain_similarity"
            ],
        )
        self.assertGreaterEqual(
            report[
                "after_excluding_selected_targets_if_validated"
            ]["fraction_ge_edge"],
            report["before_excluding_selected_targets"][
                "fraction_ge_edge"
            ],
        )


class V5eRunnerInputTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.config = yaml.safe_load(
            (V5E / "config.yaml").read_text(encoding="utf-8")
        )

    def test_identity_metadata_supports_reference_resolution(self):
        identity = _identity_metadata("CCCO")
        self.assertEqual(identity["canonical_smiles"], "CCCO")
        self.assertEqual(identity["molecular_formula"], "C3H8O")
        self.assertGreater(identity["molecular_weight_g_mol"], 50.0)
        self.assertTrue(
            identity["inchi_key"] is None
            or len(identity["inchi_key"]) > 10
        )

    def test_candidate_pool_uses_only_exploratory_full_prescreen(self):
        payload = {
            "passes": [
                {
                    "candidate_id": "explore",
                    "smiles": "CCCN",
                    "lane": "exploratory_domain_expansion",
                    "property_screen": {
                        "property_priority_score": 0.6,
                    },
                    "provenance": [{"family": "amines"}],
                },
                {
                    "candidate_id": "rankable",
                    "smiles": "CCCC",
                    "lane": "rankable",
                    "property_screen": {
                        "property_priority_score": 0.7,
                    },
                    "provenance": [{"family": "alkanes"}],
                },
            ]
        }
        pool = _candidate_pool(payload, self.config)
        self.assertEqual(len(pool), 1)
        self.assertEqual(pool[0]["candidate_id"], "explore")

    def test_default_config_is_valid(self):
        _validate_config(self.config)

    def test_invalid_weight_sum_is_rejected(self):
        config = copy.deepcopy(self.config)
        config["acquisition"]["weights"]["property_merit"] = 0.5
        with self.assertRaises(ValueError):
            _validate_config(config)

    def test_default_config_requests_substantial_calibration_set(self):
        acquisition = self.config["acquisition"]
        self.assertGreaterEqual(acquisition["target_count"], 30)
        self.assertGreater(
            acquisition["weights"]["selected_diversity"],
            0.0,
        )
        self.assertLessEqual(
            acquisition["maximum_family_fraction"],
            0.5,
        )


if __name__ == "__main__":
    unittest.main()
