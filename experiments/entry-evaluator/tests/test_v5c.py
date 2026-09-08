from __future__ import annotations

import importlib.util
from pathlib import Path
import sys
import unittest

import yaml

ENTRY = Path(__file__).resolve().parents[1]
ROOT = ENTRY.parent.parent
V5C = ROOT / "experiments" / "molecular-search-v5c"
V5B = ROOT / "experiments" / "property-predictor-v5b"
for path in (V5C, V5B):
    sys.path.insert(0, str(path))

HAS_RDKIT = importlib.util.find_spec("rdkit") is not None

from applicability import canonical_smiles  # noqa: E402
from generator import generate_candidates  # noqa: E402
from run_search import (  # noqa: E402
    _candidate_pool,
    _descriptor_check,
    conservative_score,
)


class V5cRankingTests(unittest.TestCase):
    def test_conservative_score_penalizes_uncertainty(self):
        self.assertAlmostEqual(
            conservative_score(
                3000.0,
                0.10,
                multiplier=1.0,
            ),
            3333.3333333333335,
        )
        self.assertAlmostEqual(
            conservative_score(
                3000.0,
                0.10,
                multiplier=2.0,
            ),
            3750.0,
        )

    def test_uncertainty_without_finite_upper_bound_is_rejected(self):
        with self.assertRaises(ValueError):
            conservative_score(100.0, 1.0)
        with self.assertRaises(ValueError):
            conservative_score(100.0, 0.6, multiplier=2.0)

    def test_invalid_conservative_score_inputs_rejected(self):
        for predicted, uncertainty, multiplier in (
            (0.0, 0.1, 1.0),
            (-1.0, 0.1, 1.0),
            (1.0, -0.1, 1.0),
            (1.0, 0.1, -1.0),
        ):
            with self.assertRaises(ValueError):
                conservative_score(
                    predicted,
                    uncertainty,
                    multiplier=multiplier,
                )


@unittest.skipUnless(HAS_RDKIT, "install requirements-v5b.txt")
class V5cGeneratorTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.config = yaml.safe_load(
            (V5C / "config.yaml").read_text(encoding="utf-8")
        )

    def test_generator_is_deterministic(self):
        first = generate_candidates(self.config)
        second = generate_candidates(self.config)
        self.assertEqual(first, second)
        self.assertGreater(len(first), 30)

    def test_generator_spans_expected_families(self):
        candidates = generate_candidates(self.config)
        families = {candidate.family for candidate in candidates}
        self.assertTrue(
            {
                "alkanes",
                "alkenes",
                "oxygenated",
                "aromatics",
                "cyclic_hydrocarbons",
            }.issubset(families)
        )
        operations = {
            candidate.operation
            for candidate in candidates
        }
        self.assertIn("single_methyl_branch", operations)

    def test_alkenes_are_terminal_and_ketone_range_includes_acetone(self):
        candidates = generate_candidates(self.config)
        alkenes = [
            candidate
            for candidate in candidates
            if candidate.family == "alkenes"
        ]
        self.assertTrue(alkenes)
        self.assertTrue(
            all(candidate.smiles.startswith("C=") for candidate in alkenes)
        )
        ketone_smiles = {
            canonical_smiles(candidate.smiles)
            for candidate in candidates
            if candidate.operation == "ketone_chain_and_position"
        }
        self.assertIn(canonical_smiles("CC(=O)C"), ketone_smiles)

    def test_generated_smiles_are_parseable(self):
        for candidate in generate_candidates(self.config):
            canonical = canonical_smiles(candidate.smiles)
            self.assertTrue(canonical)

    def test_candidate_pool_removes_known_references_and_duplicates(self):
        calibration = [
            {
                "candidate_id": "known-pentane",
                "smiles": "CCCCC",
                "entry_error": 0.03,
            },
            {
                "candidate_id": "known-ethanol",
                "smiles": "CCO",
                "entry_error": 0.04,
            },
            {
                "candidate_id": "known-benzene",
                "smiles": "c1ccccc1",
                "entry_error": 0.05,
            },
        ]
        candidates, stats = _candidate_pool(
            self.config,
            calibration,
        )
        smiles = {row["smiles"] for row in candidates}
        self.assertNotIn(canonical_smiles("CCCCC"), smiles)
        self.assertNotIn(canonical_smiles("CCO"), smiles)
        self.assertEqual(
            len(smiles),
            len(candidates),
        )
        self.assertGreater(stats["known_reference_removed_count"], 0)
        self.assertGreater(stats["novel_candidate_count"], 0)

    def test_descriptor_filter_accepts_neutral_co_molecule(self):
        descriptor, failure = _descriptor_check(
            canonical_smiles("CCCO"),
            self.config["candidate_constraints"],
        )
        self.assertIsNone(failure)
        self.assertIsNotNone(descriptor)
        self.assertEqual(descriptor["formal_charge"], 0)

    def test_descriptor_filter_rejects_disallowed_element(self):
        descriptor, failure = _descriptor_check(
            canonical_smiles("CCCl"),
            self.config["candidate_constraints"],
        )
        self.assertIsNone(descriptor)
        self.assertEqual(failure, "disallowed_element")


class V5cConfigTests(unittest.TestCase):
    def test_search_requires_robustness_and_domain_by_default(self):
        config = yaml.safe_load(
            (V5C / "config.yaml").read_text(encoding="utf-8")
        )
        self.assertTrue(config["search"]["require_v5c_ready"])
        self.assertTrue(config["search"]["require_in_domain"])
        self.assertTrue(
            config["search"]["require_screening_uncertainty"]
        )
        self.assertEqual(
            config["storage"]["pressure_pa"],
            101325.0,
        )


if __name__ == "__main__":
    unittest.main()
