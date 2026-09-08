from __future__ import annotations

from pathlib import Path
import copy
import sys
import unittest

import yaml

ENTRY = Path(__file__).resolve().parents[1]
ROOT = ENTRY.parent.parent
V5D = ROOT / "experiments" / "molecular-search-v5d"
V5C = ROOT / "experiments" / "molecular-search-v5c"
V5B = ROOT / "experiments" / "property-predictor-v5b"
for path in (V5D, V5C, V5B, ENTRY):
    sys.path.insert(0, str(path))

from graph_mutator import (  # noqa: E402
    classify_family,
    generate_mutations,
    mutate_once,
    structural_bucket,
)
from run_search_v5d import (  # noqa: E402
    _beam_select,
    _coarse_config,
    _generate_structural_generation,
    _structural_select,
    _validate_config,
)


class V5dMutationTests(unittest.TestCase):
    def test_mutations_are_deterministic(self):
        kwargs = {
            "generation": 2,
            "count": 20,
            "seed": "fixed-seed",
            "allowed_atomic_numbers": [6, 7, 8],
            "operators": [
                "add_terminal",
                "delete_terminal",
                "change_atom",
                "change_bond",
                "insert_atom",
                "close_ring_5_6",
                "remove_ring_bond",
            ],
        }
        first = generate_mutations("CCCO", **kwargs)
        second = generate_mutations("CCCO", **kwargs)
        self.assertEqual(first, second)
        self.assertGreater(len(first), 0)

    def test_single_mutations_return_sanitized_smiles_or_none(self):
        operators = [
            "add_terminal",
            "delete_terminal",
            "change_atom",
            "change_bond",
            "insert_atom",
            "close_ring_5_6",
            "remove_ring_bond",
        ]
        for index, operator in enumerate(operators):
            result = mutate_once(
                "CCCCCC",
                operator=operator,
                seed=100 + index,
                allowed_atomic_numbers=[6, 7, 8],
            )
            if result is not None:
                self.assertNotIn(".", result)

    def test_family_classifier_matches_search_policy(self):
        self.assertEqual(classify_family("CCCC"), "alkanes")
        self.assertEqual(classify_family("C=CCC"), "alkenes")
        self.assertEqual(classify_family("CCCO"), "oxygenated")
        self.assertEqual(classify_family("CC(=O)OC"), "esters")
        self.assertEqual(classify_family("CCC#C"), "alkynes")
        self.assertEqual(classify_family("CCCN"), "amines")
        self.assertEqual(
            classify_family("OCCC(=O)C"),
            "mixed_functional",
        )
        self.assertEqual(
            classify_family("OC1CCCCC1"),
            "mixed_functional",
        )
        self.assertEqual(
            classify_family("OCCC=C"),
            "mixed_functional",
        )
        self.assertEqual(
            classify_family("C1CCC=CC1"),
            "mixed_functional",
        )

    def test_structural_bucket_is_stable(self):
        first = structural_bucket("CCCO")
        second = structural_bucket("CCCO")
        self.assertEqual(first, second)
        self.assertEqual(first[0], "oxygenated")


class V5dSelectionTests(unittest.TestCase):
    def _rows(self):
        rows = []
        for index, (family, heavy, score) in enumerate(
            [
                ("alkanes", 5, 0.5),
                ("alkanes", 6, 0.6),
                ("oxygenated", 4, 0.8),
                ("oxygenated", 5, 0.7),
                ("amines", 4, 0.9),
                ("esters", 5, 0.65),
            ]
        ):
            rows.append(
                {
                    "candidate_id": f"c{index}",
                    "smiles": "C" * max(heavy, 2),
                    "descriptor": {
                        "heavy_atoms": heavy,
                        "atomic_numbers": [6],
                        "rings": 0,
                    },
                    "provenance": [
                        {
                            "family": family,
                            "operation": "test",
                            "parameters": {},
                            "raw_smiles": "CC",
                        }
                    ],
                    "coarse_property_screen": {
                        "property_priority_score": score,
                    },
                }
            )
        return rows

    def test_structural_select_is_deterministic_and_budgeted(self):
        rows = self._rows()
        first = _structural_select(rows, 4)
        second = _structural_select(list(reversed(rows)), 4)
        self.assertEqual(
            [row["candidate_id"] for row in first],
            [row["candidate_id"] for row in second],
        )
        self.assertEqual(len(first), 4)

    def test_beam_combines_score_and_diversity(self):
        rows = self._rows()
        beam = _beam_select(
            rows,
            width=4,
            diversity_fraction=0.5,
        )
        self.assertEqual(len(beam), 4)
        self.assertIn(
            "c4",
            {row["candidate_id"] for row in beam},
        )


class V5dGenerationTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.config = yaml.safe_load(
            (V5D / "config.yaml").read_text(encoding="utf-8")
        )

    def test_config_is_high_throughput_not_template_scale(self):
        ht = self.config["high_throughput"]
        self.assertGreaterEqual(ht["total_unique_target"], 100000)
        self.assertGreaterEqual(ht["full_prescreen_budget"], 1000)
        self.assertGreaterEqual(ht["entry_rankable_budget"], 100)
        self.assertGreaterEqual(ht["generations"], 4)

    def test_default_config_passes_validation(self):
        _validate_config(self.config)

    def test_config_validation_rejects_invalid_beam(self):
        config = copy.deepcopy(self.config)
        config["high_throughput"]["beam_width"] = 0
        with self.assertRaises(ValueError):
            _validate_config(config)

    def test_config_validation_rejects_overlapping_family_policy(self):
        config = copy.deepcopy(self.config)
        config["domain_expansion"]["exploratory_families"].append(
            "oxygenated"
        )
        with self.assertRaises(ValueError):
            _validate_config(config)

    def test_coarse_config_uses_smaller_grid(self):
        coarse = _coarse_config(self.config)
        full_states = (
            len(self.config["prescreen"]["temperature_grid_k"])
            * len(self.config["prescreen"]["pressure_grid_pa"])
        )
        coarse_states = (
            len(coarse["prescreen"]["temperature_grid_k"])
            * len(coarse["prescreen"]["pressure_grid_pa"])
        )
        self.assertLess(coarse_states, full_states)

    def test_small_generation_produces_unique_valid_candidates(self):
        config = copy.deepcopy(self.config)
        config["high_throughput"].update(
            {
                "children_per_parent": 12,
                "generation_unique_target": 20,
                "total_unique_target": 1000,
            }
        )
        seen = {"CCCO"}
        rows, report = _generate_structural_generation(
            ["CCCO"],
            generation=1,
            config=config,
            seen=seen,
        )
        smiles = [row["smiles"] for row in rows]
        self.assertEqual(len(smiles), len(set(smiles)))
        self.assertLessEqual(len(rows), 20)
        self.assertGreater(report["raw_mutation_attempt_count"], 0)
        self.assertEqual(report["new_structural_candidate_count"], len(rows))

    def test_all_configured_families_are_explicit(self):
        rankable = set(
            self.config["domain_expansion"]["rankable_families"]
        )
        exploratory = set(
            self.config["domain_expansion"]["exploratory_families"]
        )
        self.assertFalse(rankable & exploratory)
        self.assertIn("mixed_functional", exploratory)
        self.assertIn("oxygenated", rankable)


if __name__ == "__main__":
    unittest.main()
