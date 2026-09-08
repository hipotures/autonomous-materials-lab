from __future__ import annotations

from pathlib import Path
from types import SimpleNamespace
import sys
import unittest

import yaml

ENTRY = Path(__file__).resolve().parents[1]
ROOT = ENTRY.parent.parent
V5C = ROOT / "experiments" / "molecular-search-v5c"
V5B = ROOT / "experiments" / "property-predictor-v5b"
for path in (V5C, V5B, ENTRY):
    sys.path.insert(0, str(path))

from applicability import canonical_smiles  # noqa: E402
from generator_v5c2 import generate_candidates  # noqa: E402
from run_search_v5c2 import (  # noqa: E402
    _atomic_numbers,
    _calibration_atomic_numbers,
    _candidate_pool,
    _family_is_configured,
    _family_is_rankable,
)
from thermo_prescreen import (  # noqa: E402
    compare_to_reference,
    passes_property_gate,
    property_screen,
)


class _FakeProvider:
    def __init__(self, enthalpy_scale: float):
        self.enthalpy_scale = float(enthalpy_scale)

    def enthalpy_j_kg(self, temperature_k: float, pressure_pa: float) -> float:
        base = self.enthalpy_scale * (float(temperature_k) - 293.15)
        pressure_penalty = 2.0e-3 * max(float(pressure_pa) - 101325.0, 0.0)
        latent = 200000.0 if float(temperature_k) > 350.0 else 0.0
        return base - pressure_penalty + latent

    def state(self, temperature_k: float, pressure_pa: float):
        return SimpleNamespace(
            enthalpy_j_kg=self.enthalpy_j_kg(temperature_k, pressure_pa),
            density_kg_m3=800.0,
            cp_j_kg_k=self.enthalpy_scale,
            phase="liquid" if float(temperature_k) < 350.0 else "gas",
        )

    def phase(self, temperature_k: float, pressure_pa: float) -> str:
        return "liquid" if float(temperature_k) < 350.0 else "gas"

    def saturation_at_pressure(self, pressure_pa: float):
        return SimpleNamespace(
            supported=True,
            bubble_temperature_k=350.0,
            dew_temperature_k=350.0,
            failure_reason=None,
        )


class V5c2PropertyPrescreenTests(unittest.TestCase):
    def setUp(self):
        self.config = {
            "storage": {
                "temperature_k": 293.15,
                "pressure_pa": 101325.0,
            },
            "prescreen": {
                "temperature_grid_k": [350.0, 400.0, 500.0],
                "pressure_grid_pa": [101325.0, 300000.0],
                "latent_proxy_offset_k": 2.0,
                "minimum_positive_delta_h_fraction": 0.8,
                "minimum_enthalpy_ratio_vs_water": 0.25,
            },
        }

    def test_prescreen_reports_enthalpy_grid_and_latent_proxy(self):
        metrics = property_screen(_FakeProvider(2000.0), self.config)
        self.assertEqual(metrics["storage_phase"], "liquid")
        self.assertEqual(metrics["requested_grid_state_count"], 6)
        self.assertGreater(metrics["positive_delta_h_fraction"], 0.8)
        self.assertGreater(metrics["delta_h_q25_j_kg"], 0.0)
        self.assertGreater(metrics["latent_proxy_j_kg"], 0.0)

    def test_q25_includes_negative_grid_states(self):
        config = {
            "storage": self.config["storage"],
            "prescreen": {
                **self.config["prescreen"],
                "pressure_grid_pa": [101325.0, 1000000000.0],
            },
        }
        metrics = property_screen(_FakeProvider(1000.0), config)
        self.assertLess(metrics["positive_delta_h_fraction"], 0.8)
        self.assertLess(metrics["delta_h_q25_j_kg"], 0.0)

    def test_property_priority_is_water_normalized(self):
        water = property_screen(_FakeProvider(1500.0), self.config)
        candidate = property_screen(_FakeProvider(2500.0), self.config)
        targets = compare_to_reference(candidate, water)
        self.assertGreater(targets["property_priority_score"], 1.0)

        merged = {**candidate, **targets}
        passed, reason = passes_property_gate(merged, self.config)
        self.assertTrue(passed)
        self.assertIsNone(reason)


class V5c2GeneratorTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.config = yaml.safe_load(
            (V5C / "config-v5c2.yaml").read_text(encoding="utf-8")
        )

    def test_expanded_generator_contains_new_supported_motifs(self):
        candidates = generate_candidates(self.config)
        families = {candidate.family for candidate in candidates}
        self.assertTrue(
            {
                "aldehydes",
                "alkynes",
                "esters",
                "amines",
            }.issubset(families)
        )
        operations = {candidate.operation for candidate in candidates}
        self.assertIn("branched_primary_alcohol", operations)
        self.assertGreater(len(candidates), 80)

    def test_all_generated_structures_are_parseable(self):
        for candidate in generate_candidates(self.config):
            self.assertTrue(canonical_smiles(candidate.smiles))

    def test_nitrogen_is_model_space_but_not_calibration_element_space(self):
        calibration = [
            {
                "candidate_id": "ethanol",
                "smiles": "CCO",
                "entry_error": 0.04,
            },
            {
                "candidate_id": "pentane",
                "smiles": "CCCCC",
                "entry_error": 0.03,
            },
        ]
        supported = _calibration_atomic_numbers(calibration)
        self.assertEqual(supported, {6, 8})
        amine_elements = _atomic_numbers(canonical_smiles("CCCN"))
        self.assertEqual(amine_elements, {6, 7})
        self.assertFalse(amine_elements.issubset(supported))

    def test_prior_v5c1_candidates_are_removed(self):
        calibration = [
            {
                "candidate_id": "ethanol",
                "smiles": "CCO",
                "entry_error": 0.04,
            },
            {
                "candidate_id": "pentane",
                "smiles": "CCCCC",
                "entry_error": 0.03,
            },
            {
                "candidate_id": "benzene",
                "smiles": "c1ccccc1",
                "entry_error": 0.05,
            },
        ]
        prior = {canonical_smiles("CCCO")}
        candidates, stats = _candidate_pool(
            self.config,
            calibration,
            prior,
        )
        smiles = {row["smiles"] for row in candidates}
        self.assertNotIn(canonical_smiles("CCCO"), smiles)
        self.assertGreater(stats["prior_search_removed_count"], 0)

    def test_unconfigured_family_is_not_silently_accepted(self):
        candidate = {
            "provenance": [
                {
                    "family": "invented_family",
                    "operation": "test",
                    "parameters": {},
                    "raw_smiles": "CCC",
                }
            ]
        }
        self.assertFalse(
            _family_is_configured(candidate, self.config)
        )

    def test_new_functional_families_are_exploratory(self):
        candidate = {
            "provenance": [
                {
                    "family": "esters",
                    "operation": "test",
                    "parameters": {},
                    "raw_smiles": "CC(=O)OC",
                }
            ]
        }
        self.assertFalse(
            _family_is_rankable(candidate, self.config)
        )
        candidate["provenance"][0]["family"] = "oxygenated"
        self.assertTrue(
            _family_is_rankable(candidate, self.config)
        )

    def test_v5c2_config_expands_elements_but_limits_supported_rings(self):
        constraints = self.config["candidate_constraints"]
        self.assertEqual(
            constraints["allowed_atomic_numbers"],
            [6, 7, 8],
        )
        rings = self.config["generation"]["cycloalkanes"]
        self.assertEqual(rings["min_ring"], 5)
        self.assertEqual(rings["max_ring"], 6)


if __name__ == "__main__":
    unittest.main()
