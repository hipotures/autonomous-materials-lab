from __future__ import annotations

import importlib.util
from pathlib import Path
import sys
import unittest

HERE = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(HERE))

from chemistry import ChemistryLimiter
from coolant import CoolantModel
from property_provider import build_property_provider

HAS_V5B = (
    importlib.util.find_spec("feos") is not None
    and importlib.util.find_spec("rdkit") is not None
)


class V5bStaticContractTests(unittest.TestCase):
    def test_feos_factory_requires_structure_not_known_properties(self):
        with self.assertRaisesRegex(ValueError, "SMILES"):
            build_property_provider(
                {
                    "property_provider": {
                        "type": "feos",
                        "model": "gc_pcsaft_joback",
                    }
                }
            )

    def test_feos_provider_source_does_not_import_coolprop(self):
        source = (HERE / "feos_gc_provider.py").read_text(encoding="utf-8")
        self.assertNotIn("import CoolProp", source)
        self.assertNotIn("PropsSI", source)
        self.assertIn("from_json_smiles", source)
        self.assertIn("gc_pcsaft", source)
        self.assertIn(".joback(", source)

    def test_vendored_parameter_files_exist(self):
        parameter_dir = HERE / "parameters" / "v5b"
        expected = {
            "sauer2014_smarts.json",
            "rehner2023_hetero.json",
            "joback1987.json",
        }
        self.assertEqual(
            {path.name for path in parameter_dir.glob("*.json")},
            expected,
        )


@unittest.skipUnless(HAS_V5B, "install requirements-v5b.txt")
class V5bRuntimeTests(unittest.TestCase):
    def _provider(self, smiles: str, name: str):
        return build_property_provider(
            {
                "property_provider": {
                    "type": "feos",
                    "model": "gc_pcsaft_joback",
                    "name": name,
                    "smiles": smiles,
                }
            }
        )

    def test_ethanol_structure_prediction_is_physical(self):
        provider = self._provider("CCO", "ethanol")
        storage = provider.state(293.15, 101325.0)
        hot = provider.state(450.0, 25000.0)
        self.assertEqual(provider.phase(293.15, 101325.0), "liquid")
        self.assertGreater(storage.density_kg_m3, 400.0)
        self.assertLess(storage.density_kg_m3, 1200.0)
        self.assertGreater(hot.enthalpy_j_kg - storage.enthalpy_j_kg, 0.0)
        self.assertGreater(storage.cp_j_kg_k, 0.0)
        meta = provider.metadata()
        self.assertTrue(meta["structure_is_only_candidate_input"])
        self.assertFalse(meta["reference_property_backend_used"])
        self.assertEqual(meta["smiles"], "CCO")
        for item in meta["parameter_files"].values():
            self.assertEqual(len(item["sha256"]), 64)

    def test_acetone_structure_prediction_is_physical(self):
        provider = self._provider("CC(=O)C", "acetone")
        storage = provider.state(293.15, 101325.0)
        hot = provider.state(450.0, 25000.0)
        self.assertEqual(provider.phase(293.15, 101325.0), "liquid")
        self.assertGreater(storage.density_kg_m3, 400.0)
        self.assertLess(storage.density_kg_m3, 1400.0)
        self.assertGreater(hot.enthalpy_j_kg - storage.enthalpy_j_kg, 0.0)

    def test_structure_provider_runs_through_coolant_model(self):
        model = CoolantModel(
            {
                "storage_temperature_k": 293.15,
                "storage_pressure_pa": 101325.0,
                "max_exit_temperature_k": 500.0,
                "wall_to_fluid_approach_k": 25.0,
                "cooled_area_m2": 1.0,
                "property_provider": {
                    "type": "feos",
                    "model": "gc_pcsaft_joback",
                    "name": "ethanol-holdout",
                    "smiles": "CCO",
                },
            },
            ChemistryLimiter({"mode": "disabled"}),
        )
        step = model.evaluate(
            wall_temperature_k=500.0,
            surface_pressure_pa=25000.0,
            coolant_heat_flux_w_m2=1.0e6,
            area_m2=1.0,
        )
        self.assertTrue(step.feasible, step.failure_reason)
        self.assertGreater(step.usable_enthalpy_j_kg, 0.0)
        self.assertAlmostEqual(
            step.mass_flow_kg_s * step.usable_enthalpy_j_kg,
            1.0e6,
            places=5,
        )


if __name__ == "__main__":
    unittest.main()
