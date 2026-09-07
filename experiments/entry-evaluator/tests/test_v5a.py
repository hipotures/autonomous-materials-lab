from __future__ import annotations

import math
from pathlib import Path
import sys
import unittest

import CoolProp.CoolProp as CP

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from chemistry import ChemistryLimiter
from coolant import CoolantModel
from property_provider import (
    CoolPropHEOSProvider,
    build_property_provider,
)
from run_mixture_sweep import continuity_report


class PropertyProviderTests(unittest.TestCase):
    def test_legacy_pure_fluid_matches_propssi(self):
        provider = build_property_provider(
            {"coolprop_name": "Water"}
        )
        h = provider.enthalpy_j_kg(350.0, 2.0e5)
        expected = CP.PropsSI(
            "H",
            "T",
            350.0,
            "P",
            2.0e5,
            "Water",
        )
        self.assertAlmostEqual(
            h,
            expected,
            delta=abs(expected) * 1e-10 + 1e-8,
        )
        self.assertEqual(
            provider.metadata()["components"],
            ["Water"],
        )

    def test_zero_fraction_endpoint_collapses_to_true_pure_fluid(self):
        ethanol = CoolPropHEOSProvider(
            components=["Water", "Ethanol"],
            fractions=[0.0, 1.0],
        )
        water = CoolPropHEOSProvider(
            components=["Water", "Ethanol"],
            fractions=[1.0, 0.0],
        )
        self.assertEqual(
            ethanol.metadata()["components"],
            ["Ethanol"],
        )
        self.assertEqual(
            ethanol.metadata()["mole_fractions"],
            [1.0],
        )
        self.assertEqual(
            water.metadata()["components"],
            ["Water"],
        )
        self.assertEqual(
            water.metadata()["mole_fractions"],
            [1.0],
        )

    def test_stability_algorithm_is_explicit_in_metadata(self):
        provider = CoolPropHEOSProvider(
            components=["Water", "Ethanol"],
            fractions=[0.5, 0.5],
            stability_algorithm=0,
        )
        meta = provider.metadata()
        self.assertEqual(meta["stability_algorithm"], 0)
        self.assertEqual(
            meta["stability_algorithm_name"],
            "legacy_Gernert",
        )
        with self.assertRaises(ValueError):
            CoolPropHEOSProvider(
                components=["Water", "Ethanol"],
                fractions=[0.5, 0.5],
                stability_algorithm=2,
            )

    def test_invalid_fraction_vector_rejected(self):
        with self.assertRaises(ValueError):
            CoolPropHEOSProvider(
                components=["Water", "Ethanol"],
                fractions=[0.5],
            )
        with self.assertRaises(ValueError):
            CoolPropHEOSProvider(
                components=["Water", "Ethanol"],
                fractions=[-0.1, 1.1],
            )

    def test_problematic_mid_compositions_return_physical_hot_low_pressure_state(self):
        for x_water in (0.4, 0.5, 0.6):
            provider = CoolPropHEOSProvider(
                components=["Water", "Ethanol"],
                fractions=[x_water, 1.0 - x_water],
                stability_algorithm=1,
            )
            h_storage = provider.enthalpy_j_kg(
                293.15,
                101325.0,
            )

            # Exercise saturation diagnostics first.  V5a deliberately does
            # not require CoolProp's textual phase label to be "gas"; the
            # regression target is the old metastable dense-liquid root and
            # its absurdly negative enthalpy.
            for pressure_pa in (25000.0, 30000.0, 25000.0):
                saturation = provider.saturation_at_pressure(pressure_pa)
                self.assertTrue(saturation.supported)

            hot = provider.state(500.0, 25000.0)

            # At 25 kPa and 500 K the stable state must be vapor-like.  The
            # broken 8.0.0 path returned a dense liquid root with enthalpy of
            # order -1e8..-1e9 J/kg.  Use generous physical guards so this
            # test detects that wrong root without depending on phase labels.
            self.assertLess(hot.density_kg_m3, 10.0)
            self.assertGreater(
                hot.enthalpy_j_kg - h_storage,
                0.0,
            )
            self.assertLess(
                abs(hot.enthalpy_j_kg),
                1.0e7,
            )

    def test_documented_water_ethanol_bubble_point(self):
        provider = CoolPropHEOSProvider(
            components=["Water", "Ethanol"],
            fractions=[0.4, 0.6],
        )
        saturation = provider.saturation_at_pressure(101325.0)
        self.assertTrue(saturation.supported)
        self.assertIsNotNone(
            saturation.bubble_temperature_k
        )
        self.assertAlmostEqual(
            saturation.bubble_temperature_k,
            352.3522142890429,
            places=5,
        )

    def test_mixture_storage_state_and_coolant_model(self):
        coolant_config = {
            "storage_temperature_k": 293.15,
            "storage_pressure_pa": 101325.0,
            "max_exit_temperature_k": 500.0,
            "wall_to_fluid_approach_k": 25.0,
            "cooled_area_m2": 1.0,
            "injection_pressure_margin": 1.15,
            "porous_delta_p_pa": 0.0,
            "max_injection_pressure_pa": 2.0e7,
            "property_provider": {
                "type": "coolprop",
                "backend": "HEOS",
                "components": ["Water", "Ethanol"],
                "composition_basis": "mole",
                "fractions": [0.5, 0.5],
            },
        }
        model = CoolantModel(
            coolant_config,
            ChemistryLimiter({"mode": "disabled"}),
        )
        self.assertIn(
            model.properties.phase(293.15, 101325.0),
            {"liquid", "supercritical_liquid"},
        )
        step = model.evaluate(
            wall_temperature_k=450.0,
            surface_pressure_pa=2.0e5,
            coolant_heat_flux_w_m2=1.0e6,
            area_m2=1.0,
        )
        self.assertTrue(step.feasible)
        self.assertGreater(step.usable_enthalpy_j_kg, 0.0)
        self.assertGreater(step.mass_flow_kg_s, 0.0)


    def test_legacy_coolant_name_is_preserved(self):
        model = CoolantModel(
            {
                "coolprop_name": "Water",
                "storage_temperature_k": 293.15,
                "storage_pressure_pa": 101325.0,
                "max_exit_temperature_k": 500.0,
                "cooled_area_m2": 1.0,
            },
            ChemistryLimiter({"mode": "disabled"}),
        )
        self.assertEqual(model.fluid, "Water")
        self.assertEqual(
            model.properties.metadata()["identity"],
            "HEOS:Water",
        )


class MixtureValidationTests(unittest.TestCase):
    def test_endpoint_validation_only_hard_gates_pure_identity(self):
        rows = [
            {
                "fraction_component_a": 0.0,
                "provider_components": ["Ethanol"],
                "storage_liquid": True,
                "state_points_requested": 1,
                "state_points_supported": 1,
                "storage_density_kg_m3": 790.0,
                "storage_cp_j_kg_k": 2400.0,
                "storage_viscosity_pa_s": 0.001,
                "storage_conductivity_w_m_k": 0.17,
                "bubble_temperature_k": 351.0,
                "dew_temperature_k": 351.0,
                "storage_enthalpy_j_kg": 0.0,
            },
            {
                "fraction_component_a": 0.5,
                "provider_components": ["Water", "Ethanol"],
                "storage_liquid": True,
                "state_points_requested": 1,
                "state_points_supported": 1,
                "storage_density_kg_m3": 900.0,
                "storage_cp_j_kg_k": 3000.0,
                "storage_viscosity_pa_s": None,
                "storage_conductivity_w_m_k": None,
                "bubble_temperature_k": 360.0,
                "dew_temperature_k": 370.0,
                "storage_enthalpy_j_kg": 1000.0,
            },
            {
                "fraction_component_a": 1.0,
                "provider_components": ["Water"],
                "storage_liquid": True,
                "state_points_requested": 1,
                "state_points_supported": 1,
                "storage_density_kg_m3": 998.0,
                "storage_cp_j_kg_k": 4180.0,
                "storage_viscosity_pa_s": 0.001,
                "storage_conductivity_w_m_k": 0.6,
                "bubble_temperature_k": 373.0,
                "dew_temperature_k": 373.0,
                "storage_enthalpy_j_kg": 2000.0,
            },
        ]
        report = continuity_report(
            rows,
            ["Water", "Ethanol"],
        )
        self.assertTrue(report["endpoint_pure_state_pass"])
        self.assertTrue(report["all_storage_states_liquid"])
        self.assertEqual(report["composition_count"], 3)
        self.assertIsNone(
            report["metrics"]["storage_viscosity_pa_s"].get(
                "hard_reject"
            )
        )


if __name__ == "__main__":
    unittest.main()
