from __future__ import annotations

import math
from pathlib import Path
import sys
import unittest
from unittest.mock import patch

import CoolProp.CoolProp as CP

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from chemistry import ChemistryLimiter
from coolant import CoolantModel
from property_provider import CoolPropPureProvider, build_property_provider
from run_mixture_sweep import continuity_report, entry_summary, preflight_failure_reason, run_entry_task
from thermo_provider import ThermoMixtureProvider, nrtl_model


def mixture(x_water=0.5):
    return ThermoMixtureProvider(components=["Water", "Ethanol"], fractions=[x_water, 1-x_water])


class PropertyProviderTests(unittest.TestCase):
    def test_legacy_pure_fluid_and_factory_endpoints_match_propssi(self):
        for config, fluid in [({"coolprop_name": "Water"}, "Water"),
                              ({"property_provider": {"type": "thermo", "model": "NRTL",
                                "components": ["Water", "Ethanol"], "fractions": [0, 1]}}, "Ethanol"),
                              ({"property_provider": {"type": "thermo", "model": "NRTL",
                                "components": ["Water", "Ethanol"], "fractions": [1, 0]}}, "Water")]:
            provider = build_property_provider(config)
            self.assertIsInstance(provider, CoolPropPureProvider)
            expected = CP.PropsSI("H", "T", 350, "P", 2e5, fluid)
            self.assertAlmostEqual(provider.enthalpy_j_kg(350, 2e5), expected, delta=abs(expected)*1e-10)
            self.assertEqual(provider.metadata()["components"], [fluid])

    def test_heos_mixtures_and_invalid_inputs_rejected(self):
        with self.assertRaisesRegex(ValueError, "HEOS mixtures are disabled"):
            CoolPropPureProvider(components=["Water", "Ethanol"], fractions=[.5, .5])
        for fractions in ([.5], [-.1, 1.1], [math.nan, 1], [0, 0]):
            with self.assertRaises(ValueError):
                ThermoMixtureProvider(components=["Water", "Ethanol"], fractions=fractions)
        with self.assertRaises(ValueError):
            ThermoMixtureProvider(components=["Water", "Methanol"], fractions=[.5, .5])
        for t, p in ((501, 1e5), (270, 1e5), (300, 2.01e6), (300, 99), (math.nan, 1e5)):
            with self.assertRaises(ValueError):
                mixture().enthalpy_j_kg(t, p)

    def test_documented_nrtl_gamma_and_excess_enthalpy(self):
        # thermo NRTL docs / DDBST P05.01b: order ethanol, water; values are J/mol.
        model = nrtl_model(343.15, [.252, .748])
        for actual, expected in zip(model.gammas(), [1.9360516514475439, 1.1536630452052914]):
            self.assertAlmostEqual(actual, expected, places=10)
        self.assertAlmostEqual(model.HE(), 582.9648539351417, places=7)
        dT = .001
        derivative = (nrtl_model(343.15+dT, [.252, .748]).GE()
                      - nrtl_model(343.15-dT, [.252, .748]).GE())/(2*dT)
        self.assertAlmostEqual(model.HE(), model.GE()-343.15*derivative, delta=1e-5)

    def test_bubble_flash_obeys_modified_raoult_law(self):
        provider = mixture(.748)
        t = 343.15
        gammas = nrtl_model(t, [.252, .748]).gammas()
        partial = [z*g*vp(t) for z, g, vp in zip([.252, .748], gammas, provider.correlations.VaporPressures)]
        bubble = provider.flasher.flash(T=t, VF=0, zs=provider.zs)
        self.assertAlmostEqual(bubble.P, sum(partial), delta=.01)
        for y, expected in zip(bubble.gas.zs, [v/sum(partial) for v in partial]):
            self.assertAlmostEqual(y, expected, places=7)

    def test_liquid_enthalpy_includes_excess_contribution(self):
        from thermo import GibbsExcessLiquid, IdealSolution
        p = mixture(.748)
        t, pressure = 300., 101325.
        ideal = GibbsExcessLiquid(
            VaporPressures=p.correlations.VaporPressures,
            HeatCapacityGases=p.correlations.HeatCapacityGases,
            VolumeLiquids=p.correlations.VolumeLiquids,
            GibbsExcessModel=IdealSolution(T=t, xs=p.zs),
            equilibrium_basis="Psat", caloric_basis="Psat", T=t, P=pressure, zs=p.zs,
        )
        h_actual = p.enthalpy_j_kg(t, pressure)*p.molar_mass_kg_mol
        self.assertAlmostEqual(h_actual-ideal.H(), nrtl_model(t, p.zs).HE(), places=7)

    def test_component_order_and_mass_basis(self):
        original = mixture(.4)
        swapped = ThermoMixtureProvider(components=["Ethanol", "Water"], fractions=[.6, .4])
        by_mass = ThermoMixtureProvider(components=["Water", "Ethanol"],
                                       fractions=[.4*18.01528, .6*46.06844], composition_basis="mass")
        for p in (swapped, by_mass):
            for t in (293.15, 500):
                self.assertAlmostEqual(p.enthalpy_j_kg(t, 101325), original.enthalpy_j_kg(t, 101325), places=6)

    def test_problem_mid_compositions_and_call_order(self):
        for x in (.4, .5):
            provider = mixture(x)
            h_in = provider.enthalpy_j_kg(293.15, 101325)
            hot = provider.state(500, 25000)
            self.assertEqual(hot.phase, "gas")
            self.assertLess(hot.density_kg_m3, 1)
            self.assertTrue(1e6 < hot.enthalpy_j_kg-h_in < 3e6)
            for p in (25000, 30000, 101325):
                self.assertTrue(provider.saturation_at_pressure(p).supported)
            self.assertEqual(hot, provider.state(500, 25000))
            self.assertIsNone(hot.viscosity_pa_s)

    def test_single_phase_cp_is_enthalpy_derivative(self):
        p = mixture()
        for t, pressure in ((293.15, 101325), (450, 25000)):
            step = .001
            derivative = (p.enthalpy_j_kg(t+step, pressure)-p.enthalpy_j_kg(t-step, pressure))/(2*step)
            # Psat correlation derivative implementations close to 1e-4 relative.
            self.assertAlmostEqual(p.state(t, pressure).cp_j_kg_k, derivative, delta=abs(derivative)*1e-4)

    def test_full_two_phase_flash_uses_bulk_enthalpy_and_closes_balances(self):
        p = mixture()
        sat = p.saturation_at_pressure(101325)
        t = (sat.bubble_temperature_k+sat.dew_temperature_k)/2
        state = p.state(t, 101325)
        result = p._cached_flash(t, 101325)
        self.assertEqual(state.phase, "twophase")
        self.assertIsNone(state.cp_j_kg_k)
        self.assertTrue(0 < result.VF < 1)
        self.assertAlmostEqual(sum(result.betas), 1, places=12)
        for i, z in enumerate(p.zs):
            self.assertAlmostEqual(sum(b*phase.zs[i] for b, phase in zip(result.betas, result.phases)), z, places=7)
        h_molar = sum(b*phase.H() for b, phase in zip(result.betas, result.phases))
        self.assertAlmostEqual(state.enthalpy_j_kg, h_molar/p.molar_mass_kg_mol, places=6)
        self.assertAlmostEqual(state.enthalpy_j_kg, result.H_mass(), places=6)

    def test_mixture_coolant_energy_balance(self):
        model = CoolantModel({"storage_temperature_k": 293.15, "storage_pressure_pa": 101325,
                              "max_exit_temperature_k": 500, "wall_to_fluid_approach_k": 25,
                              "cooled_area_m2": 1, "property_provider": {"type": "thermo", "model": "NRTL",
                              "components": ["Water", "Ethanol"], "fractions": [.5, .5]}},
                             ChemistryLimiter({"mode": "disabled"}))
        step = model.evaluate(wall_temperature_k=525, surface_pressure_pa=25000,
                              coolant_heat_flux_w_m2=1e6, area_m2=1)
        self.assertTrue(step.feasible, step.failure_reason)
        self.assertAlmostEqual(step.mass_flow_kg_s*step.usable_enthalpy_j_kg, 1e6, places=6)
        self.assertIn("caloric validation pending", model.properties.metadata()["project_validation_status"])


class MixtureValidationTests(unittest.TestCase):
    def test_preflight_failure_is_specific_to_candidate(self):
        good = {"storage_liquid": True, "state_points_requested": 20, "state_points_supported": 20}
        bad_grid = {**good, "state_points_supported": 19}
        bad_storage = {**good, "storage_liquid": False}
        self.assertIsNone(preflight_failure_reason(good))
        self.assertIn("grid", preflight_failure_reason(bad_grid))
        self.assertIn("storage", preflight_failure_reason(bad_storage))

    def test_missing_composition_is_not_bridged(self):
        rows = [{"fraction_component_a": x, "provider_components": comps,
                 "storage_liquid": True, "state_points_requested": 1, "state_points_supported": 1,
                 "storage_density_kg_m3": density} for x, comps, density in (
                     (0, ["Ethanol"], 790), (.5, ["Water", "Ethanol"], None), (1, ["Water"], 998))]
        report = continuity_report(rows, ["Water", "Ethanol"])
        self.assertTrue(report["endpoint_pure_state_pass"])
        metric = report["metrics"]["storage_density_kg_m3"]
        self.assertEqual(metric["adjacent_changes"], [])
        self.assertEqual(len(metric["unsupported_intervals"]), 2)
        self.assertIsNone(metric["max_adjacent_relative_change"])

    def test_partial_nonfinite_and_zero_masses_cannot_win(self):
        rows = [{"candidate_id": str(x), "fraction_component_a": x,
                 "coolant_used_kg": mass, "status": status} for x, mass, status in (
                     (1, 3484.6, "terminal_velocity"), (.9, 3924.3, "terminal_velocity"),
                     (.4, 5, "coolant_failure"), (.5, 2415, "exception"),
                     (.6, math.nan, "terminal_velocity"), (.7, math.inf, "terminal_velocity"),
                     (.8, 0, "terminal_velocity"), (.3, 100, "atmospheric_exit"))]
        report = entry_summary(rows, 1)
        self.assertEqual(report["comparable_candidate_count"], 2)
        self.assertEqual(report["best_candidate"]["fraction_component_a"], 1)
        for row in rows[2:]:
            self.assertIsNone(row["score_coolant_kg"])
            self.assertIsNone(row["coolant_mass_ratio_vs_reference"])

    def test_failure_or_different_termination_does_not_get_reference_ratio(self):
        rows = [{"candidate_id": "water", "fraction_component_a": 1, "status": "exception", "coolant_used_kg": 5},
                {"candidate_id": "mix", "fraction_component_a": .5, "status": "terminal_velocity", "coolant_used_kg": 100}]
        self.assertIsNone(entry_summary(rows, 1)["best_candidate"])
        rows[0].update(status="terminal_altitude", coolant_used_kg=200)
        report = entry_summary(rows, 1)
        self.assertEqual(report["comparable_candidate_count"], 1)
        self.assertIsNone(rows[1]["coolant_mass_ratio_vs_reference"])
        with patch("run_mixture_sweep.evaluate", side_effect=ValueError("backend failed")):
            row = run_entry_task({"config": {}, "candidate_id": "bad", "fraction_component_a": .4})
        self.assertEqual(row["status"], "exception")
        self.assertIsNone(entry_summary([row], .4)["best_candidate"])


if __name__ == "__main__":
    unittest.main()
