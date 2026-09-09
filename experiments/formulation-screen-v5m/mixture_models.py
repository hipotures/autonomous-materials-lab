"""Separate cold-liquid correlations from the provisional ethanol/water flash."""
from __future__ import annotations

from importlib.metadata import version
import math
from pathlib import Path
import sys
from formulations import PROFILES, COLD_FIELDS, number

ENTRY = Path(__file__).resolve().parents[1] / "entry-evaluator"


class MixtureModels:
    def __init__(self, cp=None, nrtl_factory=None):
        if cp is None:
            if version("CoolProp") != "8.0.0":
                raise RuntimeError("use the existing V5b environment with CoolProp==8.0.0")
            import CoolProp.CoolProp as cp
        self.cp = cp
        self.nrtl_factory = nrtl_factory
        self.cache = {}

    def cold(self, formulation: dict, t: float, pressure: float, *, fit_endpoint_probe: bool = False) -> dict:
        t, pressure = number(t, "temperature"), number(pressure, "pressure")
        a, w = formulation["additive_id"], formulation["additive_mass_fraction"]
        row = {"formulation_id": formulation["formulation_id"], "temperature_k": t,
               "pressure_pa": pressure, "additive_mass_fraction": w,
               "surface_tension_n_m": None, "contact_angle_deg": None,
               "evidence_kind": "model_not_measurement", "boiling_supported": False,
               "errors": {}, "fit_endpoint_probe": fit_endpoint_probe}
        row.update({k: None for k in COLD_FIELDS})
        if a == "water":
            fluid = "HEOS::Water"
        elif a in PROFILES:
            code, tmin, tmax, xmax = PROFILES[a]
            fraction_ok = 0 < w <= xmax or (fit_endpoint_probe and w == 0)
            if not tmin <= t <= tmax or not fraction_ok:
                return {**row, "backend": code, "status": "outside_model_range"}
            fluid = f"INCOMP::{code}[{w:.12g}]"
        else:
            return {**row, "backend": None, "status": "measurement_required"}
        row["backend"] = fluid
        row["phase_basis"] = "single_phase_liquid_library_fit"
        if a == "water":
            try:
                if self.cp.PhaseSI("T", t, "P", pressure, fluid) != "liquid":
                    return {**row, "status": "not_liquid"}
            except (ValueError, RuntimeError) as exc:
                return {**row, "status": "state_failed", "errors": {"phase": str(exc)}}
        for name, output in COLD_FIELDS.items():
            try:
                value = float(self.cp.PropsSI(output, "T", t, "P", pressure, fluid))
                if not math.isfinite(value) or (name != "enthalpy_j_kg" and value <= 0):
                    raise ValueError("nonfinite_or_nonpositive_property")
                row[name] = value
            except (ValueError, RuntimeError, OverflowError) as exc:
                row["errors"][name] = str(exc)
        row["status"] = "ok" if not row["errors"] else "property_failed"
        return row

    def _nrtl(self, formulation):
        key = formulation["formulation_id"]
        if key not in self.cache:
            factory = self.nrtl_factory
            if factory is None:
                if str(ENTRY) not in sys.path:
                    sys.path.insert(0, str(ENTRY))
                from thermo_provider import ThermoMixtureProvider
                factory = ThermoMixtureProvider
            w = formulation["additive_mass_fraction"]
            self.cache[key] = factory(components=["Water", "Ethanol"], fractions=[1 - w, w],
                                      composition_basis="mass")
        return self.cache[key]

    def hot(self, formulation, inlet, outlet, pressure):
        a = formulation["additive_id"]
        row = {"formulation_id": formulation["formulation_id"], "inlet_temperature_k": inlet,
               "outlet_temperature_k": outlet, "pressure_pa": pressure, "delta_h_j_kg": None,
               "backend": None, "evidence_kind": "provisional_model_not_measurement",
               "fixed_overall_composition": True, "porous_flow_simulated": False,
               "caloric_validation_complete": False, "sequential_evaporation_simulated": False}
        if a not in {"water", "ethanol"}:
            return {**row, "status": "unsupported_phase_change_model"}
        try:
            if a == "water":
                fluid = "HEOS::Water"
                row["backend"] = fluid
                h0 = float(self.cp.PropsSI("H", "T", inlet, "P", pressure, fluid))
                h1 = float(self.cp.PropsSI("H", "T", outlet, "P", pressure, fluid))
                row["outlet_phase"] = self.cp.PhaseSI("T", outlet, "P", pressure, fluid)
            else:
                row["backend"] = "V5a_NRTL_ddbst_p05_01b"
                provider = self._nrtl(formulation)
                h0 = provider.enthalpy_j_kg(inlet, pressure)
                h1 = provider.enthalpy_j_kg(outlet, pressure)
                row["outlet_phase"] = provider.phase(outlet, pressure)
                sat = provider.saturation_at_pressure(pressure)
                row.update(bubble_temperature_k=sat.bubble_temperature_k,
                           dew_temperature_k=sat.dew_temperature_k,
                           saturation_supported=sat.supported)
            if not math.isfinite(h0) or not math.isfinite(h1):
                raise ValueError("nonfinite_enthalpy")
            row["delta_h_j_kg"] = number(h1 - h0, "same_backend_delta_h")
            row["status"] = "ok"
        except (ValueError, RuntimeError, OverflowError) as exc:
            row.update(status="model_failed", error=str(exc))
        return row
