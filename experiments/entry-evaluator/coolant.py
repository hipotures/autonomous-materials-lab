from __future__ import annotations

from dataclasses import dataclass
import math
from typing import Any

import CoolProp.CoolProp as CP

from chemistry import ChemistryLimiter


@dataclass(frozen=True)
class CoolantStep:
    exit_temperature_k: float
    usable_enthalpy_j_kg: float
    mass_flux_kg_m2_s: float
    mass_flow_kg_s: float
    required_injection_pressure_pa: float
    chemistry_source: str
    ignition_delay_s: float | None
    feasible: bool
    failure_reason: str | None


class CoolantModel:
    def __init__(
        self,
        config: dict[str, Any],
        chemistry: ChemistryLimiter,
    ):
        self.config = config
        self.chemistry = chemistry
        self.fluid = str(config["coolprop_name"])
        self.storage_temperature_k = float(
            config["storage_temperature_k"]
        )
        self.storage_pressure_pa = float(
            config["storage_pressure_pa"]
        )
        self.nominal_max_exit_temperature_k = float(
            config["max_exit_temperature_k"]
        )
        self.approach_delta_k = float(
            config.get("wall_to_fluid_approach_k", 25.0)
        )
        self.cooled_area_m2 = float(config["cooled_area_m2"])
        self.pressure_margin = float(
            config.get("injection_pressure_margin", 1.15)
        )
        self.porous_delta_p_pa = float(
            config.get("porous_delta_p_pa", 0.0)
        )
        self.max_injection_pressure_pa = float(
            config.get("max_injection_pressure_pa", 1e9)
        )
        self._cache: dict[tuple[float, float], float] = {}

        self.h_storage = float(
            CP.PropsSI(
                "H",
                "T",
                self.storage_temperature_k,
                "P",
                self.storage_pressure_pa,
                self.fluid,
            )
        )
        phase = CP.PhaseSI(
            "T",
            self.storage_temperature_k,
            "P",
            self.storage_pressure_pa,
            self.fluid,
        )
        if phase not in {"liquid", "supercritical_liquid"}:
            raise ValueError(
                "configured coolant storage state is not liquid: "
                f"{self.fluid} phase={phase}"
            )

    def _h_out(
        self,
        temperature_k: float,
        pressure_pa: float,
    ) -> float:
        key = (
            round(temperature_k, 1),
            round(math.log10(max(pressure_pa, 1.0)), 3),
        )
        if key not in self._cache:
            self._cache[key] = float(
                CP.PropsSI(
                    "H",
                    "T",
                    temperature_k,
                    "P",
                    max(pressure_pa, 100.0),
                    self.fluid,
                )
            )
        return self._cache[key]

    def evaluate(
        self,
        wall_temperature_k: float,
        surface_pressure_pa: float,
        coolant_heat_flux_w_m2: float,
    ) -> CoolantStep:
        if coolant_heat_flux_w_m2 <= 0.0:
            required_injection_pressure = (
                self.pressure_margin * surface_pressure_pa
                + self.porous_delta_p_pa
            )
            return CoolantStep(
                exit_temperature_k=self.storage_temperature_k,
                usable_enthalpy_j_kg=0.0,
                mass_flux_kg_m2_s=0.0,
                mass_flow_kg_s=0.0,
                required_injection_pressure_pa=required_injection_pressure,
                chemistry_source="inactive",
                ignition_delay_s=None,
                feasible=True,
                failure_reason=None,
            )

        chemistry_limit = self.chemistry.limit(
            surface_pressure_pa,
            self.nominal_max_exit_temperature_k,
        )
        requested_exit = min(
            chemistry_limit.max_exit_temperature_k,
            wall_temperature_k - self.approach_delta_k,
        )
        requested_exit = max(
            self.storage_temperature_k + 1.0,
            requested_exit,
        )

        required_injection_pressure = (
            self.pressure_margin * surface_pressure_pa
            + self.porous_delta_p_pa
        )
        if required_injection_pressure > self.max_injection_pressure_pa:
            return CoolantStep(
                requested_exit,
                0.0,
                0.0,
                0.0,
                required_injection_pressure,
                chemistry_limit.source,
                chemistry_limit.estimated_ignition_delay_s,
                False,
                "required injection pressure exceeds configured limit",
            )

        try:
            h_out = self._h_out(
                requested_exit,
                surface_pressure_pa,
            )
        except Exception as exc:
            return CoolantStep(
                requested_exit,
                0.0,
                0.0,
                0.0,
                required_injection_pressure,
                chemistry_limit.source,
                chemistry_limit.estimated_ignition_delay_s,
                False,
                f"CoolProp outlet state failed: {exc}",
            )

        delta_h = h_out - self.h_storage
        if delta_h <= 0.0:
            return CoolantStep(
                requested_exit,
                delta_h,
                0.0,
                0.0,
                required_injection_pressure,
                chemistry_limit.source,
                chemistry_limit.estimated_ignition_delay_s,
                False,
                "usable coolant enthalpy is non-positive",
            )

        mass_flux = max(0.0, coolant_heat_flux_w_m2) / delta_h
        mass_flow = mass_flux * self.cooled_area_m2
        return CoolantStep(
            requested_exit,
            delta_h,
            mass_flux,
            mass_flow,
            required_injection_pressure,
            chemistry_limit.source,
            chemistry_limit.estimated_ignition_delay_s,
            True,
            None,
        )
