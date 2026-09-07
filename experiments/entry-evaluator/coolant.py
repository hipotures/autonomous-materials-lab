from __future__ import annotations

from dataclasses import dataclass
import math
from typing import Any

import CoolProp.CoolProp as CP
from CoolProp import AbstractState

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
    required_ignition_delay_s: float | None = None
    ignition_margin: float | None = None


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
        self._outlet_state = AbstractState("HEOS", self.fluid)

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
        # Exact keys avoid first-visitor rounding bias between time steps/rings.
        key = (float(temperature_k), max(float(pressure_pa), 100.0))
        if key not in self._cache:
            self._outlet_state.update(CP.PT_INPUTS, key[1], key[0])
            value = float(self._outlet_state.hmass())
            if not math.isfinite(value):
                raise ValueError("non-finite outlet enthalpy")
            if len(self._cache) >= 4096:
                self._cache.clear()
            self._cache[key] = value
        return self._cache[key]

    def evaluate(
        self,
        wall_temperature_k: float,
        surface_pressure_pa: float,
        coolant_heat_flux_w_m2: float,
        *,
        area_m2: float | None = None,
    ) -> CoolantStep:
        area = self.cooled_area_m2 if area_m2 is None else float(area_m2)
        if not math.isfinite(area) or area < 0:
            raise ValueError("coolant evaluation area must be finite and nonnegative")
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
        if not chemistry_limit.feasible:
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
                chemistry_source=chemistry_limit.source,
                ignition_delay_s=None,
                feasible=False,
                failure_reason=chemistry_limit.failure_reason,
                required_ignition_delay_s=(
                    chemistry_limit.required_ignition_delay_s
                ),
                ignition_margin=None,
            )

        requested_exit = min(
            chemistry_limit.max_exit_temperature_k,
            wall_temperature_k - self.approach_delta_k,
        )
        requested_exit = max(
            self.storage_temperature_k + 1.0,
            requested_exit,
        )

        actual_ignition_delay_s = chemistry_limit.estimated_ignition_delay_s
        if chemistry_limit.source == "ignition_csv":
            try:
                actual_ignition_delay_s = self.chemistry.estimate_delay(
                    surface_pressure_pa,
                    requested_exit,
                )
            except ValueError:
                actual_ignition_delay_s = None
        required_ignition_delay_s = chemistry_limit.required_ignition_delay_s
        ignition_margin = (
            actual_ignition_delay_s / required_ignition_delay_s
            if actual_ignition_delay_s is not None
            and required_ignition_delay_s is not None
            and required_ignition_delay_s > 0.0
            else None
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
                actual_ignition_delay_s,
                False,
                "required injection pressure exceeds configured limit",
                required_ignition_delay_s,
                ignition_margin,
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
                actual_ignition_delay_s,
                False,
                f"CoolProp outlet state failed: {exc}",
                required_ignition_delay_s,
                ignition_margin,
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
                actual_ignition_delay_s,
                False,
                "usable coolant enthalpy is non-positive",
                required_ignition_delay_s,
                ignition_margin,
            )

        mass_flux = max(0.0, coolant_heat_flux_w_m2) / delta_h
        mass_flow = mass_flux * area
        return CoolantStep(
            requested_exit,
            delta_h,
            mass_flux,
            mass_flow,
            required_injection_pressure,
            chemistry_limit.source,
            actual_ignition_delay_s,
            True,
            None,
            required_ignition_delay_s,
            ignition_margin,
        )
