from __future__ import annotations

from dataclasses import dataclass, asdict
import math
from typing import Any, Protocol

import CoolProp.CoolProp as CP
from CoolProp import AbstractState


def coolprop_runtime_info() -> dict[str, str]:
    return {
        "version": CP.get_global_param_string("version"),
        "gitrevision": CP.get_global_param_string("gitrevision"),
    }


_PHASE_NAMES = {
    int(CP.iphase_liquid): "liquid",
    int(CP.iphase_gas): "gas",
    int(CP.iphase_twophase): "twophase",
    int(CP.iphase_supercritical): "supercritical",
    int(CP.iphase_supercritical_liquid): "supercritical_liquid",
    int(CP.iphase_supercritical_gas): "supercritical_gas",
    int(CP.iphase_critical_point): "critical_point",
    int(CP.iphase_unknown): "unknown",
}


@dataclass(frozen=True)
class ThermoState:
    temperature_k: float
    pressure_pa: float
    enthalpy_j_kg: float
    density_kg_m3: float
    cp_j_kg_k: float | None
    viscosity_pa_s: float | None
    conductivity_w_m_k: float | None
    phase: str


@dataclass(frozen=True)
class SaturationAtPressure:
    pressure_pa: float
    bubble_temperature_k: float | None
    dew_temperature_k: float | None
    supported: bool
    failure_reason: str | None


class PropertyProvider(Protocol):
    identity: str

    def enthalpy_j_kg(self, temperature_k: float, pressure_pa: float) -> float:
        ...

    def phase(self, temperature_k: float, pressure_pa: float) -> str:
        ...

    def state(self, temperature_k: float, pressure_pa: float) -> ThermoState:
        ...

    def saturation_at_pressure(self, pressure_pa: float) -> SaturationAtPressure:
        ...

    def metadata(self) -> dict[str, Any]:
        ...


def _validate_components(
    components: list[str],
    fractions: list[float],
) -> tuple[list[str], list[float]]:
    if not components:
        raise ValueError("property_provider.components must not be empty")
    if len(components) != len(fractions):
        raise ValueError(
            "property_provider.components and fractions must have equal length"
        )
    components = [str(c).strip() for c in components]
    if len(set(components)) != len(components):
        raise ValueError("property_provider.components must be unique")

    cleaned: list[tuple[str, float]] = []
    for component, fraction in zip(components, fractions):
        name = str(component).strip()
        value = float(fraction)
        if not name:
            raise ValueError("property provider component names must be non-empty")
        if not math.isfinite(value) or value < 0.0:
            raise ValueError(
                "property provider fractions must be finite and nonnegative"
            )
        if value > 0.0:
            cleaned.append((name, value))

    if not cleaned:
        raise ValueError("property provider fractions sum to zero")

    total = sum(value for _, value in cleaned)
    if not math.isfinite(total):
        raise ValueError("property provider fractions have non-finite sum")
    normalized = [(name, value / total) for name, value in cleaned]
    return (
        [name for name, _ in normalized],
        [value for _, value in normalized],
    )


class CoolPropPureProvider:
    """Pure-fluid HEOS only. Multicomponent HEOS is deliberately unsupported."""

    def __init__(self, *, components: list[str], fractions: list[float],
                 composition_basis: str = "mole", backend: str = "HEOS"):
        self.components, self.input_fractions = _validate_components(components, fractions)
        if len(self.components) != 1:
            raise ValueError("HEOS mixtures are disabled; select type=thermo, model=NRTL for Water/Ethanol")
        if backend.upper() != "HEOS" or composition_basis not in {"mole", "mass"}:
            raise ValueError("pure provider requires HEOS and mole or mass composition basis")
        self.composition_basis = composition_basis
        self.identity = f"HEOS:{self.components[0]}"
        self._state = AbstractState("HEOS", self.components[0])
        self._sat_state = AbstractState("HEOS", self.components[0])
        self.mole_fractions = [1.0]
        self.molar_mass_kg_mol = float(self._state.molar_mass())

    @staticmethod
    def _finite(value: float, name: str) -> float:
        value = float(value)
        if not math.isfinite(value):
            raise ValueError(f"property backend returned non-finite {name}")
        return value

    @staticmethod
    def _optional(call) -> float | None:
        try:
            value = float(call())
        except Exception:
            return None
        return value if math.isfinite(value) and value > 0 else None

    def _update_pt(self, temperature_k: float, pressure_pa: float) -> None:
        if any(not math.isfinite(v) or v <= 0 for v in (temperature_k, pressure_pa)):
            raise ValueError("property state T and P must be finite and positive")
        self._state.update(CP.PT_INPUTS, pressure_pa, temperature_k)

    def enthalpy_j_kg(self, temperature_k: float, pressure_pa: float) -> float:
        self._update_pt(temperature_k, pressure_pa)
        return self._finite(self._state.hmass(), "enthalpy")

    def phase(self, temperature_k: float, pressure_pa: float) -> str:
        self._update_pt(temperature_k, pressure_pa)
        return _PHASE_NAMES.get(int(self._state.phase()), "unknown")

    def state(self, temperature_k: float, pressure_pa: float) -> ThermoState:
        self._update_pt(temperature_k, pressure_pa)
        return ThermoState(
            temperature_k, pressure_pa, self._finite(self._state.hmass(), "enthalpy"),
            self._finite(self._state.rhomass(), "density"),
            self._optional(self._state.cpmass), self._optional(self._state.viscosity),
            self._optional(self._state.conductivity),
            _PHASE_NAMES.get(int(self._state.phase()), "unknown"),
        )

    def saturation_at_pressure(self, pressure_pa: float) -> SaturationAtPressure:
        if not math.isfinite(pressure_pa) or pressure_pa <= 0:
            raise ValueError("saturation pressure must be finite and positive")
        try:
            self._sat_state.update(CP.PQ_INPUTS, pressure_pa, 0)
            t = self._finite(self._sat_state.T(), "saturation temperature")
            return SaturationAtPressure(pressure_pa, t, t, True, None)
        except Exception as exc:
            return SaturationAtPressure(pressure_pa, None, None, False, str(exc))

    def metadata(self) -> dict[str, Any]:
        return {"provider_type": "coolprop", "backend": "HEOS", "identity": self.identity,
                "components": self.components, "composition_basis": self.composition_basis,
                "input_fractions": self.input_fractions, "mole_fractions": [1.0],
                "molar_mass_kg_mol": self.molar_mass_kg_mol,
                "interaction_model": "pure-fluid HEOS", "coolprop_runtime": coolprop_runtime_info(),
                "uncertainty_model": None, "project_validation_status": "pure-fluid reference"}


# Compatibility name for callers of the old provider; still rejects mixtures.
CoolPropHEOSProvider = CoolPropPureProvider


def build_property_provider(coolant_config: dict[str, Any]) -> PropertyProvider:
    provider = coolant_config.get("property_provider")
    if provider is None:
        fluid = coolant_config.get("coolprop_name")
        if not fluid:
            raise ValueError("coolant requires coolprop_name or property_provider")
        return CoolPropPureProvider(components=[str(fluid)], fractions=[1.0])
    components, fractions = _validate_components(
        list(provider.get("components", [])), list(provider.get("fractions", [])))
    basis = str(provider.get("composition_basis", "mole")).lower()
    kind = str(provider.get("type", "coolprop")).lower()
    if kind not in {"coolprop", "thermo"}:
        raise ValueError(f"unsupported property_provider.type: {kind}")
    if kind == "thermo" and provider.get("model", "NRTL") != "NRTL":
        raise ValueError("only NRTL is implemented for thermo mixtures")
    if len(components) == 1 or kind == "coolprop":
        return CoolPropPureProvider(components=components, fractions=fractions,
                                   composition_basis=basis, backend=provider.get("backend", "HEOS"))
    # Keep the pure-fluid runtime independent of optional thermo dependencies.
    from thermo_provider import ThermoMixtureProvider
    return ThermoMixtureProvider(components=components, fractions=fractions,
                                composition_basis=basis,
                                parameter_set=provider.get("parameter_set", "ddbst_p05_01b"))


def state_to_dict(state: ThermoState) -> dict[str, Any]:
    return asdict(state)


def saturation_to_dict(
    saturation: SaturationAtPressure,
) -> dict[str, Any]:
    return asdict(saturation)
