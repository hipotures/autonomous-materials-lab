from __future__ import annotations

from dataclasses import dataclass, asdict
import math
from typing import Any, Protocol

import CoolProp.CoolProp as CP
from CoolProp import AbstractState


V5A_COOLPROP_GITREVISION = (
    "d5b0cfb51cd9a9343284cc5af8ebd6a8bd0eecc0"
)


def coolprop_runtime_info() -> dict[str, str]:
    return {
        "version": CP.get_global_param_string("version"),
        "gitrevision": CP.get_global_param_string("gitrevision"),
    }


def require_v5a_coolprop() -> None:
    info = coolprop_runtime_info()
    if info["gitrevision"] != V5A_COOLPROP_GITREVISION:
        raise RuntimeError(
            "V5a mixture benchmark requires pinned CoolProp dev revision "
            f"{V5A_COOLPROP_GITREVISION}; runtime is "
            f"version={info['version']} gitrevision={info['gitrevision']}. "
            "Install experiments/entry-evaluator/requirements-v5a.txt."
        )


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
        if value > 1e-12:
            cleaned.append((name, value))

    if not cleaned:
        raise ValueError("property provider fractions sum to zero")

    total = sum(value for _, value in cleaned)
    normalized = [(name, value / total) for name, value in cleaned]
    return (
        [name for name, _ in normalized],
        [value for _, value in normalized],
    )


class CoolPropHEOSProvider:
    """CoolProp HEOS pure-fluid / mixture property provider.

    The same interface is used for pure fluids and mixtures. Zero-fraction
    components are removed before creating the AbstractState, which makes the
    0/100 and 100/0 sweep endpoints true pure-fluid states.
    """

    def __init__(
        self,
        *,
        components: list[str],
        fractions: list[float],
        composition_basis: str = "mole",
        backend: str = "HEOS",
        stability_algorithm: int = 1,
    ):
        components, fractions = _validate_components(components, fractions)

        self.backend = str(backend).strip().upper()
        if self.backend != "HEOS":
            raise ValueError(
                "V5a currently supports property_provider.backend=HEOS only"
            )

        self.stability_algorithm = int(stability_algorithm)
        if self.stability_algorithm not in {0, 1}:
            raise ValueError(
                "property_provider.stability_algorithm must be 0 (legacy) "
                "or 1 (Michelsen)"
            )

        self.composition_basis = str(composition_basis).strip().lower()
        if self.composition_basis not in {"mole", "mass"}:
            raise ValueError(
                "property_provider.composition_basis must be mole or mass"
            )

        self.components = components
        self.input_fractions = fractions
        self._fluid_string = "&".join(components)
        self._state = AbstractState(self.backend, self._fluid_string)
        # Keep saturation/VLE work on a separate AbstractState.  CoolProp
        # mixture state objects can retain phase-envelope/cache information
        # that changes subsequent imposed-phase PT root selection.
        self._sat_state = AbstractState(self.backend, self._fluid_string)
        self._saturation_cache: dict[
            float, SaturationAtPressure
        ] = {}

        if len(components) > 1:
            if self.composition_basis == "mole":
                self._state.set_mole_fractions(fractions)
                self._sat_state.set_mole_fractions(fractions)
            else:
                self._state.set_mass_fractions(fractions)
                self._sat_state.set_mass_fractions(fractions)

        self.mole_fractions = (
            [1.0]
            if len(self.components) == 1
            else [
                float(value)
                for value in self._state.get_mole_fractions()
            ]
        )
        self.molar_mass_kg_mol = float(self._state.molar_mass())

        if len(self.components) == 1:
            self.identity = f"{self.backend}:{self.components[0]}"
        else:
            composition = "&".join(
                f"{name}[{fraction:.12g}]"
                for name, fraction in zip(
                    self.components,
                    self.mole_fractions,
                )
            )
            self.identity = f"{self.backend}:{composition}"

    @staticmethod
    def _finite(value: float, name: str) -> float:
        value = float(value)
        if not math.isfinite(value):
            raise ValueError(f"CoolProp returned non-finite {name}")
        return value

    def _set_stability_algorithm(self) -> None:
        if len(self.components) > 1:
            CP.set_config_int(
                CP.MIXTURE_STABILITY_ALGORITHM,
                self.stability_algorithm,
            )

    def _phase_from_saturation(
        self,
        temperature_k: float,
        pressure_pa: float,
    ) -> tuple[int | None, str]:
        if len(self.components) <= 1:
            return None, "pure-fluid-auto"

        saturation = self.saturation_at_pressure(pressure_pa)
        if not saturation.supported:
            return None, "mixture-auto-no-saturation"

        bubble = saturation.bubble_temperature_k
        dew = saturation.dew_temperature_k
        if bubble is None or dew is None:
            return None, "mixture-auto-incomplete-saturation"

        low = min(bubble, dew)
        high = max(bubble, dew)
        tolerance_k = 1.0e-6

        if temperature_k < low - tolerance_k:
            return int(CP.iphase_liquid), "bubble-dew-liquid"
        if temperature_k > high + tolerance_k:
            return int(CP.iphase_gas), "bubble-dew-gas"
        return None, "bubble-dew-two-phase"

    def _update_pt(self, temperature_k: float, pressure_pa: float) -> None:
        temperature_k = float(temperature_k)
        pressure_pa = float(pressure_pa)
        if (
            not math.isfinite(temperature_k)
            or temperature_k <= 0.0
            or not math.isfinite(pressure_pa)
            or pressure_pa <= 0.0
        ):
            raise ValueError("property state T and P must be finite and positive")

        self._set_stability_algorithm()
        imposed_phase, phase_source = self._phase_from_saturation(
            temperature_k,
            pressure_pa,
        )
        if imposed_phase is None:
            self._state.unspecify_phase()
        else:
            self._state.specify_phase(imposed_phase)

        try:
            self._state.update(
                CP.PT_INPUTS,
                pressure_pa,
                temperature_k,
            )
        except Exception as exc:
            self._state.unspecify_phase()
            raise ValueError(
                "CoolProp PT flash failed "
                f"provider={self.identity} "
                f"T={temperature_k:.9g} K "
                f"P={pressure_pa:.9g} Pa "
                f"stability_algorithm={self.stability_algorithm} "
                f"phase_source={phase_source}: {exc}"
            ) from exc

        if imposed_phase is not None:
            actual_phase = int(self._state.phase())
            if actual_phase != imposed_phase:
                self._state.unspecify_phase()
                raise ValueError(
                    "CoolProp ignored imposed mixture phase "
                    f"provider={self.identity} "
                    f"T={temperature_k:.9g} K "
                    f"P={pressure_pa:.9g} Pa "
                    f"requested_phase={_PHASE_NAMES.get(imposed_phase, imposed_phase)} "
                    f"actual_phase={_PHASE_NAMES.get(actual_phase, actual_phase)} "
                    f"stability_algorithm={self.stability_algorithm} "
                    f"phase_source={phase_source}"
                )

    def enthalpy_j_kg(
        self,
        temperature_k: float,
        pressure_pa: float,
    ) -> float:
        self._update_pt(temperature_k, pressure_pa)
        return self._finite(self._state.hmass(), "enthalpy")

    def phase(
        self,
        temperature_k: float,
        pressure_pa: float,
    ) -> str:
        self._update_pt(temperature_k, pressure_pa)
        return _PHASE_NAMES.get(
            int(self._state.phase()),
            f"phase_{int(self._state.phase())}",
        )

    @staticmethod
    def _optional(call) -> float | None:
        try:
            value = float(call())
        except Exception:
            return None
        return value if math.isfinite(value) else None

    def state(
        self,
        temperature_k: float,
        pressure_pa: float,
    ) -> ThermoState:
        self._update_pt(temperature_k, pressure_pa)
        phase_code = int(self._state.phase())
        return ThermoState(
            temperature_k=float(temperature_k),
            pressure_pa=float(pressure_pa),
            enthalpy_j_kg=self._finite(
                self._state.hmass(),
                "enthalpy",
            ),
            density_kg_m3=self._finite(
                self._state.rhomass(),
                "density",
            ),
            cp_j_kg_k=self._optional(self._state.cpmass),
            viscosity_pa_s=self._optional(self._state.viscosity),
            conductivity_w_m_k=self._optional(self._state.conductivity),
            phase=_PHASE_NAMES.get(
                phase_code,
                f"phase_{phase_code}",
            ),
        )

    def saturation_at_pressure(
        self,
        pressure_pa: float,
    ) -> SaturationAtPressure:
        pressure_pa = float(pressure_pa)
        if not math.isfinite(pressure_pa) or pressure_pa <= 0.0:
            raise ValueError("saturation pressure must be finite and positive")

        cache_key = pressure_pa
        cached = self._saturation_cache.get(cache_key)
        if cached is not None:
            return cached

        bubble = dew = None
        failures: list[str] = []

        self._set_stability_algorithm()
        self._sat_state.unspecify_phase()

        try:
            self._sat_state.update(CP.PQ_INPUTS, pressure_pa, 0.0)
            bubble = self._finite(
                self._sat_state.T(),
                "bubble temperature",
            )
        except Exception as exc:
            failures.append(f"bubble: {exc}")

        try:
            self._sat_state.update(CP.PQ_INPUTS, pressure_pa, 1.0)
            dew = self._finite(
                self._sat_state.T(),
                "dew temperature",
            )
        except Exception as exc:
            failures.append(f"dew: {exc}")

        result = SaturationAtPressure(
            pressure_pa=pressure_pa,
            bubble_temperature_k=bubble,
            dew_temperature_k=dew,
            supported=bubble is not None and dew is not None,
            failure_reason=(
                None
                if not failures
                else "; ".join(failures)
            ),
        )
        if len(self._saturation_cache) >= 4096:
            self._saturation_cache.clear()
        self._saturation_cache[cache_key] = result
        return result

    def metadata(self) -> dict[str, Any]:
        return {
            "provider_type": "coolprop",
            "backend": self.backend,
            "identity": self.identity,
            "components": list(self.components),
            "composition_basis": self.composition_basis,
            "stability_algorithm": self.stability_algorithm,
            "stability_algorithm_name": (
                "Michelsen"
                if self.stability_algorithm == 1
                else "legacy_Gernert"
            ),
            "input_fractions": list(self.input_fractions),
            "mole_fractions": list(self.mole_fractions),
            "molar_mass_kg_mol": self.molar_mass_kg_mol,
            "interaction_model": (
                "CoolProp HEOS built-in mixture reducing/departure model "
                "and binary interaction parameters"
                if len(self.components) > 1
                else "pure-fluid HEOS"
            ),
            "uncertainty_model": None,
            "phase_selection": (
                "bubble/dew guarded single-phase PT on pinned CoolProp dev; "
                "dedicated saturation AbstractState; automatic PT flash "
                "inside the two-phase envelope"
                if len(self.components) > 1
                else "CoolProp pure-fluid automatic phase selection"
            ),
            "coolprop_runtime": coolprop_runtime_info(),
            "project_validation_status": "V5a screening",
        }


def build_property_provider(
    coolant_config: dict[str, Any],
) -> CoolPropHEOSProvider:
    provider = coolant_config.get("property_provider")

    if provider is None:
        fluid = coolant_config.get("coolprop_name")
        if not fluid:
            raise ValueError(
                "coolant requires coolprop_name or property_provider"
            )
        return CoolPropHEOSProvider(
            components=[str(fluid)],
            fractions=[1.0],
            composition_basis="mole",
            backend="HEOS",
            stability_algorithm=1,
        )

    provider_type = str(
        provider.get("type", "coolprop")
    ).strip().lower()
    if provider_type != "coolprop":
        raise ValueError(
            f"unsupported property_provider.type: {provider_type}"
        )

    return CoolPropHEOSProvider(
        components=[
            str(value)
            for value in provider.get("components", [])
        ],
        fractions=[
            float(value)
            for value in provider.get("fractions", [])
        ],
        composition_basis=str(
            provider.get("composition_basis", "mole")
        ),
        backend=str(provider.get("backend", "HEOS")),
        stability_algorithm=int(
            provider.get("stability_algorithm", 1)
        ),
    )


def state_to_dict(state: ThermoState) -> dict[str, Any]:
    return asdict(state)


def saturation_to_dict(
    saturation: SaturationAtPressure,
) -> dict[str, Any]:
    return asdict(saturation)
