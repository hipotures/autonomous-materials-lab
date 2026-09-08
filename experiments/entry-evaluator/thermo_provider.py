"""Water/ethanol NRTL + full vapor/liquid flash for a bounded screening benchmark."""
from __future__ import annotations

from functools import lru_cache
from importlib.metadata import version
import math

from thermo import ChemicalConstantsPackage, FlashVL, GibbsExcessLiquid, IdealGas, NRTL
from scipy.constants import R, calorie

from property_provider import ThermoState, SaturationAtPressure, _validate_components

PARAMETER_SET = "ddbst_p05_01b"
PARAMETER_SOURCE = "https://thermo.readthedocs.io/thermo.nrtl.html"
# The documented order is ETHANOL, WATER; tau_ij = B_ij / T, not B_ji / T.
TAU_B_K = [[0.0, -121.2691 * calorie / R], [1337.8574 * calorie / R, 0.0]]
ALPHA = [[0.0, 0.2974], [0.2974, 0.0]]


def nrtl_model(temperature_k: float, mole_fractions: list[float]) -> NRTL:
    return NRTL(T=temperature_k, xs=mole_fractions, tau_bs=TAU_B_K, alpha_cs=ALPHA)


class ThermoMixtureProvider:
    """VLE example reproduction, not a calorimetrically validated parameter set.

    Both endpoint enthalpies use the same flash model. Pure composition endpoints
    in the public factory route to CoolProp; no raw cross-backend h subtraction.
    """

    def __init__(self, *, components: list[str], fractions: list[float],
                 composition_basis: str = "mole", parameter_set: str = PARAMETER_SET):
        if version("thermo") != "0.6.1":
            raise RuntimeError("V5a NRTL benchmark requires thermo==0.6.1")
        components, fractions = _validate_components(components, fractions)
        if set(components) != {"Water", "Ethanol"} or len(components) != 2:
            raise ValueError("thermo NRTL provider currently supports only the Water/Ethanol binary")
        if composition_basis not in {"mole", "mass"}:
            raise ValueError("composition_basis must be mole or mass")
        if parameter_set != PARAMETER_SET:
            raise ValueError(f"unsupported NRTL parameter set: {parameter_set}")
        self.components = components
        self.input_fractions = fractions
        self.composition_basis = composition_basis
        self.constants, self.correlations = ChemicalConstantsPackage.from_IDs(["ethanol", "water"])
        ordered = [fractions[components.index(c)] for c in ("Ethanol", "Water")]
        if composition_basis == "mass":
            ordered = [w / mw for w, mw in zip(ordered, self.constants.MWs)]
        total = sum(ordered)
        self.zs = [z / total for z in ordered]
        self.mole_fractions = [self.zs[["Ethanol", "Water"].index(c)] for c in components]
        self.molar_mass_kg_mol = sum(z * mw for z, mw in zip(self.zs, self.constants.MWs)) / 1000
        self.identity = "thermo:NRTL:" + "&".join(
            f"{c}[{z:.12g}]" for c, z in zip(("Ethanol", "Water"), self.zs))
        props = self.correlations
        self.liquid = GibbsExcessLiquid(
            VaporPressures=props.VaporPressures, HeatCapacityGases=props.HeatCapacityGases,
            VolumeLiquids=props.VolumeLiquids, GibbsExcessModel=nrtl_model(298.15, self.zs),
            equilibrium_basis="Psat", caloric_basis="Psat",
        )
        self.gas = IdealGas(HeatCapacityGases=props.HeatCapacityGases)
        self.flasher = FlashVL(self.constants, props, liquid=self.liquid, gas=self.gas)
        # Per-instance bounded cache: no approximate states or phase imposition.
        self._cached_flash = lru_cache(maxsize=256)(self._flash)

    @staticmethod
    def _validate_domain(temperature_k: float | None, pressure_pa: float) -> None:
        # Operational screening bounds, NOT the published fit validity interval.
        if not math.isfinite(pressure_pa) or not 100 <= pressure_pa <= 2e6:
            raise ValueError("NRTL screening pressure domain is 100..2000000 Pa")
        if temperature_k is not None and (not math.isfinite(temperature_k) or not 273.15 <= temperature_k <= 500):
            raise ValueError("NRTL screening temperature domain is 273.15..500 K")

    def _flash(self, temperature_k: float, pressure_pa: float):
        self._validate_domain(temperature_k, pressure_pa)
        result = self.flasher.flash(T=temperature_k, P=pressure_pa, zs=self.zs)
        if result.phase not in {"L", "V", "VL"}:
            raise ValueError(f"unsupported flash phase: {result.phase}")
        if not math.isfinite(result.H()) or not math.isfinite(result.V()) or result.V() <= 0:
            raise ValueError("non-finite enthalpy or nonpositive molar volume from NRTL flash")
        betas = result.betas
        if any(not math.isfinite(b) or not 0 <= b <= 1 for b in betas) or abs(sum(betas) - 1) > 1e-8:
            raise ValueError("flash phase fractions do not close")
        for i, z in enumerate(self.zs):
            if abs(sum(b * phase.zs[i] for b, phase in zip(betas, result.phases)) - z) > 1e-7:
                raise ValueError("flash composition balance does not close")
        return result

    def enthalpy_j_kg(self, temperature_k: float, pressure_pa: float) -> float:
        # H() is J/mol; MW is kg/mol. Bulk H already includes phase fractions.
        return float(self._cached_flash(float(temperature_k), float(pressure_pa)).H() / self.molar_mass_kg_mol)

    def phase(self, temperature_k: float, pressure_pa: float) -> str:
        result = self._cached_flash(float(temperature_k), float(pressure_pa))
        return {"L": "liquid", "V": "gas", "VL": "twophase"}[result.phase]

    def state(self, temperature_k: float, pressure_pa: float) -> ThermoState:
        result = self._cached_flash(float(temperature_k), float(pressure_pa))
        cp = None
        # Avoid presenting a phase-weighted Cp as an equilibrium two-phase derivative.
        if result.phase != "VL":
            value = float(result.Cp() / self.molar_mass_kg_mol)
            if not math.isfinite(value) or value <= 0:
                raise ValueError("nonpositive single-phase heat capacity from NRTL flash")
            cp = value
        return ThermoState(float(temperature_k), float(pressure_pa),
                           self.enthalpy_j_kg(temperature_k, pressure_pa),
                           float(self.molar_mass_kg_mol / result.V()), cp, None, None,
                           self.phase(temperature_k, pressure_pa))

    def saturation_at_pressure(self, pressure_pa: float) -> SaturationAtPressure:
        self._validate_domain(None, pressure_pa)
        try:
            bubble, dew = [float(self.flasher.flash(P=pressure_pa, VF=vf, zs=self.zs).T) for vf in (0, 1)]
            self._validate_domain(bubble, pressure_pa)
            self._validate_domain(dew, pressure_pa)
            if bubble > dew + 1e-6:
                raise ValueError("bubble temperature exceeds dew temperature")
            return SaturationAtPressure(pressure_pa, bubble, dew, True, None)
        except Exception as exc:
            return SaturationAtPressure(pressure_pa, None, None, False, str(exc))

    def metadata(self) -> dict:
        return {
            "provider_type": "thermo", "backend": "FlashVL", "model": "NRTL",
            "identity": self.identity, "components": self.components,
            "composition_basis": self.composition_basis, "input_fractions": self.input_fractions,
            "mole_fractions": self.mole_fractions, "molar_mass_kg_mol": self.molar_mass_kg_mol,
            "parameter_set": PARAMETER_SET, "parameter_source": PARAMETER_SOURCE,
            "parameter_order": ["Ethanol", "Water"], "tau_b_kelvin": TAU_B_K, "alpha": ALPHA,
            "gas_model": "IdealGas", "equilibrium_basis": "Psat", "caloric_basis": "Psat",
            "enthalpy_reference": "same thermo ideal-gas reference at 298.15 K for inlet and outlet; H includes HE",
            "operational_domain": {"temperature_k": [273.15, 500], "pressure_pa": [100, 2000000]},
            "parameter_fit_validity_range": None,
            "project_validation_status": "documented VLE example reproduction; caloric validation pending",
            "transport_validation_status": "not implemented for mixtures",
            "uncertainty_model": None,
            "packages": {p: version(p) for p in ("thermo", "chemicals", "fluids", "scipy")},
            "pure_correlation_methods": {key: [x.method for x in getattr(self.correlations, key)]
                                         for key in ("VaporPressures", "HeatCapacityGases", "VolumeLiquids")},
        }
