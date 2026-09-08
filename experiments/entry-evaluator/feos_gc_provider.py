"""V5b-1 structure-derived pure-fluid provider using FeOS GC-PC-SAFT + Joback."""
from __future__ import annotations

from functools import lru_cache
import hashlib
from importlib.metadata import version
import math
from pathlib import Path
from typing import Any

from property_provider import SaturationAtPressure, ThermoState

HERE = Path(__file__).resolve().parent
PARAMETER_DIR = HERE / "parameters" / "v5b"
SMARTS_PATH = PARAMETER_DIR / "sauer2014_smarts.json"
RESIDUAL_PATH = PARAMETER_DIR / "rehner2023_hetero.json"
JOBACK_PATH = PARAMETER_DIR / "joback1987.json"

FEOS_VERSION = "0.10.1"
RDKIT_VERSION = "2026.3.6"


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


class FeosGcPcSaftProvider:
    """Predict a pure fluid from molecular structure, without reference properties.

    The residual Helmholtz model is heterosegmented GC-PC-SAFT. The ideal-gas
    heat-capacity contribution is Joback/Reid, decomposed from the same SMILES.
    CoolProp is deliberately not imported or queried by this provider.
    """

    def __init__(self, *, smiles: str, name: str | None = None):
        if version("feos") != FEOS_VERSION:
            raise RuntimeError(
                f"V5b-1 requires feos=={FEOS_VERSION}; runtime is {version('feos')}"
            )
        if version("rdkit") != RDKIT_VERSION:
            raise RuntimeError(
                f"V5b-1 requires rdkit=={RDKIT_VERSION}; runtime is {version('rdkit')}"
            )

        import feos
        import si_units as si
        from rdkit import Chem
        from rdkit.Chem import Descriptors

        mol = Chem.MolFromSmiles(str(smiles))
        if mol is None:
            raise ValueError(f"invalid SMILES: {smiles}")
        self.smiles = Chem.MolToSmiles(mol, canonical=True)
        self.name = str(name).strip() if name else self.smiles
        self.molar_mass_kg_mol = float(Descriptors.MolWt(mol)) / 1000.0
        if not math.isfinite(self.molar_mass_kg_mol) or self.molar_mass_kg_mol <= 0:
            raise ValueError("RDKit returned invalid molecular weight")

        for path in (SMARTS_PATH, RESIDUAL_PATH, JOBACK_PATH):
            if not path.is_file():
                raise FileNotFoundError(f"missing V5b parameter file: {path}")

        try:
            residual = feos.GcParameters.from_json_smiles(
                [self.smiles],
                str(SMARTS_PATH),
                str(RESIDUAL_PATH),
            )
            ideal = feos.GcParameters.from_json_smiles(
                [self.smiles],
                str(SMARTS_PATH),
                str(JOBACK_PATH),
            )
            self.eos = feos.EquationOfState.gc_pcsaft(residual).joback(ideal)
        except Exception as exc:
            raise ValueError(
                f"SMILES is outside the pinned V5b GC-PC-SAFT/Joback domain: "
                f"{self.smiles}: {exc}"
            ) from exc

        self.feos = feos
        self.si = si
        self.identity = f"feos:gc-pcsaft+joback:{self.smiles}"
        self.components = [self.name]
        self.input_fractions = [1.0]
        self.mole_fractions = [1.0]

    @staticmethod
    def _validate_tp(temperature_k: float, pressure_pa: float) -> tuple[float, float]:
        temperature_k = float(temperature_k)
        pressure_pa = float(pressure_pa)
        if (
            not math.isfinite(temperature_k)
            or temperature_k <= 0.0
            or not math.isfinite(pressure_pa)
            or pressure_pa <= 0.0
        ):
            raise ValueError("property state T and P must be finite and positive")
        return temperature_k, pressure_pa

    @lru_cache(maxsize=1024)
    def _vapor_pressure_pa(self, temperature_k: float) -> float | None:
        try:
            values = self.feos.PhaseEquilibrium.vapor_pressure(
                self.eos,
                float(temperature_k) * self.si.KELVIN,
            )
            if not values or values[0] is None:
                return None
            value = float(values[0] / self.si.PASCAL)
            return value if math.isfinite(value) and value > 0.0 else None
        except Exception:
            return None

    def _phase_hint(self, temperature_k: float, pressure_pa: float) -> str | None:
        psat = self._vapor_pressure_pa(float(temperature_k))
        if psat is None:
            return None
        tolerance = max(1.0, psat * 1.0e-6)
        if pressure_pa > psat + tolerance:
            return "liquid"
        if pressure_pa < psat - tolerance:
            return "vapor"
        return None

    @lru_cache(maxsize=4096)
    def _state(self, temperature_k: float, pressure_pa: float):
        temperature_k, pressure_pa = self._validate_tp(
            temperature_k, pressure_pa
        )
        kwargs: dict[str, Any] = {}
        hint = self._phase_hint(temperature_k, pressure_pa)
        if hint is not None:
            kwargs["density_initialization"] = hint
        try:
            return self.feos.State(
                self.eos,
                temperature=temperature_k * self.si.KELVIN,
                pressure=pressure_pa * self.si.PASCAL,
                **kwargs,
            )
        except Exception as exc:
            raise ValueError(
                f"FeOS TP state failed for {self.identity}: "
                f"T={temperature_k:g} K P={pressure_pa:g} Pa: {exc}"
            ) from exc

    def enthalpy_j_kg(self, temperature_k: float, pressure_pa: float) -> float:
        state = self._state(float(temperature_k), float(pressure_pa))
        value = float(
            state.specific_enthalpy()
            / (self.si.JOULE / self.si.KILOGRAM)
        )
        if not math.isfinite(value):
            raise ValueError("FeOS returned non-finite specific enthalpy")
        return value

    def phase(self, temperature_k: float, pressure_pa: float) -> str:
        temperature_k, pressure_pa = self._validate_tp(
            temperature_k, pressure_pa
        )
        psat = self._vapor_pressure_pa(temperature_k)
        if psat is None:
            return "supercritical"
        tolerance = max(1.0, psat * 1.0e-6)
        if pressure_pa > psat + tolerance:
            return "liquid"
        if pressure_pa < psat - tolerance:
            return "gas"
        return "twophase"

    def state(self, temperature_k: float, pressure_pa: float) -> ThermoState:
        state = self._state(float(temperature_k), float(pressure_pa))
        density = float(
            state.mass_density()
            / (self.si.KILOGRAM / self.si.METER**3)
        )
        cp = float(
            state.specific_isobaric_heat_capacity()
            / (self.si.JOULE / self.si.KILOGRAM / self.si.KELVIN)
        )
        if not math.isfinite(density) or density <= 0.0:
            raise ValueError("FeOS returned invalid mass density")
        if not math.isfinite(cp) or cp <= 0.0:
            raise ValueError("FeOS returned invalid isobaric heat capacity")
        return ThermoState(
            temperature_k=float(temperature_k),
            pressure_pa=float(pressure_pa),
            enthalpy_j_kg=self.enthalpy_j_kg(temperature_k, pressure_pa),
            density_kg_m3=density,
            cp_j_kg_k=cp,
            viscosity_pa_s=None,
            conductivity_w_m_k=None,
            phase=self.phase(temperature_k, pressure_pa),
        )

    @lru_cache(maxsize=256)
    def saturation_at_pressure(self, pressure_pa: float) -> SaturationAtPressure:
        pressure_pa = float(pressure_pa)
        if not math.isfinite(pressure_pa) or pressure_pa <= 0.0:
            raise ValueError("saturation pressure must be finite and positive")
        try:
            values = self.feos.PhaseEquilibrium.boiling_temperature(
                self.eos,
                pressure_pa * self.si.PASCAL,
            )
            if not values or values[0] is None:
                raise ValueError("no boiling temperature returned")
            temperature_k = float(values[0] / self.si.KELVIN)
            if not math.isfinite(temperature_k) or temperature_k <= 0.0:
                raise ValueError("invalid boiling temperature")
            return SaturationAtPressure(
                pressure_pa=pressure_pa,
                bubble_temperature_k=temperature_k,
                dew_temperature_k=temperature_k,
                supported=True,
                failure_reason=None,
            )
        except Exception as exc:
            return SaturationAtPressure(
                pressure_pa=pressure_pa,
                bubble_temperature_k=None,
                dew_temperature_k=None,
                supported=False,
                failure_reason=str(exc),
            )

    def metadata(self) -> dict[str, Any]:
        return {
            "provider_type": "feos",
            "backend": "GC-PC-SAFT+Joback",
            "model": "Rehner2023 heterosegmented GC-PC-SAFT + Joback1987",
            "identity": self.identity,
            "components": list(self.components),
            "smiles": self.smiles,
            "input_fractions": [1.0],
            "mole_fractions": [1.0],
            "molar_mass_kg_mol": self.molar_mass_kg_mol,
            "structure_is_only_candidate_input": True,
            "reference_property_backend_used": False,
            "parameter_files": {
                "smarts": {
                    "path": str(SMARTS_PATH.relative_to(HERE)),
                    "sha256": _sha256(SMARTS_PATH),
                    "upstream": "feos v0.10.1 parameters/pcsaft/sauer2014_smarts.json",
                },
                "residual": {
                    "path": str(RESIDUAL_PATH.relative_to(HERE)),
                    "sha256": _sha256(RESIDUAL_PATH),
                    "upstream": "feos v0.10.1 parameters/pcsaft/rehner2023_hetero.json",
                },
                "ideal_gas": {
                    "path": str(JOBACK_PATH.relative_to(HERE)),
                    "sha256": _sha256(JOBACK_PATH),
                    "upstream": "feos v0.10.1 parameters/ideal_gas/joback1987.json",
                },
            },
            "packages": {
                "feos": version("feos"),
                "rdkit": version("rdkit"),
            },
            "transport_validation_status": "not implemented in V5b-1",
            "uncertainty_model": None,
            "project_validation_status": (
                "V5b-1 structure-derived prediction; holdout validation required"
            ),
        }
