"""Isolated structure-only V5b-3 prediction worker. No reference-library imports."""
from __future__ import annotations

import argparse
from functools import lru_cache
import json
import math
from pathlib import Path
import sys
from typing import Any

ENTRY = Path(__file__).resolve().parents[1] / "entry-evaluator"


def finite_positive(value: Any) -> float:
    value = float(value)
    if not math.isfinite(value) or value <= 0:
        raise ValueError("invalid_predicted_property")
    return value


class FeosBenchmarkAdapter:
    def __init__(self, smiles: str):
        if str(ENTRY) not in sys.path:
            sys.path.insert(0, str(ENTRY))
        from feos_gc_provider import FeosGcPcSaftProvider
        self.provider = FeosGcPcSaftProvider(smiles=smiles)
        self.feos, self.si, self.eos = self.provider.feos, self.provider.si, self.provider.eos

    @lru_cache(maxsize=1)
    def critical(self):
        # Fixed initial guesses; never initialize from a reference critical point.
        errors = []
        for initial in (None, 300.0, 500.0, 700.0):
            try:
                kw = {} if initial is None else {"initial_temperature": initial * self.si.KELVIN}
                return self.feos.State.critical_point(self.eos, max_iter=100, **kw)
            except Exception as exc:
                errors.append(str(exc))
        raise ValueError("critical_point_failed: " + "; ".join(errors))

    @lru_cache(maxsize=256)
    def equilibrium(self, t: float):
        eq = self.feos.PhaseEquilibrium.pure(self.eos, t * self.si.KELVIN, max_iter=100)
        rl = finite_positive(eq.liquid.mass_density() / (self.si.KILOGRAM / self.si.METER**3))
        rv = finite_positive(eq.vapor.mass_density() / (self.si.KILOGRAM / self.si.METER**3))
        pl = finite_positive(eq.liquid.pressure() / self.si.PASCAL)
        pv = finite_positive(eq.vapor.pressure() / self.si.PASCAL)
        if abs(pl - pv) / max(pl, pv) > 1e-5:
            raise ValueError("phase_pressures_not_in_equilibrium")
        if rl <= rv * (1.0 + 1e-6):
            raise ValueError("degenerate_or_inverted_vapor_liquid_equilibrium")
        return eq

    def predict(self, request: dict[str, Any]) -> float:
        prop, t = request["property"], request["temperature_k"]
        si = self.si
        if prop == "critical_temperature_k":
            value = self.critical().temperature / si.KELVIN
        elif prop == "critical_pressure_pa":
            value = self.critical().pressure() / si.PASCAL
        elif prop == "normal_boiling_temperature_k":
            sat = self.provider.saturation_at_pressure(101325.0)
            if not sat.supported or sat.bubble_temperature_k is None:
                raise ValueError("normal_boiling_calculation_failed: " + str(sat.failure_reason))
            value = sat.bubble_temperature_k
        else:
            eq = self.equilibrium(finite_positive(t))
            if prop == "vapor_pressure_pa":
                value = eq.liquid.pressure() / si.PASCAL
            elif prop == "latent_heat_vaporization_j_kg":
                value = (eq.vapor.specific_enthalpy() - eq.liquid.specific_enthalpy()) / (si.JOULE / si.KILOGRAM)
            elif prop == "liquid_density_kg_m3":
                value = eq.liquid.mass_density() / (si.KILOGRAM / si.METER**3)
            elif prop == "liquid_cp_j_kg_k":
                value = eq.liquid.specific_isobaric_heat_capacity() / (si.JOULE / si.KILOGRAM / si.KELVIN)
            else:
                raise ValueError("unsupported_property: " + prop)
        return finite_positive(value)

    def metadata(self):
        return self.provider.metadata()


def validate_request(payload: dict[str, Any]) -> None:
    if set(payload) != {"smiles", "requests"} or not isinstance(payload["smiles"], str):
        raise ValueError("prediction_input_must_contain_only_smiles_and_requests")
    fields = {"observation_id", "property", "temperature_k", "pressure_pa", "state_basis"}
    seen = set()
    for row in payload["requests"]:
        if set(row) != fields:
            raise ValueError("reference_or_unknown_field_in_prediction_request")
        if row["observation_id"] in seen:
            raise ValueError("duplicate_prediction_request")
        seen.add(row["observation_id"])
        if row["temperature_k"] is not None:
            finite_positive(row["temperature_k"])
        if row["pressure_pa"] not in (None, 101325.0):
            raise ValueError("unsupported_request_pressure")


def run_prediction(payload: dict[str, Any], factory=FeosBenchmarkAdapter) -> dict[str, Any]:
    validate_request(payload)
    provider_error = None
    metadata = None
    try:
        adapter = factory(payload["smiles"])
        metadata = adapter.metadata()
    except Exception as exc:
        adapter = None
        provider_error = f"{type(exc).__name__}: {exc}"
    results = []
    for request in payload["requests"]:
        value, error = None, provider_error
        if adapter is not None:
            try:
                value = finite_positive(adapter.predict(request))
            except Exception as exc:
                error = f"{type(exc).__name__}: {exc}"
        results.append({"observation_id": request["observation_id"], "predicted_value": value,
                        "status": "ok" if value is not None else "prediction_failed", "error": error})
    return {"results": results, "provider_metadata": metadata,
            "reference_values_supplied": False, "parameters_fitted": False}


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    payload = json.loads(args.input.read_text(encoding="utf-8"))
    result = run_prediction(payload)
    temp = args.output.with_suffix(".tmp")
    temp.write_text(json.dumps(result, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    temp.replace(args.output)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
