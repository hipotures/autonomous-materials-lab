"""Audit existing V5e-2.1 packets without network access or predictor calls."""
from __future__ import annotations

import hashlib
import json
import math
import re
from typing import Any

TC = "critical_temperature_k"
PC = "critical_pressure_pa"
TB = "normal_boiling_temperature_k"
PS = "vapor_pressure_pa"
RHO = "liquid_density_kg_m3"
CP = "liquid_cp_j_kg_k"
HV = "latent_heat_vaporization_j_kg"
UNITS = {TC: "K", PC: "Pa", TB: "K", PS: "Pa", RHO: "kg/m^3", CP: "J/kg/K", HV: "J/kg"}
CONSTANTS = {TC, PC, TB}
ALLOW = {
    TC: {"HEOS", "IUPAC", "MATTHEWS", "CRC", "WEBBOOK"},
    PC: {"HEOS", "IUPAC", "MATTHEWS", "CRC", "WEBBOOK"},
    TB: {"HEOS", "CRC_INORG", "CRC_ORG", "COMMON_CHEMISTRY", "WEBBOOK"},
    PS: {"HEOS_FIT", "WAGNER_MCGARRY", "WAGNER_POLING", "ANTOINE_EXTENDED_POLING",
         "DIPPR_PERRY_8E", "VDI_PPDS", "COOLPROP", "ANTOINE_POLING", "VDI_TABULAR",
         "ANTOINE_WEBBOOK", "LANDOLT", "ALCOCK_ELEMENTS"},
    RHO: {"HEOS_FIT", "DIPPR_PERRY_8E", "VDI_PPDS", "COOLPROP", "VDI_TABULAR"},
    CP: {"HEOS_FIT", "COOLPROP", "VDI_TABULAR", "POLING_CONST", "CRCSTD",
         "ZABRANSKY_SPLINE", "ZABRANSKY_QUASIPOLYNOMIAL", "ZABRANSKY_SPLINE_C",
         "ZABRANSKY_QUASIPOLYNOMIAL_C", "ZABRANSKY_SPLINE_SAT",
         "ZABRANSKY_QUASIPOLYNOMIAL_SAT", "WEBBOOK_SHOMATE"},
    HV: {"HEOS_FIT", "DIPPR_PERRY_8E", "VDI_PPDS", "COOLPROP", "VDI_TABULAR",
         "CRC_HVAP_TB", "CRC_HVAP_298"},
}
DOCS = {
    TC: "https://chemicals.readthedocs.io/chemicals.critical.html",
    PC: "https://chemicals.readthedocs.io/chemicals.critical.html",
    TB: "https://chemicals.readthedocs.io/chemicals.phase_change.html",
    PS: "https://thermo.readthedocs.io/thermo.vapor_pressure.html",
    RHO: "https://thermo.readthedocs.io/thermo.volume.html",
    CP: "https://thermo.readthedocs.io/thermo.heat_capacity.html",
    HV: "https://thermo.readthedocs.io/thermo.phase_change.html",
}


def digest(value: Any) -> str:
    return hashlib.sha256(json.dumps(value, sort_keys=True, allow_nan=False,
                                     separators=(",", ":")).encode()).hexdigest()


def positive(value: Any) -> float:
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise ValueError("not_a_number")
    result = float(value)
    if not math.isfinite(result) or result <= 0:
        raise ValueError("nonpositive_or_nonfinite")
    return result


def allowed(prop: str, method: str) -> bool:
    if method in ALLOW.get(prop, set()):
        return True
    # Library interval suffixes only; no open-ended wildcard allowlist.
    base = re.sub(r"_\d+$", "", method)
    return prop == CP and base in ALLOW[CP] and base.startswith(("ZABRANSKY", "WEBBOOK_SHOMATE"))


def source_group(method: str) -> str:
    if "POLING" in method:
        return "POLING"
    if "WEBBOOK" in method:
        return "NIST"
    for prefix, group in (("HEOS", "HEOS_COOLPROP"), ("COOLPROP", "HEOS_COOLPROP"),
                          ("CRC", "CRC"), ("ZABRANSKY", "ZABRANSKY"),
                          ("DIPPR", "DIPPR_PERRY"), ("VDI", "VDI"),
                          ("WAGNER_MCGARRY", "MCGARRY")):
        if method.startswith(prefix):
            return group
    return method


def valid_cas(cas: str) -> bool:
    if not re.fullmatch(r"\d{2,7}-\d{2}-\d", cas):
        return False
    body, check = cas.replace("-", "")[:-1], int(cas[-1])
    return sum(i * int(d) for i, d in enumerate(reversed(body), 1)) % 10 == check


def structural_identity(smiles: str) -> dict[str, Any]:
    from rdkit import Chem
    from rdkit.Chem import Descriptors, rdMolDescriptors
    mol = Chem.MolFromSmiles(smiles)
    if mol is None or len(Chem.GetMolFrags(mol)) != 1:
        raise ValueError("invalid_or_multicomponent_smiles")
    key = Chem.MolToInchiKey(mol)
    if not key:
        raise ValueError("inchikey_unavailable")
    return {"smiles": Chem.MolToSmiles(mol, canonical=True, isomericSmiles=True),
            "inchi_key": key, "molecule_group": key.split("-")[0],
            "formula": rdMolDescriptors.CalcMolFormula(mol),
            "mw": float(Descriptors.MolWt(mol))}


class LocalReferenceBackend:
    """Only packaged metadata, constants and named correlations; no HTTP."""
    def __init__(self) -> None:
        from chemicals import critical, phase_change
        from chemicals.identifiers import get_pubchem_db
        from thermo import HeatCapacityLiquid, EnthalpyVaporization, VolumeLiquid, VaporPressure
        self.metadata = get_pubchem_db()
        self.constant_getters = {TC: critical.Tc, PC: critical.Pc, TB: phase_change.Tb}
        self.constant_methods = {TC: critical.Tc_methods, PC: critical.Pc_methods, TB: phase_change.Tb_methods}
        self.classes = {PS: VaporPressure, RHO: VolumeLiquid, CP: HeatCapacityLiquid, HV: EnthalpyVaporization}
        self.objects: dict[str, Any] = {}

    def identity(self, cas: str) -> dict[str, Any]:
        item = self.metadata.search_CAS(cas, autoload=True)
        # chemicals uses False for a missing lookup; some backends use None.
        # Convert only missing-record sentinels into a per-target rejection.
        if item is False or item is None:
            raise ValueError("cas_not_in_local_metadata")
        return {"cas": item.CASs, "inchi_key": item.InChI_key, "smiles": item.smiles,
                "formula": item.formula, "mw": float(item.MW)}

    def constant(self, cas: str, prop: str, method: str) -> float:
        if method not in self.constant_methods[prop](cas):
            raise ValueError("reference_method_unavailable")
        return positive(float(self.constant_getters[prop](cas, method=method)))

    def alternatives(self, cas: str, prop: str) -> list[dict[str, Any]]:
        result = []
        for method in sorted(self.constant_methods[prop](cas)):
            if allowed(prop, method):
                try:
                    result.append({"method": method, "source_group": source_group(method),
                                   "value": self.constant(cas, prop, method)})
                except (ValueError, TypeError, ArithmeticError):
                    continue
        return result

    def object(self, cas: str, prop: str, packet: dict[str, Any]):
        key = digest([cas, prop, packet.get("constants"), packet["molecular_weight_g_mol"]])
        if key not in self.objects:
            constants = packet.get("constants") or {}
            def known(name):
                row = constants.get(name) or {}
                return row.get("value") if allowed(name, str(row.get("method", ""))) else None
            kw = {"CASRN": cas, "extrapolation": None, "Tc": known(TC)}
            if prop in (PS, RHO, HV):
                kw.update(Tb=known(TB), Pc=known(PC))
            if prop in (RHO, CP):
                kw["MW"] = packet["molecular_weight_g_mol"]
            obj = self.classes[prop](**kw)
            obj.tabular_extrapolation_permitted = False
            self.objects[key] = obj
        return self.objects[key]

    def native_point(self, cas: str, prop: str, method: str, packet: dict[str, Any]):
        obj = self.object(cas, prop, packet)
        if method not in obj.all_methods:
            raise ValueError("reference_method_unavailable")
        if method == "CRC_HVAP_TB":
            return float(obj.CRC_HVAP_TB_Tb), float(obj.CRC_HVAP_TB_Hvap)
        if method == "CRC_HVAP_298":
            return 298.15, float(obj.CRC_HVAP_298)
        if prop == CP and method in {"POLING_CONST", "CRCSTD"}:
            return 298.15, float(obj.calculate(298.15, method))
        raise ValueError("not_a_native_point_method")

    def series_value(self, cas: str, prop: str, method: str, t: float, packet: dict[str, Any]):
        obj = self.object(cas, prop, packet)
        if method not in obj.all_methods:
            raise ValueError("reference_method_unavailable")
        limits = getattr(obj, "T_limits", {}).get(method)
        tabular = getattr(obj, "tabular_data", {}).get(method)
        if tabular is not None:
            limits = (min(tabular[0]), max(tabular[0]))
        if limits is None or not all(math.isfinite(float(x)) for x in limits):
            raise ValueError("reference_bounds_unknown")
        if not float(limits[0]) <= t <= float(limits[1]):
            raise ValueError("outside_reference_bounds")
        if not obj.test_method_validity(t, method):
            raise ValueError("reference_method_invalid_at_temperature")
        return positive(float(obj.calculate(t, method))), [float(x) for x in limits]


def verify_identity(target: dict[str, Any], backend: Any, tolerance: float) -> dict[str, Any]:
    local = structural_identity(str(target.get("canonical_smiles") or target.get("smiles") or ""))
    expected = str(target.get("inchi_key") or "")
    if expected != local["inchi_key"]:
        raise ValueError("target_smiles_inchikey_mismatch")
    cas = str(target.get("selected_cas_number") or "")
    if not valid_cas(cas):
        raise ValueError("missing_or_invalid_selected_cas")
    packet = target.get("packet") or {}
    if packet.get("cas_number") != cas:
        raise ValueError("packet_selected_cas_mismatch")
    meta = backend.identity(cas)
    if meta["cas"] != cas or meta["inchi_key"] != expected:
        raise ValueError("cas_identity_mismatch")
    if structural_identity(meta["smiles"])["inchi_key"] != expected:
        raise ValueError("metadata_structure_mismatch")
    if meta["formula"] != local["formula"]:
        raise ValueError("metadata_formula_mismatch")
    for value in (meta["mw"], packet.get("molecular_weight_g_mol")):
        if abs(positive(value) / local["mw"] - 1) > tolerance:
            raise ValueError("molecular_weight_mismatch")
    return {**local, "cas": cas, "identity_verified": True,
            "identity_source": "chemicals_local_metadata_exact_CAS_and_full_InChIKey",
            "metadata_independent_of_pubchem": False}


def physical_basis(prop: str, method: str) -> tuple[str, bool]:
    if prop in CONSTANTS:
        return ("normal_boiling_at_101325_pa" if prop == TB else "critical_point"), True
    if prop in (PS, HV):
        return "saturation_at_temperature", True
    # VDI tabulations explicitly use the saturation curve. Other legacy liquid
    # packets do not preserve enough pressure information for strict comparison.
    if method == "VDI_TABULAR":
        return "saturated_liquid_at_temperature", True
    return "liquid_pressure_unspecified_diagnostic_at_predicted_saturation", False


def convert(raw: float, prop: str, mw: float) -> float:
    if prop == RHO:
        return positive(mw / 1000.0 / positive(raw))
    if prop in (CP, HV):
        return positive(raw * 1000.0 / mw)
    return positive(raw)


def audit_target(target: dict[str, Any], config: dict[str, Any], backend: Any) -> dict[str, Any]:
    base = {"candidate_id": target["candidate_id"], "family": target.get("family", "unknown"),
            "legacy_reference_quality": target.get("reference_quality"),
            "rankable": False, "observations": [], "audit": []}
    try:
        identity = verify_identity(target, backend, config["identity_mw_relative_tolerance"])
    except (ValueError, TypeError, KeyError, RuntimeError) as exc:
        return {**base, "identity_verified": False, "identity_reason": str(exc)}
    base.update(identity)
    packet = target["packet"]
    mw = positive(packet["molecular_weight_g_mol"])
    for section in ("constants", "series"):
        for prop, record in sorted((packet.get(section) or {}).items()):
            if not record:
                continue
            method = str(record.get("method", ""))
            log = {"property": prop, "method": method, "accepted_points": 0, "rejected": []}
            base["audit"].append(log)
            try:
                if prop not in UNITS or record.get("unit") != UNITS[prop]:
                    raise ValueError("unknown_property_or_unit_mismatch")
                if (prop in CONSTANTS) != (section == "constants"):
                    raise ValueError("property_section_mismatch")
                if not allowed(prop, method):
                    raise ValueError("method_not_in_audited_allowlist")
                basis, strict = physical_basis(prop, method)
                provenance = {"source_group": source_group(method), "method": method,
                              "source_documentation": DOCS[prop],
                              "primary_measurement_id": None, "measurement_uncertainty": None,
                              "independent_measurement_count": None,
                              "legacy_record_sha256": digest(record),
                              "library_bounds_are_not_measurement_uncertainty": True}
                single = (prop == HV and method in {"CRC_HVAP_TB", "CRC_HVAP_298"}) or (
                    prop == CP and method in {"POLING_CONST", "CRCSTD"})
                if prop in CONSTANTS:
                    native = backend.constant(identity["cas"], prop, method)
                    if abs(native / positive(record["value"]) - 1) > config["reproduction_relative_tolerance"]:
                        raise ValueError("local_source_differs_from_frozen_packet")
                    alternatives = backend.alternatives(identity["cas"], prop)
                    conflicts = [r for r in alternatives if r["source_group"] != source_group(method)
                                 and abs(r["value"] / native - 1) > config["constant_conflict_relative_tolerance"]]
                    log["alternative_compilation_values"] = alternatives
                    if conflicts:
                        raise ValueError("constant_compilation_conflict")
                    points = [{"temperature_k": None, "value": native}]
                    kind = "compiled_reference_constant"
                elif single:
                    t, raw = backend.native_point(identity["cas"], prop, method, packet)
                    points = [{"temperature_k": positive(t), "value": convert(raw, prop, mw)}]
                    kind = "native_reference_point"
                    log["transformation"] = "replace_legacy_grid_with_native_reference_point"
                    log["legacy_grid_points_not_independent"] = len(record.get("points") or [])
                else:
                    points = list(record.get("points") or [])
                    kind = "reference_correlation_samples"
                seen = set()
                for point in points:
                    try:
                        t = positive(point["temperature_k"]) if prop not in CONSTANTS else None
                        if t in seen:
                            raise ValueError("duplicate_reference_temperature")
                        seen.add(t)
                        limits = None
                        value = positive(point["value"])
                        if prop not in CONSTANTS and not single:
                            raw, limits = backend.series_value(identity["cas"], prop, method, t, packet)
                            rebuilt = convert(raw, prop, mw)
                            if abs(rebuilt / value - 1) > config["reproduction_relative_tolerance"]:
                                raise ValueError("local_source_differs_from_frozen_packet")
                        if prop not in CONSTANTS:
                            tc_record = (packet.get("constants") or {}).get(TC) or {}
                            if allowed(TC, str(tc_record.get("method", ""))) and tc_record.get("value"):
                                if t >= positive(tc_record["value"]):
                                    raise ValueError("reference_at_or_above_critical_temperature")
                        row = {"candidate_id": target["candidate_id"], "molecule_group": identity["molecule_group"],
                               "smiles": identity["smiles"], "inchi_key": identity["inchi_key"],
                               "family": base["family"], "property": prop, "unit": UNITS[prop],
                               "temperature_k": t, "pressure_pa": 101325.0 if prop == TB else None,
                               "reference_value": value, "state_basis": basis,
                               "benchmark_eligible": strict, "evidence_kind": kind,
                               "correlation_temperature_bounds_k": limits, "provenance": provenance}
                        row["observation_id"] = digest([identity["inchi_key"], prop, t, method])[:24]
                        base["observations"].append(row)
                        log["accepted_points"] += 1
                    except (ValueError, TypeError, KeyError, ArithmeticError, RuntimeError) as exc:
                        log["rejected"].append({"temperature_k": point.get("temperature_k"), "reason": str(exc)})
            except (ValueError, TypeError, KeyError, ArithmeticError, RuntimeError, AttributeError) as exc:
                log["rejected"].append({"reason": str(exc)})
    base["source_groups"] = sorted({r["provenance"]["source_group"] for r in base["observations"]})
    base["source_group_count_is_not_independent_experiment_count"] = True
    return base
