"""Automatic parameter-backed catalog and explicit binary NRTL/Dortmund models.

Pair parameters are not independent measurements. The ideal vapor, additive
liquid volumes and log-viscosity rule are provisional engineering assumptions.
"""
from __future__ import annotations

from functools import lru_cache
import math
from typing import Any

from v5m2_core import WATER_CAS, MODELS, digest, number, mass_to_mole, summarize_flash


def water_partners(table: dict) -> list[str]:
    """Enumerate keys actually present, never a hand-authored additive list."""
    result = set()
    for key in table:
        parts = key.split()
        if len(parts) == 2 and WATER_CAS in parts:
            result.update(x for x in parts if x != WATER_CAS)
    return sorted(result, key=lambda x: (int(x.replace("-", "")), x))


def nrtl_parameters(db: Any, cas: str) -> dict:
    order = [WATER_CAS, cas]
    bs, alphas = [[0., 0.], [0., 0.]], [[0., 0.], [0., 0.]]
    for i, j in ((0, 1), (1, 0)):
        for field, matrix in (("bij", bs), ("alphaij", alphas)):
            if not db.has_ip_specific("ChemSep NRTL", [order[i], order[j]], field):
                raise ValueError(f"missing directed NRTL parameter: {field} {order[i]}->{order[j]}")
            matrix[i][j] = number(db.get_ip_specific("ChemSep NRTL", [order[i], order[j]], field), field)
    if any(not 0 < alphas[i][j] <= 1 for i, j in ((0, 1), (1, 0))):
        raise ValueError("invalid NRTL alpha")
    return {"component_order": order, "tau_bs_k": bs, "alpha_cs": alphas,
            "source": "thermo.IPDB:ChemSep NRTL", "temperature_expression": "tau_ij = bij / T",
            "fit_temperature_range_k": None, "independent_validation": False,
            "table_metadata": db.metadata["ChemSep NRTL"]}


def complete_group_interactions(groups: list[dict], subgroups: dict, interactions: dict) -> dict:
    if len(groups) != 2 or any(not g for g in groups):
        raise ValueError("missing Dortmund group assignment")
    converted = [{int(k): int(v) for k, v in g.items()} for g in groups]
    if any(v <= 0 for g in converted for v in g.values()):
        raise ValueError("invalid Dortmund group count")
    mains = set()
    for g in converted:
        for subgroup in g:
            if subgroup not in subgroups:
                raise ValueError("unknown Dortmund subgroup")
            mains.add(subgroups[subgroup].main_group_id)
    missing = [(i, j) for i in sorted(mains) for j in sorted(mains)
               if i != j and j not in interactions.get(i, {})]
    if missing:
        raise ValueError("missing directed Dortmund interactions: " + str(missing))
    coeffs = {f"{i}:{j}": [number(float(v), "UNIFAC coefficient") for v in interactions[i][j]]
              for i in sorted(mains) for j in sorted(mains) if i != j}
    return {"groups": converted, "interaction_sha256": digest(coeffs),
            "parameter_set": "DOUFIP2016", "missing_interactions": [],
            "subgroups_sha256": digest({str(k): [float(subgroups[k].R), float(subgroups[k].Q),
                                                   subgroups[k].main_group_id]
                                        for g in converted for k in g}),
            "independent_validation": False, "fit_temperature_range_k": None}


def molecular_metadata(item: Any, expected_cas: str, config: dict) -> dict:
    if item is False or item is None:
        raise ValueError("missing_local_identity")
    from rdkit import Chem
    from rdkit.Chem import Descriptors, rdMolDescriptors
    if item.CASs != expected_cas:
        raise ValueError("CAS_metadata_mismatch")
    mol = Chem.MolFromSmiles(item.smiles)
    if mol is None or len(Chem.GetMolFrags(mol)) != 1:
        raise ValueError("invalid_or_multicomponent_structure")
    if any(a.GetFormalCharge() or a.GetNumRadicalElectrons() for a in mol.GetAtoms()):
        raise ValueError("electrolytes_or_radicals_outside_molecular_model")
    if any(a.GetAtomicNum() not in config["allowed_atomic_numbers"] for a in mol.GetAtoms()):
        raise ValueError("outside_declared_element_domain")
    if not any(a.GetAtomicNum() == 6 for a in mol.GetAtoms()):
        raise ValueError("nonorganic_additive_outside_this_catalog")
    for pattern in config["excluded_smarts"]:
        query = Chem.MolFromSmarts(pattern)
        if query is None:
            raise RuntimeError("invalid configured SMARTS: " + pattern)
        if mol.HasSubstructMatch(query):
            raise ValueError("excluded_model_chemistry:" + pattern)
    key = Chem.MolToInchiKey(mol)
    if not key or key != item.InChI_key:
        raise ValueError("metadata_structure_inchikey_mismatch")
    mw = float(Descriptors.MolWt(mol))
    if not 0 < mw <= config["maximum_molecular_weight_g_mol"]:
        raise ValueError("outside_molecular_weight_domain")
    if abs(float(item.MW)/mw - 1) > .002:
        raise ValueError("metadata_molecular_weight_mismatch")
    return {"pair_id": "water-" + key, "cas_number": expected_cas,
            "name": str(item.common_name or item.iupac_name or expected_cas),
            "smiles": Chem.MolToSmiles(mol, canonical=True, isomericSmiles=True),
            "inchi_key": key, "molecule_group": key.split("-")[0], "formula": rdMolDescriptors.CalcMolFormula(mol),
            "molecular_weight_g_mol": mw, "commercial_availability_verified": False,
            "safety_review_status": "not_reviewed", "preparable_without_new_synthesis": "known_identity_not_procurement_verified"}


def discover_catalog(config: dict, additional_cas: list[str] | None = None) -> dict:
    from chemicals.identifiers import get_pubchem_db
    from chemicals import phase_change
    from thermo.interaction_parameters import IPDB
    from thermo.unifac import UNIFAC_group_assignment_DDBST, DOUFSG, DOUFIP2016
    table = IPDB.tables.get("ChemSep NRTL")
    if not isinstance(table, dict) or not table:
        raise RuntimeError("ChemSep NRTL catalog is unavailable")
    metadata = get_pubchem_db()
    catalog_cas = set(water_partners(table))
    candidates = sorted(catalog_cas | set(additional_cas or []))
    accepted, rejected = [], []
    tin = config["grid"]["inlet_temperature_k"]
    water_groups = UNIFAC_group_assignment_DDBST(WATER_CAS, "MODIFIED_UNIFAC")
    for cas in candidates:
        try:
            row = molecular_metadata(metadata.search_CAS(cas, autoload=True), cas, config["catalog"])
            # These are scope guards, not an SLE or aqueous speciation model.
            tm, tb = phase_change.Tm(cas), phase_change.Tb(cas)
            if tm is None or tb is None or not float(tm) < tin < float(tb):
                raise ValueError("pure_liquid_inlet_scope_not_established")
            row.update(pure_melting_temperature_k=float(tm), pure_boiling_temperature_k=float(tb),
                       catalog_sources=["ChemSep NRTL"] if cas in catalog_cas else ["reference_input"])
            row["models"], row["model_rejections"] = {}, {}
            try:
                row["models"][MODELS[0]] = nrtl_parameters(IPDB, cas)
            except (ValueError, KeyError, TypeError) as exc:
                row["model_rejections"][MODELS[0]] = str(exc)
            try:
                groups = UNIFAC_group_assignment_DDBST(cas, "MODIFIED_UNIFAC")
                row["models"][MODELS[1]] = complete_group_interactions([water_groups, groups], DOUFSG, DOUFIP2016)
            except (ValueError, KeyError, TypeError) as exc:
                row["model_rejections"][MODELS[1]] = str(exc)
            if not row["models"]:
                raise ValueError("no_complete_mixture_model")
            row["parameter_sha256"] = digest(row["models"])
            accepted.append(row)
        except (ValueError, KeyError, TypeError, AttributeError) as exc:
            rejected.append({"cas_number": cas, "reason": str(exc)})
    accepted.sort(key=lambda r: r["pair_id"])
    # Identical chemical identities never become extra independent candidates.
    unique = {}
    for row in accepted:
        if row["inchi_key"] in unique:
            rejected.append({"cas_number": row["cas_number"], "reason": "duplicate_full_identity"})
        else:
            unique[row["inchi_key"]] = row
    accepted = list(unique.values())
    limit = config["catalog"]["maximum_pairs"]
    selected = accepted if limit is None else accepted[:limit]
    return {"source": "local_parameter_database_not_measurement_archive", "candidates": selected,
            "rejections": rejected, "discovered_cas_count": len(candidates),
            "eligible_pair_count_before_budget": len(accepted), "budget_limited": len(selected) < len(accepted),
            "raw_parameter_table_sha256": digest(table), "handwritten_additive_list_used": False}


def bounds(obj: Any, method: str) -> tuple[float, float]:
    limits = getattr(obj, "T_limits", {}).get(method)
    data = getattr(obj, "tabular_data", {}).get(method)
    if data is not None:
        limits = [min(data[0]), max(data[0])]
    if limits is None or len(limits) != 2 or any(not math.isfinite(float(x)) for x in limits):
        raise ValueError("pure_correlation_bounds_unknown")
    return float(limits[0]), float(limits[1])


def in_bounds(obj: Any, t: float) -> None:
    lo, hi = bounds(obj, obj.method)
    if not lo <= t <= hi or not obj.test_method_validity(t, obj.method):
        raise ValueError("outside_pure_correlation_range:" + str(obj.method))


def pin_correlation(obj: Any, temperatures: list[float]) -> dict:
    available = list(getattr(obj, "all_methods", []))
    preferred = list(getattr(obj, "ranked_methods", []))
    scored = []
    for method in available:
        try:
            lo, hi = bounds(obj, method)
            score = sum(lo <= t <= hi and bool(obj.test_method_validity(t, method)) for t in temperatures)
            priority = preferred.index(method) if method in preferred else len(preferred)
            scored.append((score, -priority, method, lo, hi))
        except (ValueError, TypeError, KeyError, AttributeError):
            continue
    if not scored:
        raise ValueError("no_bounded_pure_correlation")
    _, _, method, lo, hi = max(scored)
    obj.method = method
    obj.extrapolation = None
    obj.tabular_extrapolation_permitted = False
    return {"method": method, "temperature_bounds_k": [lo, hi], "extrapolation": None,
            "evidence_kind": "library_correlation_not_independent_measurement"}


class BinaryModel:
    """One fixed activity model, ideal vapor, two possible liquid phases."""
    def __init__(self, pair: dict, model: str, config: dict):
        from thermo import ChemicalConstantsPackage, NRTL, UNIFAC, GibbsExcessLiquid, IdealGas, FlashVLN
        self.pair, self.model, self.config = pair, model, config
        if model not in pair["models"]:
            raise ValueError("mixture_model_unavailable")
        self.constants, self.props = ChemicalConstantsPackage.from_IDs([WATER_CAS, pair["cas_number"]])
        if list(self.constants.CASs) != [WATER_CAS, pair["cas_number"]]:
            raise ValueError("component_order_mismatch")
        if self.constants.InChI_Keys[1] != pair["inchi_key"]:
            raise ValueError("property_package_identity_mismatch")
        self.mw = [float(x) for x in self.constants.MWs]
        points = sorted(set([298.15, config["grid"]["inlet_temperature_k"]] + config["grid"]["outlet_temperatures_k"]))
        self.methods = {}
        for name in ("VaporPressures", "HeatCapacityGases", "VolumeLiquids"):
            self.methods[name] = [pin_correlation(o, points) for o in getattr(self.props, name)]
        self.viscosities = []
        self.viscosity_metadata = []
        for obj in self.props.ViscosityLiquids:
            try:
                self.viscosity_metadata.append(pin_correlation(obj, [config["grid"]["inlet_temperature_k"]]))
                self.viscosities.append(obj)
            except (ValueError, TypeError, AttributeError):
                self.viscosities.append(None); self.viscosity_metadata.append(None)
        params = pair["models"][model]
        if model == "chemsep_nrtl":
            self.ge = NRTL(T=298.15, xs=[.5, .5], tau_bs=params["tau_bs_k"], alpha_cs=params["alpha_cs"])
        else:
            from thermo.unifac import DOUFSG, DOUFIP2016
            groups = [{int(k): int(v) for k, v in g.items()} for g in params["groups"]]
            check = complete_group_interactions(groups, DOUFSG, DOUFIP2016)
            if check["interaction_sha256"] != params["interaction_sha256"] or check["subgroups_sha256"] != params["subgroups_sha256"]:
                raise ValueError("Dortmund_parameters_changed_since_catalog")
            self.ge = UNIFAC.from_subgroups(T=298.15, xs=[.5, .5], chemgroups=groups,
                    version=1, interaction_data=DOUFIP2016, subgroups=DOUFSG)
        self.liquid = GibbsExcessLiquid(VaporPressures=self.props.VaporPressures,
            HeatCapacityGases=self.props.HeatCapacityGases, VolumeLiquids=self.props.VolumeLiquids,
            GibbsExcessModel=self.ge, equilibrium_basis="Psat", caloric_basis="Psat")
        gas = IdealGas(HeatCapacityGases=self.props.HeatCapacityGases)
        self.flasher = FlashVLN(self.constants, self.props, liquids=[self.liquid, self.liquid], gas=gas)

    def check_temperature(self, t: float) -> None:
        for name in ("VaporPressures", "HeatCapacityGases", "VolumeLiquids"):
            for obj in getattr(self.props, name):
                in_bounds(obj, t)
        for obj in self.props.HeatCapacityGases:
            # Ideal-gas enthalpy integrates from the common 298.15 K reference.
            in_bounds(obj, 298.15)
        for tc in self.constants.Tcs:
            if tc is None or t >= tc:
                raise ValueError("pure_supercritical_component_not_supported_by_Psat_liquid_model")

    @lru_cache(maxsize=4096)
    def state(self, w: float, t: float, p: float) -> dict:
        self.check_temperature(t)
        zs = mass_to_mole(w, self.mw)
        res = self.flasher.flash(T=t, P=p, zs=zs)
        return summarize_flash(res, zs, self.mw, self.config["balance_tolerance"])

    def bubble(self, w: float, t: float) -> dict:
        self.check_temperature(t)
        xs = mass_to_mole(w, self.mw)
        ge = self.ge.to_T_xs(T=t, xs=xs)
        gammas = [number(float(v), "gamma", positive=True) for v in ge.gammas()]
        psats = [number(float(o.calculate(t, o.method)), "pure Psat", positive=True) for o in self.props.VaporPressures]
        partial = [x*g*p for x, g, p in zip(xs, gammas, psats)]
        pressure = number(sum(partial), "bubble pressure", positive=True)
        return {"bubble_pressure_pa": pressure, "water_activity": xs[0]*gammas[0],
                "gammas": gammas, "incipient_vapor_additive_mole_fraction": partial[1]/pressure,
                "excess_enthalpy_j_mol": number(float(ge.HE()), "HE"),
                "excess_heat_capacity_j_mol_k": number(float(ge.CpE()), "CpE"),
                "bubble_basis": "homogeneous_liquid_gamma_Psat_ideal_vapor_no_Poynting"}

    def viscosity_ratio(self, w: float, t: float) -> float | None:
        """Molar log-mixing proxy only: no claim of a composition-specific fit."""
        if any(o is None for o in self.viscosities):
            return None
        try:
            for obj in self.viscosities:
                in_bounds(obj, t)
            mus = [number(float(o.calculate(t, o.method)), "mu", positive=True) for o in self.viscosities]
            zs = mass_to_mole(w, self.mw)
            return math.exp(sum(z*math.log(mu) for z, mu in zip(zs, mus)))/mus[0]
        except (ValueError, TypeError, ArithmeticError):
            return None

    def evaluate(self, w: float, t: float, p: float) -> dict:
        row = {"pair_id": self.pair["pair_id"], "cas_number": self.pair["cas_number"], "name": self.pair["name"],
               "smiles": self.pair["smiles"], "model": self.model, "additive_mass_fraction": w,
               "outlet_temperature_k": t, "pressure_pa": p, "inlet_temperature_k": self.config["grid"]["inlet_temperature_k"],
               "status": "pending", "model_comparison_eligible": False, "evidence_kind": "unvalidated_model_prediction",
               "caloric_accuracy_validated": False, "surface_enhancement_predicted": False,
               "porous_flow_simulated": False, "sequential_evaporation_simulated": False,
               "calibrated_uncertainty": None, "contact_angle_deg": None, "surface_tension_n_m": None,
               "viscosity_model": "mole_log_mixing_proxy_not_mixture_data", "liquid_volume_model": "additive_pure_molar_volumes"}
        try:
            inlet = row["inlet_temperature_k"]
            before, after = self.state(w, inlet, p), self.state(w, t, p)
            row.update(inlet=before, outlet=after)
            tol = self.config["phase_fraction_tolerance"]
            if before["liquid_phase_count"] != 1 or before["vapor_mole_fraction"] > tol:
                return {**row, "status": "inlet_not_single_liquid"}
            if after["liquid_phase_count"] > 1:
                return {**row, "status": "outlet_liquid_liquid_split"}
            # The same overall composition enters and leaves; no fictitious
            # separation heat or cross-composition absolute-enthalpy subtraction.
            dh = number(after["enthalpy_j_kg"] - before["enthalpy_j_kg"], "delta h", positive=True)
            w0, w1 = self.state(0., inlet, p), self.state(0., t, p)
            dh_water_model = number(w1["enthalpy_j_kg"]-w0["enthalpy_j_kg"], "water-model delta h", positive=True)
            water = water_reference(inlet, t, p)
            if water["status"] != "ok":
                return {**row, "status": "water_reference_unavailable", "water_reference": water}
            bubble = self.bubble(w, inlet)
            water_bubble = self.bubble(0., inlet)
            ratio = dh/dh_water_model
            endpoint_error = abs(dh_water_model/water["delta_h_j_kg"]-1)
            row.update(status="ok", delta_h_j_kg=dh, same_model_water_delta_h_j_kg=dh_water_model,
                same_model_water_delta_h_ratio=ratio, heos_water_delta_h_j_kg=water["delta_h_j_kg"],
                heos_water_delta_h_ratio=dh/water["delta_h_j_kg"], endpoint_relative_error=endpoint_error,
                endpoint_within_configured_tolerance=endpoint_error <= self.config["endpoint_relative_tolerance"],
                mass_ratio_same_net_heat=1/ratio, additional_net_heat_reduction_for_mass_parity=max(0., 1-ratio),
                storage_bubble_pressure_ratio=bubble["bubble_pressure_pa"]/water_bubble["bubble_pressure_pa"],
                storage_vle=bubble, viscosity_ratio_proxy=self.viscosity_ratio(w, inlet),
                model_comparison_eligible=endpoint_error <= self.config["endpoint_relative_tolerance"],
                vapor_composition_shift=after["vapor_additive_mass_fraction"]-w if after["vapor_additive_mass_fraction"] is not None else None,
                endpoint_check_is_accuracy_validation=False)
        except (ValueError, RuntimeError, TypeError, ArithmeticError) as exc:
            row.update(status="model_state_rejected", error=f"{type(exc).__name__}:{exc}")
        return row


@lru_cache(maxsize=128)
def water_reference(inlet: float, outlet: float, pressure: float) -> dict:
    import CoolProp.CoolProp as CP
    try:
        if CP.PhaseSI("T", inlet, "P", pressure, "Water") != "liquid":
            raise ValueError("water_inlet_not_liquid")
        h0 = float(CP.PropsSI("Hmass", "T", inlet, "P", pressure, "Water"))
        h1 = float(CP.PropsSI("Hmass", "T", outlet, "P", pressure, "Water"))
        return {"status": "ok", "delta_h_j_kg": number(h1-h0, "water delta h", positive=True),
                "inlet_temperature_k": inlet, "outlet_temperature_k": outlet, "pressure_pa": pressure,
                "provider": "CoolProp HEOS Water"}
    except (ValueError, RuntimeError, ArithmeticError) as exc:
        return {"status": "water_reference_failed", "error": str(exc)}
