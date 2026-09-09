"""V5e-3 property-gap active acquisition.

Select the next small calibration batch to close property-specific V5b-3
reference deficits. The selector uses V5d thermodynamic merit, the deterministic
V5b-3 molecule split, structural distance to existing strict anchors, and
empirical source availability discovered only after exact identity resolution.
"""
from __future__ import annotations

from dataclasses import asdict, dataclass
from functools import lru_cache
import hashlib
import math
from typing import Any, Iterable

TC = "critical_temperature_k"
PC = "critical_pressure_pa"
PS = "vapor_pressure_pa"
HV = "latent_heat_vaporization_j_kg"
TARGET_PROPERTIES = (TC, PC, PS, HV)


@dataclass(frozen=True)
class GapFeature:
    candidate_id: str
    smiles: str
    family: str
    split: str
    molecule_group: str
    property_priority_score: float
    property_percentile: float
    nearest_existing_anchor_similarity: float
    global_structural_novelty: float
    gap_structural_opportunity: float
    preprobe_score: float

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)


@lru_cache(maxsize=262144)
def _fingerprint(smiles: str):
    from rdkit import Chem
    from rdkit.Chem import rdFingerprintGenerator
    molecule = Chem.MolFromSmiles(smiles)
    if molecule is None:
        raise ValueError(f"invalid SMILES: {smiles}")
    generator = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
    return generator.GetFingerprint(molecule)


def canonical_identity(smiles: str) -> dict[str, Any]:
    from rdkit import Chem
    from rdkit.Chem import Descriptors, rdMolDescriptors
    molecule = Chem.MolFromSmiles(str(smiles))
    if molecule is None or len(Chem.GetMolFrags(molecule)) != 1:
        raise ValueError(f"invalid or multicomponent SMILES: {smiles}")
    canonical = Chem.MolToSmiles(molecule, canonical=True, isomericSmiles=True)
    key = Chem.MolToInchiKey(molecule)
    if not key:
        raise ValueError(f"InChIKey unavailable: {smiles}")
    return {
        "canonical_smiles": canonical,
        "inchi_key": key,
        "molecule_group": key.split("-")[0],
        "molecular_formula": rdMolDescriptors.CalcMolFormula(molecule),
        "molecular_weight_g_mol": float(Descriptors.MolWt(molecule)),
    }


def assignment(molecule_group: str, statistics: dict[str, Any]) -> str:
    value = hashlib.sha256(
        (str(statistics["split_seed"]) + ":" + str(molecule_group)).encode()
    ).hexdigest()
    fraction = int(value[:16], 16) / 2**64
    return "evaluation" if fraction < float(statistics["evaluation_fraction"]) else "calibration"


def tanimoto(smiles_a: str, smiles_b: str) -> float:
    from rdkit import DataStructs
    return float(DataStructs.TanimotoSimilarity(_fingerprint(smiles_a), _fingerprint(smiles_b)))


def nearest_similarity(smiles: str, references: Iterable[str]) -> float:
    values = [tanimoto(smiles, ref) for ref in references]
    return max(values) if values else 0.0


def percentile_map(values: dict[str, float]) -> dict[str, float]:
    if not values:
        return {}
    ordered = sorted(values.items(), key=lambda item: (item[1], item[0]))
    if len(ordered) == 1:
        return {ordered[0][0]: 1.0}
    return {key: index / (len(ordered) - 1) for index, (key, _) in enumerate(ordered)}


def property_deficits(v5b3_summary: dict[str, Any], config: dict[str, Any]) -> dict[str, dict[str, int]]:
    targets = config["gap_targets"]
    result: dict[str, dict[str, int]] = {}
    for prop in config["priority_properties"]:
        row = (v5b3_summary.get("property_results") or {}).get(prop) or {}
        current_cal = int(row.get("calibration_group_count") or 0)
        current_eval = int(row.get("evaluation_group_count") or 0)
        target_cal = int(targets["calibration_groups"])
        target_eval = int(targets["evaluation_groups"])
        result[str(prop)] = {
            "current_calibration": current_cal,
            "current_evaluation": current_eval,
            "target_calibration": target_cal,
            "target_evaluation": target_eval,
            "calibration": max(0, target_cal - current_cal),
            "evaluation": max(0, target_eval - current_eval),
        }
    return result


def strict_property_anchors(audited_reference: dict[str, Any], properties: Iterable[str]) -> dict[str, list[str]]:
    wanted = {str(prop) for prop in properties}
    anchors: dict[str, set[str]] = {prop: set() for prop in wanted}
    for target in audited_reference.get("targets", []):
        for obs in target.get("observations", []):
            prop = str(obs.get("property") or "")
            if prop in wanted and bool(obs.get("benchmark_eligible")) and obs.get("smiles"):
                anchors[prop].add(str(obs["smiles"]))
    for obs in audited_reference.get("observations", []):
        prop = str(obs.get("property") or "")
        if prop in wanted and bool(obs.get("benchmark_eligible")) and obs.get("smiles"):
            anchors[prop].add(str(obs["smiles"]))
    return {key: sorted(value) for key, value in anchors.items()}


def _finite_positive(value: Any) -> float | None:
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        return None
    value = float(value)
    if not math.isfinite(value) or value <= 0.0:
        return None
    return value


def _family(row: dict[str, Any]) -> str:
    provenance = row.get("provenance") or []
    if provenance and isinstance(provenance[0], dict):
        return str(provenance[0].get("family") or "unknown")
    return str(row.get("family") or "unknown")


def candidate_pool(full_prescreen: dict[str, Any], *, excluded_smiles: set[str], config: dict[str, Any]) -> list[dict[str, Any]]:
    minimum_merit = float(config["candidate_gate"]["minimum_property_priority_score"])
    require_exploratory = bool(config["candidate_gate"].get("require_exploratory_lane", True))
    rows: dict[str, dict[str, Any]] = {}
    for source in full_prescreen.get("passes", []):
        if require_exploratory and source.get("lane") != "exploratory_domain_expansion":
            continue
        metrics = source.get("property_screen") or {}
        merit = _finite_positive(metrics.get("property_priority_score"))
        if merit is None or merit < minimum_merit:
            continue
        identity = canonical_identity(str(source["smiles"]))
        smiles = identity["canonical_smiles"]
        if smiles in excluded_smiles:
            continue
        row = {**source, "smiles": smiles, "_identity": identity,
               "_property_priority_score": merit, "_family": _family(source)}
        rows.setdefault(smiles, row)
    return [rows[key] for key in sorted(rows)]


def prepare_preprobe_features(candidates: list[dict[str, Any]], *, deficits: dict[str, dict[str, int]],
                              anchors: dict[str, list[str]], statistics: dict[str, Any],
                              config: dict[str, Any]) -> list[GapFeature]:
    merit_values = {str(row["candidate_id"]): float(row["_property_priority_score"]) for row in candidates}
    merit_pct = percentile_map(merit_values)
    all_anchor_smiles = sorted({smiles for values in anchors.values() for smiles in values})
    weights = config["preprobe_weights"]
    output: list[GapFeature] = []
    for row in candidates:
        identity = row["_identity"]
        split = assignment(identity["molecule_group"], statistics)
        candidate_id = str(row["candidate_id"])
        weighted_gap = 0.0
        weighted_novelty = 0.0
        for prop in config["priority_properties"]:
            remaining = int(deficits[prop][split])
            if remaining <= 0:
                continue
            target = max(1, int(deficits[prop]["target_calibration" if split == "calibration" else "target_evaluation"]))
            pressure = remaining / target
            novelty = 1.0 - nearest_similarity(row["smiles"], anchors.get(prop, []))
            weighted_gap += pressure
            weighted_novelty += pressure * novelty
        structural_opportunity = weighted_novelty / weighted_gap if weighted_gap > 0.0 else 0.0
        nearest_all = nearest_similarity(row["smiles"], all_anchor_smiles)
        global_novelty = 1.0 - nearest_all
        gap_pressure = min(1.0, weighted_gap / max(1.0, len(config["priority_properties"])))
        score = (
            float(weights["property_merit"]) * merit_pct[candidate_id]
            + float(weights["gap_pressure"]) * gap_pressure
            + float(weights["property_structural_novelty"]) * structural_opportunity
            + float(weights["global_novelty"]) * global_novelty
        )
        output.append(GapFeature(candidate_id=candidate_id, smiles=row["smiles"], family=row["_family"],
                                 split=split, molecule_group=identity["molecule_group"],
                                 property_priority_score=row["_property_priority_score"],
                                 property_percentile=merit_pct[candidate_id],
                                 nearest_existing_anchor_similarity=nearest_all,
                                 global_structural_novelty=global_novelty,
                                 gap_structural_opportunity=structural_opportunity,
                                 preprobe_score=score))
    output.sort(key=lambda row: (-row.preprobe_score, row.candidate_id))
    return output


def _minimum_distance(smiles: str, selected_smiles: list[str]) -> float:
    if not selected_smiles:
        return 1.0
    return min(1.0 - tanimoto(smiles, other) for other in selected_smiles)


def select_probe_shortlist(features: list[GapFeature], config: dict[str, Any]) -> list[dict[str, Any]]:
    probe = config["probe"]
    count = min(int(probe["candidate_count"]), len(features))
    if count <= 0:
        return []
    maximum_fraction = float(probe["maximum_family_fraction"])
    max_per_family = max(1, int(math.ceil(count * maximum_fraction)))
    diversity_weight = float(probe["selected_diversity_weight"])
    selected: list[GapFeature] = []
    selected_ids: set[str] = set()
    family_counts: dict[str, int] = {}
    while len(selected) < count:
        best: GapFeature | None = None
        best_key = None
        selected_smiles = [row.smiles for row in selected]
        for row in features:
            if row.candidate_id in selected_ids:
                continue
            if family_counts.get(row.family, 0) >= max_per_family:
                continue
            diversity = _minimum_distance(row.smiles, selected_smiles)
            dynamic = row.preprobe_score + diversity_weight * diversity
            key = (dynamic, diversity, row.gap_structural_opportunity, row.property_percentile, row.candidate_id)
            if best is None or key > best_key:
                best, best_key = row, key
        if best is None:
            break
        selected.append(best)
        selected_ids.add(best.candidate_id)
        family_counts[best.family] = family_counts.get(best.family, 0) + 1
    result = []
    prior: list[str] = []
    for index, row in enumerate(selected, start=1):
        distance = _minimum_distance(row.smiles, prior)
        result.append({**row.to_dict(), "probe_rank": index,
                       "minimum_distance_to_earlier_probe": distance,
                       "probe_score_at_selection": row.preprobe_score + diversity_weight * distance})
        prior.append(row.smiles)
    return result


def _identity_matches_candidate(candidate_identity: dict[str, Any], pubchem: dict[str, Any] | None,
                                *, mw_relative_tolerance: float) -> tuple[bool, list[str]]:
    reasons: list[str] = []
    if not pubchem:
        return False, ["pubchem_identity_not_resolved"]
    if str(pubchem.get("inchi_key") or "").upper() != str(candidate_identity["inchi_key"]).upper():
        reasons.append("inchi_key_mismatch")
    if pubchem.get("molecular_formula") and str(pubchem["molecular_formula"]) != str(candidate_identity["molecular_formula"]):
        reasons.append("molecular_formula_mismatch")
    p_mw = _finite_positive(pubchem.get("molecular_weight_g_mol"))
    c_mw = _finite_positive(candidate_identity.get("molecular_weight_g_mol"))
    if p_mw is not None and c_mw is not None and abs(p_mw / c_mw - 1.0) > float(mw_relative_tolerance):
        reasons.append("molecular_weight_mismatch")
    return not reasons, reasons


def _first_constant(backend: Any, cas: str, prop: str, allowed_fn: Any) -> tuple[str | None, float | None]:
    try:
        methods = sorted(str(x) for x in backend.constant_methods[prop](cas))
    except Exception:
        return None, None
    for method in methods:
        if not allowed_fn(prop, method):
            continue
        try:
            return method, float(backend.constant(cas, prop, method))
        except Exception:
            continue
    return None, None


def empirical_availability(*, candidate_identity: dict[str, Any], pubchem: dict[str, Any] | None,
                           backend: Any, allowed_fn: Any, structural_identity_fn: Any,
                           config: dict[str, Any]) -> dict[str, Any]:
    verified, reasons = _identity_matches_candidate(
        candidate_identity, pubchem,
        mw_relative_tolerance=float(config["identity"]["mw_relative_tolerance"]),
    )
    if not verified:
        return {"identity_verified": False, "identity_reasons": reasons,
                "selected_cas_number": None, "available_properties": [],
                "property_methods": {}, "cas_attempts": []}
    cas_numbers = sorted({str(value) for value in (pubchem or {}).get("cas_numbers", []) if value})
    attempts: list[dict[str, Any]] = []
    best: dict[str, Any] | None = None
    grid = [float(x) for x in config["property_probe"]["temperatures_k"]]
    for cas in cas_numbers:
        attempt = {"cas_number": cas, "identity_exact": False, "available_properties": [],
                   "property_methods": {}, "failure": None}
        try:
            meta = backend.identity(cas)
            if str(meta.get("inchi_key") or "").upper() != str(candidate_identity["inchi_key"]).upper():
                raise ValueError("local_cas_inchikey_mismatch")
            if structural_identity_fn(str(meta.get("smiles") or ""))["inchi_key"] != candidate_identity["inchi_key"]:
                raise ValueError("local_cas_structure_mismatch")
            attempt["identity_exact"] = True
            constants: dict[str, Any] = {}
            for prop in (TC, PC):
                method, value = _first_constant(backend, cas, prop, allowed_fn)
                if method is not None and value is not None:
                    attempt["available_properties"].append(prop)
                    attempt["property_methods"][prop] = [method]
                    constants[prop] = {"method": method, "value": value}
            tb_method, tb_value = _first_constant(backend, cas, "normal_boiling_temperature_k", allowed_fn)
            if tb_method is not None and tb_value is not None:
                constants["normal_boiling_temperature_k"] = {"method": tb_method, "value": tb_value}
            packet = {"molecular_weight_g_mol": float(candidate_identity["molecular_weight_g_mol"]),
                      "constants": constants}
            for prop in (PS, HV):
                try:
                    obj = backend.object(cas, prop, packet)
                    methods = sorted(method for method in getattr(obj, "all_methods", set()) if allowed_fn(prop, str(method)))
                except Exception:
                    methods = []
                usable_methods: list[str] = []
                for method in methods:
                    try:
                        if prop == HV and method in {"CRC_HVAP_TB", "CRC_HVAP_298"}:
                            backend.native_point(cas, prop, method, packet)
                            usable_methods.append(str(method))
                            continue
                        for temperature_k in grid:
                            tc_value = (constants.get(TC) or {}).get("value")
                            if tc_value is not None and temperature_k >= float(tc_value):
                                continue
                            try:
                                backend.series_value(cas, prop, str(method), temperature_k, packet)
                                usable_methods.append(str(method))
                                break
                            except Exception:
                                continue
                    except Exception:
                        continue
                if usable_methods:
                    attempt["available_properties"].append(prop)
                    attempt["property_methods"][prop] = sorted(set(usable_methods))
            attempt["available_properties"] = sorted(set(attempt["available_properties"]))
        except Exception as exc:
            attempt["failure"] = f"{type(exc).__name__}:{exc}"
        attempts.append(attempt)
        key = (int(attempt["identity_exact"]), len(attempt["available_properties"]), cas)
        if best is None or key > best["_sort_key"]:
            best = {**attempt, "_sort_key": key}
    if best is None or not best["identity_exact"]:
        return {"identity_verified": True, "identity_reasons": [], "selected_cas_number": None,
                "available_properties": [], "property_methods": {}, "cas_attempts": attempts}
    return {"identity_verified": True, "identity_reasons": [], "selected_cas_number": best["cas_number"],
            "available_properties": best["available_properties"], "property_methods": best["property_methods"],
            "cas_attempts": attempts}


def _deficit_gain(row: dict[str, Any], remaining: dict[str, dict[str, int]], config: dict[str, Any]) -> float:
    split = str(row["split"])
    available = set(row.get("available_properties") or [])
    weights = config["property_weights"]
    gain = 0.0
    for prop in config["priority_properties"]:
        if prop in available and int(remaining[prop][split]) > 0:
            gain += float(weights[prop])
    return gain


def _property_structural_gain(row: dict[str, Any], *, remaining: dict[str, dict[str, int]],
                              property_selected_anchors: dict[str, list[str]], config: dict[str, Any]) -> float:
    split = str(row["split"])
    available = set(row.get("available_properties") or [])
    weights = config["property_weights"]
    numerator = denominator = 0.0
    for prop in config["priority_properties"]:
        if prop not in available or int(remaining[prop][split]) <= 0:
            continue
        weight = float(weights[prop])
        numerator += weight * (1.0 - nearest_similarity(str(row["smiles"]), property_selected_anchors.get(prop, [])))
        denominator += weight
    return numerator / denominator if denominator > 0.0 else 0.0


def select_gap_targets(resolved: list[dict[str, Any]], *, deficits: dict[str, dict[str, int]],
                       anchors: dict[str, list[str]], config: dict[str, Any]) -> tuple[list[dict[str, Any]], dict[str, dict[str, int]]]:
    acquisition = config["acquisition"]
    target_count = int(acquisition["target_count"])
    max_per_family = max(1, int(math.ceil(target_count * float(acquisition["maximum_family_fraction"]))))
    weights = acquisition["weights"]
    eligible = [row for row in resolved if row.get("identity_verified") and row.get("selected_cas_number")
                and row.get("available_properties")]
    remaining = {prop: {"calibration": int(values["calibration"]), "evaluation": int(values["evaluation"])}
                 for prop, values in deficits.items()}
    selected: list[dict[str, Any]] = []
    selected_ids: set[str] = set()
    family_counts: dict[str, int] = {}
    property_selected_anchors = {prop: list(anchors.get(prop, [])) for prop in config["priority_properties"]}
    max_gain = max(1.0, sum(float(config["property_weights"][prop]) for prop in config["priority_properties"]))
    while len(selected) < target_count:
        if all(remaining[prop][split] <= 0 for prop in remaining for split in ("calibration", "evaluation")):
            break
        best: dict[str, Any] | None = None
        best_key = None
        selected_smiles = [str(row["smiles"]) for row in selected]
        for row in eligible:
            candidate_id = str(row["candidate_id"])
            if candidate_id in selected_ids:
                continue
            family = str(row["family"])
            if family_counts.get(family, 0) >= max_per_family:
                continue
            raw_gap_gain = _deficit_gain(row, remaining, config)
            if raw_gap_gain <= 0.0:
                continue
            gap_gain = raw_gap_gain / max_gain
            structural_gain = _property_structural_gain(
                row, remaining=remaining, property_selected_anchors=property_selected_anchors, config=config)
            diversity = _minimum_distance(str(row["smiles"]), selected_smiles)
            merit = float(row.get("property_percentile") or 0.0)
            dynamic = (float(weights["gap_closure"]) * gap_gain
                       + float(weights["property_structural_leverage"]) * structural_gain
                       + float(weights["property_merit"]) * merit
                       + float(weights["selected_diversity"]) * diversity)
            key = (dynamic, raw_gap_gain, structural_gain, merit, diversity, candidate_id)
            if best is None or key > best_key:
                best = {**row, "_dynamic_score": dynamic, "_raw_gap_gain": raw_gap_gain,
                        "_structural_gain": structural_gain, "_selected_diversity": diversity}
                best_key = key
        if best is None:
            break
        split = str(best["split"])
        covered: list[str] = []
        for prop in config["priority_properties"]:
            if prop in set(best.get("available_properties") or []) and remaining[prop][split] > 0:
                remaining[prop][split] -= 1
                covered.append(prop)
                property_selected_anchors[prop].append(str(best["smiles"]))
        selected_ids.add(str(best["candidate_id"]))
        family = str(best["family"])
        family_counts[family] = family_counts.get(family, 0) + 1
        selected.append({**{key: value for key, value in best.items() if not key.startswith("_")},
                         "acquisition_rank": len(selected) + 1,
                         "gap_properties_covered_at_selection": covered,
                         "gap_score_at_selection": best["_dynamic_score"],
                         "raw_weighted_gap_gain_at_selection": best["_raw_gap_gain"],
                         "property_structural_gain_at_selection": best["_structural_gain"],
                         "minimum_distance_to_earlier_selected": best["_selected_diversity"]})
    return selected, remaining


def gap_projection(deficits: dict[str, dict[str, int]], remaining: dict[str, dict[str, int]]) -> dict[str, Any]:
    output: dict[str, Any] = {}
    for prop, before in deficits.items():
        cal_gain = int(before["calibration"]) - int(remaining[prop]["calibration"])
        eval_gain = int(before["evaluation"]) - int(remaining[prop]["evaluation"])
        output[prop] = {
            "before": {"calibration": int(before["current_calibration"]), "evaluation": int(before["current_evaluation"])},
            "target": {"calibration": int(before["target_calibration"]), "evaluation": int(before["target_evaluation"])},
            "projected_selected_gain": {"calibration": cal_gain, "evaluation": eval_gain},
            "projected_after_if_all_selected_sources_audit_clean": {
                "calibration": int(before["current_calibration"]) + cal_gain,
                "evaluation": int(before["current_evaluation"]) + eval_gain,
            },
            "remaining_gap_after_selection": {"calibration": int(remaining[prop]["calibration"]),
                                              "evaluation": int(remaining[prop]["evaluation"])},
        }
    return output
