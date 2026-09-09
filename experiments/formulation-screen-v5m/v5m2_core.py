"""Deterministic design and reporting for a provisional binary-mixture grid.

This module contains no property database and never manufactures measurements.
"""
from __future__ import annotations

from collections import defaultdict
import hashlib
import json
import math
from typing import Any

VERSION = "v5m-2"
WATER_CAS = "7732-18-5"
MODELS = ("chemsep_nrtl", "unifac_dortmund")


def digest(value: Any) -> str:
    return hashlib.sha256(json.dumps(value, sort_keys=True, allow_nan=False,
                                    separators=(",", ":")).encode()).hexdigest()


def number(value: Any, label: str, *, positive: bool = False) -> float:
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise ValueError(f"{label}: expected a number")
    value = float(value)
    if not math.isfinite(value) or (positive and value <= 0):
        raise ValueError(f"{label}: nonfinite or nonpositive value")
    return value


def validate_config(c: dict) -> None:
    if not isinstance(c, dict) or c.get("version") != VERSION:
        raise ValueError("expected version v5m-2")
    if c.get("models") != list(MODELS):
        raise ValueError("this experiment requires both declared models, in order")
    catalog = c["catalog"]
    if catalog["source"] != "ChemSep NRTL":
        raise ValueError("unsupported catalog source")
    if not catalog["allowed_atomic_numbers"] or any(
        isinstance(x, bool) or not isinstance(x, int) or x <= 0
        for x in catalog["allowed_atomic_numbers"]
    ):
        raise ValueError("invalid element domain")
    number(catalog["maximum_molecular_weight_g_mol"], "maximum MW", positive=True)
    for label in ("maximum_pairs",):
        v = catalog[label]
        if v is not None and (isinstance(v, bool) or not isinstance(v, int) or v < 1):
            raise ValueError(f"invalid {label}")
    g = c["grid"]
    lo, transition, hi = [number(g[k], k, positive=True) for k in
                           ("minimum_mass_fraction", "transition_mass_fraction", "maximum_mass_fraction")]
    if not 0 < lo < transition < hi < 1:
        raise ValueError("require 0 < min < transition < max < 1")
    for k in ("log_points", "linear_points"):
        if isinstance(g[k], bool) or not isinstance(g[k], int) or not 2 <= g[k] <= 200:
            raise ValueError("grid point counts must be integers in [2, 200]")
    inlet = number(g["inlet_temperature_k"], "inlet T", positive=True)
    temps = g["outlet_temperatures_k"]
    pressures = g["pressures_pa"]
    for label, values in (("temperatures", temps), ("pressures", pressures)):
        if not isinstance(values, list) or not values:
            raise ValueError(f"empty {label}")
        nums = [number(v, label, positive=True) for v in values]
        if nums != sorted(set(nums)):
            raise ValueError(f"{label} must be sorted and unique")
    if inlet < 273.15 or any(not inlet < t <= 500 for t in temps):
        raise ValueError("operational domain: 273.15 <= inlet < outlet <= 500 K")
    if any(not 25000 <= p <= 300000 for p in pressures):
        raise ValueError("ideal-gas vapor screening limited to 25..300 kPa")
    for k in ("workers", "pair_timeout_s", "refinement_pairs", "refinement_intervals_per_pair",
              "shortlist_size"):
        v = c[k]
        if isinstance(v, bool) or not isinstance(v, int) or v < 1:
            raise ValueError(f"{k} must be a positive integer")
    if c["refinement_intervals_per_pair"] > 50:
        raise ValueError("refinement interval cap is 50")
    for k in ("phase_fraction_tolerance", "balance_tolerance"):
        if not 0 < number(c[k], k, positive=True) <= 1e-4:
            raise ValueError(f"invalid {k}")
    if not 0 < number(c["endpoint_relative_tolerance"], "endpoint tolerance") < .2:
        raise ValueError("invalid endpoint tolerance")


def mass_grid(c: dict) -> list[float]:
    g = c["grid"]
    lo, mid, hi = (g[k] for k in ("minimum_mass_fraction", "transition_mass_fraction", "maximum_mass_fraction"))
    values = [math.exp(math.log(lo) + i * math.log(mid/lo)/(g["log_points"]-1))
              for i in range(g["log_points"])]
    values += [mid + i*(hi-mid)/(g["linear_points"]-1) for i in range(g["linear_points"])]
    return sorted({float(f"{v:.12g}") for v in values})


def mass_to_mole(w: float, mw: list[float]) -> list[float]:
    w = number(w, "mass fraction")
    if not 0 <= w <= 1 or len(mw) != 2:
        raise ValueError("binary mass fraction outside [0, 1]")
    ns = [(1-w)/number(mw[0], "water MW", positive=True), w/number(mw[1], "additive MW", positive=True)]
    total = sum(ns)
    return [n/total for n in ns]


def mole_to_mass(z: list[float], mw: list[float]) -> list[float]:
    if len(z) != 2 or len(mw) != 2 or any(number(x, "mole fraction") < 0 for x in z):
        raise ValueError("invalid binary mole fractions")
    if abs(sum(z)-1) > 1e-8:
        raise ValueError("mole fractions do not sum to one")
    ms = [x*number(m, "MW", positive=True) for x, m in zip(z, mw)]
    return [v/sum(ms) for v in ms]


def summarize_flash(result: Any, zs: list[float], mw: list[float], tolerance: float) -> dict:
    """Validate phase/species balances; phase fractions returned by thermo are molar."""
    phases, betas = list(result.phases), list(result.betas)
    if len(phases) != len(betas) or not phases:
        raise ValueError("invalid phase list")
    betas = [number(b, "phase fraction") for b in betas]
    if any(b < 0 or b > 1 for b in betas) or abs(sum(betas)-1) > tolerance:
        raise ValueError("phase fraction balance failed")
    for ph in phases:
        if len(ph.zs) != len(zs) or any(number(x, "phase composition") < 0 for x in ph.zs):
            raise ValueError("invalid phase composition")
        if abs(sum(ph.zs)-1) > tolerance:
            raise ValueError("phase composition sum failed")
    error = max(abs(sum(b*ph.zs[i] for b, ph in zip(betas, phases))-zs[i]) for i in range(2))
    if error > tolerance:
        raise ValueError("species balance failed")
    M = sum(z*m for z, m in zip(zs, mw))/1000
    h = number(float(result.H())/M, "specific enthalpy")
    volume = number(float(result.V()), "molar volume", positive=True)
    liquid_count = len(result.liquids)
    gas = result.gas
    vap = number(float(result.VF), "vapor fraction")
    if not 0 <= vap <= 1:
        raise ValueError("invalid vapor fraction")
    cp = None
    if len(phases) == 1:
        cp = number(float(result.Cp())/M, "single-phase Cp", positive=True)
    return {"enthalpy_j_kg": h, "density_kg_m3": M/volume, "cp_j_kg_k": cp,
            "vapor_mole_fraction": vap, "liquid_phase_count": liquid_count,
            "phase_count": len(phases), "component_balance_error": error,
            "vapor_additive_mole_fraction": float(gas.zs[1]) if gas is not None else None,
            "vapor_additive_mass_fraction": mole_to_mass(list(gas.zs), mw)[1] if gas is not None else None,
            "liquid_compositions_mole": [list(ph.zs) for ph in result.liquids]}


def compare_models(rows: list[dict], expected_conditions: int, c: dict) -> list[dict]:
    """Only like-for-like paired predictions enter numerical tradeoffs.

    Model spread is NOT an empirical uncertainty interval. Failed states are
    retained and prevent complete-grid eligibility.
    """
    grouped = defaultdict(list)
    for row in rows:
        grouped[(row["pair_id"], row["additive_mass_fraction"])].append(row)
    out = []
    for (pair, w), values in sorted(grouped.items()):
        by_condition = defaultdict(dict)
        for row in values:
            key = (row["outlet_temperature_k"], row["pressure_pa"])
            if row["model"] in by_condition[key]:
                raise ValueError("duplicate model/condition result")
            by_condition[key][row["model"]] = row
        paired = []
        failures = defaultdict(int)
        for condition, models in sorted(by_condition.items()):
            good = (set(models) == set(MODELS) and all(
                r["status"] == "ok" and r.get("model_comparison_eligible") is True for r in models.values()))
            for r in models.values():
                if r["status"] != "ok":
                    failures[r["status"]] += 1
                elif not r.get("model_comparison_eligible"):
                    failures["model_gates_not_passed"] += 1
            if not good:
                continue
            ratios = [number(r["same_model_water_delta_h_ratio"], "enthalpy ratio", positive=True) for r in models.values()]
            paired.append({"condition": condition, "min_h_ratio": min(ratios),
                           "spread_h_ratio": abs(ratios[0]-ratios[1]),
                           "max_storage_psat_ratio": max(r["storage_bubble_pressure_ratio"] for r in models.values()),
                           "max_viscosity_ratio_proxy": max(r["viscosity_ratio_proxy"] for r in models.values())
                               if all(r.get("viscosity_ratio_proxy") is not None for r in models.values()) else None})
        complete = len(paired) == expected_conditions
        first = values[0]
        out.append({"pair_id": pair, "cas_number": first["cas_number"], "name": first["name"],
                    "smiles": first["smiles"], "additive_mass_fraction": w,
                    "paired_condition_count": len(paired), "expected_condition_count": expected_conditions,
                    "complete_grid": complete, "comparison_eligible": complete,
                    "minimum_paired_delta_h_ratio": min((p["min_h_ratio"] for p in paired), default=None),
                    "maximum_model_spread_delta_h_ratio": max((p["spread_h_ratio"] for p in paired), default=None),
                    "maximum_storage_bubble_pressure_ratio": max((p["max_storage_psat_ratio"] for p in paired), default=None),
                    "maximum_viscosity_ratio_proxy": max(p["max_viscosity_ratio_proxy"] for p in paired)
                        if paired and all(p["max_viscosity_ratio_proxy"] is not None for p in paired) else None,
                    "failure_counts": dict(failures), "calibrated_uncertainty": None,
                    "model_spread_is_uncertainty": False, "experimental_winner": False,
                    "handling_approved": False})
    return out


THERMO_OBJECTIVES = {"minimum_paired_delta_h_ratio": "max", "maximum_storage_bubble_pressure_ratio": "min"}
FLOW_OBJECTIVES = {**THERMO_OBJECTIVES, "maximum_viscosity_ratio_proxy": "min"}


def pareto(rows: list[dict], objectives: dict[str, str]) -> list[dict]:
    usable = [r for r in rows if r.get("comparison_eligible") is True and all(
        isinstance(r.get(k), (int, float)) and not isinstance(r[k], bool) and math.isfinite(r[k]) for k in objectives)]
    def dominates(a, b):
        diffs = [(a[k]-b[k])*(1 if direction == "max" else -1) for k, direction in objectives.items()]
        return all(d >= 0 for d in diffs) and any(d > 0 for d in diffs)
    return [r for r in usable if not any(dominates(s, r) for s in usable)]


def refinement_grid(base: list[float], candidate_rows: list[dict], c: dict) -> list[float]:
    """Add midpoints around non-dominated compositions and model-disagreement peaks.

    No grid-convergence claim follows from one refinement round.
    """
    front = pareto(candidate_rows, THERMO_OBJECTIVES)
    valid = [r for r in candidate_rows if r.get("comparison_eligible")]
    disagreement = sorted(valid, key=lambda r: (-r["maximum_model_spread_delta_h_ratio"], r["additive_mass_fraction"]))
    seeds = sorted(front, key=lambda r: (-r["minimum_paired_delta_h_ratio"], r["additive_mass_fraction"])) + disagreement[:2]
    intervals = []
    for r in seeds:
        w = r["additive_mass_fraction"]
        i = base.index(w)
        for j in (i-1, i):
            if 0 <= j < len(base)-1 and j not in intervals:
                intervals.append(j)
    intervals = intervals[:c["refinement_intervals_per_pair"]]
    transition = c["grid"]["transition_mass_fraction"]
    return sorted({float(f"{(math.sqrt(base[i]*base[i+1]) if base[i+1] <= transition else (base[i]+base[i+1])/2):.12g}")
                   for i in intervals} - set(base))


def shortlist(rows: list[dict], c: dict) -> list[dict]:
    """Separate improvement and disagreement queues; max one formulation per pair/queue."""
    front = pareto(rows, THERMO_OBJECTIVES)
    candidates = [("thermodynamic_tradeoff", r) for r in sorted(
        front, key=lambda r: (-r["minimum_paired_delta_h_ratio"], r["pair_id"], r["additive_mass_fraction"]))]
    candidates += [("model_disagreement", r) for r in sorted(
        [r for r in rows if r.get("comparison_eligible")],
        key=lambda r: (-r["maximum_model_spread_delta_h_ratio"], r["pair_id"], r["additive_mass_fraction"]))]
    seen, result = set(), []
    quota = max(1, c["shortlist_size"]//2)
    counts = defaultdict(int)
    for role, row in candidates:
        key = (role, row["pair_id"])
        if key in seen or counts[role] >= quota:
            continue
        seen.add(key); counts[role] += 1
        result.append({**row, "acquisition_role": role,
                       "reason": "unvalidated model tradeoff" if role == "thermodynamic_tradeoff" else "model discrepancy to resolve",
                       "requires": "identity, handling, mixture-property and matched apparatus validation"})
    return result[:c["shortlist_size"]]
