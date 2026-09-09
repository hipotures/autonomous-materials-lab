"""Deterministic adaptive studies on composition, temperature and pressure edges.

No posterior choice of operating requirements and no continuum-optimum claims.
"""
from __future__ import annotations
from collections import defaultdict
import math
from typing import Any

AXES = ("mass_fraction", "temperature_k", "pressure_pa")


def point(w: float, t: float, p: float) -> tuple[float, float, float]:
    return tuple(float(f"{float(v):.12g}") for v in (w, t, p))


def paired(rows: list[dict], models: list[str], *, phase_tolerance: float = 1e-8) -> list[dict]:
    grouped = defaultdict(dict)
    for row in rows:
        key = point(row["additive_mass_fraction"], row["outlet_temperature_k"], row["pressure_pa"])
        if row["model"] in grouped[key]:
            raise ValueError("duplicate model state")
        grouped[key][row["model"]] = row
    result = []
    for key, members in sorted(grouped.items()):
        good = (set(members) == set(models) and all(r.get("status") == "ok" and
                r.get("model_comparison_eligible") is True for r in members.values()))
        phases, ratios = {}, []
        for model in models:
            row = members.get(model, {})
            state = row.get("outlet") or {}
            phases[model] = [row.get("status", "missing"), state.get("liquid_phase_count"),
                             None if state.get("vapor_mole_fraction") is None else
                             state["vapor_mole_fraction"] > phase_tolerance]
            if good:
                r = row.get("same_model_water_delta_h_ratio")
                if isinstance(r, bool) or not isinstance(r, (int, float)) or not math.isfinite(r) or r <= 0:
                    raise ValueError("invalid eligible enthalpy ratio")
                ratios.append(r)
        result.append({"point": list(key), "eligible": good,
                       "minimum_ratio": min(ratios) if ratios else None,
                       "maximum_ratio": max(ratios) if ratios else None,
                       "spread": max(ratios)-min(ratios) if ratios else None,
                       "phase_signatures": phases, "models": members})
    return result


def triggered_edges(nodes: list[dict], policy: dict) -> list[dict]:
    """Detect adjacent collinear sampled edges, including unresolved phase changes.

    This cannot prove absence of islands between untriggered samples. Reported
    termination is only with respect to the configured sampled-edge rules.
    """
    found = {}
    margin = policy["gain_margin"]
    tolerances = policy["tolerances"]
    for axis in range(3):
        lines = defaultdict(list)
        for node in nodes:
            q = node["point"]
            lines[tuple(q[i] for i in range(3) if i != axis)].append(node)
        for line in lines.values():
            line.sort(key=lambda n: n["point"][axis])
            for a, b in zip(line, line[1:]):
                lo, hi = a["point"][axis], b["point"][axis]
                span = math.log(hi/lo) if axis == 2 else hi-lo
                tolerance = tolerances[AXES[axis]]
                if span <= tolerance:
                    continue
                reasons = []
                if a["phase_signatures"] != b["phase_signatures"]:
                    reasons.append("phase_or_feasibility_change")
                if a["eligible"] and b["eligible"]:
                    gains = [a["minimum_ratio"] - 1-margin, b["minimum_ratio"] - 1-margin]
                    if min(gains) <= 0 < max(gains):
                        reasons.append("paired_gain_boundary")
                    if max(a["spread"], b["spread"]) >= policy["model_spread_trigger"]:
                        reasons.append("model_disagreement")
                    if max(gains) > 0 and abs(a["minimum_ratio"]-b["minimum_ratio"]) >= policy["gain_change_trigger"]:
                        reasons.append("gain_region_gradient")
                if not reasons:
                    continue
                mid = list(a["point"])
                mid[axis] = math.sqrt(lo*hi) if axis == 2 or (axis == 0 and hi <= policy["log_composition_transition"]) else (lo+hi)/2
                key = point(*mid)
                if key in {tuple(n["point"]) for n in line} or key in (tuple(a["point"]), tuple(b["point"])):
                    continue
                priority = min([0 if r == "paired_gain_boundary" else 1 if r == "phase_or_feasibility_change" else 2 for r in reasons])
                edge = {"axis": AXES[axis], "left": a["point"], "right": b["point"],
                        "point": list(key), "reasons": reasons, "priority": priority,
                        "normalized_span": span/tolerance}
                if key not in found or (priority, -span/tolerance) < (found[key]["priority"], -found[key]["normalized_span"]):
                    found[key] = edge
    return sorted(found.values(), key=lambda e: (e["priority"], -e["normalized_span"], e["axis"], e["point"]))


def refinement_plan(nodes: list[dict], policy: dict, remaining: int) -> dict:
    existing = {tuple(n["point"]) for n in nodes}
    edges = [e for e in triggered_edges(nodes, policy) if tuple(e["point"]) not in existing]
    budget = min(policy["points_per_round"], remaining)
    # Round-robin axes avoid spending the entire budget on composition alone.
    by_axis = {axis: [e for e in edges if e["axis"] == axis] for axis in AXES}
    chosen, used = [], set()
    while len(chosen) < budget and any(by_axis.values()):
        for axis in AXES:
            if by_axis[axis] and len(chosen) < budget:
                e = by_axis[axis].pop(0)
                if tuple(e["point"]) not in used:
                    chosen.append(e); used.add(tuple(e["point"]))
    return {"points": [e["point"] for e in chosen], "triggers": chosen,
            "unresolved_edge_count": len(edges), "candidate_point_count": len(edges),
            "budget_limited": len(edges) > budget}


def comparisons(nodes: list[dict], policy: dict, scope: dict) -> dict:
    eligible = [n for n in nodes if n["eligible"]]
    gains = [n for n in eligible if n["minimum_ratio"] > 1 + policy["gain_margin"]]
    best = max(eligible, key=lambda n: (n["minimum_ratio"], tuple(-x for x in n["point"])), default=None)
    limits = [scope[k] for k in AXES]
    boundary = []
    if best and gains:
        boundary = [AXES[i] for i in range(3) if any(math.isclose(best["point"][i], v, abs_tol=1e-10, rel_tol=1e-10) for v in limits[i])]
    regimes = []
    grouped = defaultdict(list)
    for n in nodes:
        grouped[tuple(n["point"][1:])].append(n)
    for (t, p), group in sorted(grouped.items()):
        valid = [n for n in group if n["eligible"]]
        top = max(valid, key=lambda n: n["minimum_ratio"], default=None)
        regimes.append({"temperature_k": t, "pressure_pa": p, "sampled_compositions": len(group),
                        "eligible_compositions": len(valid),
                        "paired_gain_count": sum(n["minimum_ratio"] > 1 + policy["gain_margin"] for n in valid),
                        "best_mass_fraction": top["point"][0] if top else None,
                        "best_ratio": top["minimum_ratio"] if top else None})
    return {"sampled_point_count": len(nodes), "eligible_point_count": len(eligible),
            "paired_gain_point_count": len(gains),
            "best_point": best["point"] if best else None,
            "best_paired_ratio": best["minimum_ratio"] if best else None,
            "best_equal_net_heat_mass_ratio": 1/best["minimum_ratio"] if best else None,
            "best_touches_scope_boundary": boundary,
            "maximum_model_spread": max((n["spread"] for n in eligible), default=None),
            "regimes": regimes, "empirically_validated": False, "continuum_convergence_claimed": False}
