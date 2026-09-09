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


def reference_blockers(node: dict) -> list[str]:
    reasons = set()
    for row in node["models"].values():
        reasons.update(row.get("comparison_blockers", []))
        if row.get("status") == "water_reference_unavailable":
            reasons.add("water_control_unavailable")
        if row.get("endpoint_within_configured_tolerance") is False:
            reasons.add("water_reference_tolerance_exceeded")
    return sorted(reasons)


def edge_analysis(nodes: list[dict], policy: dict) -> dict:
    """Separate geometric resolution, reference blocks and model disagreement.

    A persistent discrepancy is evidence-limited only on a small, same-regime
    edge with two discrepant endpoints. Independent phase/gain triggers remain
    refinable. No classification establishes continuum convergence or accuracy.
    """
    found, resolved, blocked, evidence = {}, [], [], []
    margin, tolerances = policy["gain_margin"], policy["tolerances"]
    local = policy.get("disagreement_resolution", {
        "mass_fraction": .01, "temperature_k": 5.0, "pressure_pa": .08})
    stable_change = policy.get("disagreement_stability_tolerance", .01)
    for axis in range(3):
        lines = defaultdict(list)
        for node in nodes:
            q = node["point"]
            lines[tuple(q[i] for i in range(3) if i != axis)].append(node)
        for line in lines.values():
            line.sort(key=lambda n: n["point"][axis])
            existing = {tuple(n["point"]) for n in line}
            for a, b in zip(line, line[1:]):
                lo, hi = a["point"][axis], b["point"][axis]
                span = math.log(hi/lo) if axis == 2 else hi-lo
                name = AXES[axis]
                mid = list(a["point"])
                mid[axis] = math.sqrt(lo*hi) if axis == 2 or (axis == 0 and hi <= policy["log_composition_transition"]) else (lo+hi)/2
                key = point(*mid)
                edge = {"axis": name, "left": a["point"], "right": b["point"],
                        "point": list(key), "normalized_span": span/tolerances[name]}
                reasons = []
                phase_change = a["phase_signatures"] != b["phase_signatures"]
                if phase_change:
                    reasons.append("phase_or_feasibility_change")
                blocks = sorted(set(reference_blockers(a) + reference_blockers(b)))
                if blocks:
                    blocked.append({**edge, "reasons": blocks})
                if a["eligible"] and b["eligible"]:
                    gains = [a["minimum_ratio"]-1-margin, b["minimum_ratio"]-1-margin]
                    if min(gains) <= 0 < max(gains):
                        reasons.append("paired_gain_boundary")
                    discrepant = max(a["spread"], b["spread"]) >= policy["model_spread_trigger"]
                    persistent = (discrepant and not phase_change and
                        min(a["spread"], b["spread"]) >= policy["model_spread_trigger"] and
                        span <= local[name] and abs(a["spread"]-b["spread"]) <= stable_change)
                    if persistent:
                        evidence.append({**edge, "reasons": ["model_disagreement_requires_evidence"],
                                         "endpoint_spreads": [a["spread"], b["spread"]]})
                    elif discrepant:
                        reasons.append("model_disagreement")
                    if max(gains) > 0 and abs(a["minimum_ratio"]-b["minimum_ratio"]) >= policy["gain_change_trigger"]:
                        reasons.append("gain_region_gradient")
                if not reasons:
                    continue
                if span <= tolerances[name] or key in existing or key in (tuple(a["point"]), tuple(b["point"])):
                    resolved.append({**edge, "reasons": reasons,
                                     "resolution_basis": "configured_sampled_edge_tolerance" if span <= tolerances[name] else "coordinate_precision_limit"})
                    continue
                priority = min(0 if r == "paired_gain_boundary" else 1 if r == "phase_or_feasibility_change" else 2 for r in reasons)
                item = {**edge, "reasons": reasons, "priority": priority}
                if key not in found or (priority, -edge["normalized_span"]) < (found[key]["priority"], -found[key]["normalized_span"]):
                    found[key] = item
    order = lambda e: (e.get("priority", 0), -e["normalized_span"], e["axis"], e["point"])
    return {"refinable_edges": sorted(found.values(), key=order),
            "resolved_edges": sorted(resolved, key=order),
            "reference_blocked_edges": sorted(blocked, key=order),
            "evidence_limited_edges": sorted(evidence, key=order),
            "reference_blocked_point_count": sum(bool(reference_blockers(n)) for n in nodes),
            "continuum_convergence_claimed": False}


def triggered_edges(nodes: list[dict], policy: dict) -> list[dict]:
    return edge_analysis(nodes, policy)["refinable_edges"]


def refinement_plan(nodes: list[dict], policy: dict, remaining: int) -> dict:
    analysis = edge_analysis(nodes, policy)
    edges = analysis["refinable_edges"]
    budget = max(0, min(policy["points_per_round"], remaining))
    chosen, used = [], set()
    # Respect scientific priority first; balance axes within each priority.
    for priority in sorted({e["priority"] for e in edges}):
        by_axis = {axis: [e for e in edges if e["axis"] == axis and e["priority"] == priority] for axis in AXES}
        while len(chosen) < budget and any(by_axis.values()):
            for axis in AXES:
                if by_axis[axis] and len(chosen) < budget:
                    e = by_axis[axis].pop(0)
                    if tuple(e["point"]) not in used:
                        chosen.append(e); used.add(tuple(e["point"]))
    return {"points": [e["point"] for e in chosen], "triggers": chosen,
            "unresolved_edge_count": len(edges), "candidate_point_count": len(edges),
            "resolved_edge_count": len(analysis["resolved_edges"]),
            "reference_blocked_edge_count": len(analysis["reference_blocked_edges"]),
            "reference_blocked_point_count": analysis["reference_blocked_point_count"],
            "evidence_limited_edge_count": len(analysis["evidence_limited_edges"]),
            "budget_limited": len(edges) > budget, "continuum_convergence_claimed": False}


def stop_assessment(plan: dict, remaining: int, rounds_exhausted: bool) -> dict:
    """Retain simultaneous blockers rather than hiding them behind one status."""
    causes = []
    if plan["unresolved_edge_count"]:
        causes.append("point_budget" if remaining <= 0 else "max_rounds" if rounds_exhausted else "refinement_pending")
    if plan.get("reference_blocked_point_count", 0):
        causes.append("reference_blocked")
    if plan.get("evidence_limited_edge_count", 0):
        causes.append("model_disagreement_requires_evidence")
    if not causes:
        causes = ["sampled_boundaries_resolved" if plan.get("resolved_edge_count", 0) else "no_triggered_sampled_edges"]
    return {"primary_reason": causes[0], "causes": causes,
            "refinable_edge_count": plan["unresolved_edge_count"],
            "resolved_edge_count": plan.get("resolved_edge_count", 0),
            "reference_blocked_point_count": plan.get("reference_blocked_point_count", 0),
            "evidence_limited_edge_count": plan.get("evidence_limited_edge_count", 0),
            "continuum_convergence_claimed": False}


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
