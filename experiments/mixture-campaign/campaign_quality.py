"""Common quality diagnostics and source-linked, machine-readable evidence needs."""
from __future__ import annotations

from collections import Counter
import math

import campaign_design as design
from campaign_store import digest


def water_reference_audit(deps, settings):
    bundle = deps["water_reference"]
    tolerance = settings["model"]["endpoint_relative_tolerance"]
    checks = []
    for item in bundle["controls"]:
        data = item["outcome"].get("data") or {}
        error = data.get("endpoint_relative_error")
        good = data.get("status") == "ok" and isinstance(error, (int, float)) and math.isfinite(error) and error <= tolerance
        checks.append({"context_id": item["context_id"], "coordinate": item["coordinate"],
            "task": item["task"], "status": "within_tolerance_not_validated" if good else
            "water_reference_tolerance_exceeded" if data.get("status") == "ok" else data.get("status", "water_control_task_failed"),
            "endpoint_relative_error": error, "comparison_allowed": good})
    nodes = design.paired(deps["observations"]["rows"], settings["models"],
                          phase_tolerance=settings["model"]["phase_fraction_tolerance"])
    blocked = [{"point": n["point"], "reasons": design.reference_blockers(n)}
               for n in nodes if design.reference_blockers(n)]
    return {"status": "reference_blocks_present" if blocked or bundle["context_errors"] or any(not c["comparison_allowed"] for c in checks)
                      else "sampled_controls_within_tolerance_not_validated",
            "control_count": len(checks), "checks": checks, "blocked_points": blocked,
            "context_errors": bundle["context_errors"], "model_contexts": bundle["model_contexts"],
            "blocked_point_count": len(blocked), "endpoint_relative_tolerance": tolerance,
            "status_counts": dict(sorted(Counter(c["status"] for c in checks).items())),
            "accuracy_validated": False, "reference_is_measurement": False,
            "missing_diagnostics": ["independent_mixture_caloric_validation"],
            "interpretation": "Shared pure-water model controls plus exact-context checks of cached mixture denominators. No correction or tolerance relaxation."}


def _request(kind, prop, unit, conditions, priority, route, reason, *, context_id=None):
    core = {"kind": kind, "property": prop, "unit": unit, "conditions": conditions,
            "water_context_id": context_id}
    return {"id": "need-" + digest(core)[:24], **core, "priority": priority,
            "acquisition_route": route, "reason": reason, "status": "pending_external_reference",
            "reference_value": None, "sources": [], "training_overlap_known": False,
            "satisfied": False}


def evidence_needs(deps, settings):
    """Generate requests, not fake measurements or a claim of physical readiness."""
    rows = deps["observations"]["rows"]
    nodes = design.paired(rows, settings["models"],
                          phase_tolerance=settings["model"]["phase_fraction_tolerance"])
    requests = {}
    def add(r):
        if r["id"] not in requests or r["priority"] < requests[r["id"]]["priority"]:
            requests[r["id"]] = r
    water = deps["water_reference_audit"]
    for check in water["checks"]:
        if check["comparison_allowed"]:
            continue
        t, p = check["coordinate"]
        conditions = {"water_cas": "7732-18-5", "composition_basis": "mass", "additive_mass_fraction": 0.,
                      "inlet_temperature_k": settings["model"]["inlet_temperature_k"],
                      "outlet_temperature_k": t, "pressure_pa": p}
        r = _request("water_model_review", "delta_h_j_kg", "J/kg", conditions, 0,
                     "reference_backend_review", check["status"], context_id=check["context_id"])
        r.update(status="model_review_required", control_task=check["task"])
        add(r)
    # Context mismatches may occur even when the standalone control passes.
    for row in rows:
        if "water_control_context_mismatch" not in row.get("comparison_blockers", []):
            continue
        cond = {"water_cas": "7732-18-5", "additive_cas": row["cas_number"],
                "inlet_temperature_k": row["inlet_temperature_k"],
                "outlet_temperature_k": row["outlet_temperature_k"], "pressure_pa": row["pressure_pa"]}
        r = _request("water_context_review", "delta_h_j_kg", "J/kg", cond, 0,
                     "reference_backend_review", "cached_binary_water_limit_differs_from_shared_control",
                     context_id=row.get("water_context_id"))
        r["status"] = "model_review_required"
        add(r)
    for model, error in sorted(water["context_errors"].items()):
        cond = {"pair_id": rows[0]["pair_id"] if rows else None, "model": model}
        r = _request("water_context_review", "model_context", "not_applicable", cond, 0,
                     "reference_backend_review", error)
        r["status"] = "model_review_required"
        add(r)
    selected, candidates = [], {}
    def candidate(q, priority, reason):
        key = design.point(*q)
        old = candidates.get(key)
        if old is None or priority < old[0]:
            candidates[key] = (priority, reason)
    valid = [n for n in nodes if n["eligible"]]
    gains = [n for n in valid if n["minimum_ratio"] > 1 + settings["adaptive"]["gain_margin"]]
    if gains:
        best = max(gains, key=lambda n: (n["minimum_ratio"], tuple(-v for v in n["point"])))
        candidate(best["point"], 1, "best_sampled_paired_gain_requires_validation")
    diagnostic = []
    for node in nodes:
        ratios = [r.get("same_model_water_delta_h_ratio") for r in node["models"].values()]
        if design.reference_blockers(node) and len(ratios) == len(settings["models"]) and all(
                r.get("status") == "ok" for r in node["models"].values()) and all(
                isinstance(v, (int, float)) and not isinstance(v, bool) and math.isfinite(v) and v > 0 for v in ratios):
            diagnostic.append((min(ratios), node["point"]))
    if diagnostic:
        score, q = max(diagnostic, key=lambda item: (item[0], tuple(-v for v in item[1])))
        if score > 1 + settings["adaptive"]["gain_margin"]:
            candidate(q, 3, "diagnostic_only_gain_blocked_by_water_reference")
    for edge in deps["phase_boundaries"].get("evidence_limited_edges", []):
        candidate(edge["point"], 2, "persistent_model_disagreement")
    # Retain at least one representative when no gain has been sampled.
    if nodes:
        center = [(settings["scope"][axis][0]+settings["scope"][axis][1])/2 for axis in design.AXES]
        spans = [settings["scope"][axis][1]-settings["scope"][axis][0] for axis in design.AXES]
        representative = min(valid or nodes, key=lambda n: (sum(((v-c)/span)**2 for v,c,span in zip(n["point"],center,spans)), n["point"]))
        candidate(representative["point"], 3, "representative_domain_validation")
    for edge in deps["phase_boundaries"].get("brackets", []):
        if "phase_or_feasibility_change" in edge["reasons"]:
            candidate(edge["point"], 4, "unresolved_sampled_phase_or_feasibility_transition")
    limit = settings.get("evidence", {}).get("max_points_per_pair", 8)
    selected = sorted(candidates, key=lambda q: (*candidates[q], q))[:limit]
    refs = deps["references"]["records"]
    for q in selected:
        priority, reason = candidates[q]
        cond = {"water_cas": "7732-18-5", "additive_cas": rows[0]["cas_number"],
                "composition_basis": "mass", "additive_mass_fraction": q[0],
                "inlet_temperature_k": settings["model"]["inlet_temperature_k"],
                "outlet_temperature_k": q[1], "pressure_pa": q[2]}
        r = _request("mixture_caloric_reference", "delta_h_j_kg", "J/kg", cond, priority,
                     "normalized_local_measurements", reason)
        matches = [ref for ref in refs if ref.get("evidence_kind") == "measured" and
            ref.get("property") == "delta_h_j_kg" and ref.get("unit") == "J/kg" and
            ref.get("water_cas") == cond["water_cas"] and ref.get("additive_cas") == cond["additive_cas"] and
            ref.get("composition_basis") == "mass" and ref.get("inlet_temperature_k") == cond["inlet_temperature_k"] and
            all(k in ref for k in ("additive_mass_fraction", "outlet_temperature_k", "pressure_pa")) and
            design.point(ref["additive_mass_fraction"], ref["outlet_temperature_k"], ref["pressure_pa"]) == q and
            isinstance(ref.get("source"), str) and ref["source"].strip() and
            isinstance(ref.get("value"), (int, float)) and not isinstance(ref["value"], bool) and math.isfinite(ref["value"]) and ref["value"] > 0]
        if matches:
            r["status"] = "source_linked_measurement_available_not_certified"
            r["sources"] = [{"id": ref["id"], "source": ref["source"], "value": ref["value"],
                             "standard_uncertainty": ref.get("standard_uncertainty")} for ref in sorted(matches, key=lambda x: str(x["id"]))]
        add(r)
        phase_conditions = {k: v for k, v in cond.items() if k != "inlet_temperature_k"}
        r = _request("mixture_phase_reference", "equilibrium_phase_split", "phase_fractions_and_compositions",
                     phase_conditions, priority+1, "phase_measurement_adapter_required", reason)
        r["status"] = "acquisition_adapter_required"
        add(r)
    if rows:
        r = _request("hardware_validation", "matched_cooling_and_pore_transport", "protocol_required",
                     {"water_cas": "7732-18-5", "additive_cas": rows[0]["cas_number"],
                      "scope": settings["scope"], "hardware_protocol": None}, 10,
                     "physical_measurements_required", "enthalpy_and_VLE_do_not_validate_wetting_flow_fouling_or_cooling")
        r["status"] = "hardware_protocol_required"
        add(r)
    ordered = sorted(requests.values(), key=lambda r: (r["priority"], r["id"]))
    return {"schema": "mixture-campaign-evidence-needs-v1", "status": "evidence_incomplete",
            "requests": ordered, "request_count": len(ordered),
            "selected_mixture_point_count": len(selected), "candidate_mixture_point_count": len(candidates),
            "unselected_mixture_point_count": max(0, len(candidates)-len(selected)),
            "queue_budget_limited": len(candidates) > limit,
            "status_counts": dict(sorted(Counter(r["status"] for r in ordered).items())),
            "reference_values_fabricated": False, "physical_validation_complete": False,
            "rankable_promotions": 0, "automatic_online_acquisition_performed": False}
