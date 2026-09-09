"""One standard protocol for every pair, with dynamic refinement and task reuse."""
from __future__ import annotations
from collections import Counter
from copy import deepcopy
from pathlib import Path

import campaign_backend as backend
import campaign_design as design
import campaign_config as config_module
import campaign_studies as studies
import campaign_water as water
from campaign_store import Store, Artifact, digest, implementation

HERE = Path(__file__).resolve().parent


class PairRun:
    def __init__(self, pair: dict, config: dict, store: Store, driver, references: list[dict], water_audit=None):
        self.pair, self.config, self.store, self.driver = pair, config, store, driver
        self.references = references
        self.water_audit = water_audit or water.SharedWaterAudit(store, driver, config)
        self.frozen = {}
        self.states = {}
        self.rounds = []
        self.workflow_code = implementation(Path(__file__), Path(config_module.__file__), Path(water.__file__))
        for model in config["models"]:
            self.frozen[model] = store.run("model.freeze", {"pair": backend.model_pair(pair, model),
                "model_config": backend.numerical_config(config)}, driver.code,
                lambda m=model: driver.freeze(pair, m, config), environment=driver.environment)

    def failed_row(self, model, q, error):
        return {"pair_id": self.pair["pair_id"], "cas_number": self.pair["cas_number"],
                "name": self.pair["name"], "smiles": self.pair["smiles"], "model": model,
                "additive_mass_fraction": q[0], "outlet_temperature_k": q[1], "pressure_pa": q[2],
                "inlet_temperature_k": self.config["model"]["inlet_temperature_k"],
                "status": "numerical_task_failed", "error": error,
                "model_comparison_eligible": False, "evidence_kind": "no_prediction"}

    def ensure(self, coordinates):
        coordinates = sorted(set(design.point(*q) for q in coordinates))
        for model in self.config["models"]:
            frozen = self.frozen[model]
            pending = []
            for q in coordinates:
                index = (model, q)
                if index in self.states: continue
                spec = self.store.spec("model.state", {"coordinate": list(q), "model": model,
                    "pair_id": self.pair["pair_id"]}, self.driver.code,
                    dependencies={"frozen": frozen}, environment=self.driver.environment)
                saved = self.store.load(spec)
                if saved is not None:
                    self.states[index] = saved
                elif frozen.outcome["status"] == "failed":
                    self.states[index] = self.store.put(spec, error="freeze_failed: " + frozen.outcome["error"])
                else:
                    pending.append((q, spec))
            step = self.config["execution"]["batch_points"]
            for start in range(0, len(pending), step):
                chunk = pending[start:start+step]
                try:
                    rows = self.driver.states(frozen.data, [list(q) for q, _ in chunk])
                    by_point = {design.point(r["additive_mass_fraction"], r["outlet_temperature_k"], r["pressure_pa"]): r for r in rows}
                    if len(rows) != len(chunk) or set(by_point) != {q for q, _ in chunk} or any(
                            r.get("model") != model or r.get("pair_id") != self.pair["pair_id"] for r in rows):
                        raise ValueError("mismatched worker state identities")
                    for q, spec in chunk:
                        self.states[(model, q)] = self.store.put(spec, by_point[q])
                except Exception as exc:
                    for q, spec in chunk:
                        self.states[(model, q)] = self.store.put(spec, error=f"{type(exc).__name__}: {exc}")

    def bundle(self, coordinates):
        coords = sorted(set(design.point(*q) for q in coordinates))
        deps, rows = {}, []
        for model in self.config["models"]:
            for q in coords:
                a = self.states[(model, q)]
                deps[f"{model}:{digest(q)}"] = a
                rows.append(a.data if a.outcome["status"] == "complete" else self.failed_row(model, q, a.outcome["error"]))
        control = self.water_audit.bundle(self.frozen, [q[1:] for q in coords])
        deps["water_reference"] = control
        return self.store.run("observations.bundle", {"pair_id": self.pair["pair_id"], "coordinates": coords,
                    "water_audit_policy": self.config.get("water_audit", {}),
                    "endpoint_tolerance": self.config["model"]["endpoint_relative_tolerance"]},
                    self.workflow_code, lambda: {"rows": water.annotate(rows, control.data, self.config)}, dependencies=deps)

    def run(self, ordered):
        self.water_audit.bundle(self.frozen, [], preflight=True)
        basic = config_module.basic_points(self.config)
        self.ensure(basic)
        baseline = self.bundle(basic)
        coords = set(basic)
        bundle = baseline
        policy = self.config["adaptive"]
        reason = "max_rounds"
        for round_index in range(policy["max_rounds"]):
            remaining = policy["max_new_points_per_pair"] - (len(coords)-len(basic))
            plan = self.store.run("adaptive.plan", {"pair_id": self.pair["pair_id"], "policy": policy,
                "models": self.config["models"], "phase_tolerance": self.config["model"]["phase_fraction_tolerance"],
                "remaining_budget": remaining}, implementation(Path(design.__file__)),
                lambda: design.refinement_plan(design.paired(bundle.data["rows"], self.config["models"],
                    phase_tolerance=self.config["model"]["phase_fraction_tolerance"]), policy, remaining),
                dependencies={"observations": bundle})
            if plan.outcome["status"] != "complete":
                reason = "planner_failed"; break
            self.rounds.append({"round": round_index+1, "task": plan.dependency(), **plan.data})
            if not plan.data["points"]:
                reason = design.stop_assessment(plan.data, remaining, False)["primary_reason"]
                break
            new = [design.point(*q) for q in plan.data["points"]]
            self.ensure(new)
            coords.update(new)
            bundle = self.bundle(coords)
            print(f"[refine] {self.pair['cas_number']} round={round_index+1} new_points={len(new)} total={len(coords)}", flush=True)
        remaining = policy["max_new_points_per_pair"] - (len(coords)-len(basic))
        assessment = self.store.run("adaptive.assessment", {"policy": policy, "remaining_budget": remaining,
            "models": self.config["models"], "phase_tolerance": self.config["model"]["phase_fraction_tolerance"],
            "rounds_exhausted": len(self.rounds) >= policy["max_rounds"]}, implementation(Path(design.__file__)),
            lambda: design.stop_assessment(design.refinement_plan(design.paired(bundle.data["rows"], self.config["models"],
                phase_tolerance=self.config["model"]["phase_fraction_tolerance"]), policy, remaining),
                remaining, len(self.rounds) >= policy["max_rounds"]), dependencies={"observations": bundle})
        if assessment.outcome["status"] != "complete":
            reason = "planner_failed"
        elif reason != "planner_failed":
            reason = assessment.data["primary_reason"]
        refs = [r for r in self.references if r["additive_cas"] == self.pair["cas_number"]]
        ref_art = self.store.run("references.freeze", {"pair_id": self.pair["pair_id"], "records": refs},
                                 self.workflow_code, lambda: {"records": refs})
        ref_points = [design.point(r["additive_mass_fraction"], r["outlet_temperature_k"], r["pressure_pa"]) for r in refs
                      if all(self.config["scope"][axis][0] <= v <= self.config["scope"][axis][1]
                             for axis, v in zip(design.AXES, [r["additive_mass_fraction"], r["outlet_temperature_k"], r["pressure_pa"]]))
                      and r["inlet_temperature_k"] == self.config["model"]["inlet_temperature_k"]]
        self.ensure(ref_points)
        ref_predictions = self.bundle(ref_points)
        artifacts = {"observations": bundle, "baseline": baseline, "references": ref_art, "reference_predictions": ref_predictions,
                     "water_reference": self.water_audit.bundle(self.frozen, [q[1:] for q in coords], preflight=True)}
        for study in ordered:
            deps = {name: artifacts[name] for name in study.requires}
            settings = {k: deepcopy(self.config[k]) for k in study.config_keys}
            artifacts[study.name] = self.store.run("study."+study.name,
                {"pair_id": self.pair["pair_id"], "settings": settings}, study.code(),
                lambda s=study, d=deps, cfg=settings: s.run({k: deepcopy(a.data) for k, a in d.items()}, deepcopy(cfg)), dependencies=deps)
        states = bundle.data["rows"]
        regime = artifacts.get("regimes")
        r = regime.data if regime and regime.outcome["status"] == "complete" else {}
        complete = reason != "planner_failed" and all(artifacts[s.name].outcome["status"] == "complete" for s in ordered)
        failed = sum(a.outcome["status"] == "failed" or a.data.get("status") == "model_state_rejected" for a in self.states.values())
        summary = {"pair_id": self.pair["pair_id"], "name": self.pair["name"], "cas_number": self.pair["cas_number"],
                   "suite_execution_complete": complete, "numerical_failure_count": failed,
                   "required_studies": {s.name: {"version": s.version, "status": artifacts[s.name].outcome["status"],
                                                "result_hash": artifacts[s.name].result_hash} for s in ordered},
                   "sampled_point_count": len(coords), "adaptive_new_point_count": len(coords)-len(basic),
                   "adaptive_stop_reason": reason, "adaptive_round_count": len(self.rounds),
                   "adaptive_stop_causes": assessment.data["causes"] if assessment.outcome["status"] == "complete" and reason != "planner_failed" else ["planner_failed"],
                   "continuum_convergence_claimed": False,
                   "best_point": r.get("best_point"), "best_paired_ratio": r.get("best_paired_ratio"),
                   "paired_gain_point_count": r.get("paired_gain_point_count", 0),
                   "baseline_worst_ratio": r.get("baseline_worst_ratio"),
                   "baseline_complete": r.get("baseline_complete", False),
                   "best_touches_scope_boundary": r.get("best_touches_scope_boundary", []),
                   "qualification": "study_failed" if not complete else "local_model_hypothesis_requires_validation" if r.get("paired_gain_point_count", 0) else
                                    "no_sampled_paired_gain" if r.get("eligible_point_count", 0) else "no_eligible_paired_states",
                   "experimental_winner": False, "physical_validation_complete": False, "rankable": False}
        evidence = artifacts.get("reference_audit")
        summary["reference_status"] = evidence.data["status"] if evidence and evidence.outcome["status"] == "complete" else "reference_audit_unavailable"
        quality = artifacts.get("water_reference_audit")
        needs = artifacts.get("evidence_needs")
        summary["water_reference_status"] = quality.data["status"] if quality and quality.outcome["status"] == "complete" else "water_audit_failed"
        summary["reference_blocked_point_count"] = quality.data["blocked_point_count"] if quality and quality.outcome["status"] == "complete" else None
        summary["evidence_request_count"] = needs.data["request_count"] if needs and needs.outcome["status"] == "complete" else None
        summary["evidence_limited_edge_count"] = assessment.data["evidence_limited_edge_count"] if assessment.outcome["status"] == "complete" else None
        return {"summary": summary, "pair": self.pair, "baseline": baseline.data, "observations": bundle.data,
                "rounds": self.rounds, "adaptive_assessment": assessment.outcome, "studies": {s.name: artifacts[s.name].outcome for s in ordered}}
