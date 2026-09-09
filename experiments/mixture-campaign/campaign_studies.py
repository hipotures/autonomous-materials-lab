"""Versioned, dependency-declaring studies. New studies backfill every eligible pair."""
from __future__ import annotations
from collections import Counter
from dataclasses import dataclass
import importlib
from importlib.metadata import distribution
import inspect
from pathlib import Path
from typing import Callable

import campaign_design as design
import campaign_quality as quality
import campaign_water as water
from campaign_store import implementation, digest, file_hash


@dataclass(frozen=True)
class Study:
    name: str
    version: str
    requires: tuple[str, ...]
    config_keys: tuple[str, ...]
    run: Callable[[dict, dict], dict]
    source_files: tuple[Path, ...] = ()
    packages: tuple[str, ...] = ()

    def code(self):
        source = inspect.getsourcefile(self.run)
        if not source:
            raise ValueError("study requires inspectable local source: " + self.name)
        packages = {}
        for name in self.packages:
            dist = distribution(name)
            hashes = {str(p): file_hash(Path(dist.locate_file(p))) for p in (dist.files or [])
                      if Path(dist.locate_file(p)).is_file() and "__pycache__" not in str(p) and not str(p).endswith((".pyc", ".pyo"))}
            if not hashes:
                raise ValueError("cannot fingerprint study dependency: " + name)
            packages[name] = {"version": dist.version, "files_sha256": digest(hashes)}
        return {"version": self.version, "files": implementation(Path(source), *self.source_files),
                "packages": packages}


def nodes(deps, settings):
    return design.paired(deps["observations"]["rows"], settings["models"],
                         phase_tolerance=settings["model"]["phase_fraction_tolerance"])


def regimes(deps, settings):
    report = design.comparisons(nodes(deps, settings), settings["adaptive"], settings["scope"])
    baseline = design.paired(deps["baseline"]["rows"], settings["models"],
                             phase_tolerance=settings["model"]["phase_fraction_tolerance"])
    good = [n for n in baseline if n["eligible"]]
    report["baseline_complete"] = len(good) == len(baseline)
    report["baseline_worst_ratio"] = min((n["minimum_ratio"] for n in good), default=None) if len(good) == len(baseline) else None
    report["baseline_point_count"] = len(baseline)
    return report


def phase_boundaries(deps, settings):
    analysis = design.edge_analysis(nodes(deps, settings), settings["adaptive"])
    edges = analysis["refinable_edges"]
    return {**analysis, "unresolved_triggered_edge_count": len(edges), "brackets": edges,
            "brackets_are_sampled_edges_not_certified_phase_boundaries": True,
            "parameter_validity_verified": False}


def model_disagreement(deps, settings):
    worst = {}
    for n in nodes(deps, settings):
        if not n["eligible"]: continue
        key = tuple(n["point"][1:])
        if key not in worst or n["spread"] > worst[key]["spread"]: worst[key] = n
    cases = []
    for key, n in sorted(worst.items()):
        cases.append({"point": n["point"], "spread": n["spread"],
                      "phase_signatures": n["phase_signatures"],
                      "models": {m: {"ratio": r["same_model_water_delta_h_ratio"],
                                     "outlet": r.get("outlet"), "storage_vle": r.get("storage_vle")}
                                 for m, r in n["models"].items()}})
    return {"cases": cases, "spread_is_calibrated_uncertainty": False,
            "state_status_counts": dict(sorted(Counter(r["status"] for r in deps["observations"]["rows"]).items()))}


def reference_audit(deps, settings):
    """Source-declared measurements; no synthetic model labels or automatic refit."""
    references = deps["references"]["records"]
    rows = deps["reference_predictions"]["rows"]
    output = []
    for ref in references:
        key = design.point(ref["additive_mass_fraction"], ref["outlet_temperature_k"], ref["pressure_pa"])
        for model in settings["models"]:
            matches = [r for r in rows if r["model"] == model and design.point(
                r["additive_mass_fraction"], r["outlet_temperature_k"], r["pressure_pa"]) == key and r.get("inlet_temperature_k") == ref["inlet_temperature_k"]]
            pred = matches[0] if matches else {}
            value = pred.get("delta_h_j_kg") if pred.get("status") == "ok" else None
            output.append({"reference_id": ref["id"], "source": ref["source"], "model": model,
                           "point": list(key), "reference_value": ref["value"], "predicted_value": value,
                           "status": "compared" if value is not None else pred.get("status", "outside_campaign_scope"),
                           "relative_error": value/ref["value"]-1 if value is not None else None,
                           "standard_uncertainty": ref.get("standard_uncertainty")})
    return {"status": "missing_external_references" if not references else
                       "compared_not_certified" if all(r["status"] == "compared" for r in output) else "incomplete_reference_comparison",
            "comparisons": output, "external_reference_count": len(references),
            "independence_status": "source_declared_training_overlap_unknown",
            "empirical_uncertainty_calibrated": False, "physical_tps_validation": False,
            "missing_studies": ["caloric_accuracy_domain_validation", "mixture_transport_validation",
                                "pore_flow_wetting_fouling", "matched_physical_cooling_test"]}


def registry(plugins=()):
    support = (Path(design.__file__), Path(quality.__file__), Path(water.__file__))
    reg = {
        "regimes": Study("regimes", "1", ("observations", "baseline"),
                         ("models", "model", "scope", "adaptive"), regimes, support),
        "phase_boundaries": Study("phase_boundaries", "2", ("observations",),
                                  ("models", "model", "adaptive"), phase_boundaries, support),
        "model_disagreement": Study("model_disagreement", "1", ("observations",),
                                    ("models", "model"), model_disagreement, support),
        "reference_audit": Study("reference_audit", "1", ("references", "reference_predictions"),
                                 ("models",), reference_audit, support),
    }
    reg["water_reference_audit"] = Study("water_reference_audit", "1", ("observations", "water_reference"),
        ("models", "model"), quality.water_reference_audit, support)
    reg["evidence_needs"] = Study("evidence_needs", "1", ("observations", "references", "phase_boundaries", "water_reference_audit"),
        ("models", "model", "scope", "adaptive", "evidence"), quality.evidence_needs, support)
    for name in plugins:
        old = dict(reg)
        importlib.import_module(name).register_studies(reg)
        if any(reg.get(k) is not v for k, v in old.items()):
            raise ValueError("plugins may not replace existing studies")
    return reg


def ordered_studies(required, reg):
    external = {"observations", "baseline", "references", "reference_predictions", "water_reference"}
    active, done, output = set(), set(), []
    def visit(name):
        if name in external or name in done: return
        if name not in reg: raise ValueError("unknown required study: " + name)
        if name in active: raise ValueError("study dependency cycle: " + name)
        active.add(name)
        for dep in reg[name].requires: visit(dep)
        active.remove(name); done.add(name); output.append(reg[name])
    for name in required: visit(name)
    return output
