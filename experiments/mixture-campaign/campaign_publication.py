"""Immutable, lossless publication with bounded UTF-8 reports and revision diffs."""
from __future__ import annotations
from collections import Counter
import csv
import gzip
import hashlib
import io
import json
from pathlib import Path

import campaign_design as design
from campaign_store import atomic_bytes, encoded, digest, read_json, write_json


def compressed(output: Path, name: str, raw: bytes) -> dict:
    data = gzip.compress(raw, compresslevel=9, mtime=0)
    if gzip.decompress(data) != raw:
        raise ValueError("gzip round-trip failure")
    atomic_bytes(output / name, data)
    return {"file": name, "bytes": len(data), "uncompressed_bytes": len(raw),
            "sha256": hashlib.sha256(data).hexdigest(), "uncompressed_sha256": hashlib.sha256(raw).hexdigest()}


def csv_line(values) -> bytes:
    stream = io.StringIO(newline="")
    csv.writer(stream, lineterminator="\n").writerow(values)
    return stream.getvalue().encode()


def table(output: Path, name: str, rows: list[dict], limit: int) -> list[dict]:
    fields = sorted({k for row in rows for k in row}) or ["status"]
    header = csv_line(fields)
    full, block, count, records, shard = bytearray(header), bytearray(header), 0, [], 1
    def flush():
        nonlocal block, count, shard
        filename = f"report-{name}-{shard:04d}.csv"
        atomic_bytes(output / filename, bytes(block))
        records.append({"file": filename, "rows": count, "bytes": len(block),
                        "sha256": hashlib.sha256(block).hexdigest()})
        block, count, shard = bytearray(header), 0, shard+1
    for row in rows:
        line = csv_line([json.dumps(row.get(k), sort_keys=True, allow_nan=False, separators=(",", ":"))
                         if isinstance(row.get(k), (dict, list)) else row.get(k) for k in fields])
        if len(header)+len(line) > limit:
            raise ValueError("single CSV row exceeds report byte limit: " + name)
        if len(block)+len(line) > limit: flush()
        block.extend(line); full.extend(line); count += 1
    if count or not rows: flush()
    records.append(compressed(output, "table-"+name+".csv.gz", bytes(full)))
    return records


def scientific_diff(previous: dict | None, current: list[dict]) -> list[dict]:
    old = {r["pair_id"]: r for r in (previous or {}).get("pairs", [])}
    new = {r["pair_id"]: r for r in current}
    changes = []
    for key in sorted(old.keys() | new.keys()):
        a, b = old.get(key), new.get(key)
        changed = sorted(k for k in set(a or {}) | set(b or {}) if (a or {}).get(k) != (b or {}).get(k))
        changes.append({"pair_id": key, "change": "added" if a is None else "removed_or_out_of_scope" if b is None else "changed" if changed else "unchanged",
                        "changed_fields": changed, "before_qualification": a.get("qualification") if a else None,
                        "after_qualification": b.get("qualification") if b else None,
                        "before_best_ratio": a.get("best_paired_ratio") if a else None,
                        "after_best_ratio": b.get("best_paired_ratio") if b else None})
    return changes


def point_report(pair_result, config):
    rows = []
    for n in design.paired(pair_result["observations"]["rows"], config["models"],
                           phase_tolerance=config["model"]["phase_fraction_tolerance"]):
        row = {"pair_id": pair_result["pair"]["pair_id"], "mass_fraction": n["point"][0],
               "temperature_k": n["point"][1], "pressure_pa": n["point"][2],
               "paired_eligible": n["eligible"], "reference_blockers": design.reference_blockers(n), "minimum_ratio": n["minimum_ratio"], "model_spread": n["spread"],
               "heat_reduction_required_for_mass_parity": max(0., 1-n["minimum_ratio"]) if n["eligible"] else None,
               "mass_ratio_equal_net_heat": 1/n["minimum_ratio"] if n["eligible"] else None}
        for model, value in n["models"].items():
            state = value.get("outlet") or {}
            row.update({model+"_status": value["status"], model+"_ratio": value.get("same_model_water_delta_h_ratio"),
                        model+"_vapor_mole_fraction": state.get("vapor_mole_fraction"),
                        model+"_liquid_phase_count": state.get("liquid_phase_count"),
                        model+"_endpoint_error": value.get("endpoint_relative_error")})
        rows.append(row)
    return rows


def quality_artifacts(results, shared_water, config):
    """Deduplicate common needs without losing their affected-pair relationships."""
    requests, consumers = {}, set()
    for result in results:
        study = result["studies"].get("evidence_needs", {})
        for request in (study.get("data") or {}).get("requests", []):
            key = request["id"]
            if key in requests and requests[key] != request:
                raise ValueError("conflicting evidence requests for one identity: " + key)
            requests[key] = request
            consumers.add((key, result["pair"]["pair_id"]))
    items = sorted(requests.values(), key=lambda r: (r["priority"], r["id"]))
    needs = {"schema": "mixture-campaign-evidence-needs-v1", "requests": items,
             "consumers": [{"request_id": k, "pair_id": p} for k, p in sorted(consumers)],
             "reference_values_fabricated": False, "automatic_online_acquisition_performed": False}
    controls = []
    tolerance = config["model"]["endpoint_relative_tolerance"]
    for item in shared_water.get("controls", []):
        data = item["outcome"].get("data") or {}
        err = data.get("endpoint_relative_error")
        controls.append({"context_id": item["context_id"], "temperature_k": item["coordinate"][0],
            "pressure_pa": item["coordinate"][1], "status": data.get("status", "water_control_task_failed"),
            "model_delta_h_j_kg": data.get("model_delta_h_j_kg"), "heos_delta_h_j_kg": data.get("heos_delta_h_j_kg"),
            "endpoint_relative_error": err,
            "comparison_allowed": data.get("status") == "ok" and isinstance(err, (int, float)) and 0 <= err <= tolerance,
            "task_key": item["task"]["task_key"], "error": item["outcome"].get("error", data.get("error"))})
    info = {"water_context_count": len(shared_water.get("contexts", {})),
            "water_control_state_count": len(controls),
            "water_control_task_failure_count": sum(r["status"] == "water_control_task_failed" for r in controls),
            "water_control_blocked_count": sum(not r["comparison_allowed"] for r in controls),
            "reference_blocked_mixture_point_count": sum(r["summary"].get("reference_blocked_point_count") or 0 for r in results),
            "evidence_request_count": len(items), "evidence_consumer_count": len(consumers),
            "evidence_status_counts": dict(sorted(Counter(r["status"] for r in items).items())),
            "evidence_queue_budget_limited_pair_count": sum(bool((r["studies"].get("evidence_needs", {}).get("data") or {}).get("queue_budget_limited")) for r in results),
            "adaptive_stop_cause_counts": dict(sorted(Counter(c for r in results for c in r["summary"].get("adaptive_stop_causes", [])).items())),
            "continuum_convergence_claimed": False, "reference_values_fabricated": False}
    return info, needs, controls


def publish(output: Path, config: dict, results: list[dict], catalog: dict, graph: list[dict],
            events: list[dict], manifest: dict, previous: dict | None, *, shared_water: dict | None = None) -> dict:
    output.mkdir(parents=True)
    index, pairs = [], [r["summary"] for r in results]
    shared_water = shared_water or {"contexts": {}, "controls": []}
    quality, needs, controls = quality_artifacts(results, shared_water, config)
    scientific = {"schema": "mixture-campaign-results-v1", "pairs": pairs,
                  "rejections": catalog["rejections"], "required_suite": manifest["required_suite"],
                  "quality": quality, "water_controls_sha256": digest(shared_water), "evidence_needs_sha256": digest(needs)}
    scientific_hash = digest(scientific)
    summary = {**quality, "study_version": "v5m-3.1", "schema": "mixture-campaign-summary-v1", "run_id": manifest["run_id"],
               "scenario": config["scenario"], "catalog_pair_count": len(catalog["candidates"]),
               "processed_pair_count": len(pairs), "limited_run": manifest["limit_pairs"] is not None,
               "suite_execution_complete": bool(pairs) and all(r["suite_execution_complete"] for r in pairs),
               "run_status": "completed_with_recorded_statuses" if pairs else "no_eligible_pairs",
               "numerical_failure_count": sum(r["numerical_failure_count"] for r in pairs) + quality["water_control_task_failure_count"],
               "scientific_result_sha256": scientific_hash,
               "state_task_count": sum(t["spec"]["kind"] == "model.state" for t in graph),
               "executed_state_count": sum(e["kind"] == "model.state" and e["action"] == "executed" for e in events),
               "cached_state_count": sum(e["kind"] == "model.state" and e["action"] == "cache_hit" for e in events),
               "same_key_changed_outcome_count": sum(bool(e.get("changed_outcome_for_same_task")) for e in events),
               "qualification_counts": dict(sorted(Counter(r["qualification"] for r in pairs).items())),
               "adaptive_stop_counts": dict(sorted(Counter(r["adaptive_stop_reason"] for r in pairs).items())),
               "physical_validation_complete": False, "experimental_winners": 0, "rankable_promotions": 0,
               "interpretation": "Automatically executed common computational suite. Local gains are unvalidated model hypotheses. Budget stops and missing evidence remain explicit.",
               "pairs": pairs}
    changes = scientific_diff(previous, pairs)
    registry_cas = sorted(set(manifest["historical_cas"]) | {r["cas_number"] for r in catalog["candidates"]} |
                          {r["cas_number"] for r in catalog["rejections"]})
    plain = {"summary.json": summary, "manifest.json": manifest, "registry.json": {"cas_numbers": registry_cas},
             "changes.json": {"previous_run_id": manifest["previous_run_id"], "pairs": changes,
                              "identical_scientific_results": bool(previous and digest(previous) == scientific_hash)}}
    for filename, data in plain.items():
        raw = encoded(data)
        atomic_bytes(output / filename, raw)
        index.append({"file": filename, "bytes": len(raw), "sha256": hashlib.sha256(raw).hexdigest()})
        index.append(compressed(output, filename+".gz", raw))
    for filename, data in (("scientific-results.json.gz", scientific), ("task-graph.json.gz", graph),
                           ("catalog.json.gz", catalog), ("execution-events.json.gz", events),
                           ("water-controls.json.gz", shared_water), ("evidence-needs.json.gz", needs)):
        index.append(compressed(output, filename, encoded(data)))
    limit = int(config["publication"]["shard_bytes"])
    # A small cross-pair table precedes full per-pair reports. No failed points are omitted.
    brief = [{k: v for k, v in p.items() if k != "required_studies"} for p in pairs]
    index += table(output, "pair-summary", brief, limit)
    index += table(output, "changes", changes, limit)
    index += table(output, "catalog-rejections", catalog["rejections"], limit)
    counts = Counter((e["kind"], e["action"], e["status"]) for e in events)
    index += table(output, "task-audit", [{"kind": k, "action": a, "status": s, "count": n}
                                           for (k, a, s), n in sorted(counts.items())], limit)
    index += table(output, "water-controls", controls, limit)
    index += table(output, "evidence-needs", [{k: v for k, v in r.items() if k != "sources"} for r in needs["requests"]], limit)
    index += table(output, "evidence-consumers", needs["consumers"], limit)
    index += table(output, "evidence-sources", [{"request_id": r["id"], **source}
                   for r in needs["requests"] for source in r["sources"]], limit)
    pair_index = []
    for i, result in enumerate(results, 1):
        label = f"pair-{i:04d}"
        pair_index.append({"label": label, "pair_id": result["pair"]["pair_id"], "cas_number": result["pair"]["cas_number"]})
        index.append(compressed(output, f"{label}.json.gz", encoded(result)))
        index += table(output, label, point_report(result, config), limit)
        module_rows = [{"study": name, "execution_status": value["status"], "error": value.get("error"),
                        "evidence_status": (value.get("data") or {}).get("status")}
                       for name, value in sorted(result["studies"].items())]
        index += table(output, label+"-studies", module_rows, limit)
    # A shared baseline only; adaptive grids differ and are never globally ranked as equal evidence.
    from campaign_backend import legacy_import
    legacy_import()
    from v5m2_core import compare_models, pareto, THERMO_OBJECTIVES
    base_rows = [row for result in results for row in result["baseline"]["rows"]]
    comparisons = compare_models(base_rows, len(config["coarse"]["temperatures_k"])*len(config["coarse"]["pressures_pa"]), {})
    water = {"pair_id": "pure-water-reference", "name": "water", "additive_mass_fraction": 0.,
             "comparison_eligible": True, "minimum_paired_delta_h_ratio": 1.,
             "maximum_storage_bubble_pressure_ratio": 1., "maximum_viscosity_ratio_proxy": 1.,
             "evidence_kind": "normalization_identity_not_measurement"}
    index += table(output, "baseline-pareto-water", pareto(comparisons+[water], THERMO_OBJECTIVES), limit)
    write_json(output / "publication-index.json", {"schema": "mixture-campaign-publication-v1", "pairs": pair_index,
                "files": index, "gzip_connector_decode_assumed": False, "maximum_csv_shard_bytes": limit})
    return summary
