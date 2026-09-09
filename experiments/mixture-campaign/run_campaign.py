#!/usr/bin/env python3
"""Discover mixtures, execute their common suite, refine and publish in one command."""
from __future__ import annotations
import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import datetime, timezone
from importlib.metadata import version
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import uuid

import yaml
import campaign_backend as backend
import campaign_config as cc
import campaign_design as design
import campaign_publication as publication
import campaign_studies as studies
import campaign_workflow as workflow
from campaign_store import Store, digest, encoded, implementation, read_json, write_json, exclusive_lock

HERE = Path(__file__).resolve().parent


def references(path: Path | None) -> list[dict]:
    if path is None: return []
    payload = read_json(path)
    if payload.get("schema") != "mixture-campaign-measurements-v1":
        raise ValueError("expected normalized mixture-campaign-measurements-v1")
    ids, records = set(), []
    for row in payload["observations"]:
        if row["id"] in ids: raise ValueError("duplicate reference id")
        ids.add(row["id"])
        if row.get("evidence_kind") != "measured" or not isinstance(row.get("source"), str) or not row["source"].strip():
            raise ValueError("reference requires source-linked measured evidence")
        if row.get("property") != "delta_h_j_kg" or row.get("unit") != "J/kg" or row.get("composition_basis") != "mass":
            raise ValueError("initial reference adapter supports same-composition delta_h in J/kg only")
        if row.get("water_cas") != "7732-18-5" or not isinstance(row.get("additive_cas"), str):
            raise ValueError("exact binary component identifiers required")
        for k in ("additive_mass_fraction", "inlet_temperature_k", "outlet_temperature_k", "pressure_pa", "value"):
            cc.numeric(row[k], k, 0)
        if not 0 < row["additive_mass_fraction"] < 1 or row["value"] <= 0 or not 0 < row["inlet_temperature_k"] < row["outlet_temperature_k"] or row["pressure_pa"] <= 0:
            raise ValueError("invalid reference state")
        if row.get("standard_uncertainty") is not None: cc.numeric(row["standard_uncertainty"], "uncertainty", 0)
        records.append(row)
    return sorted(records, key=lambda r: str(r["id"]))


def latest_snapshot(root: Path):
    pointer = root / "latest.json"
    if not pointer.is_file(): return None, None, []
    record = read_json(pointer)
    name = record["run_id"]
    if not isinstance(name, str) or Path(name).name != name or not name.startswith("run-"):
        raise ValueError("invalid latest snapshot pointer")
    folder = root / name
    old = read_json(folder / "scientific-results.json.gz")
    if digest(old) != record["scientific_result_sha256"]:
        raise ValueError("previous scientific snapshot hash mismatch")
    return name, old, read_json(folder / "registry.json")["cas_numbers"]


def execute(config: dict, cache: Path, results: Path, *, retry_failures=False, recompute=False,
            limit_pairs=None, driver_factory=backend.Driver, reference_records=None):
    cc.validate(config)
    registry = studies.registry(config["plugins"])
    ordered = studies.ordered_studies(config["studies"], registry)
    if not {"regimes", "phase_boundaries", "model_disagreement", "reference_audit"}.issubset({s.name for s in ordered}):
        raise ValueError("the initial standard suite cannot omit a required baseline study")
    previous_id, previous, historical = latest_snapshot(results)
    run_id = "run-" + datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ") + "-" + uuid.uuid4().hex[:8]
    work = cache / "work" / run_id
    work.mkdir(parents=True)
    driver = driver_factory(work / "workers", config["execution"]["worker_timeout_s"])
    store = Store(cache, work / "events.jsonl", retry_failures=retry_failures, recompute=recompute)
    print("[catalog] discovering current and historical pairs", flush=True)
    catalog_task = store.run("catalog.discover", {"policy": config["catalog"], "additional_cas": config["additional_cas"],
                            "historical_cas": historical, "inlet": config["model"]["inlet_temperature_k"]},
                            driver.code, lambda: driver.catalog(config, historical), environment=driver.environment)
    catalog = catalog_task.data
    candidates = sorted(catalog["candidates"], key=lambda r: r["pair_id"])
    ids = [r["pair_id"] for r in candidates]
    if len(set(ids)) != len(ids): raise ValueError("duplicate catalog pair identity")
    selected = candidates[:limit_pairs] if limit_pairs else candidates
    suite = [{"name": s.name, "version": s.version, "requires": s.requires, "implementation": s.code()} for s in ordered]
    write_json(work / "plan.json", {"config": config, "run_id": run_id, "required_suite": suite,
                                  "pair_ids": [p["pair_id"] for p in selected]})
    print(f"[design] pairs={len(selected)} coarse_points_per_pair={len(cc.basic_points(config))} models={len(config['models'])} max_new_points_per_pair={config['adaptive']['max_new_points_per_pair']}", flush=True)
    completed = {}
    def pair_job(pair):
        return workflow.PairRun(pair, config, store, driver, reference_records or []).run(ordered)
    with ThreadPoolExecutor(max_workers=config["execution"]["workers"]) as pool:
        futures = {pool.submit(pair_job, p): p for p in selected}
        for future in as_completed(futures):
            pair = futures[future]
            result = future.result()
            completed[pair["pair_id"]] = result
            write_json(work / (digest(pair["pair_id"])+".json"), result)
            s = result["summary"]
            print(f"[suite {len(completed)}/{len(selected)}] {pair['name']} points={s['sampled_point_count']} gains={s['paired_gain_point_count']} stop={s['adaptive_stop_reason']} studies_ok={s['suite_execution_complete']}", flush=True)
    ordered_results = [completed[p["pair_id"]] for p in selected]
    try:
        git_commit = subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=HERE, text=True, stderr=subprocess.DEVNULL).strip()
        git_dirty = bool(subprocess.check_output(["git", "status", "--porcelain"], cwd=HERE, text=True, stderr=subprocess.DEVNULL).strip())
    except (OSError, subprocess.CalledProcessError):
        git_commit, git_dirty = None, None
    manifest = {"schema": "mixture-campaign-manifest-v1", "run_id": run_id, "previous_run_id": previous_id,
                "git_commit": git_commit, "git_dirty": git_dirty, "config": config,
                "required_suite": suite, "environment": driver.environment, "backend_implementation": driver.code,
                "driver_implementation": implementation(Path(__file__), Path(workflow.__file__), Path(cc.__file__)),
                "publication_implementation": implementation(Path(publication.__file__)),
                "references_sha256": digest(reference_records or []), "historical_cas": historical,
                "limit_pairs": limit_pairs, "retry_failures": retry_failures, "recompute": recompute,
                "reference_training_overlap_known": False, "old_v5m2_results_imported_as_cache": False}
    stage = work / "publication"
    summary = publication.publish(stage, config, ordered_results, catalog, store.graph(), store.events, manifest, previous)
    results.mkdir(parents=True, exist_ok=True)
    destination = results / run_id
    # Stage on the destination filesystem to retain atomic snapshot publication.
    staged = results / (".building-"+run_id)
    shutil.copytree(stage, staged)
    os.replace(staged, destination)
    if limit_pairs is None:
        write_json(results / "latest.json", {"run_id": run_id, "scientific_result_sha256": summary["scientific_result_sha256"]})
    print("[published] " + str(destination), flush=True)
    print(json.dumps({k: v for k, v in summary.items() if k != "pairs"}, indent=2, sort_keys=True), flush=True)
    return summary, destination


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, default=HERE / "campaign.yaml")
    parser.add_argument("--cache-dir", type=Path, default=HERE / ".campaign-cache")
    parser.add_argument("--results-dir", type=Path, default=HERE / "campaign-v5m3-results")
    parser.add_argument("--workers", type=int)
    parser.add_argument("--limit-pairs", type=int, help="Smoke subset; does not update latest full-campaign snapshot.")
    parser.add_argument("--retry-failures", action="store_true")
    parser.add_argument("--recompute", action="store_true", help="Execute tasks again; retain old immutable outcomes and report same-key differences.")
    args = parser.parse_args()
    try:
        if args.limit_pairs is not None and args.limit_pairs < 1: raise ValueError("limit-pairs must be positive")
        config = yaml.safe_load(args.config.read_text(encoding="utf-8"))
        if args.workers is not None: config["execution"]["workers"] = args.workers
        path = Path(config["references"]) if config.get("references") else None
        if path is not None and not path.is_absolute(): path = args.config.resolve().parent / path
        records = references(path)
        cache, results = args.cache_dir.resolve(), args.results_dir.resolve()
        if cache == results or cache in results.parents or results in cache.parents:
            raise ValueError("cache and publication directories must not contain each other")
        with exclusive_lock(cache), exclusive_lock(results):
            summary, _ = execute(config, cache, results, retry_failures=args.retry_failures,
                                 recompute=args.recompute, limit_pairs=args.limit_pairs, reference_records=records)
        return 1 if summary["numerical_failure_count"] or not summary["suite_execution_complete"] else 0
    except (ValueError, KeyError, TypeError, OSError, RuntimeError, ImportError) as exc:
        parser.error(str(exc))


if __name__ == "__main__":
    raise SystemExit(main())
