#!/usr/bin/env python3
"""V5m-2: automatically catalog and grid-screen water/additive binary liquids."""
from __future__ import annotations

import argparse
from collections import Counter
from concurrent.futures import ThreadPoolExecutor, as_completed
import csv
import hashlib
from importlib import import_module
from importlib.metadata import version
import json
import os
from pathlib import Path
import platform
import subprocess
import sys
from typing import Any

import yaml
from v5m2_core import (VERSION, MODELS, THERMO_OBJECTIVES, FLOW_OBJECTIVES, digest, validate_config,
                       mass_grid, compare_models, pareto, refinement_grid, shortlist)

HERE = Path(__file__).resolve().parent


def load_json(path: Path) -> Any:
    def reject(v): raise ValueError("nonfinite JSON constant: " + v)
    return json.loads(path.read_text(encoding="utf-8"), parse_constant=reject)


def write_json(path: Path, value: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(path.suffix + ".tmp")
    tmp.write_text(json.dumps(value, indent=2, sort_keys=True, allow_nan=False) + "\n", encoding="utf-8")
    tmp.replace(path)


def write_csv(path: Path, rows: list[dict]) -> None:
    fields = sorted({k for r in rows for k in r}) or ["status"]
    tmp = path.with_suffix(".csv.tmp")
    with tmp.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields); writer.writeheader()
        for r in rows:
            writer.writerow({k: json.dumps(v, sort_keys=True, allow_nan=False) if isinstance(v, (dict, list)) else v for k, v in r.items()})
    tmp.replace(path)


def file_hash(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024*1024), b""):
            h.update(block)
    return h.hexdigest()


def checkpoint(path: Path, signature: str, value: Any = None) -> Any:
    if value is not None:
        write_json(path, {"signature": signature, "sha256": digest(value), "data": value})
        return value
    if not path.exists():
        return None
    saved = load_json(path)
    if saved.get("signature") != signature or saved.get("sha256") != digest(saved.get("data")):
        raise ValueError("checkpoint incompatible or corrupt: " + str(path))
    return saved["data"]


def contract(config: dict, *, catalog_only: bool, limit_pairs: int | None) -> dict:
    pins = {"thermo": "0.6.1", "chemicals": "1.5.2", "CoolProp": "8.0.0", "rdkit": "2026.3.6"}
    packages = {k: version(k) for k in pins}
    for k, v in pins.items():
        if packages[k] != v:
            raise ValueError(f"use project environment: requires {k}=={v}, installed {packages[k]}")
    for k in ("numpy", "scipy", "pandas", "fluids", "PyYAML"):
        packages[k] = version(k)
    code = {p.name: file_hash(p) for p in [HERE / "run_predictive_grid.py", HERE / "v5m2_core.py", HERE / "v5m2_backend.py"]}
    data_hashes = {}
    # Preserve database and implementation content, not version labels alone.
    for package in ("thermo", "chemicals"):
        root = Path(import_module(package).__file__).parent
        files = {str(p.relative_to(root)): file_hash(p) for p in sorted(root.rglob("*"))
                 if p.is_file() and "__pycache__" not in p.parts and p.suffix not in {".pyc", ".pyo"}}
        data_hashes[package] = digest(files)
    return {"config": config, "packages": packages, "code_sha256": code, "library_tree_sha256": data_hashes,
            "python": platform.python_version(), "catalog_only": catalog_only, "limit_pairs": limit_pairs}


def failure_rows(pair: dict, config: dict, ws: list[float], status: str, error: str, stage: str) -> list[dict]:
    return [{"pair_id": pair["pair_id"], "cas_number": pair["cas_number"], "name": pair["name"], "smiles": pair["smiles"],
             "additive_mass_fraction": w, "outlet_temperature_k": t, "pressure_pa": p, "model": model,
             "status": status, "error": error, "model_comparison_eligible": False, "stage": stage,
             "evidence_kind": "no_prediction", "calibrated_uncertainty": None}
            for model in MODELS for w in ws for t in config["grid"]["outlet_temperatures_k"] for p in config["grid"]["pressures_pa"]]


def pair_calculation(task: dict, factory=None) -> dict:
    if factory is None:
        from v5m2_backend import BinaryModel
        factory = BinaryModel
    pair, config, ws, stage = (task[k] for k in ("pair", "config", "mass_fractions", "stage"))
    rows, meta = [], {}
    for model in MODELS:
        if model not in pair["models"]:
            rows += [r for r in failure_rows(pair, config, ws, "model_unavailable",
                      pair.get("model_rejections", {}).get(model, "no complete parameters"), stage) if r["model"] == model]
            continue
        try:
            backend = factory(pair, model, config)
        except (ValueError, TypeError, RuntimeError, AttributeError, KeyError, ArithmeticError) as exc:
            rows += [r for r in failure_rows(pair, config, ws, "model_initialization_failed", f"{type(exc).__name__}:{exc}", stage) if r["model"] == model]
            continue
        meta[model] = {"pure_methods": backend.methods, "viscosity_methods": backend.viscosity_metadata,
                       "parameter_record": pair["models"][model]}
        for i, w in enumerate(ws, 1):
            print(f"[{stage}] {pair['cas_number']} {model} composition={i}/{len(ws)} w={w:g}", flush=True)
            for t in config["grid"]["outlet_temperatures_k"]:
                for p in config["grid"]["pressures_pa"]:
                    row = backend.evaluate(w, t, p)
                    rows.append({**row, "stage": stage})
    return {"rows": rows, "model_metadata": meta}


def run_worker(task: dict, folder: Path, timeout_s: float) -> dict:
    folder.mkdir(parents=True, exist_ok=True)
    inp, out = folder / "request.json", folder / "response.json"
    write_json(inp, task)
    if out.exists(): out.unlink()
    env = {**os.environ, "OMP_NUM_THREADS": "1", "OPENBLAS_NUM_THREADS": "1", "MKL_NUM_THREADS": "1"}
    status = "worker_failed"
    try:
        with (folder / "worker.log").open("w", encoding="utf-8") as log:
            result = subprocess.run([sys.executable, str(HERE / "run_predictive_grid.py"), "--worker-input", str(inp),
                                    "--worker-output", str(out)], stdout=log, stderr=subprocess.STDOUT,
                                   timeout=timeout_s, check=False, env=env)
        if result.returncode != 0 or not out.is_file():
            raise ValueError(f"worker exit={result.returncode}; see worker.log")
        data = load_json(out)
        required = {(r["model"], r["additive_mass_fraction"], r["outlet_temperature_k"], r["pressure_pa"])
                    for r in failure_rows(task["pair"], task["config"], task["mass_fractions"], "", "", task["stage"])}
        found = [(r["model"], r["additive_mass_fraction"], r["outlet_temperature_k"], r["pressure_pa"]) for r in data["rows"]]
        if len(found) != len(set(found)) or set(found) != required or any(
            r.get("pair_id") != task["pair"]["pair_id"] or r.get("stage") != task["stage"] for r in data["rows"]):
            raise ValueError("worker returned mismatched or duplicate states")
        return data
    except subprocess.TimeoutExpired:
        status, error = "worker_timeout", "per_pair_wall_clock_timeout"
    except (OSError, ValueError, KeyError, TypeError) as exc:
        error = str(exc)
    return {"rows": failure_rows(task["pair"], task["config"], task["mass_fractions"], status, error, task["stage"]),
            "model_metadata": {}, "error": error}


def execute_stage(tasks: list[dict], output: Path, signature: str, config: dict, *, retry_failures=False, worker=run_worker) -> list[dict]:
    results, pending = {}, []
    for task in tasks:
        key = digest(task)
        path = output / "checkpoints" / (key + ".json")
        saved = checkpoint(path, signature)
        if saved is not None and retry_failures and any(r["status"] in {"worker_failed", "worker_timeout", "model_initialization_failed"} for r in saved["rows"]):
            saved = None
        if saved is not None:
            results[key] = saved
            print(f"[{task['stage']} cached] {task['pair']['name']} ({task['pair']['cas_number']})", flush=True)
        else:
            pending.append((key, task, path))
    with ThreadPoolExecutor(max_workers=config["workers"]) as pool:
        futures = {pool.submit(worker, task, output / "workers" / key[:24], config["pair_timeout_s"]): (key, task, path)
                   for key, task, path in pending}
        for future in as_completed(futures):
            key, task, path = futures[future]
            data = future.result()
            checkpoint(path, signature, data); results[key] = data
            counts = dict(Counter(r["status"] for r in data["rows"]))
            print(f"[{task['stage']} {len(results)}/{len(tasks)}] {task['pair']['name']} ({task['pair']['cas_number']}) {counts}", flush=True)
            write_json(output / "progress.json", {"signature": signature, "stage": task["stage"], "completed": len(results), "total": len(tasks)})
    return [results[digest(task)] for task in tasks]


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, default=HERE / "config-v5m2.yaml")
    parser.add_argument("--output-dir", type=Path, default=Path("predictive-v5m2-results"))
    parser.add_argument("--resume", action="store_true")
    parser.add_argument("--retry-failures", action="store_true")
    parser.add_argument("--catalog-only", action="store_true")
    parser.add_argument("--limit-pairs", type=int, help="Deterministic smoke-test subset; not a complete catalog run.")
    parser.add_argument("--workers", type=int, help="Override parallel pair workers; part of the run signature.")
    parser.add_argument("--worker-input", type=Path, help=argparse.SUPPRESS)
    parser.add_argument("--worker-output", type=Path, help=argparse.SUPPRESS)
    args = parser.parse_args()
    if args.worker_input:
        if not args.worker_output: parser.error("missing worker output")
        write_json(args.worker_output, pair_calculation(load_json(args.worker_input)))
        return 0
    output = args.output_dir.resolve()
    try:
        config = yaml.safe_load(args.config.read_text(encoding="utf-8"))
        if args.workers is not None: config["workers"] = args.workers
        validate_config(config)
        if args.limit_pairs is not None and args.limit_pairs <= 0: raise ValueError("limit-pairs must be positive")
        if args.retry_failures and not args.resume: raise ValueError("retry-failures requires resume")
        provenance = contract(config, catalog_only=args.catalog_only, limit_pairs=args.limit_pairs)
    except (OSError, ValueError, TypeError, KeyError, ImportError) as exc:
        parser.error(str(exc))
    signature = digest(provenance)
    output.mkdir(parents=True, exist_ok=True)
    import fcntl
    lock = (output / ".run.lock").open("a")
    try:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
    except BlockingIOError:
        parser.error("another process owns this output directory")
    try:
        manifest_path = output / "manifest.json"
        if args.resume:
            if not manifest_path.is_file() or load_json(manifest_path).get("signature") != signature:
                parser.error("resume requires unchanged code, inputs, configuration and libraries")
        elif any(p.name != ".run.lock" for p in output.iterdir()):
            parser.error("output directory must be empty; use --resume or a new directory")
        try:
            git_sha = subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=HERE, text=True, stderr=subprocess.DEVNULL).strip()
        except (OSError, subprocess.CalledProcessError): git_sha = None
        write_json(manifest_path, {"version": VERSION, "signature": signature, "git_commit": git_sha,
                                  "contract": provenance, "network_used": False, "ml_model_trained": False})
        from v5m2_backend import discover_catalog
        catalog_path = output / "checkpoints" / "catalog.json"
        catalog = checkpoint(catalog_path, signature)
        if catalog is None:
            print("[catalog] enumerating local water pairs and checking directed parameters", flush=True)
            catalog = checkpoint(catalog_path, signature, discover_catalog(config))
        write_json(output / "catalog.json", catalog)
        write_csv(output / "catalog.csv", catalog["candidates"])
        write_csv(output / "catalog-rejections.csv", catalog["rejections"])
        pairs = catalog["candidates"][:args.limit_pairs] if args.limit_pairs else catalog["candidates"]
        if not pairs:
            parser.error("no model-supported pairs in catalog; inspect catalog-rejections.csv")
        ws = mass_grid(config)
        condition_count = len(config["grid"]["outlet_temperatures_k"])*len(config["grid"]["pressures_pa"])
        print(f"[design] pairs={len(pairs)} compositions={len(ws)} models=2 conditions={condition_count} coarse_states={len(pairs)*len(ws)*2*condition_count}", flush=True)
        plan = {"pairs": len(pairs), "mass_fractions": ws, "conditions_per_composition_model": condition_count,
                "models": list(MODELS), "composition_basis": "mass", "catalog_only": args.catalog_only}
        write_json(output / "grid-plan.json", plan)
        if args.catalog_only:
            summary = {"study_version": VERSION, "study_complete": False, "run_status": "catalog_only", **plan}
        else:
            tasks = [{"pair": pair, "config": config, "mass_fractions": ws, "stage": "coarse"} for pair in pairs]
            results = execute_stage(tasks, output, signature, config, retry_failures=args.retry_failures)
            rows = [row for res in results for row in res["rows"]]
            coarse = compare_models(rows, condition_count, config)
            write_csv(output / "coarse-tradeoffs.csv", coarse)
            by_pair = {p["pair_id"]: p for p in pairs}
            possible = sorted([r for r in coarse if r["comparison_eligible"]], key=lambda r: (
                -r["minimum_paired_delta_h_ratio"], -r["maximum_model_spread_delta_h_ratio"], r["pair_id"], r["additive_mass_fraction"]))
            chosen = []
            # Reserve half of the pair budget for model-disagreement exploration.
            streams = [possible, sorted(possible, key=lambda r: (-r["maximum_model_spread_delta_h_ratio"], r["pair_id"]))]
            for stream in streams:
                quota = max(1, config["refinement_pairs"]//2)
                count = 0
                for row in stream:
                    if row["pair_id"] not in chosen and len(chosen) < config["refinement_pairs"]:
                        chosen.append(row["pair_id"]); count += 1
                        if count >= quota: break
            refined_tasks = []
            for pair_id in chosen:
                new_ws = refinement_grid(ws, [r for r in coarse if r["pair_id"] == pair_id], config)
                if new_ws: refined_tasks.append({"pair": by_pair[pair_id], "config": config, "mass_fractions": new_ws, "stage": "refined"})
            write_json(output / "refinement-plan.json", {"tasks": refined_tasks, "grid_convergence_claimed": False})
            refined_results = execute_stage(refined_tasks, output, signature, config, retry_failures=args.retry_failures)
            rows += [r for result in refined_results for r in result["rows"]]
            tradeoffs = compare_models(rows, condition_count, config)
            write_csv(output / "prediction-grid.csv", rows)
            write_json(output / "predictions.json", {"version": VERSION, "rows": rows,
                        "model_metadata": [r["model_metadata"] for r in results+refined_results]})
            write_csv(output / "tradeoffs.csv", tradeoffs)
            write_csv(output / "pareto-thermodynamic.csv", pareto(tradeoffs, THERMO_OBJECTIVES))
            write_csv(output / "pareto-with-flow-proxy.csv", pareto(tradeoffs, FLOW_OBJECTIVES))
            write_csv(output / "validation-shortlist.csv", shortlist(tradeoffs, config))
            changes = []
            for pair_id in chosen:
                a = [r for r in coarse if r["pair_id"] == pair_id and r["comparison_eligible"]]
                b = [r for r in tradeoffs if r["pair_id"] == pair_id and r["comparison_eligible"]]
                if a and b:
                    old = max(r["minimum_paired_delta_h_ratio"] for r in a)
                    new = max(r["minimum_paired_delta_h_ratio"] for r in b)
                    changes.append({"pair_id": pair_id, "coarse_best_delta_h_ratio": old,
                                    "refined_best_delta_h_ratio": new, "change": new-old,
                                    "grid_convergence_claimed": False})
            write_csv(output / "refinement-diagnostics.csv", changes)
            failure_counts = Counter(r["status"] for r in rows)
            eligible = [r for r in tradeoffs if r["comparison_eligible"]]
            summary = {"study_version": VERSION, "study_complete": args.limit_pairs is None,
                "run_status": "smoke_test" if args.limit_pairs else "no_complete_two_model_comparisons" if not eligible else "completed_with_recorded_statuses",
                "catalog_pair_count": len(catalog["candidates"]), "screened_pair_count": len(pairs),
                "catalog_budget_limited": catalog["budget_limited"],
                "coarse_compositions_per_pair": len(ws), "refined_pair_count": len(refined_tasks),
                "predicted_state_count": len(rows), "state_status_counts": dict(failure_counts),
                "complete_two_model_formulation_count": len(eligible),
                "paired_model_gain_hypothesis_count": sum(r["minimum_paired_delta_h_ratio"] > 1 for r in eligible),
                "pareto_thermodynamic_count": len(pareto(tradeoffs, THERMO_OBJECTIVES)),
                "pareto_with_flow_proxy_count": len(pareto(tradeoffs, FLOW_OBJECTIVES)),
                "validation_shortlist_count": len(shortlist(tradeoffs, config)),
                "ml_model_trained": False, "empirical_uncertainty_calibrated": False,
                "experimental_winners": 0, "surface_enhancement_predicted": False, "rankable_promotions": 0,
                "reference_dataset_used": None, "grid_convergence_claimed": False,
                "interpretation": "Parameter-backed NRTL/Dortmund grid hypotheses, not validated coolant gains. Missing models and failed states are not imputed; model spread is not calibrated uncertainty."}
        write_json(output / "summary.json", summary)
        write_json(output / "progress.json", {"signature": signature, "stage": summary["run_status"]})
        write_json(output / "artifact-hashes.json", {p.name: file_hash(p) for p in sorted(output.iterdir())
                   if p.is_file() and p.name not in {".run.lock", "artifact-hashes.json"}})
        print(json.dumps(summary, indent=2, sort_keys=True))
        return 1 if any(summary.get("state_status_counts", {}).get(k, 0) for k in ("worker_failed", "worker_timeout", "model_initialization_failed")) else 0
    finally:
        lock.close()


if __name__ == "__main__":
    raise SystemExit(main())
