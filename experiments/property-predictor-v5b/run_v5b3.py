#!/usr/bin/env python3
"""V5b-3: offline reference audit and isolated, molecule-grouped property benchmark."""
from __future__ import annotations

import argparse
from collections import Counter
import csv
from importlib import import_module
from importlib.metadata import PackageNotFoundError, version
import json
import math
import os
from pathlib import Path
import platform
import subprocess
import sys
import time
from typing import Any

import yaml
from v5b3_audit import LocalReferenceBackend, audit_target, digest
from v5b3_metrics import applicability, comparison_rows, summarize

HERE = Path(__file__).resolve().parent
ENTRY = HERE.parent / "entry-evaluator"


def load_json(path: Path):
    def reject_constant(text):
        raise ValueError("nonfinite_json_constant: " + text)
    return json.loads(path.read_text(encoding="utf-8"), parse_constant=reject_constant)


def write_json(path: Path, value: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temp = path.with_suffix(path.suffix + ".tmp")
    temp.write_text(json.dumps(value, indent=2, sort_keys=True, allow_nan=False) + "\n", encoding="utf-8")
    temp.replace(path)


def file_hash(path: Path) -> str:
    import hashlib
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    fields = sorted({key for row in rows for key in row}) or ["status"]
    temp = path.with_suffix(".tmp")
    with temp.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader()
        for row in rows:
            writer.writerow({k: json.dumps(v, sort_keys=True, allow_nan=False) if isinstance(v, (dict, list)) else v
                             for k, v in row.items()})
    temp.replace(path)


def validate_config(config: dict[str, Any]) -> None:
    if not isinstance(config, dict) or config.get("version") != "v5b-3":
        raise ValueError("config.version must be v5b-3")
    for key in ("identity_mw_relative_tolerance", "reproduction_relative_tolerance", "constant_conflict_relative_tolerance"):
        value = config["audit"][key]
        if isinstance(value, bool) or not math.isfinite(float(value)) or not 0 < float(value) < 1:
            raise ValueError("invalid audit tolerance: " + key)
    if not 0 < float(config["statistics"]["evaluation_fraction"]) < 1:
        raise ValueError("evaluation_fraction must be in (0, 1)")
    if not 0 < float(config["statistics"]["coverage_target"]) < 1:
        raise ValueError("coverage_target must be in (0, 1)")
    if not isinstance(config["statistics"]["split_seed"], str) or not config["statistics"]["split_seed"]:
        raise ValueError("split_seed must be a nonempty string")
    for key in ("minimum_calibration_groups", "minimum_evaluation_groups", "minimum_neighbors"):
        value = config["statistics"][key]
        if isinstance(value, bool) or not isinstance(value, int) or value < 1:
            raise ValueError("invalid group/neighbor threshold: " + key)
    edge, inner = (config["statistics"][x] for x in ("neighbor_similarity", "in_domain_similarity"))
    if not 0 <= float(edge) <= float(inner) <= 1:
        raise ValueError("invalid structural similarity thresholds")
    timeout = config["worker_timeout_s"]
    if isinstance(timeout, bool) or not math.isfinite(float(timeout)) or float(timeout) <= 0:
        raise ValueError("worker_timeout_s must be finite and positive")


def load_targets(payload: dict[str, Any]) -> list[dict[str, Any]]:
    if not isinstance(payload, dict) or payload.get("version") != "v5e-2.1":
        raise ValueError("expected V5e-2.1 empirical-reference-properties.json")
    targets = payload.get("targets")
    if not isinstance(targets, list) or not targets:
        raise ValueError("no targets in input")
    ids = set()
    for row in targets:
        if not isinstance(row, dict) or not isinstance(row.get("candidate_id"), str) or not row["candidate_id"]:
            raise ValueError("invalid candidate_id")
        if row["candidate_id"] in ids:
            raise ValueError("duplicate candidate_id")
        ids.add(row["candidate_id"])
        packet = row.get("packet")
        if packet is not None:
            if not isinstance(packet, dict):
                raise ValueError("packet must be a mapping or null")
            for key in ("constants", "series"):
                if not isinstance(packet.get(key), dict):
                    raise ValueError("invalid packet section: " + key)
                for record in packet[key].values():
                    if record is not None and not isinstance(record, dict):
                        raise ValueError("invalid property record")
    return sorted(targets, key=lambda r: r["candidate_id"])


def prediction_request(audited: dict[str, Any]) -> dict[str, Any]:
    fields = ("observation_id", "property", "temperature_k", "pressure_pa", "state_basis")
    return {"smiles": audited["smiles"],
            "requests": [{k: row[k] for k in fields} for row in audited["observations"]]}


def checkpoint(path: Path, signature: str, data=None):
    if data is not None:
        write_json(path, {"signature": signature, "data_sha256": digest(data), "data": data})
        return data
    if not path.exists():
        return None
    saved = load_json(path)
    if saved.get("signature") != signature or saved.get("data_sha256") != digest(saved["data"]):
        raise ValueError("checkpoint mismatch or corruption: " + str(path))
    return saved["data"]


def predict_subprocess(request: dict[str, Any], folder: Path, timeout: float) -> dict[str, Any]:
    folder.mkdir(parents=True, exist_ok=True)
    inp, out = folder / "request.json", folder / "worker-output.json"
    write_json(inp, request)
    if out.exists():
        out.unlink()  # Never accept stale worker output after a failed retry.
    error, status = None, None
    env = {**os.environ, "OMP_NUM_THREADS": "1", "OPENBLAS_NUM_THREADS": "1", "RAYON_NUM_THREADS": "1"}
    try:
        with (folder / "worker.log").open("w", encoding="utf-8") as log:
            result = subprocess.run([sys.executable, str(HERE / "v5b3_predict.py"), "--input", str(inp),
                                     "--output", str(out)], stdout=log, stderr=subprocess.STDOUT,
                                    timeout=timeout, check=False, env=env)
        if result.returncode != 0:
            error, status = f"worker_exit_code:{result.returncode}", "worker_failed"
        elif not out.is_file():
            error, status = "worker_did_not_write_output", "worker_failed"
        else:
            value = load_json(out)
            wanted = {r["observation_id"] for r in request["requests"]}
            found = [r["observation_id"] for r in value["results"]]
            if len(found) != len(set(found)) or set(found) != wanted:
                raise ValueError("worker_result_observation_mismatch")
            return value
    except subprocess.TimeoutExpired:
        error, status = "per_molecule_wall_clock_timeout", "worker_timeout"
    except (OSError, ValueError, KeyError, TypeError) as exc:
        error, status = str(exc), "worker_failed"
    return {"results": [{"observation_id": r["observation_id"], "predicted_value": None,
                         "status": status, "error": error} for r in request["requests"]],
            "reference_values_supplied": False, "parameters_fitted": False}


def execute(targets, config, output, signature, backend, *, audit_only=False,
            retry_failures=False, predictor=predict_subprocess):
    audited, predictions, all_observations, seen_obs = [], [], [], set()
    for i, target in enumerate(targets, 1):
        start = time.monotonic()
        name = digest(target["candidate_id"])[:24]
        path = output / "checkpoints" / (name + "-audit.json")
        row = checkpoint(path, signature)
        cached = row is not None
        if row is None:
            row = checkpoint(path, signature, audit_target(target, config["audit"], backend))
        # A repeated alias cannot increase benchmark weight.
        unique = []
        for obs in row["observations"]:
            if obs["observation_id"] not in seen_obs:
                seen_obs.add(obs["observation_id"])
                unique.append(obs)
        row = {**row, "observations": unique}
        audited.append(row)
        all_observations.extend(unique)
        strict = sum(o["benchmark_eligible"] for o in unique)
        print(f"[audit {i}/{len(targets)}] {target['candidate_id']} strict={strict} diagnostic={len(unique)-strict} "
              f"identity={row.get('identity_verified', False)} cached={cached} elapsed={time.monotonic()-start:.2f}s", flush=True)
    frozen = {"version": "v5b-3", "targets": audited, "observations": all_observations,
              "reference_values_frozen_before_prediction": True}
    write_json(output / "audited-reference.json", frozen)
    write_json(output / "prediction-requests.json", [prediction_request(r) for r in audited if r["observations"]])
    audit_rows = []
    for r in audited:
        entries = r["audit"] or [{"rejected": [{"reason": r.get("identity_reason", "no_property_packet")}]}]
        audit_rows.extend({"candidate_id": r["candidate_id"], "identity_verified": r.get("identity_verified"),
                           "identity_reason": r.get("identity_reason"), **entry} for entry in entries)
    write_csv(output / "reference-audit.csv", audit_rows)
    write_json(output / "progress.json", {"phase": "audit_complete", "target_count": len(audited),
                                         "strict_observation_count": sum(o["benchmark_eligible"] for o in all_observations)})
    if not audit_only:
        planned = [r for r in audited if r["observations"]]
        for i, row in enumerate(planned, 1):
            start = time.monotonic()
            request = prediction_request(row)
            name = digest(request)
            path = output / "checkpoints" / (name + "-prediction.json")
            pred = checkpoint(path, signature)
            if pred is not None and retry_failures and any(r["status"] != "ok" for r in pred["results"]):
                pred = None
            cached = pred is not None
            print(f"[predict {i}/{len(planned)}] {row['candidate_id']} points={len(request['requests'])} cached={cached}", flush=True)
            if pred is None:
                pred = predictor(request, output / "workers" / name[:24], config["worker_timeout_s"])
                checkpoint(path, signature, pred)
            predictions.extend(pred["results"])
            failed = sum(p["status"] != "ok" for p in pred["results"])
            print(f"  completed failures={failed} elapsed={time.monotonic()-start:.2f}s", flush=True)
            write_json(output / "progress.json", {"phase": "prediction", "completed": i, "total": len(planned)})
        write_json(output / "predictions-before-comparison.json", {"results": predictions, "reference_values_supplied": False})
    compare = comparison_rows(all_observations, predictions, config["statistics"]) if not audit_only else []
    report = summarize(compare, config["statistics"])
    write_csv(output / "property-comparison.csv", compare)
    write_csv(output / "molecule-property-errors.csv", report.pop("molecule_metrics"))
    write_json(output / "property-error-calibration.json", report)
    write_json(output / "property-applicability.json", applicability(compare, config["statistics"]) if compare else [])
    successful = [r for r in compare if r["benchmark_eligible"] and r["status"] == "ok"]
    summary = {"study_version": "v5b-3", "study_complete": not audit_only,
               "benchmark_has_comparable_results": bool(successful),
               "run_status": "audit_only" if audit_only else "no_comparable_references" if not any(o["benchmark_eligible"] for o in all_observations) else "completed_with_recorded_failures" if any(r["status"] != "ok" for r in compare) else "completed",
               "target_count": len(targets), "identity_verified_count": sum(r.get("identity_verified", False) for r in audited),
               "targets_with_strict_reference_count": sum(any(o["benchmark_eligible"] for o in r["observations"]) for r in audited),
               "strict_reference_point_count": sum(o["benchmark_eligible"] for o in all_observations),
               "diagnostic_reference_point_count": sum(not o["benchmark_eligible"] for o in all_observations),
               "strict_successful_prediction_count": len(successful),
               "strict_failed_prediction_count": sum(r["benchmark_eligible"] and r["status"] != "ok" for r in compare),
               "strict_property_coverage": dict(Counter(o["property"] for o in all_observations if o["benchmark_eligible"])),
               "prediction_status_counts": dict(Counter(r["status"] for r in compare)),
               "property_results": report["properties"], "rankable_promotions": 0,
               "predictor_parameters_changed": False, "entry_uncertainty_updated": False,
               "interpretation": "Audit and fixed-predictor benchmark only. Source labels are not independent measurements; diagnostic state mismatches do not calibrate errors."}
    write_json(output / "summary.json", summary)
    write_json(output / "progress.json", {"phase": summary["run_status"]})
    return summary


def runtime_contract(input_path, config, audit_only):
    pins = {"chemicals": "1.5.2", "thermo": "0.6.1", "rdkit": "2026.3.6"}
    if not audit_only:
        pins["feos"] = "0.10.1"
    installed = {}
    for name in ("chemicals", "thermo", "rdkit", "feos", "CoolProp", "fluids", "numpy", "scipy", "PyYAML", "si-units"):
        try:
            installed[name] = version(name)
        except PackageNotFoundError:
            installed[name] = None
    for name, wanted in pins.items():
        if installed[name] != wanted:
            raise ValueError(f"requires {name}=={wanted}; installed={installed[name]}; use the project's V5b environment")
    code = {p.name: file_hash(p) for p in [HERE / "run_v5b3.py", HERE / "v5b3_audit.py",
                                         HERE / "v5b3_predict.py", HERE / "v5b3_metrics.py"]}
    provider = ENTRY / "feos_gc_provider.py"
    if provider.exists():
        code["feos_gc_provider.py"] = file_hash(provider)
    provider_contract = ENTRY / "property_provider.py"
    if provider_contract.exists():
        code["property_provider.py"] = file_hash(provider_contract)
    if not audit_only:
        required = [provider, provider_contract] + [ENTRY / "parameters" / "v5b" / n for n in
                    ("sauer2014_smarts.json", "rehner2023_hetero.json", "joback1987.json")]
        missing = [str(p) for p in required if not p.is_file()]
        if missing:
            raise ValueError("missing predictor files: " + ", ".join(missing))
    parameters = {p.name: file_hash(p) for p in sorted((ENTRY / "parameters" / "v5b").glob("*.json"))}
    # Hash packaged source and data, not just package version strings.
    libraries = {}
    for name in ("chemicals", "thermo"):
        root = Path(import_module(name).__file__).parent
        libraries[name] = {str(p.relative_to(root)): file_hash(p) for p in sorted(root.rglob("*"))
                           if p.is_file() and p.suffix.lower() in {".py", ".json", ".tsv", ".csv", ".db", ".sqlite"}}
    return {"input_sha256": file_hash(input_path), "config": config, "code": code,
            "parameter_files": parameters, "packages": installed,
            "reference_library_files": libraries, "python": platform.python_version()}


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, default=HERE / "benchmark-v5b3.yaml")
    parser.add_argument("--input", type=Path)
    parser.add_argument("--output-dir", type=Path, default=HERE / "property-v5b3-results")
    parser.add_argument("--audit-only", action="store_true")
    parser.add_argument("--resume", action="store_true")
    parser.add_argument("--retry-failures", action="store_true")
    args = parser.parse_args()
    try:
        config = yaml.safe_load(args.config.read_text(encoding="utf-8"))
        validate_config(config)
        source = (args.input or (args.config.parent / config["input"])).resolve()
        targets = load_targets(load_json(source))
        contract = runtime_contract(source, config, args.audit_only)
        output = args.output_dir.resolve()
        if output.exists() and any(output.iterdir()) and not args.resume:
            raise ValueError("output directory is not empty; use --resume or a new directory")
        if args.retry_failures and not args.resume:
            raise ValueError("--retry-failures requires --resume")
        output.mkdir(parents=True, exist_ok=True)
        # Linux advisory lock is released on crash; it has no stale PID issue.
        import fcntl
        with (output / ".run.lock").open("w") as lock:
            fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
            signature = digest(contract)
            manifest = output / "manifest.json"
            if manifest.exists():
                prior = load_json(manifest)
                if prior.get("contract_sha256") != signature:
                    raise ValueError("resume rejected: input, configuration, code or runtime changed; use a new directory")
            else:
                write_json(manifest, {"contract_sha256": signature, "contract": contract,
                                      "input_path": str(source), "network_requests": False,
                                      "reference_state_controls_may_come_from_measurement_conditions": True})
            if args.audit_only and (output / "summary.json").exists():
                if load_json(output / "summary.json").get("study_complete"):
                    raise ValueError("audit-only cannot overwrite a completed benchmark; use a new directory")
            summary = execute(targets, config, output, signature, LocalReferenceBackend(),
                              audit_only=args.audit_only, retry_failures=args.retry_failures)
            write_json(output / "artifact-hashes.json", {str(p.relative_to(output)): file_hash(p)
                for p in sorted(output.rglob("*")) if p.is_file() and p.name not in {".run.lock", "artifact-hashes.json"}})
        print(json.dumps(summary, indent=2, sort_keys=True))
        return 0
    except KeyboardInterrupt:
        print("Interrupted. Completed checkpoints were retained; rerun with --resume.", file=sys.stderr)
        return 130
    except (OSError, ValueError, KeyError, TypeError, ImportError) as exc:
        parser.error(str(exc))
    return 2


if __name__ == "__main__":
    raise SystemExit(main())
