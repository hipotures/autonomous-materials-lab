#!/usr/bin/env python3
"""V5m-1: offline formulation screening and gated physical-test planning."""
from __future__ import annotations

import argparse
from collections import Counter
import csv
import hashlib
import json
from importlib.metadata import PackageNotFoundError, version
import math
from pathlib import Path
import sys
import time

import yaml
from formulations import (MEASUREMENT_FIELDS, build_formulations, cold_ratios, digest,
    heat_only_comparison, make_bench_plan, protocol_template, read_json, PROFILES,
    retained_additive_stress, write_json)
from mixture_models import MixtureModels, ENTRY

HERE = Path(__file__).resolve().parent


def write_csv(path, rows, fields=None):
    fields = fields or sorted({k for row in rows for k in row}) or ["status"]
    tmp = path.with_suffix(path.suffix + ".tmp")
    with tmp.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fields)
        writer.writeheader()
        for row in rows:
            writer.writerow({k: json.dumps(v, sort_keys=True, allow_nan=False) if isinstance(v, (dict, list))
                             else v for k, v in row.items()})
    tmp.replace(path)



def fingerprint(config):
    packages = {}
    for name in ("CoolProp", "thermo", "chemicals", "fluids", "numpy", "scipy", "PyYAML"):
        try:
            packages[name] = version(name)
        except PackageNotFoundError:
            packages[name] = None
    files = list(HERE.glob("*.py")) + [ENTRY / "property_provider.py", ENTRY / "thermo_provider.py"]
    hashes = {str(p.relative_to(HERE.parent)): hashlib.sha256(p.read_bytes()).hexdigest()
              for p in files if p.is_file()}
    return {"config": config, "packages": packages, "code_sha256": hashes,
            "python": sys.version, "network_required": False,
            "model_sources": {
                "cold": "https://coolprop.org/fluid_properties/Incompressibles.html",
                "ethanol_water_hot": "https://thermo.readthedocs.io/thermo.nrtl.html"},
            "model_accuracy_established_by_this_stage": False}


def checkpoint(path, signature, data=None):
    if data is not None:
        write_json(path, {"signature": signature, "data_hash": digest(data), "data": data})
        return data
    if not path.exists():
        return None
    saved = read_json(path)
    if saved.get("signature") != signature or digest(saved["data"]) != saved.get("data_hash"):
        raise ValueError(f"checkpoint mismatch or corruption: {path}")
    return saved["data"]


def execute(config, output, signature, models=None):
    formulations = build_formulations(config)
    plan_only = models is None
    cold, hot, sensible = [], [], []
    for i, f in enumerate(formulations, 1):
        started = time.monotonic()
        path = output / "checkpoints" / (f["formulation_id"] + ".json")
        saved = checkpoint(path, signature)
        cached = saved is not None
        if saved is None:
            states = [] if plan_only else [models.cold(f, t, config["pressure_pa"])
                                          for t in config["cold_temperatures_k"]]
            windows = [] if plan_only else [models.hot(f, config["cold_temperatures_k"][0], t,
                                        config["pressure_pa"]) for t in config["hot_outlet_temperatures_k"]]
            saved = checkpoint(path, signature, {"cold": states, "hot": windows})
        cold.extend(saved["cold"])
        hot.extend(saved["hot"])
        states = saved["cold"]
        dh = None
        if states and states[0]["status"] == states[-1]["status"] == "ok":
            # Same composition, backend and pressure at both endpoints only.
            delta = states[-1]["enthalpy_j_kg"] - states[0]["enthalpy_j_kg"]
            if math.isfinite(delta) and delta > 0:
                dh = delta
        sensible.append({"formulation_id": f["formulation_id"], "delta_h_j_kg": dh,
            "inlet_temperature_k": config["cold_temperatures_k"][0],
            "outlet_temperature_k": config["cold_temperatures_k"][-1],
            "basis": "cold_sensible_only_same_composition_backend",
            "is_tps_cooling_capacity": False})
        print(f"[screen {i}/{len(formulations)}] {f['additive_id']} w={f['additive_mass_fraction']:g} "
              f"cold_ok={sum(r['status']=='ok' for r in saved['cold'])}/{len(saved['cold'])} "
              f"hot_ok={sum(r['status']=='ok' for r in saved['hot'])}/{len(saved['hot'])} "
              f"cached={cached} elapsed={time.monotonic()-started:.2f}s", flush=True)
        write_json(output / "progress.json", {"signature": signature, "completed": i,
                                              "total": len(formulations), "phase": "screen"})
    water = {r["temperature_k"]: r for r in cold if r["formulation_id"] == "water"}
    ratios = [cold_ratios(r, water[r["temperature_k"]]) for r in cold]
    endpoint_rows, endpoint_points = [], {}
    if not plan_only:
        for additive in sorted({f["additive_id"] for f in formulations} & set(PROFILES)):
            key = output / "checkpoints" / (additive + "-zero-endpoint.json")
            points = checkpoint(key, signature)
            if points is None:
                probe = {"formulation_id": "endpoint-" + additive, "additive_id": additive,
                         "additive_mass_fraction": 0.0}
                points = [models.cold(probe, t, config["pressure_pa"], fit_endpoint_probe=True)
                          for t in config["cold_temperatures_k"]]
                checkpoint(key, signature, points)
            for point in points:
                endpoint_points[(additive, point["temperature_k"])] = point
                endpoint_rows.append({**cold_ratios(point, water[point["temperature_k"]]),
                    "backend": point["backend"], "status": point["status"],
                    "is_pure_water_fit_endpoint_not_an_extra_formulation": True})
    by_id = {f["formulation_id"]: f for f in formulations}
    for comparison, point in zip(ratios, cold):
        additive = by_id[point["formulation_id"]]["additive_id"]
        zero = endpoint_points.get((additive, point["temperature_k"]), {})
        for field in ("density_kg_m3", "cp_j_kg_k", "viscosity_pa_s", "conductivity_w_m_k"):
            x, y = point.get(field), zero.get(field)
            comparison[field + "_ratio_vs_own_fit_water_endpoint"] = x / y if x and y else None
        comparison["within_fit_trend_is_not_experimental_calibration"] = True
    water_h = {r["outlet_temperature_k"]: r for r in hot if r["formulation_id"] == "water"}
    for r in hot:
        ref = water_h[r["outlet_temperature_k"]]
        if r["status"] == ref["status"] == "ok":
            r.update(heat_only_comparison(r["delta_h_j_kg"], ref["delta_h_j_kg"]))
    for r in sensible:
        ref = sensible[0]["delta_h_j_kg"]
        r["mass_ratio_equal_sensible_heat"] = ref / r["delta_h_j_kg"] if ref and r["delta_h_j_kg"] else None
    plan = make_bench_plan(formulations, config)
    write_json(output / "bench-plan.json", {"version": "v5m-1", "signature": signature,
        "physical_execution_authorized": False, "runs": plan, "formulations": formulations})
    write_csv(output / "formulations.csv", formulations)
    write_csv(output / "cold-property-grid.csv", cold)
    write_csv(output / "cold-comparison.csv", ratios)
    write_csv(output / "model-endpoint-comparison.csv", endpoint_rows)
    write_csv(output / "cold-sensible-window.csv", sensible)
    write_csv(output / "enthalpy-window.csv", hot)
    write_csv(output / "concentration-stress.csv", [retained_additive_stress(f, loss)
        for f in formulations if f["additive_id"] != "water" for loss in config["water_loss_fractions"]])
    write_csv(output / "bench-plan.csv", plan)
    # Templates are not measurements and are never overwritten on resume.
    template = output / "measurements-template.csv"
    if not template.exists():
        write_csv(template, [{k: r[k] for k in ("run_id", "formulation_id", "block_id", "role", "replicate")}
                            for r in plan], MEASUREMENT_FIELDS)
    protocol = output / "protocol-template.yaml"
    if not protocol.exists():
        protocol.write_text(yaml.safe_dump(protocol_template(), sort_keys=False), encoding="utf-8")
    needs = [{"formulation_id": f["formulation_id"], "missing_system_evidence":
              "stability;wetting;surface_tension;porous_flow;residue;thermal_test",
              "suspension_protocol_required": False, "physical_execution_authorized": False}
             for f in formulations]
    write_csv(output / "measurement-needs.csv", needs)
    unexpected = sum(r["status"] in {"property_failed", "state_failed", "not_liquid"} for r in cold)
    unexpected += sum(r["status"] == "model_failed" for r in hot)
    unexpected += sum(r["status"] != "ok" for r in endpoint_rows)
    summary = {"study_version": "v5m-1", "study_complete": True,
        "run_status": "plan_only" if plan_only else "completed_with_model_failures" if unexpected else "completed",
        "formulation_count_including_water": len(formulations),
        "composition_basis": "mass", "cold_status_counts": dict(Counter(r["status"] for r in cold)),
        "hot_status_counts": dict(Counter(r["status"] for r in hot)),
        "fit_endpoint_status_counts": dict(Counter(r["status"] for r in endpoint_rows)),
        "unexpected_model_failure_count": unexpected, "physical_runs_planned": len(plan),
        "physical_results_received": 0, "experimental_winners": 0,
        "surface_enhancement_predicted": False, "entry_mass_uncertainty_updated": False,
        "interpretation": "Offline model screen and test plan, not physical validation. INCOMP is cold-liquid only; NRTL windows are provisional. Missing values do not imply zero benefit or zero risk."}
    write_json(output / "summary.json", summary)
    write_json(output / "progress.json", {"signature": signature, "phase": summary["run_status"]})
    return summary


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, default=HERE / "config.yaml")
    parser.add_argument("--output-dir", type=Path, default=Path("formulation-v5m1-results"))
    parser.add_argument("--plan-only", action="store_true")
    parser.add_argument("--resume", action="store_true")
    args = parser.parse_args()
    try:
        config = yaml.safe_load(args.config.read_text(encoding="utf-8"))
        build_formulations(config)
        models = None
        if not args.plan_only:
            # Fail before creating outputs if the pinned runtime is absent.
            for package, pin in (("CoolProp", "8.0.0"), ("thermo", "0.6.1"), ("chemicals", "1.5.2")):
                if version(package) != pin:
                    raise ValueError(f"requires {package}=={pin}; use the existing V5b environment")
            if not (ENTRY / "thermo_provider.py").is_file():
                raise ValueError("existing V5a thermo_provider.py is required")
            models = MixtureModels()
        contract = {**fingerprint(config), "plan_only": args.plan_only}
        signature = digest(contract)
        output = args.output_dir.resolve()
        if output.exists() and any(output.iterdir()) and not args.resume:
            raise ValueError("output must be empty; use --resume or a new directory")
        output.mkdir(parents=True, exist_ok=True)
        import fcntl
        with (output / ".lock").open("a") as lock:
            fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
            manifest = output / "manifest.json"
            if args.resume and (not manifest.is_file() or read_json(manifest).get("signature") != signature):
                raise ValueError("resume signature changed or manifest missing; use a new output directory")
            write_json(manifest, {"signature": signature, "contract": contract})
            summary = execute(config, output, signature, models)
        print(json.dumps(summary, indent=2, sort_keys=True))
        return 1 if summary["unexpected_model_failure_count"] else 0
    except (OSError, ValueError, RuntimeError, ImportError, KeyError, TypeError, yaml.YAMLError) as exc:
        parser.error(str(exc))


if __name__ == "__main__":
    raise SystemExit(main())
