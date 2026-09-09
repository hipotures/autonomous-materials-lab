#!/usr/bin/env python3
"""Compare fixed mixture models to externally supplied, source-linked measurements.

No values are generated as reference, no models are fitted, and no training-set
independence is inferred. This is NOT a general ThermoML XML parser.
"""
from __future__ import annotations

import argparse
from collections import defaultdict
import math
from pathlib import Path
from statistics import median

import yaml
from run_predictive_grid import load_json, write_json, write_csv, contract, file_hash
from v5m2_core import MODELS, WATER_CAS, number, validate_config, digest

SUPPORTED = {"bubble_pressure_pa": "Pa", "excess_enthalpy_j_mol": "J/mol"}


def validate_records(payload: dict) -> list[dict]:
    if not isinstance(payload, dict) or payload.get("schema") != "binary-mixture-measurements-v1":
        raise ValueError("expected binary-mixture-measurements-v1")
    rows = payload.get("observations")
    if not isinstance(rows, list) or not rows:
        raise ValueError("no measured observations supplied")
    output, seen = [], set()
    for row in rows:
        if not isinstance(row, dict) or not isinstance(row.get("id"), str) or not row["id"]:
            raise ValueError("missing observation id")
        if row["id"] in seen:
            raise ValueError("duplicate observation id")
        seen.add(row["id"])
        if row.get("evidence_kind") != "measured":
            raise ValueError("only measured evidence is accepted")
        if not isinstance(row.get("source"), str) or not row["source"].strip():
            raise ValueError("missing source/DOI")
        if row.get("water_cas") != WATER_CAS or row.get("composition_basis") != "liquid_mole":
            raise ValueError("explicit water identity and liquid_mole basis required")
        if not isinstance(row.get("additive_cas"), str) or row["additive_cas"] == WATER_CAS:
            raise ValueError("invalid additive CAS")
        prop = row.get("property")
        if SUPPORTED.get(prop) != row.get("unit") or prop not in SUPPORTED:
            raise ValueError("unsupported property or unit")
        if row.get("phase") != "homogeneous_liquid":
            raise ValueError("reference liquid phase must be explicit")
        x = number(row["liquid_additive_mole_fraction"], "x")
        if not 0 < x < 1:
            raise ValueError("reference composition outside (0, 1)")
        t = number(row["temperature_k"], "T", positive=True)
        p = number(row["pressure_pa"], "P", positive=True)
        value = number(row["value"], "reference value", positive=prop == "bubble_pressure_pa")
        uncertainty = row.get("standard_uncertainty")
        if uncertainty is not None and number(uncertainty, "uncertainty") < 0:
            raise ValueError("negative uncertainty")
        if prop == "bubble_pressure_pa" and abs(p/value-1) > 1e-8:
            raise ValueError("bubble pressure and measurement pressure must agree")
        output.append({**row, "temperature_k": t, "pressure_pa": p, "value": value})
    return sorted(output, key=lambda r: r["id"])


def compare_records(records: list[dict], predict) -> tuple[list[dict], dict]:
    output = []
    for ref in records:
        for model in MODELS:
            row = {**ref, "model": model, "prediction": None, "absolute_error": None,
                   "relative_error": None, "training_overlap_unknown": True,
                   "reference_is_independent_of_model_training": None}
            try:
                pred = number(predict(ref, model), "prediction", positive=ref["property"] == "bubble_pressure_pa")
                row.update(status="ok", prediction=pred, absolute_error=abs(pred-ref["value"]),
                           relative_error=(pred/ref["value"]-1) if ref["value"] != 0 else None)
            except (ValueError, KeyError, TypeError, RuntimeError, ArithmeticError) as exc:
                row.update(status="prediction_failed", error=str(exc))
            output.append(row)
    grouped = defaultdict(list)
    for row in output:
        grouped[(row["additive_cas"], row["property"], row["model"])].append(row)
    per_pair = []
    for (cas, prop, model), rows in sorted(grouped.items()):
        errs = [r["absolute_error"] for r in rows if r["status"] == "ok"]
        rels = [abs(r["relative_error"]) for r in rows if r["status"] == "ok" and r["relative_error"] is not None]
        per_pair.append({"additive_cas": cas, "property": prop, "model": model,
                         "attempted": len(rows), "failed": len(rows)-len(errs),
                         "median_absolute_error": median(errs) if errs else None,
                         "median_absolute_relative_error": median(rels) if rels else None,
                         "sources": sorted({r["source"] for r in rows})})
    return output, {"pair_property_results": per_pair, "models_refitted": False,
                    "calibrated_uncertainty": None, "blind_holdout_claimed": False,
                    "interpretation": "Fixed-model external-data comparison. Pairwise aggregation avoids treating temperature points as independent mixtures. Source provenance does not establish absence from model training."}


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--references", type=Path, required=True)
    parser.add_argument("--config", type=Path, default=Path(__file__).with_name("config-v5m2.yaml"))
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    try:
        records = validate_records(load_json(args.references))
        config = yaml.safe_load(args.config.read_text())
        validate_config(config)
        runtime = contract(config, catalog_only=False, limit_pairs=None)
    except (OSError, ValueError, KeyError, TypeError, ImportError) as exc:
        parser.error(str(exc))
    if args.output_dir.exists() and any(args.output_dir.iterdir()):
        parser.error("validation output must be empty")
    args.output_dir.mkdir(parents=True, exist_ok=True)
    write_json(args.output_dir / "frozen-reference.json", {"schema": "binary-mixture-measurements-v1", "observations": records})
    from v5m2_backend import discover_catalog, BinaryModel
    from v5m2_core import mole_to_mass
    catalog = discover_catalog(config, sorted({r["additive_cas"] for r in records}))
    pairs = {r["cas_number"]: r for r in catalog["candidates"]}
    backends = {}
    def predict(ref, model):
        cas = ref["additive_cas"]
        if cas not in pairs: raise ValueError("reference_pair_outside_catalog_model_scope")
        key = (cas, model)
        if key not in backends: backends[key] = BinaryModel(pairs[cas], model, config)
        backend = backends[key]
        t, p = ref["temperature_k"], ref["pressure_pa"]
        if not 273.15 <= t <= 500 or not 25000 <= p <= 300000:
            raise ValueError("reference_outside_operational_screening_domain")
        x = ref["liquid_additive_mole_fraction"]
        w = mole_to_mass([1-x, x], backend.mw)[1]
        # Check that the homogeneous reference composition is not predicted to
        # split into multiple liquid phases. A finite flash is not a validation.
        bubble = backend.bubble(w, t)
        check_p = max(p, bubble["bubble_pressure_pa"]*1.0001)
        if check_p > 300000: raise ValueError("stability_probe_outside_pressure_scope")
        state = backend.state(w, t, check_p)
        if state["liquid_phase_count"] != 1 or state["vapor_mole_fraction"] > 1e-7:
            raise ValueError("homogeneous_liquid_reference_not_reproduced")
        return bubble[ref["property"]]
    rows, summary = compare_records(records, predict)
    write_csv(args.output_dir / "reference-comparisons.csv", rows)
    write_json(args.output_dir / "summary.json", summary)
    write_json(args.output_dir / "manifest.json", {"runtime": runtime, "reference_sha256": file_hash(args.references),
               "validation_code_sha256": file_hash(Path(__file__)), "reference_schema": "binary-mixture-measurements-v1"})
    print(f"Measured records: {len(records)}; prediction failures: {sum(r['status'] != 'ok' for r in rows)}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
