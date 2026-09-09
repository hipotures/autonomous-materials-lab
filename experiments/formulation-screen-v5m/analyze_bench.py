#!/usr/bin/env python3
"""Compare declared physical measurements; never use model outputs as experiments."""
from __future__ import annotations

import argparse
from collections import Counter, defaultdict
import csv
import hashlib
import json
from pathlib import Path
from statistics import median

import yaml
from formulations import digest, number, read_json, write_json
from run_formulations import write_csv

MATCH_KEYS = ("apparatus_id", "geometry_id", "material_lot_id", "surface_spec_id", "load_program_id")
PASS_FIELDS = ("safety_review_pass", "composition_verified", "phase_stable", "cold_flow_pass",
               "fresh_coupon", "fluid_path_clean", "completed")
NUMERIC_FIELDS = ("mass_used_g", "mass_uncertainty_g", "incident_energy_j", "energy_uncertainty_j",
    "duration_s", "wall_peak_k", "wall_uncertainty_k", "feed_pressure_max_pa", "pressure_uncertainty_pa",
    "hydraulic_resistance_ratio", "residue_mg")


def validate_protocol(p):
    if p.get("version") != "v5m-1-bench":
        raise ValueError("expected v5m-1-bench protocol")
    for key in MATCH_KEYS:
        if not isinstance(p.get(key), str) or not p[key].strip():
            raise ValueError("declare and freeze protocol before tests: missing " + key)
    for key in ("wall_temperature_limit_k", "max_feed_pressure_pa", "target_incident_energy_j",
                "target_duration_s", "max_hydraulic_resistance_ratio"):
        number(p.get(key), key)
    number(p.get("max_residue_mg"), "max_residue_mg", strict=False)
    for key in ("energy_relative_tolerance", "duration_relative_tolerance",
                "max_control_mass_relative_difference", "minimum_mass_saving_fraction"):
        if number(p.get(key), key, strict=False) >= 1:
            raise ValueError(key + " must be <1")
    if not isinstance(p.get("minimum_replicates"), int) or isinstance(p["minimum_replicates"], bool) or p["minimum_replicates"] < 3:
        raise ValueError("at least three independent replicate blocks required")
    if p.get("uncertainty_semantics") != "declared_absolute_bounds_not_confidence_intervals":
        raise ValueError("declare absolute uncertainty bounds; confidence limits are not accepted as hard bounds")


def valid_measurement(row, expected, protocol):
    reasons = []
    values = {}
    for key in ("formulation_id", "block_id", "role", "replicate"):
        if str(row.get(key, "")) != str(expected[key]):
            reasons.append("plan_mismatch:" + key)
    if row.get("evidence_kind") != "measured":
        reasons.append("not_declared_physical_measurement")
    for key in MATCH_KEYS:
        if row.get(key) != protocol[key]:
            reasons.append("condition_mismatch:" + key)
    for key in ("coupon_id", "preparation_record_id", "source_data_id", "safety_review_id"):
        if not isinstance(row.get(key), str) or not row[key].strip():
            reasons.append("missing:" + key)
    for key in PASS_FIELDS:
        if str(row.get(key, "")).lower() != "true":
            reasons.append("not_passed:" + key)
    for key in NUMERIC_FIELDS:
        try:
            if isinstance(row.get(key), bool):
                raise ValueError("boolean is not a measurement")
            values[key] = number(float(row.get(key, "")), key,
                strict=key not in {"residue_mg", "mass_uncertainty_g", "energy_uncertainty_j", "wall_uncertainty_k", "pressure_uncertainty_pa"})
        except (ValueError, TypeError):
            reasons.append("missing_or_invalid:" + key)
    if len(values) == len(NUMERIC_FIELDS):
        v = values
        if v["mass_used_g"] <= v["mass_uncertainty_g"] or v["incident_energy_j"] <= v["energy_uncertainty_j"]:
            reasons.append("nonpositive_lower_bound")
        if v["wall_peak_k"] + v["wall_uncertainty_k"] > protocol["wall_temperature_limit_k"]:
            reasons.append("thermal_limit_exceeded")
        if v["feed_pressure_max_pa"] + v["pressure_uncertainty_pa"] > protocol["max_feed_pressure_pa"]:
            reasons.append("pressure_limit_exceeded")
        if abs(v["incident_energy_j"] - protocol["target_incident_energy_j"]) + v["energy_uncertainty_j"] > protocol["energy_relative_tolerance"] * protocol["target_incident_energy_j"]:
            reasons.append("incident_load_not_matched")
        if abs(v["duration_s"] / protocol["target_duration_s"] - 1) > protocol["duration_relative_tolerance"]:
            reasons.append("duration_not_matched")
        if v["hydraulic_resistance_ratio"] > protocol["max_hydraulic_resistance_ratio"]:
            reasons.append("hydraulic_resistance_limit_exceeded")
        if v["residue_mg"] > protocol["max_residue_mg"]:
            reasons.append("residue_limit_exceeded")
    return {**row, **values, "valid": not reasons, "reasons": reasons}


def analyze(plan, measurements, protocol):
    validate_protocol(protocol)
    if plan.get("version") != "v5m-1":
        raise ValueError("expected generated V5m-1 plan")
    runs = plan["runs"]
    expected = {r["run_id"]: r for r in runs}
    if not expected or len(expected) != len(runs):
        raise ValueError("empty plan or duplicate plan run ids")
    ids = [r.get("run_id") for r in measurements]
    if len(ids) != len(set(ids)) or any(i not in expected for i in ids):
        raise ValueError("unknown or duplicate measurement run id")
    coupon_counts = Counter(r.get("coupon_id") for r in measurements if r.get("coupon_id"))
    checked = {}
    for row in measurements:
        item = valid_measurement(row, expected[row["run_id"]], protocol)
        if coupon_counts.get(row.get("coupon_id"), 0) > 1:
            item["valid"] = False
            item["reasons"].append("coupon_reused_surface_confounding")
        checked[row["run_id"]] = item
    blocks = defaultdict(list)
    for r in runs:
        blocks[r["block_id"]].append(r)
    comparisons = []
    for block, members in sorted(blocks.items()):
        before = [r for r in members if r["role"] == "water_before"]
        after = [r for r in members if r["role"] == "water_after"]
        if len(before) != 1 or len(after) != 1 or any(r["formulation_id"] != "water" for r in before + after):
            raise ValueError("every block requires exactly two bracketing water controls")
        ctrl = [checked.get(r["run_id"]) for r in before + after]
        control_valid = all(c and c["valid"] for c in ctrl)
        if control_valid:
            c0, c1 = ctrl
            control_valid = abs(c0["mass_used_g"] - c1["mass_used_g"]) / min(c0["mass_used_g"], c1["mass_used_g"]) <= protocol["max_control_mass_relative_difference"]
        for r in members:
            if r["role"] != "candidate":
                continue
            measured = checked.get(r["run_id"])
            why = list(measured["reasons"]) if measured else ["measurement_missing"]
            if not control_valid:
                why.append("water_controls_missing_invalid_or_drifted")
            output = {"run_id": r["run_id"], "formulation_id": r["formulation_id"],
                "replicate": r["replicate"], "block_id": block, "valid_comparison": not why,
                "reasons": why, "upper_mass_ratio": None, "upper_specific_mass_ratio": None,
                "conservative_mass_saving_fraction": None}
            if not why:
                c = measured
                numerator = c["mass_used_g"] + c["mass_uncertainty_g"]
                lower_control_mass = min(x["mass_used_g"] - x["mass_uncertainty_g"] for x in ctrl)
                lower_control_specific = min((x["mass_used_g"] - x["mass_uncertainty_g"]) /
                                             (x["incident_energy_j"] + x["energy_uncertainty_j"]) for x in ctrl)
                mass_ratio = numerator / lower_control_mass
                specific_ratio = numerator / (c["incident_energy_j"] - c["energy_uncertainty_j"]) / lower_control_specific
                output.update(upper_mass_ratio=mass_ratio, upper_specific_mass_ratio=specific_ratio,
                    conservative_mass_saving_fraction=1 - max(mass_ratio, specific_ratio))
            comparisons.append(output)
    grouped = defaultdict(list)
    for row in comparisons:
        grouped[row["formulation_id"]].append(row)
    results = []
    for fid, rows in sorted(grouped.items()):
        valid = [r for r in rows if r["valid_comparison"]]
        enough = len({r["replicate"] for r in valid}) >= protocol["minimum_replicates"]
        worst = min((r["conservative_mass_saving_fraction"] for r in valid), default=None)
        passed = enough and len(valid) == len(rows) and worst >= protocol["minimum_mass_saving_fraction"]
        results.append({"formulation_id": fid, "planned_replicates": len(rows),
            "valid_replicates": len(valid), "worst_declared_bound_mass_saving_fraction": worst,
            "median_observed_bound_saving_fraction": median(r["conservative_mass_saving_fraction"] for r in valid) if valid else None,
            "status": "bench_promising_requires_confirmation" if passed else "incomplete_or_failed_comparisons" if not enough or len(valid) != len(rows) else "no_clear_mass_advantage",
            "flight_validated": False, "statistical_confidence_claimed": False})
    return {"study_version": "v5m-1-bench", "measurement_count": len(measurements),
        "invalid_measurement_count": sum(not r["valid"] for r in checked.values()),
        "valid_candidate_comparison_count": sum(r["valid_comparison"] for r in comparisons),
        "bench_promising_count": sum(r["status"] == "bench_promising_requires_confirmation" for r in results),
        "physical_data_authenticity_independently_verified": False,
        "interpretation": "User-declared measurements only. Absolute bounds are not statistical confidence. Specific mass uses incident exposure, not absorbed heat or latent heat. All planned replicates must be valid for a promising label.",
        "formulations": results, "comparisons": comparisons, "measurement_audit": list(checked.values())}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--plan", type=Path, required=True)
    parser.add_argument("--measurements", type=Path, required=True)
    parser.add_argument("--protocol", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, default=Path("bench-v5m1-results"))
    args = parser.parse_args()
    try:
        with args.measurements.open(newline="", encoding="utf-8") as f:
            data = list(csv.DictReader(f))
        result = analyze(read_json(args.plan), data, yaml.safe_load(args.protocol.read_text()))
        if args.output_dir.exists() and any(args.output_dir.iterdir()):
            raise ValueError("use a new output directory; measurements are never overwritten")
        args.output_dir.mkdir(parents=True, exist_ok=True)
        write_csv(args.output_dir / "measurement-audit.csv", result.pop("measurement_audit"))
        write_csv(args.output_dir / "paired-comparisons.csv", result.pop("comparisons"))
        write_json(args.output_dir / "summary.json", result)
        files = [args.plan, args.measurements, args.protocol, Path(__file__)]
        write_json(args.output_dir / "manifest.json", {"inputs_and_analysis_sha256":
            {str(p): hashlib.sha256(p.read_bytes()).hexdigest() for p in files}})
        print(json.dumps(result, indent=2, sort_keys=True))
        return 0
    except (OSError, ValueError, KeyError, TypeError, yaml.YAMLError) as exc:
        parser.error(str(exc))


if __name__ == "__main__":
    raise SystemExit(main())
