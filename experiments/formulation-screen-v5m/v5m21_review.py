"""Lossless publication and water-aware reanalysis of frozen V5m-2 results.

No property model is evaluated. Reported gains are model hypotheses. Gzip is
an archive format; bounded UTF-8 report shards support text-only connectors.
"""
from __future__ import annotations

from collections import Counter, defaultdict
import csv
import gzip
import hashlib
import io
import json
import math
from pathlib import Path
from typing import Any, Iterable

from v5m2_core import MODELS, THERMO_OBJECTIVES, FLOW_OBJECTIVES, compare_models, pareto

VERSION = "v5m-2.1"
# An exact allowlist: never include checkpoints, logs, secrets or arbitrary files.
SOURCE_FILES = (
    "artifact-hashes.json", "catalog-rejections.csv", "catalog.csv", "catalog.json",
    "coarse-tradeoffs.csv", "grid-plan.json", "manifest.json", "pareto-thermodynamic.csv",
    "pareto-with-flow-proxy.csv", "prediction-grid.csv", "predictions.json",
    "refinement-diagnostics.csv", "refinement-plan.json", "summary.json", "tradeoffs.csv",
    "validation-shortlist.csv",
)
REQUIRED_FILES = (
    "artifact-hashes.json", "catalog.json", "grid-plan.json", "manifest.json",
    "prediction-grid.csv", "refinement-plan.json", "summary.json",
)
REQUIRED_COLUMNS = (
    "pair_id", "cas_number", "name", "smiles", "model", "additive_mass_fraction",
    "outlet_temperature_k", "pressure_pa", "status", "model_comparison_eligible",
    "same_model_water_delta_h_ratio", "storage_bubble_pressure_ratio",
    "viscosity_ratio_proxy",
)
NUMERIC_COLUMNS = (
    "additive_mass_fraction", "outlet_temperature_k", "pressure_pa", "inlet_temperature_k",
    "same_model_water_delta_h_ratio", "storage_bubble_pressure_ratio", "viscosity_ratio_proxy",
    "delta_h_j_kg", "same_model_water_delta_h_j_kg", "heos_water_delta_h_j_kg",
    "heos_water_delta_h_ratio", "endpoint_relative_error",
)


def finite(value: Any, label: str, *, positive: bool = False) -> float:
    if isinstance(value, bool):
        raise ValueError(f"{label}: boolean is not numeric")
    try:
        result = float(value)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"{label}: not a number") from exc
    if not math.isfinite(result) or (positive and result <= 0):
        raise ValueError(f"{label}: nonfinite or nonpositive")
    return result


def json_loads(text: str) -> Any:
    def reject(value: str):
        raise ValueError("nonfinite JSON constant: " + value)
    def unique(items):
        result = {}
        for k, v in items:
            if k in result:
                raise ValueError("duplicate JSON key: " + k)
            result[k] = v
        return result
    return json.loads(text, parse_constant=reject, object_pairs_hook=unique)


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for data in iter(lambda: stream.read(1024 * 1024), b""):
            h.update(data)
    return h.hexdigest()


def write_json(path: Path, data: Any) -> None:
    tmp = path.with_name(path.name + ".tmp")
    tmp.write_text(json.dumps(data, indent=2, sort_keys=True, allow_nan=False) + "\n", encoding="utf-8")
    tmp.replace(path)


def load_json(path: Path) -> Any:
    return json_loads(path.read_text(encoding="utf-8"))


def compress_file(source: Path, destination: Path, expected_sha256: str | None = None) -> dict:
    """Stream exact bytes; omit timestamp/name, verify decompression, never unlink input."""
    if source.is_symlink() or not source.is_file() or source.resolve() == destination.resolve():
        raise ValueError("compression requires distinct regular input and output files")
    if destination.exists() or destination.is_symlink():
        raise ValueError("compressed output already exists: " + str(destination))
    tmp = destination.with_name(destination.name + ".tmp")
    if tmp.exists() or tmp.is_symlink():
        raise ValueError("compression temporary file already exists")
    digest = hashlib.sha256()
    size = 0
    try:
        with source.open("rb") as inp, tmp.open("xb") as out:
            with gzip.GzipFile(filename="", mode="wb", fileobj=out, mtime=0, compresslevel=6) as zipped:
                for block in iter(lambda: inp.read(1024 * 1024), b""):
                    digest.update(block); size += len(block); zipped.write(block)
        if expected_sha256 is not None and digest.hexdigest() != expected_sha256:
            raise ValueError("source changed before/during compression: " + source.name)
        check = hashlib.sha256()
        with gzip.open(tmp, "rb") as stream:
            for block in iter(lambda: stream.read(1024 * 1024), b""):
                check.update(block)
        if digest.digest() != check.digest():
            raise ValueError("gzip round-trip verification failed")
        tmp.replace(destination)
    finally:
        if tmp.exists():
            tmp.unlink()
    return {"file": destination.name, "original_bytes": size,
            "compressed_bytes": destination.stat().st_size,
            "original_sha256": digest.hexdigest(), "compressed_sha256": sha256(destination),
            "round_trip_verified": True}


def csv_text(rows: Iterable[dict], fields: list[str]) -> str:
    stream = io.StringIO(newline="")
    writer = csv.DictWriter(stream, fieldnames=fields, lineterminator="\n")
    writer.writeheader()
    for row in rows:
        writer.writerow({k: json.dumps(v, sort_keys=True, allow_nan=False) if isinstance(v, (list, dict)) else v
                         for k, v in row.items()})
    return stream.getvalue()


def write_table(output: Path, label: str, rows: list[dict], *, text: bool,
                max_bytes: int = 48 * 1024) -> list[dict]:
    """Archive every table and optionally publish all rows in bounded text shards."""
    if not label.replace("-", "").isalnum():
        raise ValueError("invalid table label")
    fields = sorted({k for row in rows for k in row}) or ["status"]
    temp = output / ("table-" + label + ".csv.tmp")
    temp.write_text(csv_text(rows, fields), encoding="utf-8", newline="")
    try:
        archived = compress_file(temp, output / ("table-" + label + ".csv.gz"))
    finally:
        temp.unlink(missing_ok=True)
    entries = [{**archived, "table": label, "rows": len(rows), "kind": "full_table_gzip"}]
    if not text:
        return entries
    header = csv_text([], fields)
    chunk = header; count = 0; part = 1
    def flush():
        name = f"report-{label}-{part:04d}.csv"
        path = output / name
        path.write_text(chunk, encoding="utf-8", newline="")
        return {"file": name, "table": label, "kind": "utf8_report", "rows": count,
                "bytes": path.stat().st_size, "sha256": sha256(path)}
    for row in rows:
        line = csv_text([row], fields)[len(header):]
        if len((header + line).encode()) > max_bytes:
            raise ValueError("single report row exceeds byte limit: " + label)
        if len((chunk + line).encode()) > max_bytes and count:
            entries.append(flush()); part += 1; chunk = header; count = 0
        chunk += line; count += 1
    entries.append(flush())
    return entries


def read_grid(path: Path) -> list[dict]:
    result = []
    csv.field_size_limit(16 * 1024 * 1024)
    with path.open("r", encoding="utf-8-sig", newline="") as stream:
        reader = csv.DictReader(stream)
        header = reader.fieldnames or []
        if len(header) != len(set(header)) or set(REQUIRED_COLUMNS) - set(header):
            raise ValueError("invalid prediction-grid header")
        for n, raw in enumerate(reader, 2):
            if None in raw or any(v is None for v in raw.values()):
                raise ValueError(f"malformed CSV record {n}")
            row = {k: raw[k] for k in ("pair_id", "cas_number", "name", "smiles", "model", "status")}
            flag = raw["model_comparison_eligible"].lower()
            if flag not in {"true", "false"}:
                raise ValueError(f"invalid eligibility flag at record {n}")
            row["model_comparison_eligible"] = flag == "true"
            for key in NUMERIC_COLUMNS:
                val = raw.get(key, "")
                row[key] = None if val == "" else finite(val, f"{key} at record {n}")
            for key in ("inlet", "outlet", "storage_vle"):
                value = raw.get(key)
                row[key] = json_loads(value) if value else None
                if row[key] is not None and not isinstance(row[key], dict):
                    raise ValueError("invalid state object at record " + str(n))
            for key in ("outlet_temperature_k", "pressure_pa"):
                finite(row[key], key, positive=True)
            w = finite(row["additive_mass_fraction"], "mass fraction")
            if not 0 < w < 1:
                raise ValueError("input grid must contain additive fractions in (0, 1)")
            if row["model"] not in MODELS or not row["pair_id"] or not row["status"]:
                raise ValueError("invalid grid identity/model/status")
            if row["status"] == "ok":
                for key in ("same_model_water_delta_h_ratio", "storage_bubble_pressure_ratio"):
                    finite(row[key], key, positive=True)
                if row["viscosity_ratio_proxy"] is not None:
                    finite(row["viscosity_ratio_proxy"], "viscosity proxy", positive=True)
            result.append(row)
    return result


def validate_source(source: Path) -> tuple[dict, dict, list[dict], dict[str, str]]:
    """Verify original hashes, model/config contract and every planned state key."""
    for name in REQUIRED_FILES:
        p = source / name
        if p.is_symlink() or not p.is_file():
            raise ValueError("missing regular source artifact: " + name)
    expected_hashes = load_json(source / "artifact-hashes.json")
    actual = {}
    for name in SOURCE_FILES:
        p = source / name
        if p.is_symlink():
            raise ValueError("source artifact is a symlink: " + name)
        if p.is_file():
            actual[name] = sha256(p)
            if name != "artifact-hashes.json" and expected_hashes.get(name) != actual[name]:
                raise ValueError("source artifact hash mismatch: " + name)
    summary = load_json(source / "summary.json")
    if summary.get("study_version") != "v5m-2" or summary.get("study_complete") is not True:
        raise ValueError("requires a completed, non-smoke V5m-2 run")
    manifest = load_json(source / "manifest.json")
    contract = manifest["contract"]
    from v5m2_core import digest
    if manifest.get("version") != "v5m-2" or manifest.get("signature") != digest(contract):
        raise ValueError("source manifest signature mismatch")
    config = contract["config"]
    if config.get("models") != list(MODELS):
        raise ValueError("unexpected source model list")
    grid = load_json(source / "grid-plan.json")
    catalog = load_json(source / "catalog.json")["candidates"]
    if grid.get("models") != list(MODELS) or grid.get("composition_basis") != "mass":
        raise ValueError("unsupported source grid")
    pairs = {p["pair_id"]: p for p in catalog}
    if len(pairs) != len(catalog) or len(pairs) != grid["pairs"]:
        raise ValueError("catalog/grid pair mismatch")
    planned = {pair: set(grid["mass_fractions"]) for pair in pairs}
    for task in load_json(source / "refinement-plan.json")["tasks"]:
        pair = task["pair"]["pair_id"]
        if pair not in pairs or set(task["mass_fractions"]) & planned[pair]:
            raise ValueError("invalid/duplicate refinement task")
        planned[pair].update(task["mass_fractions"])
    conditions = [(float(t), float(p)) for t in config["grid"]["outlet_temperatures_k"]
                  for p in config["grid"]["pressures_pa"]]
    if grid["conditions_per_composition_model"] != len(conditions):
        raise ValueError("source condition-count mismatch")
    wanted = {(pair, float(w), t, p, model) for pair, ws in planned.items() for w in ws
              for t, p in conditions for model in MODELS}
    rows = read_grid(source / "prediction-grid.csv")
    found = set()
    for row in rows:
        key = (row["pair_id"], row["additive_mass_fraction"], row["outlet_temperature_k"],
               row["pressure_pa"], row["model"])
        if key in found:
            raise ValueError("duplicate prediction state")
        found.add(key)
        pair = pairs.get(row["pair_id"])
        if pair is None or any(row[k] != pair[k] for k in ("cas_number", "name", "smiles")):
            raise ValueError("prediction/catalog identity mismatch")
    if found != wanted:
        raise ValueError(f"grid plan mismatch: missing={len(wanted-found)} unexpected={len(found-wanted)}")
    if len(rows) != summary["predicted_state_count"] or dict(Counter(r["status"] for r in rows)) != summary["state_status_counts"]:
        raise ValueError("source summary/grid mismatch")
    return summary, config, rows, actual


def water_row(conditions: int) -> dict:
    return {"pair_id": "pure-water-reference", "cas_number": "7732-18-5", "name": "water",
            "smiles": "O", "additive_mass_fraction": 0., "comparison_eligible": True,
            "complete_grid": True, "paired_condition_count": conditions, "expected_condition_count": conditions,
            "minimum_paired_delta_h_ratio": 1., "maximum_storage_bubble_pressure_ratio": 1.,
            "maximum_viscosity_ratio_proxy": 1., "maximum_model_spread_delta_h_ratio": 0.,
            "reference_basis": "normalization_identity_not_new_measurement", "failure_counts": {},
            "handling_approved": False, "experimental_winner": False}


def water_dominates(row: dict) -> bool:
    h = row["minimum_paired_delta_h_ratio"]
    p = row["maximum_storage_bubble_pressure_ratio"]
    return row["comparison_eligible"] and h <= 1 and p >= 1 and (h < 1 or p > 1)


def annotate(row: dict, gain_margin: float) -> dict:
    if not row["comparison_eligible"]:
        return {**row, "classification": "incomplete_comparison", "water_dominated_thermodynamic": None}
    ratio = row["minimum_paired_delta_h_ratio"]
    return {**row, "water_dominated_thermodynamic": water_dominates(row),
            "enthalpy_penalty_fraction": max(0., 1-ratio),
            "mass_ratio_equal_net_heat": 1/ratio,
            "required_net_heat_reduction_for_mass_parity": max(0., 1-ratio),
            "model_gain_above_margin": ratio > 1+gain_margin,
            "classification": "model_gain_hypothesis" if ratio > 1+gain_margin else
                "below_reporting_margin" if ratio > 1 else "no_paired_enthalpy_gain",
            "margin_is_calibrated_uncertainty": False, "rankable": False}


def calculate_review(rows: list[dict], config: dict, gain_margin: float = 0.001) -> dict:
    gain_margin = finite(gain_margin, "gain margin")
    if not 0 <= gain_margin < 1:
        raise ValueError("gain margin outside [0, 1)")
    conditions = [(float(t), float(p)) for t in config["grid"]["outlet_temperatures_k"]
                  for p in config["grid"]["pressures_pa"]]
    global_rows = [annotate(r, gain_margin) for r in compare_models(rows, len(conditions), config)]
    base = water_row(len(conditions))
    global_front = pareto(global_rows + [base], THERMO_OBJECTIVES)
    flow_front = pareto(global_rows + [base], FLOW_OBJECTIVES)
    by_regime, per_model, paired_conditions, regime_summary, disagreements = [], [], [], [], []
    for t, pressure in conditions:
        subset = [r for r in rows if (r["outlet_temperature_k"], r["pressure_pa"]) == (t, pressure)]
        paired = [annotate(r, gain_margin) for r in compare_models(subset, 1, config)]
        eligible = [r for r in paired if r["comparison_eligible"]]
        front = pareto(paired + [water_row(1)], THERMO_OBJECTIVES)
        front_keys = {(r["pair_id"], r["additive_mass_fraction"]) for r in front}
        regime = {"outlet_temperature_k": t, "pressure_pa": pressure}
        paired_conditions.extend({**r, **regime, "pareto_with_water":
                                  (r["pair_id"], r["additive_mass_fraction"]) in front_keys} for r in paired)
        grouped = defaultdict(list)
        for row in paired:
            grouped[row["pair_id"]].append(row)
        for pair, values in sorted(grouped.items()):
            good = [r for r in values if r["comparison_eligible"]]
            best = max(good, key=lambda r: (r["minimum_paired_delta_h_ratio"], -r["additive_mass_fraction"])) if good else None
            if good:
                discrepant = max(good, key=lambda r: (r["maximum_model_spread_delta_h_ratio"], -r["additive_mass_fraction"]))
                w = discrepant["additive_mass_fraction"]
                comparison = {**regime, "pair_id": pair, "name": values[0]["name"],
                              "additive_mass_fraction": w,
                              "spread_delta_h_ratio": discrepant["maximum_model_spread_delta_h_ratio"]}
                for model in MODELS:
                    raw = next(r for r in subset if r["pair_id"] == pair and r["additive_mass_fraction"] == w and r["model"] == model)
                    state = raw.get("outlet") or {}
                    comparison[model + "_delta_h_ratio"] = raw["same_model_water_delta_h_ratio"]
                    comparison[model + "_vapor_mole_fraction"] = state.get("vapor_mole_fraction")
                    comparison[model + "_liquid_phase_count"] = state.get("liquid_phase_count")
                    comparison[model + "_inlet_excess_enthalpy_j_mol"] = (raw.get("storage_vle") or {}).get("excess_enthalpy_j_mol")
                disagreements.append(comparison)
            by_regime.append({**regime, "pair_id": pair, "name": values[0]["name"], "cas_number": values[0]["cas_number"],
                "composition_count": len(values), "paired_composition_count": len(good),
                "failed_or_unpaired_composition_count": len(values)-len(good),
                "water_dominated_composition_count": sum(water_dominates(r) for r in good),
                "paired_gain_above_margin_count": sum(r["minimum_paired_delta_h_ratio"] > 1+gain_margin for r in good),
                "best_mass_fraction": best["additive_mass_fraction"] if best else None,
                "best_paired_delta_h_ratio": best["minimum_paired_delta_h_ratio"] if best else None,
                "best_mass_ratio_equal_net_heat": best["mass_ratio_equal_net_heat"] if best else None,
                "required_net_heat_reduction_at_best": best["required_net_heat_reduction_for_mass_parity"] if best else None,
                "max_model_spread_ratio": max((r["maximum_model_spread_delta_h_ratio"] for r in good), default=None)})
        for model in MODELS:
            groups = defaultdict(list)
            for row in subset:
                if row["model"] == model:
                    groups[row["pair_id"]].append(row)
            for pair, values in sorted(groups.items()):
                good = [r for r in values if r["status"] == "ok" and r["model_comparison_eligible"]]
                best = max(good, key=lambda r: (r["same_model_water_delta_h_ratio"], -r["additive_mass_fraction"])) if good else None
                out = (best.get("outlet") or {}) if best else {}
                per_model.append({**regime, "pair_id": pair, "name": values[0]["name"], "model": model,
                    "eligible_composition_count": len(good), "status_counts": dict(Counter(r["status"] for r in values)),
                    "model_gain_above_margin_count": sum(r["same_model_water_delta_h_ratio"] > 1+gain_margin for r in good),
                    "best_mass_fraction": best["additive_mass_fraction"] if best else None,
                    "best_delta_h_ratio": best["same_model_water_delta_h_ratio"] if best else None,
                    "best_outlet_vapor_mole_fraction": out.get("vapor_mole_fraction"),
                    "best_outlet_liquid_phase_count": out.get("liquid_phase_count")})
        regime_summary.append({**regime, "eligible_formulation_count": len(eligible),
            "paired_gain_any_positive_count": sum(r["minimum_paired_delta_h_ratio"] > 1 for r in eligible),
            "paired_gain_above_margin_count": sum(r["minimum_paired_delta_h_ratio"] > 1+gain_margin for r in eligible),
            "water_dominated_count": sum(water_dominates(r) for r in eligible),
            "pareto_count_including_water": len(front),
            "water_on_pareto": any(r["pair_id"] == base["pair_id"] for r in front)})
    eligible = [r for r in global_rows if r["comparison_eligible"]]
    summary = {"study_version": VERSION, "study_complete": True, "run_kind": "frozen_grid_reanalysis",
        "source_state_count": len(rows), "complete_two_model_formulation_count": len(eligible),
        "water_dominated_global_thermodynamic_count": sum(water_dominates(r) for r in eligible),
        "global_pareto_count_including_water": len(global_front),
        "water_on_global_pareto": any(r["pair_id"] == base["pair_id"] for r in global_front),
        "global_pareto_additive_count": sum(r["pair_id"] != base["pair_id"] for r in global_front),
        "global_paired_gain_any_positive_count": sum(r["minimum_paired_delta_h_ratio"] > 1 for r in eligible),
        "global_paired_gain_above_margin_count": sum(r["minimum_paired_delta_h_ratio"] > 1+gain_margin for r in eligible),
        "gain_reporting_margin_fraction": gain_margin, "margin_is_calibrated_uncertainty": False,
        "condition_results": regime_summary, "new_thermodynamic_states_computed": 0,
        "experimental_winners": 0, "rankable_promotions": 0, "empirical_uncertainty_calibrated": False,
        "interpretation": "Retrospective regime diagnostics, not a new validation set. Pressure remains a descriptive Pareto objective, not an apparatus constraint. Viscosity is the legacy proxy. Break-even heat reduction is required, not achieved."}
    return {"summary": summary, "global-formulations": global_rows, "paired-conditions": paired_conditions,
            "regime-summary": by_regime, "model-regimes": per_model, "disagreements": disagreements,
            "pareto-water": global_front, "pareto-water-flow-proxy": flow_front}
