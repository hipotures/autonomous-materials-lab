"""Known-component formulations and an explicitly non-executable bench plan."""
from __future__ import annotations

from collections import Counter
import hashlib
import json
import math
from pathlib import Path
import random
from typing import Any

# Bounds and mass-fraction convention: CoolProp Incompressibles documentation.
# These are library fit bounds, not validated boiling or hardware safety limits.
PROFILES = {
    "ethanol": ("MEA", 173.15, 313.15, 0.60),
    "propylene_glycol": ("MPG", 173.15, 373.15, 0.60),
    "glycerol": ("MGL", 173.15, 313.15, 0.60),
    "sodium_chloride": ("MNA", 173.15, 313.15, 0.23),
}
EXPECTED_CAS = {"ethanol": "64-17-5", "propylene_glycol": "57-55-6",
                "glycerol": "56-81-5", "sodium_chloride": "7647-14-5",
                "sds": "151-21-3", "silica_suspension": "7631-86-9"}
COLD_FIELDS = {"density_kg_m3": "D", "cp_j_kg_k": "C", "viscosity_pa_s": "V",
               "conductivity_w_m_k": "L", "enthalpy_j_kg": "H"}
MEASUREMENT_FIELDS = [
    "run_id", "formulation_id", "block_id", "role", "replicate", "evidence_kind",
    "apparatus_id", "geometry_id", "material_lot_id", "surface_spec_id", "coupon_id",
    "load_program_id", "preparation_record_id", "source_data_id", "safety_review_id",
    "safety_review_pass", "composition_verified", "phase_stable", "cold_flow_pass",
    "fresh_coupon", "fluid_path_clean", "completed", "mass_used_g", "mass_uncertainty_g",
    "incident_energy_j", "energy_uncertainty_j", "duration_s", "wall_peak_k",
    "wall_uncertainty_k", "feed_pressure_max_pa", "pressure_uncertainty_pa",
    "hydraulic_resistance_ratio", "residue_mg", "notes",
]


def digest(value: Any) -> str:
    return hashlib.sha256(json.dumps(value, sort_keys=True, allow_nan=False,
        separators=(",", ":")).encode()).hexdigest()


def number(value: Any, name: str, *, minimum: float = 0.0, strict: bool = True) -> float:
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise ValueError(f"{name}: expected a number")
    value = float(value)
    if not math.isfinite(value) or (value <= minimum if strict else value < minimum):
        raise ValueError(f"{name}: out of range or nonfinite")
    return value


def validate_config(config: dict) -> None:
    if config.get("version") != "v5m-1":
        raise ValueError("expected config version v5m-1")
    number(config["batch_mass_g"], "batch_mass_g")
    ts = config["cold_temperatures_k"]
    if len(ts) < 2 or ts != sorted(set(ts)):
        raise ValueError("cold_temperatures_k must contain at least two increasing temperatures")
    for t in ts:
        number(t, "cold_temperature")
    if not 101325.0 <= number(config["pressure_pa"], "pressure_pa") <= 2e6:
        raise ValueError("initial screen pressure must be 101325..2000000 Pa")
    hot = config["hot_outlet_temperatures_k"]
    if not hot or hot != sorted(set(hot)):
        raise ValueError("hot temperatures must be nonempty, sorted and unique")
    for t in config["hot_outlet_temperatures_k"]:
        if number(t, "hot_temperature") <= ts[0] or t > 500.0:
            raise ValueError("hot temperatures must exceed inlet temperature and be <=500 K")
    for x in config["water_loss_fractions"]:
        if number(x, "water_loss_fraction", strict=False) >= 1:
            raise ValueError("water_loss_fraction must be <1")
    for k in ("replicates", "block_size"):
        if isinstance(config[k], bool) or not isinstance(config[k], int) or config[k] < 1:
            raise ValueError(f"{k} must be a positive integer")
    if not isinstance(config["random_seed"], str) or not config["random_seed"]:
        raise ValueError("random_seed must be nonempty")
    ids = set()
    for additive in config["additives"]:
        key = additive["id"]
        if not isinstance(key, str) or not key or key == "water" or key in ids:
            raise ValueError("invalid or duplicate additive id")
        ids.add(key)
        if key in EXPECTED_CAS and additive.get("cas") != EXPECTED_CAS[key]:
            raise ValueError("catalog identity/backend mismatch: " + key)
        if not isinstance(additive["enabled"], bool):
            raise ValueError("enabled must be boolean")
        if additive["kind"] not in {"miscible_liquid", "dissolved_solid", "surfactant", "suspension"}:
            raise ValueError("unknown formulation kind")
        fractions = additive["mass_fractions"]
        if additive["enabled"] and not fractions:
            raise ValueError("enabled additive has no mass fractions")
        if len(set(fractions)) != len(fractions):
            raise ValueError("duplicate concentrations")
        for w in fractions:
            if number(w, "additive_mass_fraction") >= 1:
                raise ValueError("water-rich fractions must be <1")
        if additive["kind"] == "suspension" and additive["enabled"]:
            raise ValueError("suspensions require a separate particle/porous-flow protocol; disabled in V5m-1")


def build_formulations(config: dict) -> list[dict]:
    validate_config(config)
    mass = config["batch_mass_g"]
    rows = [{"formulation_id": "water", "additive_id": "water", "name": "Water reference",
             "cas": "7732-18-5", "kind": "pure_reference", "additive_mass_fraction": 0.0,
             "water_mass_fraction": 1.0, "batch_mass_g": mass, "additive_mass_g": 0.0,
             "water_mass_g": mass, "physical_execution_authorized": False}]
    for a in sorted(config["additives"], key=lambda x: x["id"]):
        if not a["enabled"]:
            continue
        for w in sorted(a["mass_fractions"]):
            rows.append({"formulation_id": "v5m-" + digest([a["id"], w])[:12],
                "additive_id": a["id"], "name": a["name"], "cas": a["cas"], "kind": a["kind"],
                "additive_mass_fraction": float(w), "water_mass_fraction": 1.0 - w,
                "batch_mass_g": mass, "additive_mass_g": mass * w, "water_mass_g": mass * (1 - w),
                "physical_execution_authorized": False})
    return rows


def retained_additive_stress(formulation: dict, water_loss: float) -> dict:
    """Mass balance only; assumes additive retention, not actual volatility/solubility."""
    w = number(formulation["additive_mass_fraction"], "mass_fraction", strict=False)
    loss = number(water_loss, "water_loss", strict=False)
    if loss >= 1 or w >= 1:
        raise ValueError("fractions must be <1")
    residual = w + (1 - w) * (1 - loss)
    return {"formulation_id": formulation["formulation_id"], "initial_mass_fraction": w,
            "fraction_of_initial_water_removed": loss, "remaining_total_mass_fraction": residual,
            "retained_additive_mass_fraction": w / residual,
            "assumption": "all_additive_retained_water_only_removed",
            "actual_evaporation_predicted": False, "precipitation_predicted": False}


def cold_ratios(row: dict, water: dict) -> dict:
    result = {"formulation_id": row["formulation_id"], "temperature_k": row["temperature_k"],
              "pressure_pa": row["pressure_pa"], "basis": "same_geometry_single_phase_Darcy_only"}
    for name in ("density_kg_m3", "cp_j_kg_k", "viscosity_pa_s", "conductivity_w_m_k"):
        x, y = row.get(name), water.get(name)
        result[name + "_ratio_vs_water"] = x / y if x and y else None
    mu, rho = result["viscosity_pa_s_ratio_vs_water"], result["density_kg_m3_ratio_vs_water"]
    result["pressure_drop_ratio_equal_volume_flux"] = mu
    result["pressure_drop_ratio_equal_mass_flux"] = mu / rho if mu and rho else None
    result["boiling_or_clogging_model"] = False
    return result


def heat_only_comparison(delta_h: float, water_delta_h: float) -> dict:
    dh = number(delta_h, "delta_h")
    ref = number(water_delta_h, "water_delta_h")
    ratio = dh / ref
    return {"heat_uptake_ratio_vs_water": ratio, "mass_ratio_equal_net_heat": ref / dh,
            "net_heat_ratio_for_mass_parity": ratio,
            "additional_net_heat_reduction_needed_for_mass_parity": max(0.0, 1 - ratio),
            "real_surface_heat_reduction_predicted": False, "system_winner": False}


def make_bench_plan(formulations: list[dict], config: dict) -> list[dict]:
    candidates = sorted(r["formulation_id"] for r in formulations if r["formulation_id"] != "water")
    rows = []
    for replicate in range(1, config["replicates"] + 1):
        shuffled = candidates[:]
        random.Random(f"{config['random_seed']}:{replicate}").shuffle(shuffled)
        for start in range(0, len(shuffled), config["block_size"]):
            block = f"r{replicate:02d}-b{start // config['block_size'] + 1:02d}"
            group = [("water", "water_before")] + [(c, "candidate") for c in shuffled[start:start + config["block_size"]]] + [("water", "water_after")]
            for order, (fid, role) in enumerate(group, 1):
                rows.append({"run_id": f"{block}-{order:02d}", "formulation_id": fid,
                    "block_id": block, "replicate": replicate, "role": role, "order_in_block": order,
                    "physical_execution_authorized": False,
                    "prerequisites": "review;composition;stability;cold_flow;fresh_equivalent_coupon"})
    return rows


def protocol_template() -> dict:
    return {"version": "v5m-1-bench", "scope": "matched_coupon_tests_not_flight_validation",
            "apparatus_id": None, "geometry_id": None, "material_lot_id": None,
            "surface_spec_id": None, "load_program_id": None,
            "wall_temperature_limit_k": None, "max_feed_pressure_pa": None,
            "target_incident_energy_j": None, "target_duration_s": None,
            "max_hydraulic_resistance_ratio": None, "max_residue_mg": None,
            "energy_relative_tolerance": 0.03, "duration_relative_tolerance": 0.01,
            "max_control_mass_relative_difference": 0.10, "minimum_replicates": 3,
            "minimum_mass_saving_fraction": 0.05,
            "uncertainty_semantics": "declared_absolute_bounds_not_confidence_intervals"}


def write_json(path: Path, data: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(path.suffix + ".tmp")
    tmp.write_text(json.dumps(data, indent=2, sort_keys=True, allow_nan=False) + "\n", encoding="utf-8")
    tmp.replace(path)


def read_json(path: Path):
    def reject(token):
        raise ValueError("nonfinite JSON: " + token)
    return json.loads(path.read_text(encoding="utf-8"), parse_constant=reject)
