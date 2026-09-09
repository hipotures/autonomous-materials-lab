"""Campaign input validation and reproducible basic design (no property calls)."""
from __future__ import annotations
import math
from campaign_design import point, AXES


def numeric(value, name, lo=None, hi=None):
    if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value):
        raise ValueError("finite number required: " + name)
    if lo is not None and value < lo or hi is not None and value > hi:
        raise ValueError("outside allowed range: " + name)
    return float(value)


def validate(c):
    if not isinstance(c, dict) or c.get("schema") != "mixture-campaign-v1":
        raise ValueError("expected mixture-campaign-v1")
    if c["models"] != ["chemsep_nrtl", "unifac_dortmund"]:
        raise ValueError("current adapter requires the two declared binary models")
    for axis in AXES:
        values = c["scope"][axis]
        if not isinstance(values, list) or len(values) != 2:
            raise ValueError("scope requires two bounds: " + axis)
        lo, hi = [numeric(x, axis) for x in values]
        if lo >= hi: raise ValueError("scope bounds must increase")
    wlo, whi = c["scope"]["mass_fraction"]
    if not 0 < wlo < c["coarse"]["log_transition"] < whi < 1:
        raise ValueError("mass-fraction grid must be strictly inside (0,1)")
    tin = numeric(c["model"]["inlet_temperature_k"], "inlet", 273.15)
    tlo, thi = c["scope"]["temperature_k"]
    plo, phi = c["scope"]["pressure_pa"]
    if not tin < tlo < thi <= 500 or not 25000 <= plo < phi <= 300000:
        raise ValueError("outside declared legacy-model screening bounds")
    for key, axis in (("temperatures_k", "temperature_k"), ("pressures_pa", "pressure_pa")):
        values = c["coarse"][key]
        for v in values: numeric(v, key, *c["scope"][axis])
        if values != sorted(set(values)) or len(values) < 2 or [values[0], values[-1]] != c["scope"][axis]:
            raise ValueError("coarse axes must be sorted, unique and include scope endpoints")
    temps = c["model"]["correlation_temperatures_k"]
    if not temps or temps != sorted(set(temps)):
        raise ValueError("fixed correlation temperatures must be sorted and unique")
    for t in temps: numeric(t, "correlation T", tlo, thi)
    for section, names in (("coarse", ["log_points", "linear_points"]),
                            ("execution", ["workers", "batch_points"]),
                            ("adaptive", ["max_rounds", "max_new_points_per_pair", "points_per_round"])):
        for k in names:
            v = c[section][k]
            if isinstance(v, bool) or not isinstance(v, int) or v < (2 if section == "coarse" else 1):
                raise ValueError("positive integer required: " + k)
    numeric(c["execution"]["worker_timeout_s"], "timeout", 0.01)
    numeric(c["publication"]["shard_bytes"], "shard bytes", 2048, 49152)
    for k in ("gain_margin", "model_spread_trigger", "gain_change_trigger"):
        numeric(c["adaptive"][k], k, 0, 1)
    for k in AXES: numeric(c["adaptive"]["tolerances"][k], k, 1e-12)
    for k in ("balance_tolerance", "phase_fraction_tolerance"):
        numeric(c["model"][k], k, 1e-14, 1e-4)
    numeric(c["model"]["endpoint_relative_tolerance"], "endpoint tolerance", 0, .2)
    if not isinstance(c["studies"], list) or not c["studies"] or len(set(c["studies"])) != len(c["studies"]):
        raise ValueError("required studies must be unique and nonempty")
    if not isinstance(c["additional_cas"], list) or any(not isinstance(v, str) for v in c["additional_cas"]):
        raise ValueError("additional_cas must be a list of identifiers")
    if not isinstance(c["plugins"], list): raise ValueError("plugins must be an explicit list")
    c.setdefault("water_audit", {})
    c.setdefault("evidence", {})
    for key, default in (("temperature_points", 21), ("pressure_points", 5)):
        value = c["water_audit"].setdefault(key, default)
        if isinstance(value, bool) or not isinstance(value, int) or not 2 <= value <= 101:
            raise ValueError("water audit axis requires 2..101 points: " + key)
    numeric(c["water_audit"].setdefault("consistency_relative_tolerance", 1e-7), "water consistency", 1e-12, 1e-3)
    n = c["evidence"].setdefault("max_points_per_pair", 8)
    if isinstance(n, bool) or not isinstance(n, int) or not 1 <= n <= 100:
        raise ValueError("evidence max_points_per_pair requires 1..100")
    local = c["adaptive"].setdefault("disagreement_resolution", {"mass_fraction": .01, "temperature_k": 5., "pressure_pa": .08})
    for axis in AXES:
        numeric(local[axis], "disagreement resolution " + axis, c["adaptive"]["tolerances"][axis])
    numeric(c["adaptive"].setdefault("disagreement_stability_tolerance", .01), "disagreement stability", 0, 1)
    return c


def basic_points(c):
    lo, hi = c["scope"]["mass_fraction"]
    grid = c["coarse"]
    mid = grid["log_transition"]
    ws = [math.exp(math.log(lo)+i*math.log(mid/lo)/(grid["log_points"]-1)) for i in range(grid["log_points"])]
    ws += [mid+i*(hi-mid)/(grid["linear_points"]-1) for i in range(grid["linear_points"])]
    return sorted({point(w, t, p) for w in ws for t in grid["temperatures_k"] for p in grid["pressures_pa"]})
