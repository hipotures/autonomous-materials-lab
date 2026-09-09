"""Shared pure-water controls for the existing Psat caloric model.

Controls are model-to-model checks, never independent mixture measurements.
The legacy mixture backend and its task signature deliberately remain unchanged.
"""
from __future__ import annotations

from copy import deepcopy
import math
import os
from pathlib import Path
import subprocess
import sys
import threading
import uuid

import campaign_backend as backend
from campaign_store import digest, encoded, implementation, read_json, write_json

PROPERTIES = ("VaporPressures", "HeatCapacityGases", "VolumeLiquids")
DEFAULT_POLICY = {"temperature_points": 21, "pressure_points": 5,
                  "consistency_relative_tolerance": 1e-7}


def coordinate(t, p):
    return tuple(float(f"{float(x):.12g}") for x in (t, p))


def context(frozen):
    """Project only dependencies of the pure-water limit; verify reuse per row."""
    grid = frozen["config"]["grid"]
    methods = frozen["selected_correlations"]["pure_methods"]
    return {"schema": "water-control-context-v1", "water_cas": "7732-18-5",
            "inlet_temperature_k": grid["inlet_temperature_k"],
            "selection_temperatures_k": sorted(set([298.15, grid["inlet_temperature_k"]]
                                                   + grid["outlet_temperatures_k"])),
            "pure_methods": {k: deepcopy(methods[k][0]) for k in PROPERTIES},
            "equilibrium_basis": "Psat", "caloric_basis": "Psat", "gas": "IdealGas",
            "reference_provider": "CoolProp HEOS Water"}


def control_grid(config):
    policy = {**DEFAULT_POLICY, **config.get("water_audit", {})}
    t0, t1 = config["scope"]["temperature_k"]
    p0, p1 = config["scope"]["pressure_pa"]
    ts = [t0 + i*(t1-t0)/(policy["temperature_points"]-1)
          for i in range(policy["temperature_points"])] + config["coarse"]["temperatures_k"]
    ps = [p0*math.exp(i*math.log(p1/p0)/(policy["pressure_points"]-1))
          for i in range(policy["pressure_points"])] + config["coarse"]["pressures_pa"]
    return sorted({coordinate(t, p) for t in ts for p in ps})


def calculate(request):
    """Evaluate the pure limit with the same bounded correlations and caloric basis."""
    from thermo import ChemicalConstantsPackage, GibbsExcessLiquid, IdealGas, FlashPureVLS
    import CoolProp.CoolProp as CP
    legacy = backend.legacy_import()
    ctx = request["context"]
    constants, props = ChemicalConstantsPackage.from_IDs([ctx["water_cas"]])
    for name in PROPERTIES:
        selected = legacy.pin_correlation(getattr(props, name)[0], ctx["selection_temperatures_k"])
        if selected != ctx["pure_methods"][name]:
            raise ValueError("water_control_correlation_mismatch:" + name)
    liquid = GibbsExcessLiquid(VaporPressures=props.VaporPressures,
        HeatCapacityGases=props.HeatCapacityGases, VolumeLiquids=props.VolumeLiquids,
        equilibrium_basis=ctx["equilibrium_basis"], caloric_basis=ctx["caloric_basis"])
    gas = IdealGas(HeatCapacityGases=props.HeatCapacityGases)
    flasher = FlashPureVLS(constants, props, gas=gas, liquids=[liquid], solids=[])
    factor = 1000.0/float(constants.MWs[0])
    tin = ctx["inlet_temperature_k"]
    output = []
    for t, p in request["coordinates"]:
        row = {"coordinate": [t, p], "inlet_temperature_k": tin,
               "evidence_kind": "reference_model_comparison_not_measurement"}
        try:
            for temp in (298.15, tin, t):
                for name in PROPERTIES:
                    legacy.in_bounds(getattr(props, name)[0], temp)
            phases, states, checks = [], [], []
            for temp in (tin, t):
                psat = float(props.VaporPressures[0](temp))
                heos_psat = float(CP.PropsSI("P", "T", temp, "Q", 0, "HEOS::Water"))
                if min(abs(p/psat-1), abs(p/heos_psat-1)) <= 1e-6:
                    raise ValueError("water_saturation_ambiguous_at_TP")
                state = flasher.flash(T=temp, P=p)
                phase = "gas" if state.gas is not None else "liquid"
                heos_phase = CP.PhaseSI("T", temp, "P", p, "HEOS::Water")
                states.append(float(state.H())*factor)
                phases.append((phase, heos_phase))
                checks.append({"temperature_k": temp, "model_phase": phase,
                               "heos_phase": heos_phase, "model_psat_pa": psat,
                               "heos_psat_pa": heos_psat})
            model_dh = states[1]-states[0]
            heos = legacy.water_reference(tin, t, p)
            if heos["status"] != "ok":
                raise ValueError("water_HEOS_unavailable:" + heos.get("error", "unknown"))
            if model_dh <= 0 or not math.isfinite(model_dh):
                raise ValueError("invalid_water_control_delta_h")
            row.update(status="ok" if all(a == b for a, b in phases)
                       and phases[0] == ("liquid", "liquid") else "water_phase_mismatch",
                       model_delta_h_j_kg=model_dh, heos_delta_h_j_kg=heos["delta_h_j_kg"],
                       endpoint_relative_error=abs(model_dh/heos["delta_h_j_kg"]-1),
                       endpoint_checks=checks)
            # Diagnostic errors do not erase a valid delta-h control.
            try:
                l = liquid.to(T=t, P=p, zs=[1.0])
                g = gas.to(T=t, P=p, zs=[1.0])
                row["caloric_diagnostics"] = {
                    "model_saturated_latent_heat_j_kg": (g.H()-l.H())*factor,
                    "heos_saturated_latent_heat_j_kg": float(CP.PropsSI("Hmass", "T", t, "Q", 1, "HEOS::Water")
                        - CP.PropsSI("Hmass", "T", t, "Q", 0, "HEOS::Water")),
                    "model_liquid_cp_j_kg_k": l.Cp()*factor,
                    "heos_saturated_liquid_cp_j_kg_k": float(CP.PropsSI("Cpmass", "T", t, "Q", 0, "HEOS::Water")),
                    "cp_comparison_basis": "model_liquid_at_P_vs_HEOS_saturated_liquid_diagnostic_only"}
            except (ValueError, TypeError, ArithmeticError, RuntimeError) as exc:
                row["caloric_diagnostics"] = {"status": "unavailable", "error": str(exc)}
        except (ValueError, TypeError, ArithmeticError, RuntimeError) as exc:
            row.update(status="water_control_unavailable", error=f"{type(exc).__name__}:{exc}")
        encoded(row)
        output.append(row)
    return {"rows": output}


class AuditDriver(backend.Driver):
    """Additional isolated worker, without invalidating legacy mixture states."""
    def water_states(self, ctx, coordinates):
        folder = self.work / ("water-" + uuid.uuid4().hex)
        folder.mkdir(parents=True)
        inp, out = folder/"request.json", folder/"response.json"
        write_json(inp, {"context": ctx, "coordinates": coordinates})
        env = {**os.environ, "OMP_NUM_THREADS": "1", "OPENBLAS_NUM_THREADS": "1", "MKL_NUM_THREADS": "1"}
        with (folder/"worker.log").open("w", encoding="utf-8") as log:
            try:
                proc = subprocess.run([sys.executable, str(Path(__file__).resolve()), str(inp), str(out)],
                    stdout=log, stderr=subprocess.STDOUT, env=env, timeout=self.timeout, check=False)
            except subprocess.TimeoutExpired as exc:
                raise RuntimeError("water_worker_timeout; log=" + str(folder/"worker.log")) from exc
        if proc.returncode or not out.is_file():
            raise RuntimeError("water_worker_failed; log=" + str(folder/"worker.log"))
        return read_json(out)["rows"]


class SharedWaterAudit:
    def __init__(self, store, driver, config):
        self.store, self.driver, self.config = store, driver, config
        self.code = {"water_control": implementation(Path(__file__)), "legacy_backend": driver.code}
        self.lock = threading.RLock()
        self.controls, self.contexts = {}, {}

    def ensure(self, ctx, coordinates):
        ctx_id = digest(ctx)
        coordinates = sorted({coordinate(*q) for q in coordinates})
        # Single flight across pairs and models, including --recompute.
        with self.lock:
            self.contexts[ctx_id] = ctx
            pending = []
            for q in coordinates:
                key = (ctx_id, q)
                if key in self.controls:
                    continue
                spec = self.store.spec("reference.water_state", {"context": ctx, "coordinate": q},
                                       self.code, environment=self.driver.environment)
                saved = self.store.load(spec)
                if saved is None:
                    pending.append((q, spec))
                else:
                    self.controls[key] = saved
            size = self.config["execution"]["batch_points"]
            for start in range(0, len(pending), size):
                batch = pending[start:start+size]
                try:
                    rows = self.driver.water_states(ctx, [list(q) for q, _ in batch])
                    by_q = {coordinate(*r["coordinate"]): r for r in rows}
                    if len(rows) != len(batch) or set(by_q) != {q for q, _ in batch}:
                        raise ValueError("missing_duplicate_or_mismatched_water_controls")
                    encoded(rows)
                    for q, spec in batch:
                        self.controls[(ctx_id, q)] = self.store.put(spec, by_q[q])
                except Exception as exc:
                    for q, spec in batch:
                        self.controls[(ctx_id, q)] = self.store.put(spec, error=f"{type(exc).__name__}:{exc}")
            if pending:
                print(f"[water] context={ctx_id[:12]} new_controls={len(pending)} requested={len(coordinates)}", flush=True)
        return ctx_id

    def bundle(self, frozen, coordinates, *, preflight=False):
        qs = sorted({coordinate(*q) for q in coordinates} | (set(control_grid(self.config)) if preflight else set()))
        deps, model_contexts, errors = {}, {}, {}
        for model, artifact in sorted(frozen.items()):
            try:
                ctx = context(artifact.data)
            except (ValueError, KeyError, TypeError) as exc:
                errors[model] = f"water_context_unavailable:{exc}"
                continue
            # Cache-integrity errors must propagate, not become scientific gaps.
            ctx_id = self.ensure(ctx, qs)
            model_contexts[model] = ctx_id
            for q in qs:
                a = self.controls[(ctx_id, q)]
                deps[a.key] = a
        with self.lock:
            controls = [{"context_id": cid, "coordinate": list(q), "task": a.dependency(),
                         "outcome": a.outcome}
                        for (cid, q), a in sorted(self.controls.items())
                        if cid in model_contexts.values() and q in qs]
        return self.store.run("reference.water_bundle", {"model_contexts": model_contexts,
            "context_errors": errors, "coordinates": qs}, self.code,
            lambda: {"model_contexts": model_contexts, "context_errors": errors, "controls": controls},
            dependencies=deps)

    def export(self):
        return {"schema": "mixture-campaign-water-controls-v1", "contexts": deepcopy(self.contexts),
                "controls": [{"context_id": cid, "coordinate": list(q), "task": a.dependency(),
                              "outcome": a.outcome} for (cid, q), a in sorted(self.controls.items())],
                "evidence_kind": "reference_model_comparison_not_measurement"}


def annotate(rows, controls, config):
    """Gate comparisons without changing raw cached mixture solver outputs."""
    policy = {**DEFAULT_POLICY, **config.get("water_audit", {})}
    tolerance = config["model"]["endpoint_relative_tolerance"]
    index = {(c["context_id"], coordinate(*c["coordinate"])): c for c in controls["controls"]}
    output = []
    for raw in rows:
        row = deepcopy(raw)
        cid = controls["model_contexts"].get(row["model"])
        q = coordinate(row["outlet_temperature_k"], row["pressure_pa"])
        control = index.get((cid, q))
        data = (control or {}).get("outcome", {}).get("data") or {}
        reasons = []
        if data.get("status") != "ok":
            reasons.append(data.get("status", "water_control_unavailable"))
        error = data.get("endpoint_relative_error")
        if isinstance(error, bool) or not isinstance(error, (int, float)) or not math.isfinite(error) or not 0 <= error <= tolerance:
            reasons.append("water_reference_tolerance_exceeded" if error is not None else "water_reference_error_unknown")
        if row.get("status") == "ok":
            for field, expected in (("same_model_water_delta_h_j_kg", data.get("model_delta_h_j_kg")),
                                    ("heos_water_delta_h_j_kg", data.get("heos_delta_h_j_kg"))):
                actual = row.get(field)
                if isinstance(actual, bool) or isinstance(expected, bool) or not isinstance(actual, (int, float)) or not isinstance(expected, (int, float)) or not math.isfinite(actual) or not math.isfinite(expected) or actual <= 0 or expected <= 0 or not math.isclose(
                        actual, expected, rel_tol=policy["consistency_relative_tolerance"], abs_tol=1e-6):
                    reasons.append("water_control_context_mismatch")
                    break
        row["raw_model_comparison_eligible"] = raw.get("model_comparison_eligible") is True
        row["comparison_blockers"] = sorted(set(reasons))
        row["water_context_id"] = cid
        row["water_control_task"] = control.get("task") if control else None
        row["model_comparison_eligible"] = row["raw_model_comparison_eligible"] and not reasons
        output.append(row)
    return output


if __name__ == "__main__":
    if len(sys.argv) != 3:
        raise SystemExit("internal water worker requires input and output paths")
    write_json(Path(sys.argv[2]), calculate(read_json(Path(sys.argv[1]))))
