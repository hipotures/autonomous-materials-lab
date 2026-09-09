"""Isolated adapter for the existing binary-water models; no new fitted physics."""
from __future__ import annotations

from importlib import import_module
from importlib.metadata import version
import json
import os
from pathlib import Path
import platform
import subprocess
import sys
import uuid

from campaign_store import implementation, digest, file_hash, read_json, write_json, encoded

HERE = Path(__file__).resolve().parent
LEGACY = HERE.parent / "formulation-screen-v5m"
PINNED = {"thermo": "0.6.1", "chemicals": "1.5.2", "CoolProp": "8.0.0", "rdkit": "2026.3.6"}
MODELS = ["chemsep_nrtl", "unifac_dortmund"]


def legacy_import():
    if str(LEGACY) not in sys.path:
        sys.path.insert(0, str(LEGACY))
    return import_module("v5m2_backend")


def code_signature() -> dict:
    return implementation(Path(__file__), LEGACY / "v5m2_backend.py", LEGACY / "v5m2_core.py")


def runtime() -> dict:
    packages = {name: version(name) for name in PINNED}
    for name, required in PINNED.items():
        if packages[name] != required:
            raise ValueError(f"requires existing project environment: {name}=={required}; found {packages[name]}")
    libraries = {}
    for name in ("thermo", "chemicals", "CoolProp", "rdkit", "numpy", "scipy", "fluids", "pandas"):
        packages[name] = version(name)
        root = Path(import_module(name).__file__).resolve().parent
        files = {str(p.relative_to(root)): file_hash(p) for p in sorted(root.rglob("*"))
                 if p.is_file() and "__pycache__" not in p.parts and p.suffix not in {".pyc", ".pyo"}}
        libraries[name] = digest(files)
    return {"python": platform.python_version(), "machine": platform.machine(),
            "system": platform.system(), "packages": packages, "library_tree_sha256": libraries}


def numerical_config(config: dict) -> dict:
    m = config["model"]
    # Fixed for the whole campaign. Adaptive requests NEVER influence selection.
    temps = sorted(set(m["correlation_temperatures_k"] + config["scope"]["temperature_k"]))
    return {"grid": {"inlet_temperature_k": m["inlet_temperature_k"], "outlet_temperatures_k": temps},
            **{k: m[k] for k in ("phase_fraction_tolerance", "balance_tolerance", "endpoint_relative_tolerance")}}


def model_pair(pair: dict, model: str) -> dict:
    # Pair-specific projection: discovering another compound cannot invalidate this one.
    return {**{k: pair[k] for k in ("pair_id", "cas_number", "name", "smiles", "inchi_key")},
            "models": {model: pair["models"][model]} if model in pair["models"] else {}}


def metadata(backend) -> dict:
    return {"pure_methods": backend.methods, "viscosity_methods": backend.viscosity_metadata}


def worker(request: dict) -> dict:
    legacy = legacy_import()
    if request["operation"] == "catalog":
        return legacy.discover_catalog(request["config"], additional_cas=request["additional_cas"])
    if request["operation"] == "freeze":
        obj = legacy.BinaryModel(request["pair"], request["model"], request["config"])
        return {"pair": request["pair"], "model": request["model"], "config": request["config"],
                "selected_correlations": metadata(obj), "evidence_kind": "unvalidated_model_prediction"}
    if request["operation"] != "states":
        raise ValueError("unknown numerical operation")
    frozen = request["frozen"]
    obj = legacy.BinaryModel(frozen["pair"], frozen["model"], frozen["config"])
    if metadata(obj) != frozen["selected_correlations"]:
        raise ValueError("frozen correlations changed; refusing to mix numerical models")
    rows = []
    for w, t, p in request["points"]:
        row = obj.evaluate(w, t, p)
        rows.append(row)
    return {"rows": rows}


class Driver:
    """A subprocess per bounded batch; wall-clock timeout really kills the child."""
    def __init__(self, work: Path, timeout: float):
        self.work, self.timeout = work, timeout
        self.code = code_signature()
        self.environment = runtime()

    def call(self, request: dict) -> dict:
        folder = self.work / uuid.uuid4().hex
        folder.mkdir(parents=True)
        inp, out = folder / "request.json", folder / "response.json"
        write_json(inp, request)
        env = {**os.environ, "OMP_NUM_THREADS": "1", "OPENBLAS_NUM_THREADS": "1", "MKL_NUM_THREADS": "1"}
        with (folder / "worker.log").open("w", encoding="utf-8") as log:
            try:
                proc = subprocess.run([sys.executable, str(Path(__file__).resolve()), str(inp), str(out)],
                                      stdout=log, stderr=subprocess.STDOUT, env=env,
                                      timeout=self.timeout, check=False)
            except subprocess.TimeoutExpired as exc:
                raise RuntimeError(f"worker_timeout; log={folder / 'worker.log'}") from exc
        if proc.returncode != 0 or not out.is_file():
            raise RuntimeError(f"worker_failed exit={proc.returncode}; log={folder / 'worker.log'}")
        return read_json(out)

    def catalog(self, config: dict, historical: list[str]) -> dict:
        return self.call({"operation": "catalog", "config": {"catalog": config["catalog"],
                    "grid": {"inlet_temperature_k": config["model"]["inlet_temperature_k"]}},
                    "additional_cas": sorted(set(historical + config["additional_cas"]))})

    def freeze(self, pair: dict, model: str, config: dict) -> dict:
        return self.call({"operation": "freeze", "pair": model_pair(pair, model),
                          "model": model, "config": numerical_config(config)})

    def states(self, frozen: dict, points: list[list[float]]) -> list[dict]:
        output = self.call({"operation": "states", "frozen": frozen, "points": points})
        rows = output["rows"]
        expected = {tuple(p) for p in points}
        found = [(r["additive_mass_fraction"], r["outlet_temperature_k"], r["pressure_pa"]) for r in rows]
        if len(found) != len(expected) or set(found) != expected or any(
                r["model"] != frozen["model"] or r["pair_id"] != frozen["pair"]["pair_id"] for r in rows):
            raise ValueError("worker returned missing, duplicate or mismatched states")
        encoded(rows)
        return sorted(rows, key=lambda r: (r["additive_mass_fraction"], r["outlet_temperature_k"], r["pressure_pa"]))


if __name__ == "__main__":
    if len(sys.argv) != 3:
        raise SystemExit("internal worker requires input and output files")
    write_json(Path(sys.argv[2]), worker(read_json(Path(sys.argv[1]))))
