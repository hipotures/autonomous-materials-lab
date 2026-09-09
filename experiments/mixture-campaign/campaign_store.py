"""Content-addressed task graph, immutable outcomes and an append-only run journal."""
from __future__ import annotations

from contextlib import contextmanager
from dataclasses import dataclass
import gzip
import hashlib
import json
import os
from pathlib import Path
import threading
import time
import uuid
from typing import Any, Callable

CACHE_ABI = "mixture-campaign-task-v1"


def encoded(value: Any) -> bytes:
    return (json.dumps(value, sort_keys=True, separators=(",", ":"), allow_nan=False) + "\n").encode()


def digest(value: Any) -> str:
    return hashlib.sha256(encoded(value)).hexdigest()


def file_hash(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def read_json(path: Path) -> Any:
    def reject(value):
        raise ValueError("nonfinite JSON: " + value)
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt", encoding="utf-8") as f:
        return json.load(f, parse_constant=reject)


def atomic_bytes(path: Path, data: bytes) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temp = path.with_name(path.name + "." + uuid.uuid4().hex + ".tmp")
    try:
        with temp.open("xb") as f:
            f.write(data); f.flush(); os.fsync(f.fileno())
        os.replace(temp, path)
    finally:
        temp.unlink(missing_ok=True)


def write_json(path: Path, value: Any) -> None:
    atomic_bytes(path, encoded(value))


def implementation(*paths: Path) -> dict[str, str]:
    return {str(p.name): file_hash(p) for p in paths}


@dataclass(frozen=True)
class Artifact:
    key: str
    result_hash: str
    outcome: dict
    spec: dict

    @property
    def data(self) -> Any:
        if self.outcome["status"] != "complete":
            raise ValueError("failed dependency: " + self.key)
        return self.outcome["data"]

    def dependency(self) -> dict:
        return {"task_key": self.key, "result_hash": self.result_hash}


class Store:
    """One externally locked campaign may use this store from several threads.

    Refs are mutable cache hints. Blobs, attempts and published snapshots are not.
    Downstream keys contain BOTH dependency task keys and outcome hashes.
    """
    def __init__(self, root: Path, journal: Path, *, retry_failures: bool = False,
                 recompute: bool = False):
        self.root, self.journal = root, journal
        self.retry_failures, self.recompute = retry_failures, recompute
        self._lock = threading.RLock()
        self.used: dict[str, Artifact] = {}
        self.events: list[dict] = []

    def spec(self, kind: str, inputs: Any, code: dict, *, dependencies: dict[str, Artifact] | None = None,
             environment: dict | None = None) -> dict:
        return json.loads(encoded({"abi": CACHE_ABI, "kind": kind, "inputs": inputs,
                "implementation": code, "store_implementation": file_hash(Path(__file__)),
                "environment": environment or {},
                "dependencies": {k: v.dependency() for k, v in sorted((dependencies or {}).items())}}))

    def event(self, event: dict) -> None:
        with self._lock:
            self.events.append(event)
            self.journal.parent.mkdir(parents=True, exist_ok=True)
            with self.journal.open("ab") as f:
                f.write(encoded(event))

    def load(self, spec: dict) -> Artifact | None:
        key = digest(spec)
        with self._lock:
            if key in self.used:
                return self.used[key]
            if self.recompute:
                return None
            ref = self.root / "refs" / (key + ".json")
            if not ref.exists():
                return None
            pointer = read_json(ref)
            sha = pointer.get("result_hash", "")
            if not isinstance(sha, str) or len(sha) != 64 or any(c not in "0123456789abcdef" for c in sha):
                raise ValueError("corrupt cache reference: " + key)
            payload = read_json(self.root / "objects" / (sha + ".json.gz"))
            if digest(payload) != sha or payload.get("spec") != spec:
                raise ValueError("cache integrity failure: " + key)
            if payload["outcome"]["status"] == "failed" and self.retry_failures:
                return None
            artifact = Artifact(key, sha, payload["outcome"], spec)
            self.used[key] = artifact
            self.event({"task_key": key, "kind": spec["kind"], "action": "cache_hit",
                        "result_hash": sha, "status": artifact.outcome["status"]})
            return artifact

    def put(self, spec: dict, data: Any = None, *, error: str | None = None) -> Artifact:
        key = digest(spec)
        outcome = {"status": "failed", "error": error} if error is not None else {"status": "complete", "data": data}
        payload = json.loads(encoded({"spec": spec, "outcome": outcome}))
        outcome = payload["outcome"]
        sha = digest(payload)
        artifact = Artifact(key, sha, outcome, spec)
        with self._lock:
            target = self.root / "objects" / (sha + ".json.gz")
            if not target.exists():
                atomic_bytes(target, gzip.compress(encoded(payload), compresslevel=6, mtime=0))
            elif read_json(target) != payload:
                raise ValueError("corrupt immutable cache object: " + sha)
            old_ref = self.root / "refs" / (key + ".json")
            old = read_json(old_ref) if old_ref.exists() else None
            changed = bool(old and old.get("result_hash") != sha)
            write_json(old_ref, {"task_key": key, "result_hash": sha})
            self.used[key] = artifact
            self.event({"task_key": key, "kind": spec["kind"], "action": "executed",
                        "result_hash": sha, "status": outcome["status"],
                        "changed_outcome_for_same_task": changed,
                        "previous_result_hash": old.get("result_hash") if old else None})
        return artifact

    def run(self, kind: str, inputs: Any, code: dict, fn: Callable[[], Any], *,
            dependencies: dict[str, Artifact] | None = None, environment: dict | None = None) -> Artifact:
        spec = self.spec(kind, inputs, code, dependencies=dependencies, environment=environment)
        saved = self.load(spec)
        if saved is not None:
            return saved
        # Exceptions remain explicit task failures; interruption is NOT swallowed.
        try:
            value = fn()
            encoded(value)
        except Exception as exc:
            return self.put(spec, error=f"{type(exc).__name__}: {exc}")
        return self.put(spec, value)

    def graph(self) -> list[dict]:
        return [{"task_key": a.key, "result_hash": a.result_hash, "spec": a.spec,
                 "status": a.outcome["status"]} for a in sorted(self.used.values(), key=lambda a: a.key)]


@contextmanager
def exclusive_lock(root: Path):
    import fcntl
    root.mkdir(parents=True, exist_ok=True)
    with (root / ".campaign.lock").open("a") as f:
        try:
            fcntl.flock(f, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError as exc:
            raise RuntimeError("another campaign owns this cache/output lock") from exc
        yield
