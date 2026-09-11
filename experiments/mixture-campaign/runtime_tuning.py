"""Runtime-only scheduling heuristics for mixture campaigns.

These choices affect concurrency and batching only. They are deliberately kept
out of numerical task signatures so changing worker counts does not invalidate
scientific cache entries.
"""
from __future__ import annotations

import os
from typing import Mapping, Sequence

MIB = 1024 ** 2
GIB = 1024 ** 3
DEFAULT_MAX_WORKERS = 32
DEFAULT_WORKER_MEMORY_MIB = 768


def available_cpu_count() -> int:
    """Return CPUs available to this process, respecting Linux affinity."""
    try:
        return max(1, len(os.sched_getaffinity(0)))
    except (AttributeError, OSError):
        return max(1, os.cpu_count() or 1)


def available_memory_bytes() -> int | None:
    """Best-effort estimate of currently available physical memory."""
    try:
        pages = int(os.sysconf("SC_AVPHYS_PAGES"))
        page_size = int(os.sysconf("SC_PAGE_SIZE"))
        value = pages * page_size
        return value if value > 0 else None
    except (AttributeError, OSError, ValueError, TypeError):
        return None


def _positive_int(value: str | int | None, name: str) -> int | None:
    if value is None:
        return None
    try:
        number = int(value)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"{name} must be a positive integer") from exc
    if number < 1:
        raise ValueError(f"{name} must be a positive integer")
    return number


def recommend_workers(*, cpus: int | None = None, memory_bytes: int | None = None,
                      max_workers: int | None = None, worker_memory_mib: int | None = None,
                      environ: Mapping[str, str] | None = None) -> int:
    """Choose conservative process-level concurrency for single-threaded workers."""
    env = os.environ if environ is None else environ
    cpus = available_cpu_count() if cpus is None else max(1, int(cpus))
    memory_bytes = available_memory_bytes() if memory_bytes is None else memory_bytes
    max_workers = (_positive_int(env.get("AML_MAX_WORKERS"), "AML_MAX_WORKERS")
                   if max_workers is None else _positive_int(max_workers, "max_workers"))
    if max_workers is None:
        max_workers = DEFAULT_MAX_WORKERS
    worker_memory_mib = (_positive_int(env.get("AML_WORKER_MEMORY_MIB"), "AML_WORKER_MEMORY_MIB")
                         if worker_memory_mib is None else _positive_int(worker_memory_mib, "worker_memory_mib"))
    if worker_memory_mib is None:
        worker_memory_mib = DEFAULT_WORKER_MEMORY_MIB

    reserve_cpus = 2 if cpus >= 8 else 1 if cpus >= 3 else 0
    cpu_budget = max(1, cpus - reserve_cpus)
    memory_budget = max_workers
    if memory_bytes is not None and memory_bytes > 0:
        reserve = min(4 * GIB, max(GIB, memory_bytes // 8))
        usable = max(0, memory_bytes - reserve)
        memory_budget = max(1, usable // (worker_memory_mib * MIB))
    return max(1, min(cpu_budget, int(memory_budget), max_workers))


def recommend_batch_pairs(workers: int) -> int:
    """Checkpoint batch size, not scientific search scope."""
    workers = _positive_int(workers, "workers") or 1
    return max(50, min(1024, workers * 16))


def _has_option(argv: Sequence[str], option: str) -> bool:
    return option in argv or any(value.startswith(option + "=") for value in argv)


def tune_mission_argv(argv: Sequence[str], *, cpus: int | None = None,
                      memory_bytes: int | None = None,
                      environ: Mapping[str, str] | None = None) -> tuple[list[str], dict]:
    """Inject execution-only defaults while preserving explicit CLI choices."""
    env = os.environ if environ is None else environ
    result = list(argv)
    injected = {}

    if not _has_option(result, "--workers"):
        explicit = _positive_int(env.get("AML_WORKERS"), "AML_WORKERS")
        workers = explicit or recommend_workers(cpus=cpus, memory_bytes=memory_bytes, environ=env)
        result += ["--workers", str(workers)]
        injected["workers"] = workers
    else:
        workers = None

    if not _has_option(result, "--batch-pairs") and not _has_option(result, "--limit-pairs"):
        explicit_batch = _positive_int(env.get("AML_BATCH_PAIRS"), "AML_BATCH_PAIRS")
        if explicit_batch is not None:
            batch_pairs = explicit_batch
        else:
            if workers is None:
                # Explicit CLI worker count is intentionally not reparsed here; a stable
                # moderate batch still removes the old 50-pair synchronization barrier.
                batch_pairs = 256
            else:
                batch_pairs = recommend_batch_pairs(workers)
        result += ["--batch-pairs", str(batch_pairs)]
        injected["batch_pairs"] = batch_pairs

    info = {
        "available_cpus": available_cpu_count() if cpus is None else int(cpus),
        "available_memory_bytes": available_memory_bytes() if memory_bytes is None else memory_bytes,
        "injected": injected,
    }
    return result, info
