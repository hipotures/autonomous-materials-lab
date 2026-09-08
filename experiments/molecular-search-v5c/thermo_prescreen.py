"""Thermodynamic property pre-screen for V5c-2."""
from __future__ import annotations

import math
import statistics
from typing import Any


def _finite(value: Any) -> float | None:
    if not isinstance(value, (int, float)):
        return None
    value = float(value)
    return value if math.isfinite(value) else None


def _quantile(values: list[float], fraction: float) -> float | None:
    finite = sorted(float(value) for value in values if math.isfinite(float(value)))
    if not finite:
        return None
    if not 0.0 <= fraction <= 1.0:
        raise ValueError("quantile fraction must be in [0, 1]")
    if len(finite) == 1:
        return finite[0]
    position = fraction * (len(finite) - 1)
    low = int(math.floor(position))
    high = int(math.ceil(position))
    if low == high:
        return finite[low]
    weight = position - low
    return finite[low] * (1.0 - weight) + finite[high] * weight


def property_screen(
    provider,
    config: dict[str, Any],
) -> dict[str, Any]:
    """Evaluate a cheap thermodynamic grid before the entry trajectory."""
    storage = config["storage"]
    screen = config["prescreen"]
    storage_t = float(storage["temperature_k"])
    storage_p = float(storage["pressure_pa"])
    h_storage = float(provider.enthalpy_j_kg(storage_t, storage_p))
    storage_state = provider.state(storage_t, storage_p)
    storage_phase = provider.phase(storage_t, storage_p)

    temperatures = [
        float(value)
        for value in screen["temperature_grid_k"]
    ]
    pressures = [
        float(value)
        for value in screen["pressure_grid_pa"]
    ]
    if not temperatures or not pressures:
        raise ValueError("prescreen temperature/pressure grids must be non-empty")

    grid: list[dict[str, Any]] = []
    failures: list[dict[str, Any]] = []
    delta_h_values: list[float] = []
    cp_values: list[float] = []

    for temperature_k in temperatures:
        for pressure_pa in pressures:
            try:
                state = provider.state(temperature_k, pressure_pa)
                delta_h = float(state.enthalpy_j_kg) - h_storage
                cp = _finite(state.cp_j_kg_k)
                row = {
                    "temperature_k": temperature_k,
                    "pressure_pa": pressure_pa,
                    "phase": state.phase,
                    "delta_h_j_kg": delta_h,
                    "cp_j_kg_k": cp,
                    "density_kg_m3": state.density_kg_m3,
                }
                grid.append(row)
                if math.isfinite(delta_h):
                    delta_h_values.append(delta_h)
                if cp is not None and cp > 0.0:
                    cp_values.append(cp)
            except Exception as exc:
                failures.append(
                    {
                        "temperature_k": temperature_k,
                        "pressure_pa": pressure_pa,
                        "failure_reason": str(exc),
                    }
                )

    positive_delta_h = [
        value
        for value in delta_h_values
        if value > 0.0
    ]
    requested_count = len(temperatures) * len(pressures)
    positive_fraction = (
        len(positive_delta_h) / requested_count
        if requested_count
        else 0.0
    )

    saturation = provider.saturation_at_pressure(storage_p)
    boiling_temperature_k = (
        saturation.bubble_temperature_k
        if saturation.supported
        else None
    )
    latent_proxy_j_kg = None
    latent_offset = float(screen.get("latent_proxy_offset_k", 2.0))
    if (
        saturation.supported
        and isinstance(boiling_temperature_k, (int, float))
        and float(boiling_temperature_k) > latent_offset
    ):
        try:
            t_sat = float(boiling_temperature_k)
            h_below = provider.enthalpy_j_kg(
                t_sat - latent_offset,
                storage_p,
            )
            h_above = provider.enthalpy_j_kg(
                t_sat + latent_offset,
                storage_p,
            )
            latent_proxy_j_kg = float(h_above) - float(h_below)
            if (
                not math.isfinite(latent_proxy_j_kg)
                or latent_proxy_j_kg <= 0.0
            ):
                latent_proxy_j_kg = None
        except Exception:
            latent_proxy_j_kg = None

    return {
        "storage_phase": storage_phase,
        "storage_enthalpy_j_kg": h_storage,
        "storage_density_kg_m3": storage_state.density_kg_m3,
        "storage_cp_j_kg_k": storage_state.cp_j_kg_k,
        "boiling_temperature_k": boiling_temperature_k,
        "saturation_supported": bool(saturation.supported),
        "saturation_failure_reason": saturation.failure_reason,
        "latent_proxy_j_kg": latent_proxy_j_kg,
        "requested_grid_state_count": requested_count,
        "successful_grid_state_count": len(grid),
        "failed_grid_state_count": len(failures),
        "positive_delta_h_state_count": len(positive_delta_h),
        "positive_delta_h_fraction": positive_fraction,
        "delta_h_min_j_kg": min(delta_h_values) if delta_h_values else None,
        "delta_h_q25_j_kg": _quantile(delta_h_values, 0.25),
        "delta_h_median_j_kg": (
            statistics.median(delta_h_values)
            if delta_h_values
            else None
        ),
        "delta_h_max_j_kg": max(delta_h_values) if delta_h_values else None,
        "cp_median_j_kg_k": (
            statistics.median(cp_values)
            if cp_values
            else None
        ),
        "grid": grid,
        "failures": failures,
    }


def compare_to_reference(
    candidate: dict[str, Any],
    reference: dict[str, Any],
) -> dict[str, Any]:
    """Add water-normalized property targets used to prioritize full entry runs."""
    candidate_median = _finite(candidate.get("delta_h_median_j_kg"))
    reference_median = _finite(reference.get("delta_h_median_j_kg"))
    candidate_q25 = _finite(candidate.get("delta_h_q25_j_kg"))
    reference_q25 = _finite(reference.get("delta_h_q25_j_kg"))
    candidate_latent = _finite(candidate.get("latent_proxy_j_kg"))
    reference_latent = _finite(reference.get("latent_proxy_j_kg"))

    median_ratio = (
        candidate_median / reference_median
        if candidate_median is not None
        and reference_median is not None
        and reference_median > 0.0
        else None
    )
    q25_ratio = (
        candidate_q25 / reference_q25
        if candidate_q25 is not None
        and reference_q25 is not None
        and reference_q25 > 0.0
        else None
    )
    latent_ratio = (
        candidate_latent / reference_latent
        if candidate_latent is not None
        and reference_latent is not None
        and reference_latent > 0.0
        else None
    )

    # Primary target follows the entry mass physics: maximize robust usable
    # enthalpy. q25 is preferred over the median when available because it
    # penalizes pressure/temperature pockets with weak heat uptake.
    priority = q25_ratio if q25_ratio is not None else median_ratio

    return {
        "delta_h_median_ratio_vs_water": median_ratio,
        "delta_h_q25_ratio_vs_water": q25_ratio,
        "latent_proxy_ratio_vs_water": latent_ratio,
        "property_priority_score": priority,
    }


def passes_property_gate(
    metrics: dict[str, Any],
    config: dict[str, Any],
) -> tuple[bool, str | None]:
    screen = config["prescreen"]
    positive_fraction = float(metrics.get("positive_delta_h_fraction") or 0.0)
    if positive_fraction < float(screen["minimum_positive_delta_h_fraction"]):
        return False, "insufficient_positive_delta_h_fraction"

    merit = metrics.get("delta_h_q25_ratio_vs_water")
    if not isinstance(merit, (int, float)):
        merit = metrics.get("delta_h_median_ratio_vs_water")
    if not isinstance(merit, (int, float)):
        return False, "missing_enthalpy_merit"
    if float(merit) < float(screen["minimum_enthalpy_ratio_vs_water"]):
        return False, "enthalpy_merit_below_threshold"

    return True, None
