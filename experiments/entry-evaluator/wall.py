from __future__ import annotations

from dataclasses import dataclass

SIGMA = 5.670374419e-8


@dataclass
class WallState:
    temperature_k: float


@dataclass(frozen=True)
class WallStep:
    radiation_out_w_m2: float
    coolant_required_w_m2: float
    net_after_cooling_w_m2: float
    next_temperature_k: float


def step_wall(
    wall: WallState,
    external_heat_flux_w_m2: float,
    environment_temperature_k: float,
    dt: float,
    config: dict,
) -> WallStep:
    emissivity = float(config.get("emissivity", 0.85))
    areal_heat_capacity = float(config["areal_heat_capacity_j_m2_k"])
    setpoint = float(config["temperature_setpoint_k"])
    backface_loss = float(config.get("backface_loss_w_m2", 0.0))

    radiation_out = emissivity * SIGMA * max(
        0.0,
        wall.temperature_k**4 - environment_temperature_k**4,
    )
    no_cooling_net = external_heat_flux_w_m2 - radiation_out - backface_loss

    allowable_storage_rate = (
        areal_heat_capacity
        * max(0.0, setpoint - wall.temperature_k)
        / max(dt, 1e-9)
    )
    coolant_required = max(0.0, no_cooling_net - allowable_storage_rate)
    net_after = no_cooling_net - coolant_required
    next_temperature = wall.temperature_k + dt * net_after / areal_heat_capacity

    return WallStep(
        radiation_out_w_m2=radiation_out,
        coolant_required_w_m2=coolant_required,
        net_after_cooling_w_m2=net_after,
        next_temperature_k=next_temperature,
    )
