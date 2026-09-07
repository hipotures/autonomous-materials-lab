from __future__ import annotations

from dataclasses import dataclass
from math import cos, sin
from typing import Callable

from atmosphere import AtmosphereState

EARTH_RADIUS_M = 6_371_000.0
EARTH_MU_M3_S2 = 3.986004418e14
STANDARD_G = 9.80665


@dataclass(frozen=True)
class TrajectoryState:
    altitude_m: float
    velocity_m_s: float
    flight_path_angle_rad: float
    downrange_m: float


@dataclass(frozen=True)
class AeroState:
    dynamic_pressure_pa: float
    drag_n: float
    lift_n: float
    surface_pressure_pa: float
    acceleration_g: float


def aero_state(
    state: TrajectoryState,
    atmosphere: AtmosphereState,
    mass_kg: float,
    vehicle: dict,
) -> AeroState:
    q = 0.5 * atmosphere.density_kg_m3 * state.velocity_m_s**2
    area = float(vehicle["reference_area_m2"])
    cd = float(vehicle["drag_coefficient"])
    lift_to_drag = float(vehicle.get("lift_to_drag", 0.0))
    drag = q * cd * area
    lift = drag * lift_to_drag
    cp_stag = float(vehicle.get("stagnation_pressure_coefficient", 2.0))
    surface_pressure = atmosphere.pressure_pa + cp_stag * q
    accel = (drag / max(mass_kg, 1.0)) / STANDARD_G
    return AeroState(q, drag, lift, surface_pressure, accel)


def derivatives(
    state: TrajectoryState,
    mass_kg: float,
    vehicle: dict,
    atmosphere_fn: Callable[[float], AtmosphereState],
) -> tuple[float, float, float, float]:
    atm = atmosphere_fn(state.altitude_m)
    aero = aero_state(state, atm, mass_kg, vehicle)
    r = EARTH_RADIUS_M + state.altitude_m
    g = EARTH_MU_M3_S2 / r**2
    v = max(state.velocity_m_s, 1.0)
    gamma = state.flight_path_angle_rad
    bank = float(vehicle.get("bank_angle_rad", 0.0))

    dh = v * sin(gamma)
    dv = -(aero.drag_n / mass_kg) - g * sin(gamma)
    dgamma = (
        aero.lift_n * cos(bank) / (mass_kg * v)
        + (v / r - g / v) * cos(gamma)
    )
    ds = v * cos(gamma) * EARTH_RADIUS_M / r
    return dh, dv, dgamma, ds


def rk4_step(
    state: TrajectoryState,
    dt: float,
    mass_kg: float,
    vehicle: dict,
    atmosphere_fn: Callable[[float], AtmosphereState],
) -> TrajectoryState:
    def add(
        s: TrajectoryState,
        k: tuple[float, float, float, float],
        scale: float,
    ) -> TrajectoryState:
        return TrajectoryState(
            s.altitude_m + scale * k[0],
            max(0.0, s.velocity_m_s + scale * k[1]),
            s.flight_path_angle_rad + scale * k[2],
            s.downrange_m + scale * k[3],
        )

    k1 = derivatives(state, mass_kg, vehicle, atmosphere_fn)
    k2 = derivatives(add(state, k1, 0.5 * dt), mass_kg, vehicle, atmosphere_fn)
    k3 = derivatives(add(state, k2, 0.5 * dt), mass_kg, vehicle, atmosphere_fn)
    k4 = derivatives(add(state, k3, dt), mass_kg, vehicle, atmosphere_fn)

    return TrajectoryState(
        altitude_m=state.altitude_m
        + dt * (k1[0] + 2 * k2[0] + 2 * k3[0] + k4[0]) / 6.0,
        velocity_m_s=max(
            0.0,
            state.velocity_m_s
            + dt * (k1[1] + 2 * k2[1] + 2 * k3[1] + k4[1]) / 6.0,
        ),
        flight_path_angle_rad=state.flight_path_angle_rad
        + dt * (k1[2] + 2 * k2[2] + 2 * k3[2] + k4[2]) / 6.0,
        downrange_m=state.downrange_m
        + dt * (k1[3] + 2 * k2[3] + 2 * k3[3] + k4[3]) / 6.0,
    )
