from __future__ import annotations

from dataclasses import dataclass
import math

from atmosphere import AtmosphereState

SUTTON_GRAVES_EARTH = 1.74153e-4

TAUBER_SUTTON_FV = [
    (0.0, 0.0),
    (9000.0, 1.5),
    (9250.0, 4.3),
    (9500.0, 9.7),
    (9750.0, 19.5),
    (10000.0, 35.0),
    (10250.0, 55.0),
    (10500.0, 81.0),
    (10750.0, 115.0),
    (11000.0, 151.0),
    (11500.0, 238.0),
    (12000.0, 359.0),
    (12500.0, 495.0),
    (13000.0, 660.0),
    (13500.0, 850.0),
    (14000.0, 1065.0),
    (14500.0, 1313.0),
    (15000.0, 1550.0),
    (15500.0, 1780.0),
    (16000.0, 2040.0),
]


@dataclass(frozen=True)
class HeatingState:
    convective_w_m2: float
    radiative_w_m2: float
    total_external_w_m2: float
    continuum_factor: float
    radiation_nominal_validity: bool


def _linear_table(value: float, table: list[tuple[float, float]]) -> float:
    if value <= table[0][0]:
        return table[0][1]
    if value >= table[-1][0]:
        return table[-1][1]
    for (x0, y0), (x1, y1) in zip(table, table[1:]):
        if x0 <= value <= x1:
            f = (value - x0) / (x1 - x0)
            return y0 + f * (y1 - y0)
    raise RuntimeError("table interpolation failed")


def sutton_graves(
    atmosphere: AtmosphereState,
    velocity_m_s: float,
    nose_radius_m: float,
) -> float:
    if atmosphere.density_kg_m3 <= 0 or velocity_m_s <= 0 or nose_radius_m <= 0:
        return 0.0
    q = (
        SUTTON_GRAVES_EARTH
        * math.sqrt(atmosphere.density_kg_m3 / nose_radius_m)
        * velocity_m_s**3
    )
    return max(0.0, q) * atmosphere.continuum_factor


def tauber_sutton(
    atmosphere: AtmosphereState,
    velocity_m_s: float,
    nose_radius_m: float,
) -> tuple[float, bool]:
    """Earth stagnation-point radiative heating correlation.

    q = 4.736e8 * Rn^a * rho^1.22 * f(V) [W/m2]
    a = 1.072e6 * V^-1.88 * rho^-0.325

    Nominal Earth validity:
    V=10-16 km/s, rho=6.66e-5..6.31e-4 kg/m3, Rn=0.3..3 m.
    """
    rho = atmosphere.density_kg_m3
    if rho <= 0 or velocity_m_s < 9000.0 or nose_radius_m <= 0:
        return 0.0, False

    fv = _linear_table(min(velocity_m_s, 16000.0), TAUBER_SUTTON_FV)
    a = 1.072e6 * velocity_m_s ** (-1.88) * rho ** (-0.325)
    a = min(1.0, max(0.0, a))
    q = 4.736e8 * nose_radius_m**a * rho**1.22 * fv
    q *= atmosphere.continuum_factor

    valid = (
        10000.0 <= velocity_m_s <= 16000.0
        and 6.66e-5 <= rho <= 6.31e-4
        and 0.3 <= nose_radius_m <= 3.0
        and atmosphere.continuum_factor >= 0.999
    )
    return max(0.0, q), valid


def total_heating(
    atmosphere: AtmosphereState,
    velocity_m_s: float,
    nose_radius_m: float,
    config: dict,
) -> HeatingState:
    conv = sutton_graves(atmosphere, velocity_m_s, nose_radius_m)
    if not bool(config.get("radiative_enabled", True)):
        rad, valid = 0.0, False
    else:
        rad, valid = tauber_sutton(atmosphere, velocity_m_s, nose_radius_m)

    conv *= float(config.get("convective_scale", 1.0))
    rad *= float(config.get("radiative_scale", 1.0))
    return HeatingState(
        conv,
        rad,
        conv + rad,
        atmosphere.continuum_factor,
        valid,
    )
