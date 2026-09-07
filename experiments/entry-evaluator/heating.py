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

# Brandis & Johnston, AIAA 2014-2374. The equations are reproduced in
# later NASA/AIAA engineering work and open technical literature.
BJ_RHO_MIN_KG_M3 = 1.0e-5
BJ_RHO_MAX_KG_M3 = 5.0e-3
BJ_RN_MIN_M = 0.2
BJ_RN_MAX_M = 10.0
BJ_CONV_V_MIN_M_S = 3000.0
BJ_RADIATIVE_V_MIN_M_S = 9500.0
BJ_V_MAX_M_S = 17000.0
W_CM2_TO_W_M2 = 1.0e4


@dataclass(frozen=True)
class HeatingState:
    convective_w_m2: float
    radiative_w_m2: float
    total_external_w_m2: float
    continuum_factor: float
    radiation_nominal_validity: bool
    backend: str = "legacy"
    convective_nominal_validity: bool = False


def normalize_backend(value: object) -> str:
    backend = str(value if value is not None else "legacy").strip().lower()
    if backend in {"legacy", "sutton_graves_tauber_sutton"}:
        return "legacy"
    if backend in {"brandis_johnston", "brandis_johnston_2014"}:
        return "brandis_johnston_2014"
    raise ValueError(f"unsupported heating.backend: {backend}")


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


def _scale(config: dict, key: str) -> float:
    value = float(config.get(key, 1.0))
    if not math.isfinite(value) or value < 0.0:
        raise ValueError(f"heating.{key} must be finite and nonnegative")
    return value


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


def _brandis_johnston_common_validity(
    atmosphere: AtmosphereState,
    nose_radius_m: float,
) -> bool:
    return (
        BJ_RHO_MIN_KG_M3 <= atmosphere.density_kg_m3 <= BJ_RHO_MAX_KG_M3
        and BJ_RN_MIN_M <= nose_radius_m <= BJ_RN_MAX_M
        and atmosphere.continuum_factor >= 0.999
    )


def brandis_johnston_convective(
    atmosphere: AtmosphereState,
    velocity_m_s: float,
    nose_radius_m: float,
) -> tuple[float, bool]:
    """Brandis-Johnston Earth stagnation convective heat flux.

    Inputs are SI. The published correlation returns W/cm2; this function
    converts to W/m2.

    3.0 <= V < 9.5 km/s:
      q = 7.455e-9 rho^0.4705 V^3.089 Rn^-0.52

    9.5 <= V < 17 km/s:
      q = 1.270e-6 rho^0.4678 V^2.524 Rn^-0.52
    """
    rho = atmosphere.density_kg_m3
    if (
        rho <= 0.0
        or nose_radius_m <= 0.0
        or not BJ_CONV_V_MIN_M_S <= velocity_m_s < BJ_V_MAX_M_S
    ):
        return 0.0, False

    if velocity_m_s < BJ_RADIATIVE_V_MIN_M_S:
        q_w_cm2 = (
            7.455e-9
            * rho**0.4705
            * velocity_m_s**3.089
            * nose_radius_m**-0.52
        )
    else:
        q_w_cm2 = (
            1.270e-6
            * rho**0.4678
            * velocity_m_s**2.524
            * nose_radius_m**-0.52
        )

    valid = _brandis_johnston_common_validity(atmosphere, nose_radius_m)
    q = max(0.0, q_w_cm2) * W_CM2_TO_W_M2 * atmosphere.continuum_factor
    return q, valid


def _brandis_johnston_a_max(nose_radius_m: float) -> float | None:
    if 0.0 < nose_radius_m <= 0.5:
        return 0.61
    if nose_radius_m <= 2.0:
        return 1.23
    if nose_radius_m <= 10.0:
        return 0.49
    return None


def brandis_johnston_radiative(
    atmosphere: AtmosphereState,
    velocity_m_s: float,
    nose_radius_m: float,
) -> tuple[float, bool]:
    """Brandis-Johnston Earth stagnation radiative heat flux.

    Inputs are SI. The published correlation returns W/cm2; this function
    converts to W/m2. No value is invented outside the published 9.5-17 km/s
    velocity branch or the defined radius-cap branches.
    """
    rho = atmosphere.density_kg_m3
    a_max = _brandis_johnston_a_max(nose_radius_m)
    if (
        rho <= 0.0
        or a_max is None
        or not BJ_RADIATIVE_V_MIN_M_S <= velocity_m_s < BJ_V_MAX_M_S
    ):
        return 0.0, False

    a = min(
        3.175e6 * velocity_m_s**-1.80 * rho**-0.1575,
        a_max,
    )
    f_v = -53.26 + 6555.0 / (
        1.0 + (16000.0 / velocity_m_s) ** 8.25
    )
    q_w_cm2 = 3.416e4 * nose_radius_m**a * rho**1.261 * f_v
    valid = _brandis_johnston_common_validity(atmosphere, nose_radius_m)
    q = max(0.0, q_w_cm2) * W_CM2_TO_W_M2 * atmosphere.continuum_factor
    return q, valid


def total_heating(
    atmosphere: AtmosphereState,
    velocity_m_s: float,
    nose_radius_m: float,
    config: dict,
) -> HeatingState:
    backend = normalize_backend(config.get("backend", "legacy"))

    if backend == "legacy":
        conv = sutton_graves(atmosphere, velocity_m_s, nose_radius_m)
        # V0/V1 did not encode an explicit Sutton-Graves validity envelope.
        conv_valid = False
        if not bool(config.get("radiative_enabled", True)):
            rad, rad_valid = 0.0, False
        else:
            rad, rad_valid = tauber_sutton(
                atmosphere,
                velocity_m_s,
                nose_radius_m,
            )
    else:
        conv, conv_valid = brandis_johnston_convective(
            atmosphere,
            velocity_m_s,
            nose_radius_m,
        )
        if not bool(config.get("radiative_enabled", True)):
            rad, rad_valid = 0.0, False
        else:
            rad, rad_valid = brandis_johnston_radiative(
                atmosphere,
                velocity_m_s,
                nose_radius_m,
            )

    conv *= _scale(config, "convective_scale")
    rad *= _scale(config, "radiative_scale")
    return HeatingState(
        conv,
        rad,
        conv + rad,
        atmosphere.continuum_factor,
        rad_valid,
        backend,
        conv_valid,
    )
