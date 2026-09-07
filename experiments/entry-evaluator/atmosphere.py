from __future__ import annotations

from dataclasses import dataclass
from datetime import datetime, timezone
from math import pi, sqrt
from typing import Any

import numpy as np
import pymsis

BOLTZMANN = 1.380649e-23
DEFAULT_COLLISION_DIAMETER_M = 3.7e-10


@dataclass(frozen=True)
class AtmosphereState:
    altitude_m: float
    density_kg_m3: float
    temperature_k: float
    pressure_pa: float
    number_density_m3: float
    mean_free_path_m: float
    knudsen: float
    continuum_factor: float
    n2_m3: float
    o2_m3: float
    o_m3: float


class AtmosphereModel:
    """Precomputed NRLMSIS 2.1 atmosphere interpolated by altitude.

    The trajectory evaluator uses a fixed date/location and fixed F10.7/Ap values
    during one entry. This keeps the atmosphere deterministic and avoids any
    runtime download of space-weather data.
    """

    def __init__(self, config: dict[str, Any], characteristic_length_m: float):
        self.config = config
        self.characteristic_length_m = characteristic_length_m
        self.collision_diameter_m = float(
            config.get("collision_diameter_m", DEFAULT_COLLISION_DIAMETER_M)
        )

        min_alt_km = float(config.get("min_altitude_km", 0.0))
        max_alt_km = float(config.get("max_altitude_km", 220.0))
        step_km = float(config.get("grid_step_km", 0.5))
        if step_km <= 0:
            raise ValueError("atmosphere.grid_step_km must be positive")

        self.altitudes_km = np.arange(min_alt_km, max_alt_km + 0.5 * step_km, step_km)
        n = self.altitudes_km.size

        date = np.datetime64(
            self._parse_datetime(config.get("date", "2026-09-07T12:00:00Z"))
        )
        dates = np.full(n, date, dtype="datetime64[s]")
        lons = np.full(n, float(config.get("longitude_deg", 0.0)))
        lats = np.full(n, float(config.get("latitude_deg", 0.0)))
        f107 = float(config.get("f107", 150.0))
        f107a = float(config.get("f107a", 150.0))
        ap = float(config.get("ap", 7.0))
        f107s = np.full(n, f107)
        f107as = np.full(n, f107a)
        aps = np.full((n, 7), ap)

        raw = pymsis.calculate(
            dates,
            lons,
            lats,
            self.altitudes_km,
            f107s=f107s,
            f107as=f107as,
            aps=aps,
            version=2.1,
        )
        data = np.asarray(raw).reshape(n, -1)

        self.density = np.nan_to_num(
            data[:, pymsis.Variable.MASS_DENSITY].astype(float),
            nan=1e-16,
            posinf=1e-16,
            neginf=1e-16,
        )
        self.temperature = np.nan_to_num(
            data[:, pymsis.Variable.TEMPERATURE].astype(float),
            nan=200.0,
            posinf=200.0,
            neginf=200.0,
        )

        def species(variable: pymsis.Variable) -> np.ndarray:
            return np.nan_to_num(
                data[:, variable].astype(float),
                nan=0.0,
                posinf=0.0,
                neginf=0.0,
            )

        self.n2 = species(pymsis.Variable.N2)
        self.o2 = species(pymsis.Variable.O2)
        self.o = species(pymsis.Variable.O)
        self.he = species(pymsis.Variable.HE)
        self.h = species(pymsis.Variable.H)
        self.ar = species(pymsis.Variable.AR)
        self.n = species(pymsis.Variable.N)
        self.no = species(pymsis.Variable.NO)

        species_arrays = [
            self.n2,
            self.o2,
            self.o,
            self.he,
            self.h,
            self.ar,
            self.n,
            self.no,
        ]
        self.number_density = np.sum(
            np.vstack(species_arrays),
            axis=0,
        )
        self.pressure = self.number_density * BOLTZMANN * self.temperature

    @staticmethod
    def _parse_datetime(value: str) -> str:
        text = str(value).strip()
        if text.endswith("Z"):
            text = text[:-1]
        datetime.fromisoformat(text).replace(tzinfo=timezone.utc)
        return text

    @staticmethod
    def _continuum_factor(knudsen: float) -> float:
        """Smooth screening transition from continuum to rarefied flow.

        This is not a rarefied-flow solution. It only prevents continuum
        heating correlations from being applied at clearly invalid Knudsen
        numbers in the V0 evaluator.
        """
        if knudsen <= 0.01:
            return 1.0
        if knudsen >= 0.1:
            return 0.0
        return -np.log10(knudsen) - 1.0

    def sample(self, altitude_m: float) -> AtmosphereState:
        altitude_km = np.clip(
            altitude_m / 1000.0,
            self.altitudes_km[0],
            self.altitudes_km[-1],
        )

        def interp(values: np.ndarray) -> float:
            return float(np.interp(altitude_km, self.altitudes_km, values))

        density = max(interp(self.density), 1e-16)
        temperature = max(interp(self.temperature), 1.0)
        pressure = max(interp(self.pressure), 1e-12)
        number_density = max(interp(self.number_density), 1.0)

        mean_free_path = 1.0 / (
            sqrt(2.0) * pi * self.collision_diameter_m**2 * number_density
        )
        kn = mean_free_path / max(self.characteristic_length_m, 1e-6)

        return AtmosphereState(
            altitude_m=altitude_m,
            density_kg_m3=density,
            temperature_k=temperature,
            pressure_pa=pressure,
            number_density_m3=number_density,
            mean_free_path_m=mean_free_path,
            knudsen=kn,
            continuum_factor=self._continuum_factor(kn),
            n2_m3=interp(self.n2),
            o2_m3=interp(self.o2),
            o_m3=interp(self.o),
        )
