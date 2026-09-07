from __future__ import annotations

import csv
import math
from dataclasses import dataclass
from pathlib import Path


_TABLE_CACHE: dict[Path, list[dict]] = {}


@dataclass(frozen=True)
class ChemistryLimit:
    max_exit_temperature_k: float
    source: str
    estimated_ignition_delay_s: float | None
    required_ignition_delay_s: float | None = None
    lookup_pressure_bar: float | None = None
    lookup_phi: float | None = None
    feasible: bool = True
    failure_reason: str | None = None


class ChemistryLimiter:
    """Optional H2 ignition-delay lookup with a fixed-limit fallback.

    ignition_csv mode interpolates log10(tau) in temperature and log-pressure.
    It does not extrapolate beyond the table pressure range. Equivalence ratio
    must exist explicitly in the table unless allow_nearest_phi=true.
    """

    def __init__(self, config: dict):
        self.mode = str(config.get("mode", "fixed_limit"))
        self.fixed_limit_k = float(
            config.get("fixed_max_exit_temperature_k", 750.0)
        )
        self.phi = float(config.get("phi", 4.0))
        self.residence_time_s = float(config.get("residence_time_s", 1e-3))
        self.safety_factor = float(config.get("ignition_safety_factor", 10.0))
        self.allow_nearest_phi = bool(config.get("allow_nearest_phi", False))
        self.table: list[dict] = []
        self._limit_cache: dict[tuple[float, float], ChemistryLimit] = {}

        for name, value in (
            ("phi", self.phi),
            ("residence_time_s", self.residence_time_s),
            ("ignition_safety_factor", self.safety_factor),
        ):
            if not math.isfinite(value) or value <= 0.0:
                raise ValueError(
                    f"chemistry.{name} must be finite and positive"
                )

        if self.mode == "ignition_csv":
            path = Path(str(config.get("ignition_csv", ""))).expanduser()
            if not path.exists():
                raise FileNotFoundError(
                    f"chemistry ignition CSV not found: {path}"
                )
            resolved = path.resolve()
            if resolved not in _TABLE_CACHE:
                _TABLE_CACHE[resolved] = self._load(resolved)
            self.table = _TABLE_CACHE[resolved]
            self.curves = {
                (phi, pressure): sorted(
                    [
                        row
                        for row in self.table
                        if row["phi"] == phi
                        and row["pressure_bar"] == pressure
                    ],
                    key=lambda row: row["temperature_k"],
                )
                for phi in sorted({r["phi"] for r in self.table})
                for pressure in sorted({
                    r["pressure_bar"]
                    for r in self.table
                    if r["phi"] == phi
                })
            }
            self.available_phis = sorted({r["phi"] for r in self.table})
            self.pressures_by_phi = {
                phi: sorted({
                    row["pressure_bar"]
                    for row in self.table
                    if row["phi"] == phi
                })
                for phi in self.available_phis
            }
            self.temperature_min_k = min(
                r["temperature_k"] for r in self.table
            )
            self.temperature_max_k = max(
                r["temperature_k"] for r in self.table
            )
        elif self.mode not in {"fixed_limit", "disabled"}:
            raise ValueError(f"unsupported chemistry.mode: {self.mode}")

    @property
    def required_ignition_delay_s(self) -> float | None:
        if self.mode != "ignition_csv":
            return None
        return self.residence_time_s * self.safety_factor

    @staticmethod
    def _load(path: Path) -> list[dict]:
        rows: list[dict] = []
        with path.open("r", newline="", encoding="utf-8") as handle:
            for row in csv.DictReader(handle):
                status = row.get("status", "")
                max_time = float(row.get("max_time_s") or 0.0)
                tau_text = row.get("tau_dTdt_s", "")
                tau = (
                    float(tau_text)
                    if tau_text not in {"", None}
                    else None
                )
                # No ignition within the horizon means tau is larger than the
                # horizon. Using the horizon itself is conservative here.
                if tau is None and status == "no_ignition_within_limit":
                    tau = max_time
                if (
                    tau is None
                    or tau <= 0.0
                    or not math.isfinite(tau)
                ):
                    continue
                rows.append(
                    {
                        "temperature_k": float(
                            row["temperature_initial_K"]
                        ),
                        "pressure_bar": float(row["pressure_bar"]),
                        "phi": float(row["phi"]),
                        "tau_s": tau,
                    }
                )
        if not rows:
            raise ValueError("ignition CSV contains no usable rows")
        return rows

    def _select_phi(self) -> float:
        exact = [
            value
            for value in self.available_phis
            if math.isclose(
                value,
                self.phi,
                rel_tol=0.0,
                abs_tol=1e-12,
            )
        ]
        if exact:
            return exact[0]
        if self.allow_nearest_phi:
            return min(
                self.available_phis,
                key=lambda x: abs(x - self.phi),
            )
        raise ValueError(
            f"chemistry phi={self.phi:g} not present in ignition table; "
            f"available={self.available_phis}"
        )

    def _curve(
        self,
        pressure_bar: float,
        phi: float,
    ) -> list[dict]:
        return self.curves.get((phi, pressure_bar), [])

    @staticmethod
    def _log_interp_tau(
        curve: list[dict],
        temperature_k: float,
    ) -> float | None:
        if not curve:
            return None
        if not (
            curve[0]["temperature_k"]
            <= temperature_k
            <= curve[-1]["temperature_k"]
        ):
            return None

        for row in curve:
            if math.isclose(
                temperature_k,
                row["temperature_k"],
                rel_tol=0.0,
                abs_tol=1e-12,
            ):
                return row["tau_s"]

        for a, b in zip(curve, curve[1:]):
            if (
                a["temperature_k"]
                <= temperature_k
                <= b["temperature_k"]
            ):
                x = (
                    temperature_k - a["temperature_k"]
                ) / (
                    b["temperature_k"] - a["temperature_k"]
                )
                la = math.log10(a["tau_s"])
                lb = math.log10(b["tau_s"])
                return 10 ** (la + x * (lb - la))
        return None

    def _pressure_bracket(
        self,
        pressure_bar: float,
        phi: float,
    ) -> tuple[float, float]:
        pressures = self.pressures_by_phi.get(phi, [])
        if not pressures:
            raise ValueError(
                f"ignition table contains no curve for phi={phi:g}"
            )
        if (
            pressure_bar < pressures[0]
            or pressure_bar > pressures[-1]
        ):
            raise ValueError(
                f"surface pressure {pressure_bar:g} bar outside ignition "
                f"table range {pressures[0]:g}..{pressures[-1]:g} bar "
                f"for phi={phi:g}"
            )
        if math.isclose(
            pressure_bar,
            pressures[-1],
            rel_tol=0.0,
            abs_tol=1e-12,
        ):
            return pressures[-1], pressures[-1]
        for p0, p1 in zip(pressures, pressures[1:]):
            if p0 <= pressure_bar <= p1:
                return p0, p1
        raise RuntimeError("pressure bracketing failed")

    def estimate_delay(
        self,
        pressure_pa: float,
        temperature_k: float,
    ) -> float | None:
        if self.mode != "ignition_csv":
            return None

        pressure_bar = max(float(pressure_pa), 1.0) / 1e5
        phi = self._select_phi()
        p0, p1 = self._pressure_bracket(
            pressure_bar,
            phi,
        )
        tau0 = self._log_interp_tau(
            self._curve(p0, phi),
            temperature_k,
        )
        if tau0 is None:
            return None
        if p0 == p1:
            return tau0

        tau1 = self._log_interp_tau(
            self._curve(p1, phi),
            temperature_k,
        )
        if tau1 is None:
            return None

        x = math.log(pressure_bar / p0) / math.log(p1 / p0)
        return 10 ** (
            math.log10(tau0)
            + x * (
                math.log10(tau1)
                - math.log10(tau0)
            )
        )

    def _temperature_grid(
        self,
        pressure_bar: float,
        phi: float,
        requested_max_k: float,
    ) -> list[float]:
        p0, p1 = self._pressure_bracket(
            pressure_bar,
            phi,
        )
        temperatures = sorted(
            {
                row["temperature_k"]
                for row in self.table
                if row["phi"] == phi
                and row["pressure_bar"] in {p0, p1}
            }
        )
        upper = min(
            requested_max_k,
            max(temperatures),
        )
        grid = [
            temperature_k
            for temperature_k in temperatures
            if temperature_k <= upper
        ]
        if grid and upper > grid[-1]:
            grid.append(upper)
        return grid

    @staticmethod
    def _crossing_temperature(
        t0: float,
        tau0: float,
        t1: float,
        tau1: float,
        required_tau: float,
    ) -> float | None:
        y0 = math.log10(tau0)
        y1 = math.log10(tau1)
        target = math.log10(required_tau)
        if math.isclose(y0, y1):
            return None
        x = (target - y0) / (y1 - y0)
        if not 0.0 <= x <= 1.0:
            return None
        return t0 + x * (t1 - t0)

    def limit(
        self,
        pressure_pa: float,
        requested_max_k: float,
    ) -> ChemistryLimit:
        if self.mode == "disabled":
            return ChemistryLimit(
                requested_max_k,
                "disabled",
                None,
            )
        if self.mode == "fixed_limit":
            return ChemistryLimit(
                min(requested_max_k, self.fixed_limit_k),
                "fixed_limit",
                None,
            )

        cache_key = (
            float(pressure_pa),
            float(requested_max_k),
        )
        if cache_key in self._limit_cache:
            return self._limit_cache[cache_key]

        pressure_bar = max(float(pressure_pa), 1.0) / 1e5
        required_tau = (
            self.residence_time_s
            * self.safety_factor
        )
        try:
            phi = self._select_phi()
            grid = self._temperature_grid(
                pressure_bar,
                phi,
                requested_max_k,
            )
        except ValueError as exc:
            result = ChemistryLimit(
                requested_max_k,
                "ignition_csv",
                None,
                required_tau,
                pressure_bar,
                None,
                False,
                str(exc),
            )
            self._store_cache(cache_key, result)
            return result

        if not grid:
            result = ChemistryLimit(
                requested_max_k,
                "ignition_csv",
                None,
                required_tau,
                pressure_bar,
                phi,
                False,
                "requested temperature does not overlap ignition table range",
            )
            self._store_cache(cache_key, result)
            return result

        samples: list[tuple[float, float]] = []
        for temperature_k in grid:
            try:
                tau = self.estimate_delay(
                    pressure_pa,
                    temperature_k,
                )
            except ValueError as exc:
                result = ChemistryLimit(
                    requested_max_k,
                    "ignition_csv",
                    None,
                    required_tau,
                    pressure_bar,
                    phi,
                    False,
                    str(exc),
                )
                self._store_cache(cache_key, result)
                return result
            if tau is not None:
                samples.append(
                    (temperature_k, tau)
                )

        safe_candidates = [
            temperature_k
            for temperature_k, tau in samples
            if tau >= required_tau
        ]
        for (
            (t0, tau0),
            (t1, tau1),
        ) in zip(samples, samples[1:]):
            if (
                (tau0 - required_tau)
                * (tau1 - required_tau)
                < 0.0
            ):
                crossing = self._crossing_temperature(
                    t0,
                    tau0,
                    t1,
                    tau1,
                    required_tau,
                )
                if crossing is not None:
                    safe_candidates.append(crossing)

        if not safe_candidates:
            result = ChemistryLimit(
                min(grid),
                "ignition_csv",
                None,
                required_tau,
                pressure_bar,
                phi,
                False,
                (
                    "no temperature in ignition table satisfies "
                    f"required delay {required_tau:g} s"
                ),
            )
            self._store_cache(cache_key, result)
            return result

        safe_temperature = min(
            requested_max_k,
            max(safe_candidates),
            self.temperature_max_k,
        )
        try:
            tau = self.estimate_delay(
                pressure_pa,
                safe_temperature,
            )
        except ValueError:
            tau = None

        result = ChemistryLimit(
            safe_temperature,
            "ignition_csv",
            tau,
            required_tau,
            pressure_bar,
            phi,
            True,
            None,
        )
        self._store_cache(
            cache_key,
            result,
        )
        return result

    def _store_cache(
        self,
        key: tuple[float, float],
        value: ChemistryLimit,
    ) -> None:
        if len(self._limit_cache) >= 4096:
            self._limit_cache.clear()
        self._limit_cache[key] = value
