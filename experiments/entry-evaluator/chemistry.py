from __future__ import annotations

import csv
import math
from dataclasses import dataclass
from pathlib import Path


@dataclass(frozen=True)
class ChemistryLimit:
    max_exit_temperature_k: float
    source: str
    estimated_ignition_delay_s: float | None


class ChemistryLimiter:
    """Optional H2 ignition-delay lookup with a fixed-limit fallback.

    The CSV format is the output produced by
    experiments/hydrogen-ignition-delay/ignition_delay.py.
    No reacting-flow calculation is performed inside a trajectory run.
    """

    def __init__(self, config: dict):
        self.mode = str(config.get("mode", "fixed_limit"))
        self.fixed_limit_k = float(
            config.get("fixed_max_exit_temperature_k", 750.0)
        )
        self.phi = float(config.get("phi", 4.0))
        self.residence_time_s = float(config.get("residence_time_s", 1e-3))
        self.safety_factor = float(config.get("ignition_safety_factor", 10.0))
        self.table: list[dict] = []

        if self.mode == "ignition_csv":
            path = Path(str(config.get("ignition_csv", ""))).expanduser()
            if not path.exists():
                raise FileNotFoundError(
                    f"chemistry ignition CSV not found: {path}"
                )
            self.table = self._load(path)
        elif self.mode not in {"fixed_limit", "disabled"}:
            raise ValueError(f"unsupported chemistry.mode: {self.mode}")

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
                if tau is None and status == "no_ignition_within_limit":
                    tau = max_time
                if tau is None or tau <= 0:
                    continue
                rows.append(
                    {
                        "temperature_k": float(row["temperature_initial_K"]),
                        "pressure_bar": float(row["pressure_bar"]),
                        "phi": float(row["phi"]),
                        "tau_s": tau,
                    }
                )
        if not rows:
            raise ValueError("ignition CSV contains no usable rows")
        return rows

    def _nearest_curve(self, pressure_bar: float) -> list[dict]:
        pressures = sorted({r["pressure_bar"] for r in self.table})
        phis = sorted({r["phi"] for r in self.table})
        p = min(
            pressures,
            key=lambda x: abs(
                math.log(max(x, 1e-9) / max(pressure_bar, 1e-9))
            ),
        )
        phi = min(phis, key=lambda x: abs(x - self.phi))
        curve = [
            r
            for r in self.table
            if r["pressure_bar"] == p and r["phi"] == phi
        ]
        return sorted(curve, key=lambda r: r["temperature_k"])

    @staticmethod
    def _log_interp_tau(
        curve: list[dict],
        temperature_k: float,
    ) -> float | None:
        if not curve:
            return None
        if temperature_k <= curve[0]["temperature_k"]:
            return curve[0]["tau_s"]
        if temperature_k >= curve[-1]["temperature_k"]:
            return curve[-1]["tau_s"]
        for a, b in zip(curve, curve[1:]):
            if a["temperature_k"] <= temperature_k <= b["temperature_k"]:
                x = (
                    temperature_k - a["temperature_k"]
                ) / (
                    b["temperature_k"] - a["temperature_k"]
                )
                la = math.log10(a["tau_s"])
                lb = math.log10(b["tau_s"])
                return 10 ** (la + x * (lb - la))
        return None

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

        curve = self._nearest_curve(max(pressure_pa, 1.0) / 1e5)
        required_tau = self.residence_time_s * self.safety_factor
        candidates = sorted({r["temperature_k"] for r in curve})
        safe = [
            t
            for t in candidates
            if t <= requested_max_k
            and (self._log_interp_tau(curve, t) or 0.0) >= required_tau
        ]
        t = max(safe) if safe else min(candidates)
        tau = self._log_interp_tau(curve, t)
        return ChemistryLimit(t, "ignition_csv", tau)
