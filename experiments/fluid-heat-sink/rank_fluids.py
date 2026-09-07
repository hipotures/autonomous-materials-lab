#!/usr/bin/env python3
"""Rank known CoolProp fluids by idealized heat absorption before vaporization.

The model is intentionally simple:

    q_total = h_sat_vapor(P) - h_initial(T0, P)

Only fluids that are liquid at the initial state and have a saturation state at
the selected pressure are included.

This script is a screening experiment, not an atmospheric-entry simulation.
"""

from __future__ import annotations

import argparse
import csv
import math
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

import CoolProp
import CoolProp.CoolProp as CP


@dataclass(frozen=True)
class FluidResult:
    fluid: str
    boiling_temperature_k: float
    density_initial_kg_m3: float
    sensible_mj_kg: float
    latent_mj_kg: float
    total_mj_kg: float
    volumetric_mj_l: float
    mass_ratio_vs_water: float


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Rank CoolProp fluids by idealized heat absorbed from an initial "
            "liquid state to saturated vapor at a fixed pressure."
        )
    )
    parser.add_argument(
        "--temperature-k",
        type=float,
        default=293.15,
        help="Initial fluid temperature in kelvin (default: 293.15 K).",
    )
    parser.add_argument(
        "--pressure-pa",
        type=float,
        default=101325.0,
        help="Reference pressure in pascals (default: 101325 Pa).",
    )
    parser.add_argument(
        "--top",
        type=int,
        default=25,
        help="Number of rows to print. Use 0 to print all (default: 25).",
    )
    parser.add_argument(
        "--csv",
        type=Path,
        default=None,
        help="Optional CSV output path.",
    )
    return parser.parse_args()


def all_coolprop_fluids() -> list[str]:
    raw = CP.get_global_param_string("fluids_list")
    return sorted({name.strip() for name in raw.split(",") if name.strip()})


def finite_positive(value: float) -> bool:
    return math.isfinite(value) and value > 0.0


def evaluate_fluid(fluid: str, temperature_k: float, pressure_pa: float) -> dict[str, float] | None:
    """Return thermodynamic screening values or None if the fluid is unsuitable.

    A fluid is considered eligible only when:
    - a saturation temperature exists at the selected pressure;
    - the initial temperature is below that saturation temperature;
    - the initial state can be evaluated as a liquid-like single phase;
    - all required enthalpy and density values are finite.
    """

    try:
        critical_pressure = CP.PropsSI("Pcrit", fluid)
        if not finite_positive(critical_pressure) or pressure_pa >= critical_pressure:
            return None

        boiling_temperature_k = CP.PropsSI(
            "T", "P", pressure_pa, "Q", 0.0, fluid
        )
        if not math.isfinite(boiling_temperature_k):
            return None

        # Require a meaningful subcooled-liquid initial state.
        if temperature_k >= boiling_temperature_k - 1e-6:
            return None

        phase = CP.PhaseSI("T", temperature_k, "P", pressure_pa, fluid)
        if phase not in {"liquid", "supercritical_liquid"}:
            return None

        h_initial = CP.PropsSI("H", "T", temperature_k, "P", pressure_pa, fluid)
        h_sat_liquid = CP.PropsSI("H", "P", pressure_pa, "Q", 0.0, fluid)
        h_sat_vapor = CP.PropsSI("H", "P", pressure_pa, "Q", 1.0, fluid)
        density_initial = CP.PropsSI(
            "D", "T", temperature_k, "P", pressure_pa, fluid
        )

        values = [
            h_initial,
            h_sat_liquid,
            h_sat_vapor,
            density_initial,
        ]
        if not all(math.isfinite(v) for v in values):
            return None

        sensible = h_sat_liquid - h_initial
        latent = h_sat_vapor - h_sat_liquid
        total = h_sat_vapor - h_initial

        if sensible < 0.0 or latent <= 0.0 or total <= 0.0 or density_initial <= 0.0:
            return None

        return {
            "boiling_temperature_k": boiling_temperature_k,
            "density_initial_kg_m3": density_initial,
            "sensible_j_kg": sensible,
            "latent_j_kg": latent,
            "total_j_kg": total,
        }
    except Exception:
        # CoolProp contains many fluids with different validity ranges. A failed
        # property lookup is treated as "not eligible for this simple screen",
        # not as a scientific rejection of the fluid.
        return None


def build_results(
    fluids: Iterable[str],
    temperature_k: float,
    pressure_pa: float,
) -> list[FluidResult]:
    water = evaluate_fluid("Water", temperature_k, pressure_pa)
    if water is None:
        raise RuntimeError(
            "Water is not a valid reference for the selected temperature/pressure."
        )

    q_water = water["total_j_kg"]
    results: list[FluidResult] = []

    for fluid in fluids:
        values = evaluate_fluid(fluid, temperature_k, pressure_pa)
        if values is None:
            continue

        total_j_kg = values["total_j_kg"]
        density = values["density_initial_kg_m3"]

        results.append(
            FluidResult(
                fluid=fluid,
                boiling_temperature_k=values["boiling_temperature_k"],
                density_initial_kg_m3=density,
                sensible_mj_kg=values["sensible_j_kg"] / 1e6,
                latent_mj_kg=values["latent_j_kg"] / 1e6,
                total_mj_kg=total_j_kg / 1e6,
                volumetric_mj_l=(total_j_kg * density) / 1e9,
                mass_ratio_vs_water=q_water / total_j_kg,
            )
        )

    return sorted(results, key=lambda item: item.total_mj_kg, reverse=True)


def print_table(results: list[FluidResult], top: int) -> None:
    shown = results if top == 0 else results[:top]

    headers = (
        "rank",
        "fluid",
        "boil_C",
        "rho_kg_m3",
        "sensible_MJ_kg",
        "latent_MJ_kg",
        "total_MJ_kg",
        "MJ_L",
        "mass_vs_water",
    )

    rows = []
    for index, item in enumerate(shown, start=1):
        rows.append(
            (
                str(index),
                item.fluid,
                f"{item.boiling_temperature_k - 273.15:.2f}",
                f"{item.density_initial_kg_m3:.1f}",
                f"{item.sensible_mj_kg:.4f}",
                f"{item.latent_mj_kg:.4f}",
                f"{item.total_mj_kg:.4f}",
                f"{item.volumetric_mj_l:.4f}",
                f"{item.mass_ratio_vs_water:.4f}",
            )
        )

    widths = [
        max(len(headers[i]), *(len(row[i]) for row in rows))
        for i in range(len(headers))
    ]

    def format_row(row: tuple[str, ...]) -> str:
        return "  ".join(value.ljust(widths[i]) for i, value in enumerate(row))

    print(format_row(headers))
    print(format_row(tuple("-" * width for width in widths)))
    for row in rows:
        print(format_row(row))


def write_csv(path: Path, results: list[FluidResult]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(
            [
                "rank",
                "fluid",
                "boiling_temperature_k",
                "density_initial_kg_m3",
                "sensible_mj_kg",
                "latent_mj_kg",
                "total_mj_kg",
                "volumetric_mj_l",
                "mass_ratio_vs_water",
            ]
        )
        for rank, item in enumerate(results, start=1):
            writer.writerow(
                [
                    rank,
                    item.fluid,
                    item.boiling_temperature_k,
                    item.density_initial_kg_m3,
                    item.sensible_mj_kg,
                    item.latent_mj_kg,
                    item.total_mj_kg,
                    item.volumetric_mj_l,
                    item.mass_ratio_vs_water,
                ]
            )


def main() -> None:
    args = parse_args()

    if args.temperature_k <= 0.0:
        raise SystemExit("--temperature-k must be > 0")
    if args.pressure_pa <= 0.0:
        raise SystemExit("--pressure-pa must be > 0")
    if args.top < 0:
        raise SystemExit("--top must be >= 0")

    results = build_results(
        all_coolprop_fluids(),
        temperature_k=args.temperature_k,
        pressure_pa=args.pressure_pa,
    )

    print(
        f"CoolProp {CoolProp.__version__} | "
        f"T0={args.temperature_k:.2f} K | "
        f"P={args.pressure_pa:.0f} Pa | "
        f"eligible fluids={len(results)}"
    )
    print()
    print_table(results, args.top)

    water_rank = next(
        (index for index, item in enumerate(results, start=1) if item.fluid == "Water"),
        None,
    )
    if water_rank is not None:
        water = results[water_rank - 1]
        print()
        print(
            f"Water reference: rank={water_rank}, "
            f"q_total={water.total_mj_kg:.4f} MJ/kg, "
            f"volumetric={water.volumetric_mj_l:.4f} MJ/L"
        )

    if args.csv is not None:
        write_csv(args.csv, results)
        print(f"CSV written to: {args.csv}")


if __name__ == "__main__":
    main()
