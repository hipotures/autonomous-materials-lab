#!/usr/bin/env python3
"""Equilibrium sanity check for hot reactive working fluids mixed with air.

This script deliberately separates:

1. idealized thermal enthalpy uptake from a configured storage state; and
2. theoretical chemical energy release after mixing the hot gas with air.

It does not model ignition delay, flame propagation, wall heat flux, mixing,
hypersonic flow, or reaction location.
"""

from __future__ import annotations

import argparse
import csv
import math
from dataclasses import dataclass
from pathlib import Path

import cantera as ct
import CoolProp
import CoolProp.CoolProp as CP


AIR = "O2:1,N2:3.76"
DEFAULT_MECHANISM = "gri30.yaml"


@dataclass(frozen=True)
class Fuel:
    label: str
    cantera_fuel: str
    coolprop_name: str
    storage_temperature_k: float
    storage_pressure_pa: float
    notes: str


@dataclass(frozen=True)
class CaseResult:
    fuel: str
    temperature_k: float
    pressure_pa: float
    phi: float
    initial_fuel_mass_fraction: float
    thermal_q_mj_kg_fuel: float | None
    equilibrium_temperature_hp_k: float | None
    chemical_release_tp_mj_kg_fuel: float | None
    chemical_to_thermal_ratio: float | None
    x_h2o_tp: float | None
    x_o2_tp: float | None
    x_h2_tp: float | None
    x_co2_tp: float | None
    x_co_tp: float | None
    x_no_tp: float | None
    status: str


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Compare thermal enthalpy uptake with equilibrium chemical "
            "oxidation potential for hot fuels mixed with air."
        )
    )
    parser.add_argument(
        "--fuels",
        type=Path,
        default=Path(__file__).with_name("fuels.csv"),
        help="Fuel/storage-state CSV file.",
    )
    parser.add_argument(
        "--temperatures-k",
        type=float,
        nargs="+",
        default=[500.0, 750.0, 1000.0, 1250.0, 1500.0],
        help="Initial hot-gas temperatures in kelvin.",
    )
    parser.add_argument(
        "--phi",
        type=float,
        nargs="+",
        default=[0.5, 1.0, 2.0, 4.0],
        help="Equivalence ratios (default: 0.5 1 2 4).",
    )
    parser.add_argument(
        "--pressure-pa",
        type=float,
        default=ct.one_atm,
        help="Hot fuel/air mixture pressure in pascals.",
    )
    parser.add_argument(
        "--mechanism",
        default=DEFAULT_MECHANISM,
        help=f"Cantera mechanism (default: {DEFAULT_MECHANISM}).",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=None,
        help="Optional CSV output path.",
    )
    return parser.parse_args()


def finite(value: float | None) -> bool:
    return value is not None and math.isfinite(value)


def load_fuels(path: Path) -> list[Fuel]:
    fuels: list[Fuel] = []

    with path.open("r", newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        required = {
            "label",
            "cantera_fuel",
            "coolprop_name",
            "T_storage_K",
            "P_storage_Pa",
            "notes",
        }
        missing = required - set(reader.fieldnames or [])
        if missing:
            raise ValueError(
                "Missing required CSV columns: "
                + ", ".join(sorted(missing))
            )

        for row_number, row in enumerate(reader, start=2):
            try:
                temperature_k = float(row["T_storage_K"])
                pressure_pa = float(row["P_storage_Pa"])
            except (TypeError, ValueError) as exc:
                raise ValueError(
                    f"Invalid storage state on row {row_number}"
                ) from exc

            if temperature_k <= 0.0 or pressure_pa <= 0.0:
                raise ValueError(
                    f"Storage T/P must be positive on row {row_number}"
                )

            fuels.append(
                Fuel(
                    label=row["label"].strip(),
                    cantera_fuel=row["cantera_fuel"].strip(),
                    coolprop_name=row["coolprop_name"].strip(),
                    storage_temperature_k=temperature_k,
                    storage_pressure_pa=pressure_pa,
                    notes=row["notes"].strip(),
                )
            )

    return fuels


def coolprop_enthalpy(
    fluid: str,
    temperature_k: float,
    pressure_pa: float,
) -> float | None:
    try:
        value = float(
            CP.PropsSI("H", "T", temperature_k, "P", pressure_pa, fluid)
        )
        return value if math.isfinite(value) else None
    except Exception:
        return None


def thermal_q_from_storage(
    fuel: Fuel,
    hot_temperature_k: float,
    hot_pressure_pa: float,
) -> float | None:
    h_storage = coolprop_enthalpy(
        fuel.coolprop_name,
        fuel.storage_temperature_k,
        fuel.storage_pressure_pa,
    )
    h_hot = coolprop_enthalpy(
        fuel.coolprop_name,
        hot_temperature_k,
        hot_pressure_pa,
    )

    if h_storage is None or h_hot is None:
        return None

    q_j_kg = h_hot - h_storage
    if q_j_kg <= 0.0:
        return None

    return q_j_kg / 1e6


def species_mole_fraction(gas: ct.Solution, species: str) -> float | None:
    if species not in gas.species_names:
        return None
    return float(gas[species].X[0])


def fuel_species_name(cantera_fuel: str) -> str:
    """Return the first configured fuel species name.

    The current experiment uses one-species fuels such as H2:1.
    """

    first = cantera_fuel.split(",")[0].strip()
    return first.split(":")[0].strip()


def evaluate_case(
    fuel: Fuel,
    temperature_k: float,
    pressure_pa: float,
    phi: float,
    mechanism: str,
) -> CaseResult:
    gas = ct.Solution(mechanism)

    fuel_species = fuel_species_name(fuel.cantera_fuel)
    if fuel_species not in gas.species_names:
        return CaseResult(
            fuel=fuel.label,
            temperature_k=temperature_k,
            pressure_pa=pressure_pa,
            phi=phi,
            initial_fuel_mass_fraction=math.nan,
            thermal_q_mj_kg_fuel=thermal_q_from_storage(
                fuel, temperature_k, pressure_pa
            ),
            equilibrium_temperature_hp_k=None,
            chemical_release_tp_mj_kg_fuel=None,
            chemical_to_thermal_ratio=None,
            x_h2o_tp=None,
            x_o2_tp=None,
            x_h2_tp=None,
            x_co2_tp=None,
            x_co_tp=None,
            x_no_tp=None,
            status=f"fuel species {fuel_species!r} missing from mechanism",
        )

    try:
        gas.TP = temperature_k, pressure_pa
        gas.set_equivalence_ratio(
            phi,
            fuel=fuel.cantera_fuel,
            oxidizer=AIR,
        )

        initial_x = gas.X
        initial_enthalpy_mass = float(gas.enthalpy_mass)
        fuel_index = gas.species_index(fuel_species)
        fuel_mass_fraction = float(gas.Y[fuel_index])

        if fuel_mass_fraction <= 0.0:
            raise ValueError("initial fuel mass fraction is zero")

        # Equilibrium at fixed T,P: theoretical chemical enthalpy decrease at
        # the same thermodynamic state.
        gas.equilibrate("TP")
        equilibrium_enthalpy_mass = float(gas.enthalpy_mass)

        q_chemical_mixture_j_kg = (
            initial_enthalpy_mass - equilibrium_enthalpy_mass
        )
        q_chemical_fuel_mj_kg = (
            q_chemical_mixture_j_kg / fuel_mass_fraction / 1e6
        )

        x_h2o = species_mole_fraction(gas, "H2O")
        x_o2 = species_mole_fraction(gas, "O2")
        x_h2 = species_mole_fraction(gas, "H2")
        x_co2 = species_mole_fraction(gas, "CO2")
        x_co = species_mole_fraction(gas, "CO")
        x_no = species_mole_fraction(gas, "NO")

        # Reset the exact initial mixture and compute adiabatic HP equilibrium.
        gas.TPX = temperature_k, pressure_pa, initial_x
        gas.equilibrate("HP")
        equilibrium_temperature_hp_k = float(gas.T)

        thermal_q = thermal_q_from_storage(
            fuel,
            temperature_k,
            pressure_pa,
        )

        ratio = None
        if (
            thermal_q is not None
            and thermal_q > 0.0
            and q_chemical_fuel_mj_kg > 0.0
        ):
            ratio = q_chemical_fuel_mj_kg / thermal_q

        return CaseResult(
            fuel=fuel.label,
            temperature_k=temperature_k,
            pressure_pa=pressure_pa,
            phi=phi,
            initial_fuel_mass_fraction=fuel_mass_fraction,
            thermal_q_mj_kg_fuel=thermal_q,
            equilibrium_temperature_hp_k=equilibrium_temperature_hp_k,
            chemical_release_tp_mj_kg_fuel=q_chemical_fuel_mj_kg,
            chemical_to_thermal_ratio=ratio,
            x_h2o_tp=x_h2o,
            x_o2_tp=x_o2,
            x_h2_tp=x_h2,
            x_co2_tp=x_co2,
            x_co_tp=x_co,
            x_no_tp=x_no,
            status="ok",
        )
    except Exception as exc:
        return CaseResult(
            fuel=fuel.label,
            temperature_k=temperature_k,
            pressure_pa=pressure_pa,
            phi=phi,
            initial_fuel_mass_fraction=math.nan,
            thermal_q_mj_kg_fuel=thermal_q_from_storage(
                fuel, temperature_k, pressure_pa
            ),
            equilibrium_temperature_hp_k=None,
            chemical_release_tp_mj_kg_fuel=None,
            chemical_to_thermal_ratio=None,
            x_h2o_tp=None,
            x_o2_tp=None,
            x_h2_tp=None,
            x_co2_tp=None,
            x_co_tp=None,
            x_no_tp=None,
            status=f"failed: {exc}",
        )


def fmt(value: float | None, digits: int = 3) -> str:
    if value is None or not math.isfinite(value):
        return ""
    return f"{value:.{digits}f}"


def print_results(results: list[CaseResult]) -> None:
    headers = [
        "fuel",
        "T_K",
        "phi",
        "Y_fuel",
        "q_thermal",
        "q_chem",
        "chem/thermal",
        "T_eq_HP",
        "X_H2O",
        "X_O2",
        "X_H2",
        "X_CO2",
        "X_CO",
        "X_NO",
        "status",
    ]

    rows: list[list[str]] = []

    for result in results:
        rows.append(
            [
                result.fuel,
                f"{result.temperature_k:.0f}",
                f"{result.phi:.2f}",
                fmt(result.initial_fuel_mass_fraction, 4),
                fmt(result.thermal_q_mj_kg_fuel),
                fmt(result.chemical_release_tp_mj_kg_fuel),
                fmt(result.chemical_to_thermal_ratio),
                fmt(result.equilibrium_temperature_hp_k, 1),
                fmt(result.x_h2o_tp, 4),
                fmt(result.x_o2_tp, 4),
                fmt(result.x_h2_tp, 4),
                fmt(result.x_co2_tp, 4),
                fmt(result.x_co_tp, 4),
                fmt(result.x_no_tp, 5),
                result.status,
            ]
        )

    widths = [
        max(len(headers[index]), *(len(row[index]) for row in rows))
        for index in range(len(headers))
    ]

    def render(row: list[str]) -> str:
        return "  ".join(
            value.ljust(widths[index])
            for index, value in enumerate(row)
        )

    print(render(headers))
    print(render(["-" * width for width in widths]))
    for row in rows:
        print(render(row))


def write_csv(path: Path, results: list[CaseResult]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)

    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(
            [
                "fuel",
                "temperature_K",
                "pressure_Pa",
                "phi",
                "initial_fuel_mass_fraction",
                "thermal_q_from_storage_MJ_kg_fuel",
                "chemical_release_TP_MJ_kg_fuel",
                "chemical_to_thermal_ratio",
                "equilibrium_temperature_HP_K",
                "equilibrium_TP_X_H2O",
                "equilibrium_TP_X_O2",
                "equilibrium_TP_X_H2",
                "equilibrium_TP_X_CO2",
                "equilibrium_TP_X_CO",
                "equilibrium_TP_X_NO",
                "status",
            ]
        )

        for result in results:
            writer.writerow(
                [
                    result.fuel,
                    result.temperature_k,
                    result.pressure_pa,
                    result.phi,
                    result.initial_fuel_mass_fraction,
                    result.thermal_q_mj_kg_fuel,
                    result.chemical_release_tp_mj_kg_fuel,
                    result.chemical_to_thermal_ratio,
                    result.equilibrium_temperature_hp_k,
                    result.x_h2o_tp,
                    result.x_o2_tp,
                    result.x_h2_tp,
                    result.x_co2_tp,
                    result.x_co_tp,
                    result.x_no_tp,
                    result.status,
                ]
            )


def main() -> None:
    args = parse_args()

    if args.pressure_pa <= 0.0:
        raise SystemExit("--pressure-pa must be positive")
    if any(value <= 0.0 for value in args.temperatures_k):
        raise SystemExit("All temperatures must be positive")
    if any(value <= 0.0 for value in args.phi):
        raise SystemExit("All equivalence ratios must be positive")

    fuels = load_fuels(args.fuels)
    if not fuels:
        raise SystemExit("No fuels configured")

    results = [
        evaluate_case(
            fuel=fuel,
            temperature_k=temperature_k,
            pressure_pa=args.pressure_pa,
            phi=phi,
            mechanism=args.mechanism,
        )
        for fuel in fuels
        for temperature_k in args.temperatures_k
        for phi in args.phi
    ]

    print(
        f"Cantera {ct.__version__} | "
        f"CoolProp {CoolProp.__version__} | "
        f"mechanism={args.mechanism} | "
        f"pressure={args.pressure_pa:.0f} Pa | "
        f"cases={len(results)}"
    )
    print()
    print_results(results)

    if args.output is not None:
        write_csv(args.output, results)
        print()
        print(f"CSV written to: {args.output}")


if __name__ == "__main__":
    main()
