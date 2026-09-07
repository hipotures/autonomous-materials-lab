#!/usr/bin/env python3
"""Sweep homogeneous H2/air ignition delay with Cantera 3.2.

This is a zero-dimensional finite-rate kinetics sanity check. It does not model
mixing, wall heat flux, boundary layers, shocks, or hypersonic flow.
"""

from __future__ import annotations

import argparse
import csv
import math
from dataclasses import dataclass
from pathlib import Path

import cantera as ct


AIR = "O2:1,N2:3.76"
HYDROGEN = "H2:1"
DEFAULT_MECHANISM = "h2o2.yaml"


@dataclass(frozen=True)
class Result:
    temperature_initial_k: float
    pressure_bar: float
    phi: float
    reactor_type: str
    max_time_s: float
    ignition_rise_k: float
    tau_dtdt_s: float | None
    tau_oh_peak_s: float | None
    temperature_max_k: float
    temperature_rise_k: float
    max_dtdt_k_s: float | None
    oh_peak: float | None
    final_time_s: float
    final_temperature_k: float
    status: str
    steps_recorded: int


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Sweep homogeneous H2/air ignition delay with Cantera."
    )
    parser.add_argument(
        "--temperatures-k",
        nargs="+",
        type=float,
        default=[500.0, 750.0, 1000.0, 1250.0, 1500.0],
        help="Initial temperatures in kelvin.",
    )
    parser.add_argument(
        "--pressures-bar",
        nargs="+",
        type=float,
        default=[0.1, 0.3, 1.0, 3.0],
        help="Initial pressures in bar.",
    )
    parser.add_argument(
        "--phi",
        nargs="+",
        type=float,
        default=[0.5, 1.0, 2.0, 4.0, 8.0, 16.0],
        help="Hydrogen/air equivalence ratios.",
    )
    parser.add_argument(
        "--mechanism",
        default=DEFAULT_MECHANISM,
        help=f"Cantera mechanism (default: {DEFAULT_MECHANISM}).",
    )
    parser.add_argument(
        "--reactor",
        choices=["constant-pressure", "constant-volume"],
        default="constant-pressure",
        help="Zero-dimensional adiabatic reactor model.",
    )
    parser.add_argument(
        "--max-time-s",
        type=float,
        default=10.0,
        help="Maximum integration time per case (default: 10 s).",
    )
    parser.add_argument(
        "--ignition-rise-k",
        type=float,
        default=200.0,
        help=(
            "Minimum temperature rise required to classify substantial "
            "ignition (default: 200 K)."
        ),
    )
    parser.add_argument(
        "--advance-limit-k",
        type=float,
        default=5.0,
        help=(
            "Maximum temperature change per recorded integration advance "
            "(default: 5 K)."
        ),
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("ignition_delay.csv"),
        help="CSV output path (default: ignition_delay.csv).",
    )
    return parser.parse_args()


def validate_args(args: argparse.Namespace) -> None:
    if any(value <= 0.0 for value in args.temperatures_k):
        raise SystemExit("All temperatures must be positive.")
    if any(value <= 0.0 for value in args.pressures_bar):
        raise SystemExit("All pressures must be positive.")
    if any(value <= 0.0 for value in args.phi):
        raise SystemExit("All equivalence ratios must be positive.")
    if args.max_time_s <= 0.0:
        raise SystemExit("--max-time-s must be positive.")
    if args.ignition_rise_k <= 0.0:
        raise SystemExit("--ignition-rise-k must be positive.")
    if args.advance_limit_k <= 0.0:
        raise SystemExit("--advance-limit-k must be positive.")


def make_reactor(
    gas: ct.Solution,
    reactor_type: str,
) -> ct.Reactor:
    if reactor_type == "constant-pressure":
        return ct.IdealGasConstPressureReactor(
            gas,
            energy="on",
            clone=True,
        )

    return ct.IdealGasReactor(
        gas,
        energy="on",
        clone=True,
    )


def derivative_peak(
    times: list[float],
    temperatures: list[float],
) -> tuple[float | None, float | None]:
    best_rate: float | None = None
    best_time: float | None = None

    for index in range(1, len(times)):
        dt = times[index] - times[index - 1]
        if dt <= 0.0:
            continue

        rate = (
            temperatures[index] - temperatures[index - 1]
        ) / dt

        if not math.isfinite(rate):
            continue

        if best_rate is None or rate > best_rate:
            best_rate = rate
            best_time = 0.5 * (
                times[index] + times[index - 1]
            )

    return best_time, best_rate


def evaluate_case(
    mechanism: str,
    temperature_k: float,
    pressure_bar: float,
    phi: float,
    reactor_type: str,
    max_time_s: float,
    ignition_rise_k: float,
    advance_limit_k: float,
) -> Result:
    pressure_pa = pressure_bar * 1e5

    gas = ct.Solution(mechanism)
    gas.TP = temperature_k, pressure_pa
    gas.set_equivalence_ratio(
        phi,
        fuel=HYDROGEN,
        oxidizer=AIR,
    )

    if "OH" not in gas.species_names:
        return Result(
            temperature_initial_k=temperature_k,
            pressure_bar=pressure_bar,
            phi=phi,
            reactor_type=reactor_type,
            max_time_s=max_time_s,
            ignition_rise_k=ignition_rise_k,
            tau_dtdt_s=None,
            tau_oh_peak_s=None,
            temperature_max_k=temperature_k,
            temperature_rise_k=0.0,
            max_dtdt_k_s=None,
            oh_peak=None,
            final_time_s=0.0,
            final_temperature_k=temperature_k,
            status="OH_missing_from_mechanism",
            steps_recorded=0,
        )

    reactor = make_reactor(gas, reactor_type)
    reactor.set_advance_limit("temperature", advance_limit_k)
    network = ct.ReactorNet([reactor])

    times = [0.0]
    temperatures = [float(reactor.T)]
    oh = [float(reactor.phase["OH"].X[0])]

    try:
        while network.time < max_time_s:
            previous_time = network.time
            network.advance(max_time_s)

            current_time = float(network.time)
            current_temperature = float(reactor.T)
            current_oh = float(reactor.phase["OH"].X[0])

            if current_time <= previous_time:
                break

            times.append(current_time)
            temperatures.append(current_temperature)
            oh.append(current_oh)

            # Once a large post-ignition temperature rise has been sampled,
            # continue until dT/dt has clearly passed its maximum. The advance
            # limit keeps the rapid event resolved without forcing tiny steps
            # during a non-reacting induction period.
            if (
                current_temperature - temperature_k
                >= max(ignition_rise_k + 500.0, 700.0)
                and len(times) >= 5
            ):
                recent_rates = []
                for index in range(max(1, len(times) - 4), len(times)):
                    dt = times[index] - times[index - 1]
                    if dt > 0.0:
                        recent_rates.append(
                            (
                                temperatures[index]
                                - temperatures[index - 1]
                            )
                            / dt
                        )

                if (
                    len(recent_rates) >= 3
                    and recent_rates[-1] < max(recent_rates[:-1])
                ):
                    break

    except Exception as exc:
        return Result(
            temperature_initial_k=temperature_k,
            pressure_bar=pressure_bar,
            phi=phi,
            reactor_type=reactor_type,
            max_time_s=max_time_s,
            ignition_rise_k=ignition_rise_k,
            tau_dtdt_s=None,
            tau_oh_peak_s=None,
            temperature_max_k=max(temperatures),
            temperature_rise_k=max(temperatures) - temperature_k,
            max_dtdt_k_s=None,
            oh_peak=max(oh),
            final_time_s=times[-1],
            final_temperature_k=temperatures[-1],
            status=f"integration_failed: {exc}",
            steps_recorded=len(times),
        )

    temperature_max = max(temperatures)
    temperature_rise = temperature_max - temperature_k
    ignited = temperature_rise >= ignition_rise_k

    tau_dtdt, max_dtdt = derivative_peak(times, temperatures)

    oh_peak = max(oh)
    oh_peak_index = oh.index(oh_peak)
    tau_oh_peak = times[oh_peak_index]

    if not ignited:
        tau_dtdt = None
        tau_oh_peak = None
        status = "no_ignition_within_limit"
    else:
        status = "ignited"

    return Result(
        temperature_initial_k=temperature_k,
        pressure_bar=pressure_bar,
        phi=phi,
        reactor_type=reactor_type,
        max_time_s=max_time_s,
        ignition_rise_k=ignition_rise_k,
        tau_dtdt_s=tau_dtdt,
        tau_oh_peak_s=tau_oh_peak,
        temperature_max_k=temperature_max,
        temperature_rise_k=temperature_rise,
        max_dtdt_k_s=max_dtdt if ignited else None,
        oh_peak=oh_peak,
        final_time_s=times[-1],
        final_temperature_k=temperatures[-1],
        status=status,
        steps_recorded=len(times),
    )


def fmt(value: float | None, digits: int = 4) -> str:
    if value is None or not math.isfinite(value):
        return ""
    return f"{value:.{digits}g}"


def print_results(results: list[Result]) -> None:
    headers = [
        "T0_K",
        "P_bar",
        "phi",
        "tau_s",
        "tau_ms",
        "tau_us",
        "Tmax_K",
        "dT_K",
        "OH_peak",
        "status",
    ]

    rows: list[list[str]] = []

    for result in results:
        tau = result.tau_dtdt_s
        rows.append(
            [
                f"{result.temperature_initial_k:.0f}",
                f"{result.pressure_bar:g}",
                f"{result.phi:g}",
                fmt(tau),
                fmt(None if tau is None else tau * 1e3),
                fmt(None if tau is None else tau * 1e6),
                f"{result.temperature_max_k:.1f}",
                f"{result.temperature_rise_k:.1f}",
                fmt(result.oh_peak, 3),
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


def write_csv(path: Path, results: list[Result]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)

    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(
            [
                "temperature_initial_K",
                "pressure_bar",
                "phi",
                "reactor_type",
                "max_time_s",
                "ignition_rise_K",
                "tau_dTdt_s",
                "tau_dTdt_ms",
                "tau_dTdt_us",
                "tau_OH_peak_s",
                "temperature_max_K",
                "temperature_rise_K",
                "max_dTdt_K_s",
                "OH_peak_mole_fraction",
                "final_time_s",
                "final_temperature_K",
                "status",
                "steps_recorded",
            ]
        )

        for result in results:
            tau = result.tau_dtdt_s
            writer.writerow(
                [
                    result.temperature_initial_k,
                    result.pressure_bar,
                    result.phi,
                    result.reactor_type,
                    result.max_time_s,
                    result.ignition_rise_k,
                    tau,
                    None if tau is None else tau * 1e3,
                    None if tau is None else tau * 1e6,
                    result.tau_oh_peak_s,
                    result.temperature_max_k,
                    result.temperature_rise_k,
                    result.max_dtdt_k_s,
                    result.oh_peak,
                    result.final_time_s,
                    result.final_temperature_k,
                    result.status,
                    result.steps_recorded,
                ]
            )


def main() -> None:
    args = parse_args()
    validate_args(args)

    results = [
        evaluate_case(
            mechanism=args.mechanism,
            temperature_k=temperature_k,
            pressure_bar=pressure_bar,
            phi=phi,
            reactor_type=args.reactor,
            max_time_s=args.max_time_s,
            ignition_rise_k=args.ignition_rise_k,
            advance_limit_k=args.advance_limit_k,
        )
        for temperature_k in args.temperatures_k
        for pressure_bar in args.pressures_bar
        for phi in args.phi
    ]

    print(
        f"Cantera {ct.__version__} | "
        f"mechanism={args.mechanism} | "
        f"reactor={args.reactor} | "
        f"cases={len(results)} | "
        f"max_time={args.max_time_s:g} s | "
        f"ignition_rise={args.ignition_rise_k:g} K"
    )
    print()
    print_results(results)

    write_csv(args.output, results)
    print()
    print(f"CSV written to: {args.output}")


if __name__ == "__main__":
    main()
