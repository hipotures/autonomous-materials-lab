#!/usr/bin/env python3
"""Compare candidate liquids from their own initial storage states.

For every configured fluid the script computes idealized enthalpy uptake from
the initial liquid state (T0, P0) to:

1. saturated vapor at P0;
2. one or more target temperatures at P0.

This is a thermodynamic screening calculation only. It does not model chemical
decomposition, aerodynamics, nozzles, trajectories, tank mass, or insulation.
"""

from __future__ import annotations

import argparse
import csv
import math
from dataclasses import dataclass
from pathlib import Path

import CoolProp
import CoolProp.CoolProp as CP


@dataclass(frozen=True)
class Candidate:
    label: str
    coolprop_name: str
    temperature_k: float
    pressure_pa: float
    notes: str


@dataclass
class Result:
    candidate: Candidate
    phase_initial: str
    density_kg_m3: float
    boiling_temperature_k: float | None
    q_to_vapor_mj_kg: float | None
    q_to_vapor_mj_l: float | None
    q_targets_mj_kg: dict[float, float | None]
    q_targets_mj_l: dict[float, float | None]
    errors: list[str]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Compare liquid working-fluid candidates from candidate-specific "
            "initial storage states."
        )
    )
    parser.add_argument(
        "--fluids",
        type=Path,
        default=Path(__file__).with_name("fluids.csv"),
        help="Candidate CSV file.",
    )
    parser.add_argument(
        "--targets-k",
        type=float,
        nargs="*",
        default=[500.0, 1000.0],
        help="Final temperatures in kelvin (default: 500 1000).",
    )
    parser.add_argument(
        "--water-label",
        default="Water",
        help="Label used as the normalization reference (default: Water).",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=None,
        help="Optional output CSV path.",
    )
    return parser.parse_args()


def finite(value: float) -> bool:
    return math.isfinite(value)


def safe_props(*args: object) -> float | None:
    try:
        value = float(CP.PropsSI(*args))
        return value if finite(value) else None
    except Exception:
        return None


def load_candidates(path: Path) -> list[Candidate]:
    candidates: list[Candidate] = []

    with path.open("r", newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        required = {"label", "coolprop_name", "T0_K", "P0_Pa", "notes"}
        missing = required - set(reader.fieldnames or [])
        if missing:
            raise ValueError(
                f"Missing required CSV columns: {', '.join(sorted(missing))}"
            )

        for row_number, row in enumerate(reader, start=2):
            try:
                temperature_k = float(row["T0_K"])
                pressure_pa = float(row["P0_Pa"])
            except (TypeError, ValueError) as exc:
                raise ValueError(
                    f"Invalid numeric state on CSV row {row_number}"
                ) from exc

            if temperature_k <= 0.0 or pressure_pa <= 0.0:
                raise ValueError(
                    f"T0_K and P0_Pa must be positive on CSV row {row_number}"
                )

            candidates.append(
                Candidate(
                    label=row["label"].strip(),
                    coolprop_name=row["coolprop_name"].strip(),
                    temperature_k=temperature_k,
                    pressure_pa=pressure_pa,
                    notes=row["notes"].strip(),
                )
            )

    return candidates


def evaluate(candidate: Candidate, target_temperatures_k: list[float]) -> Result:
    errors: list[str] = []

    try:
        phase = CP.PhaseSI(
            "T",
            candidate.temperature_k,
            "P",
            candidate.pressure_pa,
            candidate.coolprop_name,
        )
    except Exception as exc:
        return Result(
            candidate=candidate,
            phase_initial="unsupported",
            density_kg_m3=math.nan,
            boiling_temperature_k=None,
            q_to_vapor_mj_kg=None,
            q_to_vapor_mj_l=None,
            q_targets_mj_kg={t: None for t in target_temperatures_k},
            q_targets_mj_l={t: None for t in target_temperatures_k},
            errors=[f"initial state unsupported: {exc}"],
        )

    if phase not in {"liquid", "supercritical_liquid"}:
        return Result(
            candidate=candidate,
            phase_initial=phase,
            density_kg_m3=math.nan,
            boiling_temperature_k=None,
            q_to_vapor_mj_kg=None,
            q_to_vapor_mj_l=None,
            q_targets_mj_kg={t: None for t in target_temperatures_k},
            q_targets_mj_l={t: None for t in target_temperatures_k},
            errors=[f"initial phase is {phase}, not liquid"],
        )

    h_initial = safe_props(
        "H",
        "T",
        candidate.temperature_k,
        "P",
        candidate.pressure_pa,
        candidate.coolprop_name,
    )
    density = safe_props(
        "D",
        "T",
        candidate.temperature_k,
        "P",
        candidate.pressure_pa,
        candidate.coolprop_name,
    )

    if h_initial is None:
        errors.append("initial enthalpy unavailable")
    if density is None or density <= 0.0:
        errors.append("initial density unavailable")

    boiling_temperature_k: float | None = None
    q_to_vapor_mj_kg: float | None = None
    q_to_vapor_mj_l: float | None = None

    p_crit = safe_props("Pcrit", candidate.coolprop_name)
    p_triple = safe_props("ptriple", candidate.coolprop_name)

    saturation_allowed = (
        p_crit is not None
        and candidate.pressure_pa < p_crit
        and (p_triple is None or candidate.pressure_pa >= p_triple)
    )

    if saturation_allowed and h_initial is not None:
        boiling_temperature_k = safe_props(
            "T",
            "P",
            candidate.pressure_pa,
            "Q",
            0.0,
            candidate.coolprop_name,
        )
        h_sat_vapor = safe_props(
            "H",
            "P",
            candidate.pressure_pa,
            "Q",
            1.0,
            candidate.coolprop_name,
        )

        if h_sat_vapor is not None:
            q_j_kg = h_sat_vapor - h_initial
            if q_j_kg > 0.0:
                q_to_vapor_mj_kg = q_j_kg / 1e6
                if density is not None and density > 0.0:
                    q_to_vapor_mj_l = q_j_kg * density / 1e9
            else:
                errors.append("saturated-vapor enthalpy is not above initial enthalpy")
        else:
            errors.append("saturated-vapor state unavailable")
    else:
        errors.append("saturation calculation unavailable at candidate pressure")

    q_targets_mj_kg: dict[float, float | None] = {}
    q_targets_mj_l: dict[float, float | None] = {}

    for target_k in target_temperatures_k:
        q_targets_mj_kg[target_k] = None
        q_targets_mj_l[target_k] = None

        if target_k <= candidate.temperature_k:
            errors.append(
                f"target {target_k:g} K is not above initial temperature"
            )
            continue

        if h_initial is None:
            continue

        h_target = safe_props(
            "H",
            "T",
            target_k,
            "P",
            candidate.pressure_pa,
            candidate.coolprop_name,
        )

        if h_target is None:
            errors.append(f"{target_k:g} K state unavailable")
            continue

        q_j_kg = h_target - h_initial
        if q_j_kg <= 0.0:
            errors.append(f"{target_k:g} K enthalpy rise is non-positive")
            continue

        q_targets_mj_kg[target_k] = q_j_kg / 1e6
        if density is not None and density > 0.0:
            q_targets_mj_l[target_k] = q_j_kg * density / 1e9

    return Result(
        candidate=candidate,
        phase_initial=phase,
        density_kg_m3=density if density is not None else math.nan,
        boiling_temperature_k=boiling_temperature_k,
        q_to_vapor_mj_kg=q_to_vapor_mj_kg,
        q_to_vapor_mj_l=q_to_vapor_mj_l,
        q_targets_mj_kg=q_targets_mj_kg,
        q_targets_mj_l=q_targets_mj_l,
        errors=errors,
    )


def ratio(reference: float | None, candidate: float | None) -> float | None:
    if (
        reference is None
        or candidate is None
        or reference <= 0.0
        or candidate <= 0.0
    ):
        return None
    return reference / candidate


def fmt(value: float | None, digits: int = 4) -> str:
    return "" if value is None or not math.isfinite(value) else f"{value:.{digits}f}"


def water_reference(
    results: list[Result],
    water_label: str,
) -> Result:
    for result in results:
        if result.candidate.label == water_label:
            return result
    raise ValueError(f"Reference label not found: {water_label!r}")


def print_results(
    results: list[Result],
    reference: Result,
    targets_k: list[float],
) -> None:
    headers = [
        "rank",
        "fluid",
        "T0_K",
        "P0_bar",
        "phase",
        "rho",
        "boil_C",
        "q_vap",
        "vs_H2O",
    ]
    for target in targets_k:
        headers.extend([f"q_{target:g}K", f"vs_H2O_{target:g}K"])

    rows: list[list[str]] = []

    for rank_number, result in enumerate(results, start=1):
        row = [
            str(rank_number),
            result.candidate.label,
            f"{result.candidate.temperature_k:.2f}",
            f"{result.candidate.pressure_pa / 1e5:.3f}",
            result.phase_initial,
            fmt(result.density_kg_m3, 1),
            (
                ""
                if result.boiling_temperature_k is None
                else f"{result.boiling_temperature_k - 273.15:.2f}"
            ),
            fmt(result.q_to_vapor_mj_kg),
            fmt(
                ratio(
                    reference.q_to_vapor_mj_kg,
                    result.q_to_vapor_mj_kg,
                )
            ),
        ]

        for target in targets_k:
            q_value = result.q_targets_mj_kg[target]
            row.extend(
                [
                    fmt(q_value),
                    fmt(
                        ratio(
                            reference.q_targets_mj_kg[target],
                            q_value,
                        )
                    ),
                ]
            )

        rows.append(row)

    widths = [
        max(len(headers[i]), *(len(row[i]) for row in rows))
        for i in range(len(headers))
    ]

    def format_row(row: list[str]) -> str:
        return "  ".join(
            value.ljust(widths[index]) for index, value in enumerate(row)
        )

    print(format_row(headers))
    print(format_row(["-" * width for width in widths]))
    for row in rows:
        print(format_row(row))


def write_csv(
    path: Path,
    results: list[Result],
    reference: Result,
    targets_k: list[float],
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)

    headers = [
        "label",
        "coolprop_name",
        "T0_K",
        "P0_Pa",
        "phase_initial",
        "density_initial_kg_m3",
        "boiling_temperature_K",
        "q_to_vapor_MJ_kg",
        "q_to_vapor_MJ_L",
        "mass_ratio_vs_water_to_vapor",
    ]
    for target in targets_k:
        headers.extend(
            [
                f"q_to_{target:g}K_MJ_kg",
                f"q_to_{target:g}K_MJ_L",
                f"mass_ratio_vs_water_to_{target:g}K",
            ]
        )
    headers.extend(["notes", "errors"])

    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(headers)

        for result in results:
            row: list[object] = [
                result.candidate.label,
                result.candidate.coolprop_name,
                result.candidate.temperature_k,
                result.candidate.pressure_pa,
                result.phase_initial,
                result.density_kg_m3,
                result.boiling_temperature_k,
                result.q_to_vapor_mj_kg,
                result.q_to_vapor_mj_l,
                ratio(
                    reference.q_to_vapor_mj_kg,
                    result.q_to_vapor_mj_kg,
                ),
            ]

            for target in targets_k:
                q_value = result.q_targets_mj_kg[target]
                row.extend(
                    [
                        q_value,
                        result.q_targets_mj_l[target],
                        ratio(reference.q_targets_mj_kg[target], q_value),
                    ]
                )

            row.extend(
                [
                    result.candidate.notes,
                    " | ".join(result.errors),
                ]
            )
            writer.writerow(row)


def main() -> None:
    args = parse_args()

    if any(target <= 0.0 for target in args.targets_k):
        raise SystemExit("All --targets-k values must be positive.")

    candidates = load_candidates(args.fluids)
    if not candidates:
        raise SystemExit("No candidate fluids found.")

    results = [
        evaluate(candidate, args.targets_k)
        for candidate in candidates
    ]

    reference = water_reference(results, args.water_label)

    results.sort(
        key=lambda item: (
            item.q_to_vapor_mj_kg is not None,
            item.q_to_vapor_mj_kg if item.q_to_vapor_mj_kg is not None else -math.inf,
        ),
        reverse=True,
    )

    print(
        f"CoolProp {CoolProp.__version__} | "
        f"candidates={len(results)} | "
        f"reference={reference.candidate.label}"
    )
    print()
    print_results(results, reference, args.targets_k)

    problems = [
        result
        for result in results
        if result.errors
    ]
    if problems:
        print()
        print("Notes / unsupported states:")
        for result in problems:
            print(
                f"- {result.candidate.label}: "
                + "; ".join(result.errors)
            )

    if args.output is not None:
        write_csv(args.output, results, reference, args.targets_k)
        print()
        print(f"CSV written to: {args.output}")


if __name__ == "__main__":
    main()
