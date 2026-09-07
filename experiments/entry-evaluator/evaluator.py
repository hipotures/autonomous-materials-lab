from __future__ import annotations

import csv
import json
from math import degrees, radians
from pathlib import Path
from typing import Any

from atmosphere import AtmosphereModel
from chemistry import ChemistryLimiter
from coolant import CoolantModel
from heating import total_heating
from trajectory import TrajectoryState, aero_state, rk4_step
from wall import WallState, step_wall


def _initial_state(config: dict[str, Any]) -> TrajectoryState:
    entry = config["entry"]
    return TrajectoryState(
        altitude_m=float(entry["altitude_km"]) * 1000.0,
        velocity_m_s=float(entry["velocity_km_s"]) * 1000.0,
        flight_path_angle_rad=radians(
            float(entry["flight_path_angle_deg"])
        ),
        downrange_m=0.0,
    )


def evaluate(
    config: dict[str, Any],
    history_path: Path | None = None,
) -> dict[str, Any]:
    vehicle = config["vehicle"]
    numerics = config["numerics"]
    terminal = config["terminal"]
    wall_cfg = config["wall"]
    coolant_cfg = config["coolant"]

    atmosphere = AtmosphereModel(
        config["atmosphere"],
        float(vehicle["nose_radius_m"]),
    )
    chemistry = ChemistryLimiter(config.get("chemistry", {}))
    coolant = CoolantModel(coolant_cfg, chemistry)

    state = _initial_state(config)
    wall = WallState(float(wall_cfg["initial_temperature_k"]))
    dt = float(numerics.get("dt_s", 0.05))
    max_time = float(numerics.get("max_time_s", 2000.0))
    initial_mass = float(vehicle["initial_mass_kg"])
    dry_mass = float(
        vehicle.get(
            "dry_mass_kg",
            initial_mass
            - float(coolant_cfg.get("available_mass_kg", 0.0)),
        )
    )
    available_coolant = float(
        coolant_cfg.get("available_mass_kg", 1e12)
    )
    mass = initial_mass
    coolant_used = 0.0
    time_s = 0.0
    heat_load_j_m2 = 0.0
    rad_valid_steps = 0
    heating_steps = 0
    failure_reason: str | None = None
    status = "running"

    peaks = {
        "heat_flux_w_m2": 0.0,
        "convective_heat_flux_w_m2": 0.0,
        "radiative_heat_flux_w_m2": 0.0,
        "dynamic_pressure_pa": 0.0,
        "deceleration_g": 0.0,
        "wall_temperature_k": wall.temperature_k,
        "coolant_flow_kg_s": 0.0,
        "injection_pressure_pa": 0.0,
    }
    minimum_ignition_delay_s: float | None = None
    history: list[dict[str, Any]] = []

    while time_s < max_time:
        if (
            state.altitude_m
            <= float(terminal.get("altitude_km", 20.0)) * 1000.0
        ):
            status = "terminal_altitude"
            break
        if (
            state.velocity_m_s
            <= float(terminal.get("velocity_km_s", 1.0)) * 1000.0
        ):
            status = "terminal_velocity"
            break
        if state.altitude_m < 0.0:
            status = "ground"
            break

        atm = atmosphere.sample(state.altitude_m)
        aero = aero_state(state, atm, mass, vehicle)
        heating = total_heating(
            atm,
            state.velocity_m_s,
            float(vehicle["nose_radius_m"]),
            config.get("heating", {}),
        )
        if heating.total_external_w_m2 > 0:
            heating_steps += 1
            if heating.radiation_nominal_validity:
                rad_valid_steps += 1

        provisional_wall = step_wall(
            wall,
            heating.total_external_w_m2,
            atm.temperature_k,
            dt,
            wall_cfg,
        )
        coolant_step = coolant.evaluate(
            wall.temperature_k,
            aero.surface_pressure_pa,
            provisional_wall.coolant_required_w_m2,
        )
        if (
            not coolant_step.feasible
            and provisional_wall.coolant_required_w_m2 > 0
        ):
            status = "failed"
            failure_reason = coolant_step.failure_reason
            break

        mdot = coolant_step.mass_flow_kg_s
        dm = mdot * dt
        if coolant_used + dm > available_coolant:
            status = "failed"
            failure_reason = "configured coolant mass exhausted"
            break

        coolant_used += dm
        mass = max(dry_mass, initial_mass - coolant_used)
        wall.temperature_k = provisional_wall.next_temperature_k
        heat_load_j_m2 += heating.total_external_w_m2 * dt

        if wall.temperature_k > float(
            wall_cfg["maximum_temperature_k"]
        ):
            status = "failed"
            failure_reason = "wall maximum temperature exceeded"
            break

        peaks["heat_flux_w_m2"] = max(
            peaks["heat_flux_w_m2"],
            heating.total_external_w_m2,
        )
        peaks["convective_heat_flux_w_m2"] = max(
            peaks["convective_heat_flux_w_m2"],
            heating.convective_w_m2,
        )
        peaks["radiative_heat_flux_w_m2"] = max(
            peaks["radiative_heat_flux_w_m2"],
            heating.radiative_w_m2,
        )
        peaks["dynamic_pressure_pa"] = max(
            peaks["dynamic_pressure_pa"],
            aero.dynamic_pressure_pa,
        )
        peaks["deceleration_g"] = max(
            peaks["deceleration_g"],
            aero.acceleration_g,
        )
        peaks["wall_temperature_k"] = max(
            peaks["wall_temperature_k"],
            wall.temperature_k,
        )
        peaks["coolant_flow_kg_s"] = max(
            peaks["coolant_flow_kg_s"],
            mdot,
        )
        peaks["injection_pressure_pa"] = max(
            peaks["injection_pressure_pa"],
            coolant_step.required_injection_pressure_pa,
        )

        if coolant_step.ignition_delay_s is not None:
            minimum_ignition_delay_s = (
                coolant_step.ignition_delay_s
                if minimum_ignition_delay_s is None
                else min(
                    minimum_ignition_delay_s,
                    coolant_step.ignition_delay_s,
                )
            )

        history.append(
            {
                "time_s": time_s,
                "altitude_m": state.altitude_m,
                "velocity_m_s": state.velocity_m_s,
                "flight_path_angle_deg": degrees(
                    state.flight_path_angle_rad
                ),
                "downrange_m": state.downrange_m,
                "mass_kg": mass,
                "rho_kg_m3": atm.density_kg_m3,
                "atmosphere_temperature_k": atm.temperature_k,
                "freestream_pressure_pa": atm.pressure_pa,
                "surface_pressure_pa": aero.surface_pressure_pa,
                "knudsen": atm.knudsen,
                "continuum_factor": atm.continuum_factor,
                "dynamic_pressure_pa": aero.dynamic_pressure_pa,
                "deceleration_g": aero.acceleration_g,
                "convective_heat_flux_w_m2": heating.convective_w_m2,
                "radiative_heat_flux_w_m2": heating.radiative_w_m2,
                "radiation_nominal_validity": int(
                    heating.radiation_nominal_validity
                ),
                "total_heat_flux_w_m2": heating.total_external_w_m2,
                "wall_temperature_k": wall.temperature_k,
                "wall_radiation_out_w_m2": (
                    provisional_wall.radiation_out_w_m2
                ),
                "coolant_heat_flux_w_m2": (
                    provisional_wall.coolant_required_w_m2
                ),
                "coolant_exit_temperature_k": (
                    coolant_step.exit_temperature_k
                ),
                "coolant_usable_enthalpy_j_kg": (
                    coolant_step.usable_enthalpy_j_kg
                ),
                "coolant_mass_flow_kg_s": mdot,
                "coolant_used_kg": coolant_used,
                "required_injection_pressure_pa": (
                    coolant_step.required_injection_pressure_pa
                ),
                "ignition_delay_s": coolant_step.ignition_delay_s,
            }
        )

        state = rk4_step(
            state,
            dt,
            mass,
            vehicle,
            atmosphere.sample,
        )
        time_s += dt

    if status == "running":
        status = "max_time"

    summary = {
        "case_name": config.get("name", "entry-case"),
        "status": status,
        "failure_reason": failure_reason,
        "fluid": coolant_cfg["coolprop_name"],
        "time_s": time_s,
        "final_altitude_km": state.altitude_m / 1000.0,
        "final_velocity_km_s": state.velocity_m_s / 1000.0,
        "final_flight_path_angle_deg": degrees(
            state.flight_path_angle_rad
        ),
        "downrange_km": state.downrange_m / 1000.0,
        "coolant_used_kg": coolant_used,
        "coolant_remaining_kg": max(
            0.0,
            available_coolant - coolant_used,
        ),
        "heat_load_mj_m2": heat_load_j_m2 / 1e6,
        "minimum_ignition_delay_s": minimum_ignition_delay_s,
        "radiative_correlation_valid_fraction": (
            rad_valid_steps / heating_steps
            if heating_steps
            else 0.0
        ),
        **{
            f"peak_{key}": value
            for key, value in peaks.items()
        },
    }

    if history_path is not None:
        history_path.parent.mkdir(
            parents=True,
            exist_ok=True,
        )
        if history:
            with history_path.open(
                "w",
                newline="",
                encoding="utf-8",
            ) as handle:
                writer = csv.DictWriter(
                    handle,
                    fieldnames=list(history[0].keys()),
                )
                writer.writeheader()
                writer.writerows(history)

    return summary


def write_summary(
    path: Path,
    summary: dict[str, Any],
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        json.dumps(summary, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
