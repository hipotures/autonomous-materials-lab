from __future__ import annotations

import csv
import json
from math import degrees, radians, isfinite
from pathlib import Path
from typing import Any

from atmosphere import AtmosphereModel
from chemistry import ChemistryLimiter
from coolant import CoolantModel
from heating import total_heating
from trajectory import TrajectoryState, aero_state, rk4_step
from surface import Forebody


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
    surface = Forebody(config.get("surface", {}), float(coolant_cfg["cooled_area_m2"]), wall_cfg)
    wall = surface.stagnation
    dt = float(numerics.get("dt_s", 0.05))
    max_time = float(numerics.get("max_time_s", 2000.0))
    if not isfinite(dt) or dt <= 0 or not isfinite(max_time) or max_time <= 0:
        raise ValueError("numerics dt_s and max_time_s must be finite and positive")
    initial_mass = float(vehicle["initial_mass_kg"])
    available_raw = coolant_cfg.get("available_mass_kg")
    available_coolant = (
        None
        if available_raw is None
        else float(available_raw)
    )
    dry_mass_raw = vehicle.get("dry_mass_kg")
    if dry_mass_raw is not None:
        dry_mass = float(dry_mass_raw)
    elif available_coolant is not None:
        dry_mass = initial_mass - available_coolant
    else:
        dry_mass = initial_mass
    couple_mass = bool(
        vehicle.get("couple_coolant_mass_to_trajectory", False)
    )
    mass = initial_mass
    coolant_used = 0.0
    time_s = 0.0
    heat_load_j_m2 = 0.0
    convective_heat_load_j_m2 = 0.0
    radiative_heat_load_j_m2 = 0.0
    radiative_valid_heat_load_j_m2 = 0.0
    incident_energy_j = coolant_energy_j = 0.0
    vehicle_convective_energy_j = 0.0
    vehicle_radiative_energy_j = 0.0
    vehicle_radiative_valid_energy_j = 0.0
    energy_residual_abs_j = 0.0
    rarefied_time_s = 0.0
    initial_altitude_m = state.altitude_m
    minimum_altitude_m = state.altitude_m
    minimum_knudsen = float("inf")
    has_descended = False
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
        "vehicle_heating_power_w": 0.0,
        "vehicle_coolant_power_w": 0.0,
        "coolant_mass_flux_kg_m2_s": 0.0,
    }
    minimum_ignition_delay_s: float | None = None
    history: list[dict[str, Any]] = []

    while time_s < max_time:
        step_dt = min(dt, max_time - time_s)
        minimum_altitude_m = min(minimum_altitude_m, state.altitude_m)
        if state.altitude_m < initial_altitude_m - 1000.0:
            has_descended = True
        if (
            has_descended
            and state.altitude_m >= initial_altitude_m
            and state.flight_path_angle_rad > 0.0
        ):
            status = "atmospheric_exit"
            break
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
        minimum_knudsen = min(minimum_knudsen, atm.knudsen)
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

        surface_step = surface.evaluate(
            heating, atm.temperature_k, aero.surface_pressure_pa, step_dt, coolant,
        )
        provisional_wall = surface_step.stagnation_wall
        coolant_step = surface_step.stagnation_coolant
        if surface_step.failure_reason is not None:
            status = "failed"
            failure_reason = surface_step.failure_reason
            break

        mdot = surface_step.mass_flow_kg_s
        dm = mdot * step_dt
        if (
            available_coolant is not None
            and coolant_used + dm > available_coolant
        ):
            status = "failed"
            failure_reason = "configured coolant mass exhausted"
            break

        coolant_used += dm
        if couple_mass:
            mass = max(dry_mass, initial_mass - coolant_used)
        else:
            mass = initial_mass
        surface.commit(surface_step)
        heat_load_j_m2 += heating.total_external_w_m2 * step_dt
        convective_heat_load_j_m2 += heating.convective_w_m2 * step_dt
        radiative_step_j_m2 = heating.radiative_w_m2 * step_dt
        radiative_heat_load_j_m2 += radiative_step_j_m2
        if heating.radiation_nominal_validity:
            radiative_valid_heat_load_j_m2 += radiative_step_j_m2

        incident_energy_j += surface_step.incident_power_w * step_dt
        coolant_energy_j += surface_step.coolant_power_w * step_dt
        surface_meta = surface.metadata()
        convective_vehicle_step_j = (
            heating.convective_w_m2
            * surface_meta["surface_convective_area_factor"]
            * surface_meta["cooled_surface_area_m2"]
            * step_dt
        )
        radiative_vehicle_step_j = (
            heating.radiative_w_m2
            * surface_meta["surface_radiative_area_factor"]
            * surface_meta["cooled_surface_area_m2"]
            * step_dt
        )
        vehicle_convective_energy_j += convective_vehicle_step_j
        vehicle_radiative_energy_j += radiative_vehicle_step_j
        if heating.radiation_nominal_validity:
            vehicle_radiative_valid_energy_j += radiative_vehicle_step_j
        energy_residual_abs_j += abs(surface_step.energy_residual_w) * step_dt
        if atm.knudsen >= 0.1:
            rarefied_time_s += step_dt

        if surface_step.peak_temperature_k > float(
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
            surface_step.peak_temperature_k,
        )
        peaks["coolant_flow_kg_s"] = max(
            peaks["coolant_flow_kg_s"],
            mdot,
        )
        peaks["injection_pressure_pa"] = max(
            peaks["injection_pressure_pa"],
            coolant_step.required_injection_pressure_pa,
        )

        peaks["vehicle_heating_power_w"] = max(peaks["vehicle_heating_power_w"], surface_step.incident_power_w)
        peaks["vehicle_coolant_power_w"] = max(peaks["vehicle_coolant_power_w"], surface_step.coolant_power_w)
        peaks["coolant_mass_flux_kg_m2_s"] = max(peaks["coolant_mass_flux_kg_m2_s"], surface_step.peak_mass_flux_kg_m2_s)

        if coolant_step.ignition_delay_s is not None:
            minimum_ignition_delay_s = (
                coolant_step.ignition_delay_s
                if minimum_ignition_delay_s is None
                else min(
                    minimum_ignition_delay_s,
                    coolant_step.ignition_delay_s,
                )
            )

        if history_path is not None:
            history.append(
                {
                    "time_s": time_s,
                    "step_end_time_s": time_s + step_dt,
                    "vehicle_heating_power_w": surface_step.incident_power_w,
                    "vehicle_coolant_power_w": surface_step.coolant_power_w,
                    "vehicle_reradiation_power_w": surface_step.radiation_out_power_w,
                    "vehicle_storage_power_w": surface_step.storage_power_w,
                    "wall_energy_residual_w": surface_step.energy_residual_w,
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
            step_dt,
            mass,
            vehicle,
            atmosphere.sample,
        )
        time_s += step_dt
        minimum_altitude_m = min(minimum_altitude_m, state.altitude_m)

    if status == "running":
        status = "max_time"

    summary = {
        "evaluator_version": "v1",
        "dt_s": dt,
        **surface.metadata(),
        "vehicle_incident_heat_mj": incident_energy_j / 1e6,
        "vehicle_convective_incident_heat_mj": vehicle_convective_energy_j / 1e6,
        "vehicle_radiative_incident_heat_mj": vehicle_radiative_energy_j / 1e6,
        "vehicle_radiative_valid_heat_mj": vehicle_radiative_valid_energy_j / 1e6,
        "vehicle_radiative_energy_valid_fraction": (
            vehicle_radiative_valid_energy_j / vehicle_radiative_energy_j
            if vehicle_radiative_energy_j > 0.0
            else None
        ),
        "vehicle_radiative_energy_fraction": (
            vehicle_radiative_energy_j / incident_energy_j
            if incident_energy_j > 0.0
            else None
        ),
        "vehicle_coolant_heat_mj": coolant_energy_j / 1e6,
        "wall_energy_residual_abs_j": energy_residual_abs_j,
        "wall_energy_relative_residual": energy_residual_abs_j / max(incident_energy_j, 1.0),
        "rarefied_heating_disabled_time_s": rarefied_time_s,
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
        "coolant_remaining_kg": (
            None
            if available_coolant is None
            else max(0.0, available_coolant - coolant_used)
        ),
        "trajectory_mass_coupled_to_coolant": couple_mass,
        "minimum_altitude_km": minimum_altitude_m / 1000.0,
        "minimum_knudsen": (
            None
            if minimum_knudsen == float("inf")
            else minimum_knudsen
        ),
        "heat_load_mj_m2": heat_load_j_m2 / 1e6,
        "convective_heat_load_mj_m2": convective_heat_load_j_m2 / 1e6,
        "radiative_heat_load_mj_m2": radiative_heat_load_j_m2 / 1e6,
        "radiative_valid_heat_load_mj_m2": radiative_valid_heat_load_j_m2 / 1e6,
        "radiative_energy_valid_fraction": (
            radiative_valid_heat_load_j_m2 / radiative_heat_load_j_m2
            if radiative_heat_load_j_m2 > 0.0
            else None
        ),
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
