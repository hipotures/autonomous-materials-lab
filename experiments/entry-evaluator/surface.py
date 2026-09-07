"""Area-conserving angular forebody surrogate; not a calibrated geometry model."""
from __future__ import annotations

from dataclasses import dataclass
from math import cos, isfinite, radians

from coolant import CoolantModel, CoolantStep
from heating import HeatingState
from wall import WallState, WallStep, step_wall


@dataclass(frozen=True)
class SurfaceZone:
    area_m2: float
    convective_factor: float
    radiative_factor: float


def build_zones(config: dict, area_m2: float) -> list[SurfaceZone]:
    if not isfinite(area_m2) or area_m2 <= 0:
        raise ValueError("coolant.cooled_area_m2 must be finite and positive")
    model = config.get("model", "uniform")
    if model == "uniform":
        return [SurfaceZone(area_m2, 1.0, 1.0)]
    if model != "cosine_cap":
        raise ValueError(f"unsupported surface.model: {model}")
    count = config.get("zones", 32)
    if isinstance(count, bool) or not isinstance(count, int) or not 1 <= count <= 4096:
        raise ValueError("surface.zones must be an integer in [1, 4096]")
    theta = float(config.get("theta_max_deg", 90.0))
    if not isfinite(theta) or not 0 < theta <= 90:
        raise ValueError("surface.theta_max_deg must be in (0, 90]")
    exponents = [float(config.get(key, 1.0)) for key in
                 ("convective_exponent", "radiative_exponent")]
    if any(not isfinite(n) or n < 0 or n > 20 for n in exponents):
        raise ValueError("surface exponents must be finite and in [0, 20]")
    mu_min = max(0.0, cos(radians(theta)))
    width = (1.0 - mu_min) / count
    if width <= 1e-14:
        raise ValueError("surface angular span is too small to resolve")
    zones = []
    for i in range(count):
        lo, hi = mu_min + i * width, mu_min + (i + 1) * width
        factors = [(hi ** (n + 1) - lo ** (n + 1)) / ((n + 1) * width)
                   for n in exponents]
        zones.append(SurfaceZone(area_m2 / count, *factors))
    return zones


@dataclass(frozen=True)
class SurfaceStep:
    stagnation_wall: WallStep
    stagnation_coolant: CoolantStep
    zone_walls: list[WallStep]
    incident_power_w: float
    coolant_power_w: float
    radiation_out_power_w: float
    storage_power_w: float
    mass_flow_kg_s: float
    peak_mass_flux_kg_m2_s: float
    peak_temperature_k: float
    energy_residual_w: float
    failure_reason: str | None


class Forebody:
    def __init__(self, config: dict, area_m2: float, wall_config: dict):
        self.config = config
        self.zones = build_zones(config, area_m2)
        self.area_m2 = area_m2
        self.wall_config = wall_config
        initial = float(wall_config["initial_temperature_k"])
        self.stagnation = WallState(initial)
        self.walls = [WallState(initial) for _ in self.zones]

    def evaluate(self, heating: HeatingState, environment_k: float,
                 surface_pressure_pa: float, dt: float, coolant: CoolantModel) -> SurfaceStep:
        # A zero-area stagnation probe checks local feasibility without adding mass.
        stag = step_wall(self.stagnation, heating.total_external_w_m2,
                         environment_k, dt, self.wall_config)
        stag_coolant = coolant.evaluate(self.stagnation.temperature_k, surface_pressure_pa,
                                       stag.coolant_required_w_m2, area_m2=0.0)
        failure = stag_coolant.failure_reason if not stag_coolant.feasible else None
        wall_steps = []
        incident = cooling = radiation = storage = mdot = 0.0
        max_flux = stag_coolant.mass_flux_kg_m2_s
        # V1 holds V0's stagnation pressure approximation over all rings.
        # No pressure distribution or lateral conduction is introduced here.
        for zone, wall in zip(self.zones, self.walls):
            q = (heating.convective_w_m2 * zone.convective_factor
                 + heating.radiative_w_m2 * zone.radiative_factor)
            step = step_wall(wall, q, environment_k, dt, self.wall_config)
            fluid = coolant.evaluate(wall.temperature_k, surface_pressure_pa,
                                     step.coolant_required_w_m2, area_m2=zone.area_m2)
            if not fluid.feasible and failure is None:
                failure = fluid.failure_reason
            incident += q * zone.area_m2
            cooling += step.coolant_required_w_m2 * zone.area_m2
            radiation += step.radiation_out_w_m2 * zone.area_m2
            storage += (step.next_temperature_k - wall.temperature_k) * float(
                self.wall_config["areal_heat_capacity_j_m2_k"]) * zone.area_m2 / dt
            mdot += fluid.mass_flow_kg_s
            max_flux = max(max_flux, fluid.mass_flux_kg_m2_s)
            wall_steps.append(step)
        backface = float(self.wall_config.get("backface_loss_w_m2", 0.0)) * self.area_m2
        peak = max(stag.next_temperature_k, *(w.next_temperature_k for w in wall_steps))
        return SurfaceStep(stag, stag_coolant, wall_steps, incident, cooling, radiation,
                           storage, mdot, max_flux, peak,
                           incident - cooling - radiation - backface - storage, failure)

    def commit(self, result: SurfaceStep) -> None:
        self.stagnation.temperature_k = result.stagnation_wall.next_temperature_k
        for wall, step in zip(self.walls, result.zone_walls):
            wall.temperature_k = step.next_temperature_k

    def metadata(self) -> dict:
        return {
            "surface_model": self.config.get("model", "uniform"),
            "surface_zones": len(self.zones),
            "cooled_surface_area_m2": self.area_m2,
            "surface_convective_area_factor": sum(
                z.area_m2 * z.convective_factor for z in self.zones) / self.area_m2,
            "surface_radiative_area_factor": sum(
                z.area_m2 * z.radiative_factor for z in self.zones) / self.area_m2,
        }
