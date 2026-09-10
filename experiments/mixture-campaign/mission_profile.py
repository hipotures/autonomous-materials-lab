"""Shared, fixed-trajectory heat-demand replay using the existing entry modules.

An unlimited heat-removal controller creates DEMANDS, not a successful fluid
simulation. Candidate delivery/enthalpy feasibility is checked subsequently.
No fluid-dependent aerodynamic, chemical, ablative or boundary-layer gain is
assumed. The original entry evaluator and its property providers are unchanged.
"""
from __future__ import annotations
from copy import deepcopy
from math import radians, degrees
from pathlib import Path
import sys

ENTRY = Path(__file__).resolve().parent.parent / 'entry-evaluator'


def entry_configuration(base, campaign):
    result = deepcopy(base)
    for key, value in campaign['mission'].get('overrides', {}).items():
        parts = key.split('.')
        target = result
        for part in parts[:-1]:
            if part not in target or not isinstance(target[part], dict):
                raise ValueError('unknown entry override: '+key)
            target = target[part]
        if parts[-1] not in target:
            raise ValueError('unknown entry override: '+key)
        target[parts[-1]] = value
    result['vehicle']['couple_coolant_mass_to_trajectory'] = False
    result['chemistry'] = {'mode': 'disabled'}
    result['coolant'] = {**campaign['fluid'], 'available_mass_kg': None}
    return result


class DemandCollector:
    def __init__(self, config):
        self.config = config
        self.rows = []
        self.time_s, self.dt_s = 0., 0.

    def evaluate(self, wall_temperature_k, surface_pressure_pa, coolant_heat_flux_w_m2, *, area_m2=None):
        from coolant import CoolantStep
        f = self.config
        temperature = max(f['storage_temperature_k']+1.,
                          min(f['max_exit_temperature_k'], wall_temperature_k-f['wall_to_fluid_approach_k']))
        if coolant_heat_flux_w_m2 > 0:
            # The 100 Pa evaluation floor is inherited explicitly from CoolantModel.
            self.rows.append({'time_s': self.time_s, 'dt_s': self.dt_s,
                'temperature_k': temperature, 'pressure_pa': max(float(surface_pressure_pa), 100.),
                'surface_pressure_pa': float(surface_pressure_pa), 'wall_temperature_k': wall_temperature_k,
                'heat_flux_w_m2': float(coolant_heat_flux_w_m2), 'area_m2': float(area_m2),
                'heat_j': float(coolant_heat_flux_w_m2)*float(area_m2)*self.dt_s})
        return CoolantStep(temperature, 0., 0., 0., 0., 'demand_only_no_fluid_model', None, True, None)


def build_profile(entry, campaign):
    if str(ENTRY) not in sys.path:
        sys.path.insert(0, str(ENTRY))
    from atmosphere import AtmosphereModel
    from trajectory import TrajectoryState, aero_state, rk4_step
    from heating import total_heating
    from surface import Forebody
    cfg = entry_configuration(entry, campaign)
    vehicle, wall = cfg['vehicle'], cfg['wall']
    dt, maximum = float(cfg['numerics']['dt_s']), float(cfg['numerics']['max_time_s'])
    if dt <= 0 or maximum <= 0:
        raise ValueError('invalid mission integration limits')
    atmosphere = AtmosphereModel(cfg['atmosphere'], float(vehicle['nose_radius_m']))
    state = TrajectoryState(float(cfg['entry']['altitude_km'])*1000.,
        float(cfg['entry']['velocity_km_s'])*1000., radians(cfg['entry']['flight_path_angle_deg']), 0.)
    start_altitude, descended = state.altitude_m, False
    surface = Forebody(cfg.get('surface', {}), cfg['coolant']['cooled_area_m2'], wall)
    collector = DemandCollector(cfg['coolant'])
    t, status, reason = 0., 'max_time', None
    incident = cooling = residual = conv = rad = conv_valid = rad_valid = 0.
    peaks = {'wall_temperature_k': wall['initial_temperature_k'], 'external_heat_flux_w_m2': 0.,
             'surface_pressure_pa': 0., 'deceleration_g': 0.}
    history = []
    while t < maximum:
        step_dt = min(dt, maximum-t)
        if state.altitude_m < start_altitude-1000:
            descended = True
        if descended and state.altitude_m >= start_altitude and state.flight_path_angle_rad > 0:
            status = 'atmospheric_exit'; break
        if state.altitude_m <= cfg['terminal']['altitude_km']*1000:
            status = 'terminal_altitude'; break
        if state.velocity_m_s <= cfg['terminal']['velocity_km_s']*1000:
            status = 'terminal_velocity'; break
        if state.altitude_m < 0:
            status = 'ground'; break
        atm = atmosphere.sample(state.altitude_m)
        aero = aero_state(state, atm, vehicle['initial_mass_kg'], vehicle)
        heat = total_heating(atm, state.velocity_m_s, vehicle['nose_radius_m'], cfg['heating'])
        collector.time_s, collector.dt_s = t, step_dt
        result = surface.evaluate(heat, atm.temperature_k, aero.surface_pressure_pa, step_dt, collector)
        meta = surface.metadata()
        dc = heat.convective_w_m2*meta['surface_convective_area_factor']*meta['cooled_surface_area_m2']*step_dt
        dr = heat.radiative_w_m2*meta['surface_radiative_area_factor']*meta['cooled_surface_area_m2']*step_dt
        conv += dc; rad += dr
        conv_valid += dc*bool(heat.convective_nominal_validity)
        rad_valid += dr*bool(heat.radiation_nominal_validity)
        incident += result.incident_power_w*step_dt
        cooling += result.coolant_power_w*step_dt
        residual += abs(result.energy_residual_w)*step_dt
        peaks['wall_temperature_k'] = max(peaks['wall_temperature_k'], result.peak_temperature_k)
        peaks['external_heat_flux_w_m2'] = max(peaks['external_heat_flux_w_m2'], heat.total_external_w_m2)
        peaks['surface_pressure_pa'] = max(peaks['surface_pressure_pa'], aero.surface_pressure_pa)
        peaks['deceleration_g'] = max(peaks['deceleration_g'], aero.acceleration_g)
        history.append({'time_s': t, 'dt_s': step_dt, 'altitude_m': state.altitude_m,
            'velocity_m_s': state.velocity_m_s, 'surface_pressure_pa': aero.surface_pressure_pa,
            'external_heat_flux_w_m2': heat.total_external_w_m2, 'wall_temperature_k': result.peak_temperature_k,
            'required_cooling_power_w': result.coolant_power_w})
        surface.commit(result)
        if result.peak_temperature_k > wall['maximum_temperature_k']:
            status, reason = 'failed', 'wall_temperature_limit'; break
        state = rk4_step(state, step_dt, vehicle['initial_mass_kg'], vehicle, atmosphere.sample)
        t += step_dt
    total = sum(r['heat_j'] for r in collector.rows)
    if abs(total-cooling) > max(1., abs(cooling)*1e-8):
        raise ValueError('recorded zone demand does not close the surface heat balance')
    return {'schema': 'entry-demand-profile-v1', 'entry_config': cfg, 'demands': collector.rows,
        'history': history, 'summary': {'status': status, 'failure_reason': reason, 'time_s': t,
        'terminal_altitude_km': state.altitude_m/1000, 'terminal_velocity_km_s': state.velocity_m_s/1000,
        'terminal_angle_deg': degrees(state.flight_path_angle_rad), 'required_heat_j': total,
        'incident_heat_j': incident, 'energy_relative_residual': residual/max(incident, 1.),
        'convective_energy_valid_fraction': conv_valid/conv if conv else None,
        'radiative_energy_valid_fraction': rad_valid/rad if rad else None,
        'demand_rows': len(collector.rows), 'peaks': peaks,
        'coolant_temperature_range_k': [min((r['temperature_k'] for r in collector.rows), default=None),
                                        max((r['temperature_k'] for r in collector.rows), default=None)],
        'active_pressure_range_pa': [min((r['pressure_pa'] for r in collector.rows), default=None),
                                     max((r['pressure_pa'] for r in collector.rows), default=None)],
        'fluid_independent_profile': True, 'unlimited_heat_removal_assumed_for_profile': True,
        'coolant_mass_kg': None, 'ablative_wall_simulated': False, 'injection_heat_reduction_assumed': 0.,
        'frozen_trajectory_mass_kg': vehicle['initial_mass_kg'], 'hot_boundary_layer_temperature_k': None,
        'profile_complete': status in {'terminal_altitude', 'terminal_velocity'} and total > 0}}
