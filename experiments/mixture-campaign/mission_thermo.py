"""Bounded liquid flash plus an explicitly provisional frozen-species gas branch.

No NRTL/UNIFAC liquid correlation is extrapolated to 900 or 3000 K. Gas sensible
enthalpy uses bounded ideal-gas Cp correlations. Continuing a flash-confirmed
vapor branch to a higher temperature assumes no reentrant condensation and no
chemical reaction; that assumption is recorded, not certified as real cooling.
"""
from __future__ import annotations
from functools import lru_cache
import math
from campaign_backend import legacy_import
from mission_core import finite


def heos_states(nodes, fluid):
    import CoolProp.CoolProp as CP
    tin, pin = fluid['storage_temperature_k'], fluid['storage_pressure_pa']
    if CP.PhaseSI('T', tin, 'P', pin, 'Water') != 'liquid':
        raise ValueError('HEOS water storage is not liquid')
    h0 = CP.PropsSI('Hmass', 'T', tin, 'P', pin, 'Water')
    rho = CP.PropsSI('Dmass', 'T', tin, 'P', pin, 'Water')
    output = []
    for n in nodes:
        r = {'temperature_k': n['temperature_k'], 'pressure_pa': n['pressure_pa']}
        try:
            h = CP.PropsSI('Hmass', 'T', r['temperature_k'], 'P', r['pressure_pa'], 'Water')
            r.update(status='ok', delta_h_j_kg=finite(float(h-h0), 'HEOS delta h'),
                storage_density_kg_m3=float(rho), provider='CoolProp HEOS Water')
        except (ValueError, RuntimeError, ArithmeticError) as exc:
            r.update(status='unsupported_water_state', reason=str(exc))
        output.append(r)
    return output


class MissionModel:
    def __init__(self, pair, model, config):
        from thermo import IdealGas, GibbsExcessLiquid, FlashVLN
        self.legacy, self.config = legacy_import(), config
        self.pair, self.model = pair, model
        f, th = config['fluid'], config['thermodynamics']
        # Select liquid correlations over a fixed envelope, never over a shortlist.
        cfg = {'grid': {'inlet_temperature_k': f['storage_temperature_k'],
                       'outlet_temperatures_k': th['liquid_selection_temperatures_k']},
               'phase_fraction_tolerance': 1e-8, 'balance_tolerance': 1e-7,
               'endpoint_relative_tolerance': th['endpoint_relative_tolerance']}
        self.base = self.legacy.BinaryModel(pair, model, cfg)
        self.base.methods['HeatCapacityGases'] = [self.legacy.pin_correlation(o,
            [298.15, f['storage_temperature_k'], f['max_exit_temperature_k']]) for o in self.base.props.HeatCapacityGases]
        for obj in self.base.props.HeatCapacityGases:
            for t in (298.15, f['storage_temperature_k'], f['max_exit_temperature_k']):
                self.legacy.in_bounds(obj, t)
        self.base.liquid = GibbsExcessLiquid(VaporPressures=self.base.props.VaporPressures,
            HeatCapacityGases=self.base.props.HeatCapacityGases, VolumeLiquids=self.base.props.VolumeLiquids,
            GibbsExcessModel=self.base.ge, equilibrium_basis='Psat', caloric_basis='Psat')
        self.gas = IdealGas(HeatCapacityGases=self.base.props.HeatCapacityGases)
        self.base.flasher = FlashVLN(self.base.constants, self.base.props,
            liquids=[self.base.liquid, self.base.liquid], gas=self.gas)
        self.water_liquid = GibbsExcessLiquid(VaporPressures=[self.base.props.VaporPressures[0]],
            HeatCapacityGases=[self.base.props.HeatCapacityGases[0]], VolumeLiquids=[self.base.props.VolumeLiquids[0]],
            equilibrium_basis='Psat', caloric_basis='Psat')
        self.water_gas = IdealGas(HeatCapacityGases=[self.base.props.HeatCapacityGases[0]])
        self.certificates = {}
        self.metadata = {'pure_methods': self.base.methods, 'model': model,
            'caloric_basis': 'Psat', 'gas_basis': 'bounded_ideal_gas_fixed_species',
            'gas_phase_continuation_assumption': 'no_reentrant_condensation_no_reaction',
            'chemistry_validated': False, 'caloric_accuracy_validated': False,
            'mixture_parameter_sha256': pair['parameter_sha256']}

    def certify_vapor(self, w, p):
        if p > self.config['thermodynamics']['maximum_flash_pressure_pa']:
            return
        known = self.certificates.get(w)
        if known and p <= known['pressure_pa']:
            return
        for t in sorted(self.config['thermodynamics']['bridge_temperatures_k'], reverse=True):
            try:
                state = self.base.state(w, t, p)
                if state['liquid_phase_count'] == 0 and state['vapor_mole_fraction'] >= 1-1e-8:
                    self.certificates[w] = {'temperature_k': t, 'pressure_pa': p,
                                           'source': 'bounded_nonreactive_flash'}
                    return
            except (ValueError, RuntimeError, ArithmeticError, TypeError):
                pass

    def gas_state(self, w, t, p, bridge):
        from v5m2_core import mass_to_mole
        for obj in self.base.props.HeatCapacityGases:
            self.legacy.in_bounds(obj, t)
        zs = mass_to_mole(w, self.base.mw)
        molar_mass = sum(z*m for z, m in zip(zs, self.base.mw))/1000
        phase = self.gas.to_TP_zs(T=t, P=p, zs=zs)
        h = float(phase.H())/molar_mass
        if not math.isfinite(h):
            raise ValueError('nonfinite gas enthalpy')
        return {'enthalpy_j_kg': h, 'vapor_mole_fraction': 1., 'liquid_phase_count': 0,
            'gas_branch': 'frozen_species_continuation', 'phase_bridge': bridge}

    @lru_cache(maxsize=65536)
    def state(self, w, t, p):
        if t > self.config['fluid']['max_exit_temperature_k'] or p > self.config['thermodynamics']['maximum_flash_pressure_pa']:
            raise ValueError('outside_declared_mission_thermo_domain')
        bridge = self.certificates.get(w)
        if bridge and t >= bridge['temperature_k'] and p <= bridge['pressure_pa']:
            return self.gas_state(w, t, p, bridge)
        try:
            return {**self.base.state(w, t, p), 'gas_branch': 'bounded_flash'}
        except (ValueError, RuntimeError, ArithmeticError, TypeError):
            self.certify_vapor(w, p)
            bridge = self.certificates.get(w)
            if bridge and t >= bridge['temperature_k'] and p <= bridge['pressure_pa']:
                return self.gas_state(w, t, p, bridge)
            raise ValueError('no_bounded_flash_or_vapor_bridge_for_requested_state')

    @lru_cache(maxsize=32768)
    def water_enthalpy(self, t, p):
        obj = self.base.props.VaporPressures[0]
        tc = self.base.constants.Tcs[0]
        self.legacy.in_bounds(self.base.props.HeatCapacityGases[0], t)
        if t < tc:
            self.legacy.in_bounds(obj, t)
            psat = float(obj.calculate(t, obj.method))
            if abs(p/psat-1) < 1e-7:
                raise ValueError('pure_water_saturation_quality_ambiguous')
            liquid = p > psat
        else:
            liquid = False
        if liquid:
            self.legacy.in_bounds(self.base.props.VolumeLiquids[0], t)
        phase = (self.water_liquid if liquid else self.water_gas).to_TP_zs(T=t, P=p, zs=[1.])
        return float(phase.H())/(self.base.mw[0]/1000)

    def evaluate(self, w, nodes):
        f = self.config['fluid']
        output = []
        try:
            before = self.base.state(w, f['storage_temperature_k'], f['storage_pressure_pa'])
            if before['liquid_phase_count'] != 1 or before['vapor_mole_fraction'] > 1e-8:
                raise ValueError('storage_not_single_liquid')
            self.certify_vapor(w, max(n['pressure_pa'] for n in nodes))
            water_before = self.water_enthalpy(f['storage_temperature_k'], f['storage_pressure_pa'])
            water = {(r['temperature_k'], r['pressure_pa']): r for r in heos_states(nodes, f)}
        except (ValueError, RuntimeError, ArithmeticError, TypeError) as exc:
            return [{'temperature_k': n['temperature_k'], 'pressure_pa': n['pressure_pa'],
                     'status': 'storage_or_reference_unavailable', 'reason': str(exc)} for n in nodes]
        for node in nodes:
            t, p = node['temperature_k'], node['pressure_pa']
            r = {'temperature_k': t, 'pressure_pa': p, 'model': self.model,
                 'mass_fraction': w, 'chemistry_validated': False}
            try:
                after = self.state(w, t, p)
                if after['liquid_phase_count'] > 1:
                    raise ValueError('outlet_liquid_liquid_split')
                dh = finite(after['enthalpy_j_kg']-before['enthalpy_j_kg'], 'mixture delta h')
                ref = water[(t, p)]
                if ref['status'] != 'ok':
                    raise ValueError('HEOS_reference_unavailable')
                dhw = self.water_enthalpy(t, p)-water_before
                error = abs(dhw/ref['delta_h_j_kg']-1)
                r.update(delta_h_j_kg=dh, storage_density_kg_m3=before['density_kg_m3'],
                    pure_water_model_delta_h_j_kg=dhw, heos_water_delta_h_j_kg=ref['delta_h_j_kg'],
                    water_endpoint_relative_error=error, vapor_mole_fraction=after['vapor_mole_fraction'],
                    gas_branch=after['gas_branch'], phase_bridge=after.get('phase_bridge'))
                if error > self.config['thermodynamics']['endpoint_relative_tolerance']:
                    raise ValueError('water_caloric_reference_gate')
                r['status'] = 'ok'
            except (ValueError, RuntimeError, ArithmeticError, TypeError) as exc:
                r.update(status='unsupported_or_blocked_state', reason=str(exc))
            output.append(r)
        return output
