"""Model-independent design, quadrature, complete-profile scoring and selection."""
from __future__ import annotations
from collections import Counter, defaultdict
import math
from campaign_store import digest

SCHEMA = 'mixture-mission-campaign-v1'


def finite(value, name, minimum=0.0, strict=True):
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise ValueError(f'{name} must be a finite number')
    if not math.isfinite(value) or (value <= minimum if strict else value < minimum):
        raise ValueError(f'{name} outside its permitted range')
    return float(value)


def validate(c):
    if c.get('schema') != SCHEMA:
        raise ValueError('unsupported mission campaign schema')
    for section, keys in {'execution': ['workers', 'timeout_s', 'batch_compositions'],
            'screening': ['coarse_nodes', 'refined_nodes', 'detail_pairs', 'per_pair', 'exact_candidates'],
            'composition': ['log_points', 'linear_points']}.items():
        for key in keys:
            n = c[section][key]
            if isinstance(n, bool) or not isinstance(n, int) or n < 1:
                raise ValueError(f'{section}.{key} must be a positive integer')
    d = c['composition']
    if d['basis'] != 'mass' or not 0 < d['minimum'] < d['transition'] < d['maximum'] < 1:
        raise ValueError('invalid mass-composition design')
    if min(d['log_points'], d['linear_points']) < 2:
        raise ValueError('composition axes need at least two points')
    maximum = c['catalog']['maximum_pairs']
    if maximum is not None and (isinstance(maximum, bool) or not isinstance(maximum, int) or maximum < 1):
        raise ValueError('maximum_pairs must be null or a positive integer')
    if c['screening']['refined_nodes'] < c['screening']['coarse_nodes']:
        raise ValueError('refined_nodes must not reduce resolution')
    for key in ('relative_mass_tolerance', 'gain_margin'):
        finite(c['screening'][key], key, strict=False)
    f = c['fluid']
    for key in ('storage_temperature_k', 'storage_pressure_pa', 'max_exit_temperature_k',
                'wall_to_fluid_approach_k', 'max_injection_pressure_pa', 'injection_pressure_margin',
                'pump_efficiency', 'cooled_area_m2'):
        finite(f[key], key)
    finite(f['porous_delta_p_pa'], 'porous pressure loss', strict=False)
    if f['pump_efficiency'] > 1 or f['max_exit_temperature_k'] <= f['storage_temperature_k']:
        raise ValueError('invalid pump efficiency or exit temperature')
    if f.get('maximum_mass_flux_kg_m2_s') is not None:
        finite(f['maximum_mass_flux_kg_m2_s'], 'maximum_mass_flux')
    for key in ('maximum_flash_pressure_pa', 'endpoint_relative_tolerance'):
        finite(c['thermodynamics'][key], key)
    if not isinstance(c['thermodynamics']['bridge_temperatures_k'], list) or not c['thermodynamics']['bridge_temperatures_k']:
        raise ValueError('a fixed gas-bridge temperature protocol is required')
    for t in c['thermodynamics']['bridge_temperatures_k']:
        finite(t, 'bridge temperature')
    if not c['mission']['entry_config'] or c['publication']['shard_bytes'] < 4096:
        raise ValueError('missing mission or undersized report shards')


def compositions(d):
    a, b = d['log_points'], d['linear_points']
    lo, mid, hi = d['minimum'], d['transition'], d['maximum']
    values = [math.exp(math.log(lo)+(math.log(mid)-math.log(lo))*i/(a-1)) for i in range(a)]
    values += [mid+(hi-mid)*i/(b-1) for i in range(b)]
    return sorted({round(v, 12) for v in values})


def select_catalog(rows, maximum, seed):
    # No preference for NRTL membership, predicted merit, CAS ordering or names.
    ordered = sorted(rows, key=lambda r: (digest([seed, r['inchi_key']]), r['pair_id']))
    return sorted(ordered if maximum is None else ordered[:maximum], key=lambda r: r['pair_id'])


def coord(r):
    return float(r['temperature_k']), float(r['pressure_pa'])


def exact_nodes(profile):
    groups = {}
    for r in profile['demands']:
        q = coord(r)
        g = groups.setdefault(q, {'temperature_k': q[0], 'pressure_pa': q[1], 'heat_j': 0., 'peak_flux_w_m2': 0.})
        g['heat_j'] += finite(r['heat_j'], 'heat', strict=False)
        g['peak_flux_w_m2'] = max(g['peak_flux_w_m2'], r['heat_flux_w_m2'])
    return [groups[q] for q in sorted(groups)]


def quadrature(profile, count):
    """Heat-preserving medoid bins. Passing a bin does not validate unsampled T/P."""
    rows = exact_nodes(profile)
    total = math.fsum(r['heat_j'] for r in rows)
    if not rows or total <= 0:
        raise ValueError('profile contains no positive cooling demand')
    if count >= len(rows):
        return {'nodes': rows, 'exact': True, 'total_heat_j': total, 'source_nodes': len(rows)}
    bins = defaultdict(list)
    accumulated = 0.
    for r in rows:
        index = min(count-1, int((accumulated+0.5*r['heat_j'])/total*count))
        bins[index].append(r)
        accumulated += r['heat_j']
    output = []
    for chunk in bins.values():
        energy = math.fsum(r['heat_j'] for r in chunk)
        running, chosen = 0., chunk[-1]
        for r in chunk:
            running += r['heat_j']
            if running >= energy/2:
                chosen = r; break
        output.append({**chosen, 'heat_j': energy,
            'peak_flux_w_m2': max(r['peak_flux_w_m2'] for r in chunk),
            'source_point_count': len(chunk), 'source_temperature_range_k': [min(r['temperature_k'] for r in chunk), max(r['temperature_k'] for r in chunk)],
            'source_pressure_range_pa': [min(r['pressure_pa'] for r in chunk), max(r['pressure_pa'] for r in chunk)]})
    return {'nodes': output, 'exact': False, 'total_heat_j': total, 'source_nodes': len(rows)}


def injection_pressure(p, fluid):
    return fluid['injection_pressure_margin']*p + fluid['porous_delta_p_pa']


def score(nodes, states, fluid, *, profile=None):
    """Integrate delivered heat; pump work is not credited as wall heat uptake."""
    expected = {coord(n) for n in nodes}
    if len(states) != len(expected) or {coord(s) for s in states} != expected:
        raise ValueError('missing, duplicate or mismatched thermodynamic states')
    by_q = {coord(s): s for s in states}
    mass, covered, pump_j, failures, effective = 0., 0., 0., [], {}
    total = math.fsum(n['heat_j'] for n in nodes)
    if total <= 0:
        raise ValueError('zero heat is not a valid mission comparison')
    max_flux, max_pressure = 0., 0.
    for node in nodes:
        q, state = coord(node), by_q[coord(node)]
        try:
            if state.get('status') != 'ok':
                raise ValueError(state.get('reason', state.get('status', 'missing_state')))
            dh = finite(state['delta_h_j_kg'], 'enthalpy window')
            rho = finite(state['storage_density_kg_m3'], 'storage density')
            required = injection_pressure(q[1], fluid)
            if required > fluid['max_injection_pressure_pa']:
                raise ValueError('injection_pressure_limit')
            pump = max(0., required-fluid['storage_pressure_pa'])/(rho*fluid['pump_efficiency'])
            usable = finite(dh-pump, 'net wall-heat capacity after pump work')
            flux = node['peak_flux_w_m2']/usable
            limit = fluid.get('maximum_mass_flux_kg_m2_s')
            if limit is not None and flux > limit:
                raise ValueError('mass_flux_limit')
            dm = node['heat_j']/usable
            mass += dm; covered += node['heat_j']; pump_j += dm*pump
            max_flux = max(max_flux, flux); max_pressure = max(max_pressure, required)
            effective[q] = usable
        except (ValueError, KeyError, TypeError, ArithmeticError) as exc:
            failures.append({'temperature_k': q[0], 'pressure_pa': q[1], 'heat_j': node['heat_j'], 'reason': str(exc)})
    peak_flow = None
    if profile is not None:
        # Includes zero-area stagnation probes: they add no mass but can fail delivery.
        flows = defaultdict(float)
        for row in profile['demands']:
            q = coord(row)
            if q not in effective:
                continue
            flux = row['heat_flux_w_m2']/effective[q]
            max_flux = max(max_flux, flux)
            if fluid.get('maximum_mass_flux_kg_m2_s') is not None and flux > fluid['maximum_mass_flux_kg_m2_s']:
                failures.append({'temperature_k': q[0], 'pressure_pa': q[1], 'heat_j': row['heat_j'], 'reason': 'exact_mass_flux_limit'})
            flows[row['time_s']] += row['heat_j']/row['dt_s']/effective[q]
        peak_flow = max(flows.values(), default=0.)
    complete = not failures and len(effective) == len(expected)
    return {'status': 'complete_model_estimate' if complete else 'incomplete_profile',
        'mass_kg': mass if complete else None, 'partial_mass_diagnostic_kg': mass,
        'heat_coverage_fraction': covered/total, 'failed_state_count': len(failures), 'failures': failures,
        'peak_mass_flux_kg_m2_s': max_flux if complete else None,
        'peak_mass_flow_kg_s': peak_flow if complete else None,
        'maximum_injection_pressure_pa': max_pressure if complete else None,
        'pump_work_j': pump_j if complete else None,
        'rankable': False, 'experimental_winner': False, 'chemistry_validated': False,
        'mass_score_basis': 'fixed-load-nonreactive-demand-replay'}


def combine(pair, w, model_results, water_mass):
    good = {m: r for m, r in model_results.items() if r.get('status') == 'complete_model_estimate'}
    expected = sorted(pair['models'])
    complete = bool(expected) and set(good) == set(expected) and water_mass is not None and water_mass > 0
    worst = max((r['mass_kg'] for r in good.values()), default=None) if complete else None
    return {'candidate_id': digest([pair['pair_id'], w])[:24], 'pair_id': pair['pair_id'],
        'cas_number': pair['cas_number'], 'name': pair['name'], 'mass_fraction': w,
        'available_model_count': len(expected), 'successful_model_count': len(good),
        'status': 'complete_model_estimate' if complete else 'incomplete_models_or_reference',
        'worst_model_mass_kg': worst, 'mass_ratio_vs_heos_water': worst/water_mass if complete else None,
        'model_relative_spread': ((max(r['mass_kg'] for r in good.values())/min(r['mass_kg'] for r in good.values()))-1) if len(good)>1 else None,
        'evidence_tier': 'two_model_exploratory' if len(expected)>1 else 'single_model_exploratory',
        'heat_coverage_min': min((r.get('heat_coverage_fraction', 0.) for r in model_results.values()), default=0.),
        'rankable': False, 'experimental_winner': False}


def shortlist(rows, pair_budget, per_pair):
    """Every supported pair is screened. Detail budget combines merit and spread."""
    usable = [r for r in rows if r['status'] == 'complete_model_estimate']
    grouped = defaultdict(list)
    for r in usable:
        grouped[r['pair_id']].append(r)
    best = sorted(grouped, key=lambda p: (min(r['mass_ratio_vs_heos_water'] for r in grouped[p]), p))
    uncertain = sorted(grouped, key=lambda p: (-max(r.get('model_relative_spread') or 0. for r in grouped[p]), p))
    pairs = []
    for i in range(max(len(best), len(uncertain))):
        for order in (best, uncertain):
            if i < len(order) and order[i] not in pairs and len(pairs)<pair_budget:
                pairs.append(order[i])
    selected = []
    for p in pairs:
        candidates = sorted(grouped[p], key=lambda r: (r['mass_ratio_vs_heos_water'], r['mass_fraction']))
        selected.extend(candidates[:per_pair])
    return sorted(selected, key=lambda r: (r['pair_id'], r['mass_fraction']))


def refinement_fractions(all_fractions, winners):
    values = set()
    for w in winners:
        i = all_fractions.index(w)
        for j in (i-1, i+1):
            if 0 <= j < len(all_fractions):
                values.add(round((w+all_fractions[j])/2, 12))
        values.add(w)
    return sorted(values)
