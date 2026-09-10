"""Common V5m-4.1 studies and selection. Predictions are never measurements.

Failure categories route further work; they do not establish physical infeasibility.
The zero-additive limit is evaluated by the same binary implementation, not HEOS.
"""
from __future__ import annotations

from collections import Counter, defaultdict
from copy import deepcopy
import importlib
import inspect
import math
from pathlib import Path
from typing import Callable

from campaign_store import digest, file_hash

Study = Callable[[dict, dict], dict]


def classify(reason: str) -> str:
    text = str(reason).lower()
    for category, tokens in (
        ('reference_blocked', ('water_caloric_reference_gate', 'heos_reference', 'water_control',
                               'pure_water_saturation', 'reference_unavailable')),
        ('storage_phase_model_rejection', ('storage_not_single_liquid',)),
        ('outlet_phase_model_rejection', ('outlet_liquid_liquid_split',)),
        ('gas_bridge_or_flash_unsupported', ('no_bounded_flash_or_vapor_bridge',)),
        ('property_domain_or_parameters_missing', ('outside_', 'correlation', 'parameter', 'no_complete',
                                                   'missing dortmund', 'unknown dortmund')),
        ('delivery_limit', ('injection_pressure_limit', 'mass_flux_limit')),
        ('nonpositive_or_invalid_enthalpy', ('enthalpy window', 'net wall-heat capacity')),
        ('numerical_task_failed', ('worker_', 'timeout', 'model_freeze_failed', 'numerical_task_failed')),
    ):
        if any(token in text for token in tokens):
            return category
    return 'unclassified_requires_diagnosis'


ROUTES = {
    'reference_blocked': 'reference_model_review',
    'storage_phase_model_rejection': 'storage_phase_validation',
    'outlet_phase_model_rejection': 'phase_equilibrium_validation',
    'gas_bridge_or_flash_unsupported': 'phase_and_hot_gas_model_extension',
    'property_domain_or_parameters_missing': 'property_data_or_model_acquisition',
    'delivery_limit': 'hardware_limits_review',
    'nonpositive_or_invalid_enthalpy': 'caloric_model_review',
    'numerical_task_failed': 'retry_or_solver_diagnosis',
    'unclassified_requires_diagnosis': 'preserved_error_review',
}


def annotate(records: list[dict], raw: list[dict], water_mass: float | None) -> tuple[list[dict], list[dict]]:
    """Keep HEOS scores intact and add matched-model effects without a correction."""
    controls = {}
    scores = {}
    for item in raw:
        scores[(item['mass_fraction'], item['model'])] = item['score']
        if item['mass_fraction'] == 0:
            score = item['score']
            mass = score.get('mass_kg') if score.get('status') == 'complete_model_estimate' else None
            controls[item['model']] = {
                'model': item['model'], 'status': score['status'], 'mass_kg': mass,
                'mass_ratio_vs_heos_water': mass / water_mass if mass and water_mass else None,
                'control_kind': 'same_binary_implementation_zero_additive',
                'source_is_measurement': False,
            }
    output = []
    for source in records:
        if source['mass_fraction'] == 0:
            continue
        row = deepcopy(source)
        effects = {}
        for model, zero in controls.items():
            score = scores.get((row['mass_fraction'], model), {})
            value = score.get('mass_kg') if score.get('status') == 'complete_model_estimate' else None
            effects[model] = value / zero['mass_kg'] if value and zero['mass_kg'] else None
        complete = (row['status'] == 'complete_model_estimate'
                    and len(effects) == row['available_model_count']
                    and all(v is not None and math.isfinite(v) and v > 0 for v in effects.values()))
        row.update(same_model_mass_ratios=effects,
                   worst_same_model_mass_ratio=max(effects.values()) if complete else None,
                   zero_control_status='complete' if complete else 'incomplete',
                   zero_limit_does_not_validate_mixture=True)
        output.append(row)
    return output, [controls[m] for m in sorted(controls)]


def coverage(data: dict, config: dict) -> dict:
    records = data['records']
    complete = sum(r['status'] == 'complete_model_estimate' for r in records)
    return {'formulation_count': len(records), 'complete_formulation_count': complete,
            'unresolved_formulation_count': len(records)-complete,
            'controlled_formulation_count': sum(r['zero_control_status'] == 'complete' for r in records),
            'all_formulations_evaluated': True, 'physical_validation_complete': False}


def failure_diagnosis(data: dict, config: dict) -> dict:
    """One primary category per failed model-state; retain exact original errors."""
    counts, examples, affected = Counter(), {}, defaultdict(set)
    failed_tasks = 0
    for item in data['raw']:
        if item['mass_fraction'] == 0:
            continue  # Control failures have a separate audit; do not inflate formulation counts.
        outcome = item['states']
        if outcome['status'] != 'complete':
            failed_tasks += 1
            errors = [{'reason': outcome.get('error', 'numerical_task_failed')}]
        else:
            errors = item['score'].get('failures', [])
        for error in errors:
            reason = str(error.get('reason', 'missing_failure_reason'))
            category = classify(reason)
            if outcome['status'] != 'complete' and category == 'unclassified_requires_diagnosis':
                category = 'numerical_task_failed'
            counts[category] += 1
            affected[category].add(item['mass_fraction'])
            examples.setdefault(category, {'reason': reason, 'model': item['model'],
                'mass_fraction': item['mass_fraction'], 'temperature_k': error.get('temperature_k'),
                'pressure_pa': error.get('pressure_pa')})
    return {'failed_task_count': failed_tasks, 'categories': [
        {'category': k, 'failed_model_state_or_task_count': counts[k],
         'affected_formulation_count': len(affected[k]), 'example': examples[k],
         'next_action': ROUTES[k], 'empirical_rejection': False}
        for k in sorted(counts)], 'unknown_errors_preserved': True,
        'counts_are_not_independent_experiments': True}


def zero_limit(data: dict, config: dict) -> dict:
    return {'controls': data['controls'], 'all_controls_complete': bool(data['controls']) and all(
        c['status'] == 'complete_model_estimate' for c in data['controls']),
        'heos_scores_modified': False, 'empirical_uncertainty_calibrated': False}


def evidence_needs(data: dict, config: dict) -> dict:
    diagnoses = failure_diagnosis(data, config)
    requests = []
    for row in diagnoses['categories']:
        requests.append({'category': row['category'], 'route': row['next_action'],
            'example_conditions': row['example'], 'status': 'pending', 'satisfied': False,
            'source': None, 'reference_value': None})
    for kind in ('mixture_phase_and_caloric_validation', 'high_temperature_chemistry',
                 'porous_delivery_and_residue', 'injection_heat_transfer', 'procurement_and_safety'):
        requests.append({'category': kind, 'route': 'independent_data_or_measurements',
                         'status': 'pending', 'satisfied': False, 'source': None, 'reference_value': None})
    return {'requests': requests, 'automatic_online_acquisition_performed': False}


def delivery_requirements(data: dict, config: dict) -> dict:
    """Report already computed flow/pressure demands, not a fabricated pore model."""
    values = []
    for item in data['raw']:
        score = item['score']
        if item['mass_fraction'] == 0 or score.get('status') != 'complete_model_estimate':
            continue
        values.append({'mass_fraction': item['mass_fraction'], 'model': item['model'],
                      **{k: score.get(k) for k in ('peak_mass_flux_kg_m2_s', 'peak_mass_flow_kg_s',
                                                   'maximum_injection_pressure_pa', 'pump_work_j')}})
    return {'requirements': values, 'hardware_feasibility_validated': False,
            'configured_mass_flux_limit': config['fluid'].get('maximum_mass_flux_kg_m2_s'),
            'pore_permeability_inferred': False, 'heat_transfer_coefficient_inferred': False}


def load_studies(modules: list[str] = ()) -> dict[str, Study]:
    studies = {'coverage': coverage, 'failure_diagnosis': failure_diagnosis,
               'zero_limit': zero_limit, 'evidence_needs': evidence_needs,
               'delivery_requirements': delivery_requirements}
    for name in modules:
        extra = importlib.import_module(name).STUDIES
        if not isinstance(extra, dict) or set(extra) & set(studies):
            raise ValueError('study plugins must add unique study names')
        for key, fn in extra.items():
            if not isinstance(key, str) or not key.isidentifier() or not callable(fn):
                raise ValueError('invalid study plugin')
        studies.update(extra)
    return studies


def study_signature(studies: dict[str, Study]) -> dict:
    result = {}
    for name, fn in sorted(studies.items()):
        source = inspect.getsourcefile(fn)
        if not source:
            raise ValueError('study must have inspectable source: ' + name)
        # Plugins may declare additional data/code dependencies via DEPENDENCY_FILES.
        module = inspect.getmodule(fn)
        files = [Path(source), *[Path(p) for p in getattr(module, 'DEPENDENCY_FILES', ())]]
        result[name] = {'module': fn.__module__, 'function': fn.__name__,
                        'files': {str(i)+':'+p.name: file_hash(p) for i, p in enumerate(files)}}
    return result


def merit(row: dict) -> tuple:
    value = row.get('worst_same_model_mass_ratio')
    return (value is None, value if value is not None else math.inf,
            row.get('mass_ratio_vs_heos_water') or math.inf, row['pair_id'], row['mass_fraction'])


def compact(data: dict, studies: dict, config: dict) -> dict:
    rows = data['records']
    good = sorted((r for r in rows if r['status'] == 'complete_model_estimate'), key=merit)
    uncertain = sorted((r for r in rows if r.get('model_relative_spread') is not None),
                       key=lambda r: (-r['model_relative_spread'], r['mass_fraction']))
    unresolved = sorted((r for r in rows if r['status'] != 'complete_model_estimate'),
                        key=lambda r: (-r.get('heat_coverage_min', 0), r['mass_fraction']))
    keep = config['screening']['per_pair']
    diagnosis = studies['failure_diagnosis']
    return {'pair_id': data['pair']['pair_id'], 'cas_number': data['pair']['cas_number'],
            'name': data['pair']['name'], 'stage': data['stage'],
            **studies['coverage'], 'top_records': good[:keep],
            'uncertain_record': uncertain[0] if uncertain else None,
            'unresolved_record': unresolved[0] if unresolved else None,
            'failure_categories': diagnosis['categories'],
            'failed_task_count': diagnosis['failed_task_count'],
            'required_studies': sorted(studies),
            'study_failures': [k for k, v in studies.items() if v.get('execution_status') == 'failed'],
            'controls': studies['zero_limit']['controls'],
            'rankable': False, 'experimental_winner': False}


def detail_selection(briefs: list[dict], history: dict, config: dict) -> dict[str, list[float]]:
    """Cumulative selection, plus mandatory historical hypotheses and bounded evidence work."""
    by_id = {b['pair_id']: b for b in briefs}
    selected: dict[str, set] = defaultdict(set)
    budget = config['screening']['detail_pairs']
    best = sorted((b for b in briefs if b['top_records']), key=lambda b: merit(b['top_records'][0]))
    for b in best[:budget]:
        selected[b['pair_id']].update(r['mass_fraction'] for r in b['top_records'])
    uncertain = sorted((b for b in briefs if b['uncertain_record'] and
                        (b['uncertain_record'].get('model_relative_spread') or 0) > 0),
                       key=lambda b: (-b['uncertain_record']['model_relative_spread'], b['pair_id']))
    for b in uncertain[:config['continuation']['evidence_pairs']]:
        selected[b['pair_id']].add(b['uncertain_record']['mass_fraction'])
    # A few incomplete cases receive identical follow-up rules, not favorable imputation.
    unresolved = sorted((b for b in briefs if b['unresolved_record']),
                        key=lambda b: (-b['unresolved_record'].get('heat_coverage_min', 0), b['pair_id']))
    for b in unresolved[:config['continuation']['evidence_pairs']]:
        selected[b['pair_id']].add(b['unresolved_record']['mass_fraction'])
    lo, hi = config['composition']['minimum'], config['composition']['maximum']
    for p, h in history.items():
        if p not in by_id or not h['hypothesis']:
            continue
        values = [w for w in h['fractions'] if lo <= w <= hi]
        selected[p].update(values)
    return {p: sorted(ws) for p, ws in sorted(selected.items()) if ws}


def qualify_exact(row: dict, refined: dict | None, config: dict) -> dict:
    result = deepcopy(row)
    a, b = row.get('worst_model_mass_kg'), (refined or {}).get('worst_model_mass_kg')
    change = a/b-1 if a is not None and b else None
    stable = change is not None and abs(change) <= config['screening']['relative_mass_tolerance']
    gain = config['screening']['gain_margin']
    complete = row['status'] == 'complete_model_estimate' and row['zero_control_status'] == 'complete'
    heos, same = row.get('mass_ratio_vs_heos_water'), row.get('worst_same_model_mass_ratio')
    hypothesis = complete and stable and heos < 1-gain and same < 1-gain
    result.update(relative_change_refined_to_exact=change, within_quadrature_tolerance=stable,
        qualification=('controlled_fixed_load_hypothesis_requires_validation' if hypothesis else
                       'numerical_resolution_unresolved' if complete and not stable else
                       'incomplete_models_or_controls' if not complete else 'no_controlled_mission_gain'),
        rankable=False, experimental_winner=False)
    return result
