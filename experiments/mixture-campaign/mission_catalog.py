"""Enumerate known local identities BEFORE checking pair-model availability."""
from __future__ import annotations
from collections import Counter
import math
from campaign_backend import legacy_import
from campaign_store import digest
from mission_core import select_catalog


def discover(config):
    from chemicals.identifiers import get_pubchem_db
    from chemicals import phase_change
    from thermo.interaction_parameters import IPDB
    from thermo.unifac import UNIFAC_group_assignment_DDBST, DOUFSG, DOUFIP2016
    legacy = legacy_import()
    db = get_pubchem_db()
    db.finish_loading()
    items = list(db.CAS_index.values())
    items.sort(key=lambda x: str(x.CASs))
    policy = config['catalog']
    tin = config['fluid']['storage_temperature_k']
    water_groups = UNIFAC_group_assignment_DDBST('7732-18-5', 'MODIFIED_UNIFAC')
    accepted, rejected, identities = [], [], set()
    for i, item in enumerate(items, 1):
        cas = str(item.CASs)
        if i % 10000 == 0:
            print(f'[catalog inventory] scanned={i}/{len(items)} model_backed={len(accepted)}', flush=True)
        if cas == '7732-18-5':
            continue
        try:
            if not 0 < float(item.MW) <= policy['maximum_molecular_weight_g_mol']:
                raise ValueError('outside_molecular_weight_domain')
            # Reuse exact structure checks; this does NOT impose NRTL membership.
            row = legacy.molecular_metadata(item, cas, policy)
            if row['inchi_key'] in identities:
                raise ValueError('duplicate_chemical_identity')
            identities.add(row['inchi_key'])
            tm, tb = phase_change.Tm(cas), phase_change.Tb(cas)
            if tm is None or tb is None or not math.isfinite(float(tm)) or not math.isfinite(float(tb)):
                raise ValueError('pure_storage_phase_not_established')
            if not float(tm) < tin < float(tb):
                raise ValueError('outside_neutral_liquid_additive_scope')
            row.update(pure_melting_temperature_k=float(tm), pure_boiling_temperature_k=float(tb),
                phase_constants_are_scope_guards_not_measurements=True,
                catalog_sources=['chemicals-local-identity-inventory'], models={}, model_rejections={})
            try:
                row['models']['chemsep_nrtl'] = legacy.nrtl_parameters(IPDB, cas)
            except (ValueError, KeyError, TypeError) as exc:
                row['model_rejections']['chemsep_nrtl'] = str(exc)
            try:
                groups = UNIFAC_group_assignment_DDBST(cas, 'MODIFIED_UNIFAC')
                row['models']['unifac_dortmund'] = legacy.complete_group_interactions([water_groups, groups], DOUFSG, DOUFIP2016)
            except (ValueError, KeyError, TypeError) as exc:
                row['model_rejections']['unifac_dortmund'] = str(exc)
            if not row['models']:
                raise ValueError('no_complete_supported_activity_model')
            row['parameter_sha256'] = digest(row['models'])
            accepted.append(row)
        except (ValueError, KeyError, TypeError, AttributeError) as exc:
            rejected.append({'cas_number': cas, 'reason': str(exc)})
    selected = select_catalog(accepted, policy['maximum_pairs'], policy['selection_seed'])
    return {'schema': 'broad-water-catalog-v1', 'candidates': selected, 'eligible_inventory': accepted,
        'rejections': rejected, 'scanned_identity_count': len(items),
        'eligible_pair_count_before_budget': len(accepted), 'selected_pair_count': len(selected),
        'budget_limited': len(selected)<len(accepted), 'source': 'local_identity_inventory_not_NRTL_membership',
        'selection': 'stable_hash_budget_no_merit_or_NRTL_preference',
        'model_coverage_counts': dict(Counter(m for r in selected for m in r['models'])),
        'rejection_counts': dict(Counter(r['reason'].split(':')[0] for r in rejected)),
        'commercial_availability_verified': False, 'electrolytes_and_suspensions_supported': False,
        'property_database_is_not_independent_validation': True}
