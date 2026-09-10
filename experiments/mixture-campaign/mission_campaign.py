"""One reproducible broad catalog -> mission screening -> exact replay campaign."""
from __future__ import annotations
import argparse
from collections import Counter, defaultdict
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import datetime, timezone
import json
import os
from pathlib import Path
import shutil
import subprocess
import uuid
import yaml

import mission_core as core
from campaign_store import Store, digest, encoded, implementation, read_json, write_json, exclusive_lock
from campaign_publication import compressed, table

HERE = Path(__file__).resolve().parent


def require(artifact):
    if artifact.outcome['status'] != 'complete':
        raise RuntimeError(artifact.outcome['error'])
    return artifact.data


class Campaign:
    def __init__(self, config, entry, store, driver):
        self.config, self.entry, self.store, self.driver = config, entry, store, driver
        self.analysis_code = implementation(Path(core.__file__))
        self.records, self.raw = [], []
        self.catalog = self.profile = None

    def prepare(self, limit_pairs):
        c, d, s = self.config, self.driver, self.store
        self.catalog = s.run('mission.catalog', {'catalog': c['catalog'], 'storage_temperature_k': c['fluid']['storage_temperature_k']},
            d.catalog_code, lambda: d.call('catalog', config=c), environment=d.environment)
        catalog = require(self.catalog)
        self.pairs = catalog['candidates'][:limit_pairs] if limit_pairs else catalog['candidates']
        self.by_pair = {p['pair_id']: p for p in self.pairs}
        print(f"[catalog] identities={catalog['scanned_identity_count']} eligible={catalog['eligible_pair_count_before_budget']} selected={len(self.pairs)} coverage={catalog['model_coverage_counts']}", flush=True)
        self.profile = s.run('mission.profile', {'entry': self.entry, 'mission': c['mission'], 'fluid': c['fluid']},
            d.profile_code, lambda: d.call('profile', entry=self.entry, config=c), environment=d.environment)
        profile = require(self.profile)
        if not profile['summary']['profile_complete']:
            raise ValueError('mission demand profile incomplete: '+profile['summary']['status'])
        print('[mission] '+json.dumps(profile['summary'], sort_keys=True), flush=True)
        self.nodes = {}
        self.water = {}
        for name, count in [('coarse', c['screening']['coarse_nodes']), ('refined', c['screening']['refined_nodes']), ('exact', None)]:
            self.nodes[name] = s.run('mission.quadrature', {'count': count}, self.analysis_code,
                lambda n=count: core.quadrature(profile, n if n is not None else len(profile['demands'])),
                dependencies={'profile': self.profile})
            ns = require(self.nodes[name])['nodes']
            water_states = s.run('mission.water_states', {'fluid': c['fluid']}, d.thermo_code,
                lambda nn=ns: d.call('water', nodes=nn, fluid=c['fluid']),
                dependencies={'nodes': self.nodes[name]}, environment=d.environment)
            self.water[name] = s.run('mission.water_score', {'fluid': c['fluid'], 'exact': name == 'exact'}, self.analysis_code,
                lambda a=water_states, nn=ns, stage=name: core.score(nn, require(a)['states'], c['fluid'],
                        profile=profile if stage == 'exact' else None),
                dependencies={'states': water_states, 'nodes': self.nodes[name], 'profile': self.profile})
            wr = require(self.water[name])
            print(f"[water {name}] nodes={len(ns)} mass_kg={wr['mass_kg']} coverage={wr['heat_coverage_fraction']:.6f}", flush=True)
        return catalog

    def evaluate_pair(self, pair, fractions, stage):
        c, d, s = self.config, self.driver, self.store
        ns_art = self.nodes[stage]
        nodes = require(ns_art)['nodes']
        by_w = {w: {} for w in fractions}
        raw = []
        for model in sorted(pair['models']):
            freeze = s.run('mission.model.freeze', {'pair': {**pair, 'models': {model: pair['models'][model]}},
                'model': model, 'fluid': c['fluid'], 'thermodynamics': c['thermodynamics']}, d.thermo_code,
                lambda m=model: d.call('freeze', pair=pair, model=m, config=c), environment=d.environment)
            pending, states = [], {}
            for w in fractions:
                spec = s.spec('mission.model.states', {'mass_fraction': w, 'model': model, 'pair_id': pair['pair_id']},
                    d.thermo_code, dependencies={'model': freeze, 'nodes': ns_art}, environment=d.environment)
                saved = s.load(spec)
                if saved is not None:
                    states[w] = saved
                elif freeze.outcome['status'] != 'complete':
                    states[w] = s.put(spec, error='model_freeze_failed: '+freeze.outcome['error'])
                else:
                    pending.append((w, spec))
            batch = c['execution']['batch_compositions']
            for i in range(0, len(pending), batch):
                chunk = pending[i:i+batch]
                try:
                    result = d.call('states', pair=pair, model=model, config=c, frozen=freeze.data,
                                    nodes=nodes, fractions=[w for w, _ in chunk])
                    rows = result['results']
                    indexed = {r['mass_fraction']: r for r in rows}
                    if len(indexed) != len(rows) or set(indexed) != {w for w, _ in chunk}:
                        raise ValueError('worker returned mismatched compositions')
                    for w, spec in chunk:
                        returned = indexed[w]['states']
                        if len(returned) != len(nodes) or {core.coord(r) for r in returned} != {core.coord(n) for n in nodes}:
                            raise ValueError('worker returned mismatched mission states')
                    for w, spec in chunk:
                        states[w] = s.put(spec, indexed[w])
                except Exception as exc:
                    for w, spec in chunk:
                        states[w] = s.put(spec, error=f'{type(exc).__name__}: {exc}')
                print(f"[batch {stage}] {pair['cas_number']} {model} {min(i+batch,len(pending))}/{len(pending)}", flush=True)
            for w in fractions:
                a = states[w]
                def evaluate(a=a):
                    if a.outcome['status'] != 'complete':
                        return {'status': 'numerical_task_failed', 'mass_kg': None, 'heat_coverage_fraction': 0.,
                                'failed_state_count': len(nodes), 'error': a.outcome['error']}
                    return core.score(nodes, a.data['states'], c['fluid'],
                                      profile=self.profile.data if stage == 'exact' else None)
                score = s.run('mission.score', {'fluid': c['fluid'], 'exact': stage == 'exact'}, self.analysis_code,
                    evaluate, dependencies={'states': a, 'nodes': ns_art, 'profile': self.profile})
                by_w[w][model] = require(score)
                raw.append({'pair_id': pair['pair_id'], 'mass_fraction': w, 'model': model, 'stage': stage,
                            'score': score.data, 'states': a.outcome, 'states_task': a.dependency(),
                            'model_metadata': freeze.outcome})
        records = []
        for w in fractions:
            r = core.combine(pair, w, by_w[w], require(self.water[stage])['mass_kg'])
            r.update(stage=stage, profile_hash=self.profile.result_hash,
                     quadrature_node_count=len(nodes), full_profile_evaluated=stage == 'exact')
            records.append(r)
        return records, raw

    def stage(self, design, name):
        rows, raw = [], []
        with ThreadPoolExecutor(max_workers=self.config['execution']['workers']) as pool:
            futures = {pool.submit(self.evaluate_pair, self.by_pair[p], sorted(set(ws)), name): p for p, ws in sorted(design.items())}
            for n, future in enumerate(as_completed(futures), 1):
                rr, dd = future.result(); rows.extend(rr); raw.extend(dd)
                print(f'[screen {name} {n}/{len(futures)}] pair={futures[future]} complete={sum(r["status"]=="complete_model_estimate" for r in rr)}/{len(rr)}', flush=True)
        rows.sort(key=lambda r: (r['pair_id'], r['mass_fraction']))
        self.records.extend(rows); self.raw.extend(raw)
        return rows

    def run(self, limit_pairs=None):
        self.prepare(limit_pairs)
        fractions = core.compositions(self.config['composition'])
        print(f'[design] distinct_pairs={len(self.pairs)} compositions={len(fractions)} initial_formulations={len(self.pairs)*len(fractions)}', flush=True)
        coarse = self.stage({p['pair_id']: fractions for p in self.pairs}, 'coarse')
        policy = self.config['screening']
        selected = core.shortlist(coarse, policy['detail_pairs'], policy['per_pair'])
        groups = defaultdict(list)
        for row in selected:
            groups[row['pair_id']].append(row['mass_fraction'])
        detail_design = {p: core.refinement_fractions(fractions, ws) for p, ws in groups.items()}
        refined = self.stage(detail_design, 'refined') if detail_design else []
        exact = core.shortlist(refined, policy['exact_candidates'], 1)
        exact_design = {r['pair_id']: [r['mass_fraction']] for r in exact}
        final = self.stage(exact_design, 'exact') if exact_design else []
        previous = {(r['pair_id'], r['mass_fraction']): r for r in refined}
        for r in final:
            a = previous[(r['pair_id'], r['mass_fraction'])]
            ratio = r['worst_model_mass_kg']/a['worst_model_mass_kg']-1 if r['worst_model_mass_kg'] is not None and a['worst_model_mass_kg'] else None
            r['relative_change_refined_to_exact'] = ratio
            r['within_quadrature_tolerance'] = ratio is not None and abs(ratio)<=policy['relative_mass_tolerance']
            r['qualification'] = ('fixed_load_model_hypothesis_requires_chemistry_and_transport_validation'
                if r['status']=='complete_model_estimate' and r['mass_ratio_vs_heos_water']<1-policy['gain_margin']
                else 'no_complete_mission_gain_hypothesis')
        return final


def publish(output, campaign, final, manifest, previous):
    s, c = campaign.store, campaign.config
    catalog, profile = campaign.catalog.data, campaign.profile.data
    events = s.events
    summary = {'study_version': 'v5m-4', 'schema': 'broad-mixture-mission-summary-v1',
        'run_id': manifest['run_id'], 'catalog_identity_count': catalog['scanned_identity_count'],
        'eligible_pairs_before_budget': catalog['eligible_pair_count_before_budget'], 'catalog_budget_limited': catalog['budget_limited'],
        'screened_pair_count': len(campaign.pairs), 'composition_count': len(core.compositions(c['composition'])),
        'initial_formulation_count': len(campaign.pairs)*len(core.compositions(c['composition'])),
        'evaluations_by_stage': dict(Counter(r['stage'] for r in campaign.records)),
        'status_counts_by_stage': {stage: dict(Counter(r['status'] for r in campaign.records if r['stage']==stage)) for stage in ['coarse','refined','exact']},
        'model_coverage_counts': catalog['model_coverage_counts'], 'exact_replay_candidates': len(final),
        'complete_exact_replay_count': sum(r['status']=='complete_model_estimate' for r in final),
        'mission_gain_hypothesis_count': sum(r['qualification'].startswith('fixed_load_model_hypothesis') for r in final),
        'water_reference_mass_kg': campaign.water['exact'].data['mass_kg'], 'mission_profile': profile['summary'],
        'cached_task_count': sum(e['action']=='cache_hit' for e in events),
        'executed_task_count': sum(e['action']=='executed' for e in events),
        'failed_task_count': sum(e['status']=='failed' for e in events),
        'limited_run': manifest['limit_pairs'] is not None,
        'experimental_winners': 0, 'rankable_promotions': 0, 'ml_model_trained': False,
        'physical_chemistry_validated': False, 'injection_flowfield_simulated': False,
        'continuum_convergence_claimed': False, 'exhaustive_catalog_claimed': False,
        'interpretation': 'Broad local-identity search followed by fixed-entry-load mass screening and exact discrete-profile replay. Gas continuation is frozen-species and provisional. This is not a coupled reactive TPS simulation.'}
    scientific = {'summary': {k:v for k,v in summary.items() if k not in ['run_id','cached_task_count','executed_task_count']},
                  'finalists': final, 'records': campaign.records}
    summary['scientific_result_sha256'] = digest(scientific)
    index = []
    output.mkdir(parents=True)
    for name, data in [('summary',summary), ('manifest',manifest)]:
        write_json(output/(name+'.json'), data)
        index.append(compressed(output, name+'.json.gz', encoded(data)))
    for name, data in [('scientific-results',scientific), ('catalog',catalog), ('mission-profile',profile),
                       ('task-graph',s.graph()), ('execution-events',s.events)]:
        index.append(compressed(output, name+'.json.gz', encoded(data)))
    catalog_rows = [{k:r.get(k) for k in ('pair_id','cas_number','name','smiles','inchi_key')} |
                    {'models': sorted(r['models']), 'procurement_verified': False} for r in catalog['eligible_inventory']]
    for name, rows in [('catalog',catalog_rows), ('catalog-rejections',[{'reason':k,'count':v} for k,v in catalog['rejection_counts'].items()]),
                       ('ranking',campaign.records), ('finalists',final)]:
        index.extend(table(output, name, rows, c['publication']['shard_bytes']))
    old = {r['candidate_id']:r for r in (previous or {}).get('finalists',[])}
    new = {r['candidate_id']:r for r in final}
    changes = [{'candidate_id': k, 'before': old.get(k), 'after': new.get(k)} for k in sorted(old.keys()|new.keys()) if old.get(k)!=new.get(k)]
    index.append(compressed(output, 'changes.json.gz', encoded(changes)))
    index.extend(table(output, 'changes', changes, c['publication']['shard_bytes']))
    needs = [{'pair_id': p['pair_id'], 'cas_number': p['cas_number'],
              'needed': ['mixture_caloric_and_phase_validation','thermal_decomposition_and_air_reactions',
                         'porous_delivery_and_residue','injection_heat_transfer_coupling','procurement_and_hardware_safety'],
              'status': 'not_satisfied', 'automatically_acquired': False} for p in campaign.pairs]
    index.extend(table(output, 'evidence-needs', needs, c['publication']['shard_bytes']))
    for i, p in enumerate(campaign.pairs, 1):
        data = [r for r in campaign.raw if r['pair_id']==p['pair_id']]
        index.append(compressed(output, f'pair-{i:04d}.json.gz', encoded({'pair':p,'evaluations':data})))
    write_json(output/'publication-index.json', {'files':index, 'interpretation':'gzip data plus bounded text reports'})
    return summary


def execute(config, entry, cache, results, *, limit_pairs=None, retry_failures=False, recompute=False, driver_factory=None):
    core.validate(config)
    if driver_factory is None:
        from mission_driver import Driver
        driver_factory = Driver
    previous = None
    if (results/'latest.json').exists():
        pointer = read_json(results/'latest.json')
        name = pointer['run_id']
        if Path(name).name != name or not name.startswith('run-'):
            raise ValueError('invalid historical snapshot path')
        previous = read_json(results/name/'scientific-results.json.gz')
        if digest(previous)!=pointer['scientific_result_sha256']:
            raise ValueError('historical snapshot hash mismatch')
    run_id = 'run-'+datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%SZ')+'-'+uuid.uuid4().hex[:8]
    work = cache/'work'/run_id
    work.mkdir(parents=True)
    store = Store(cache, work/'events.jsonl', retry_failures=retry_failures, recompute=recompute)
    driver = driver_factory(work/'workers',config['execution']['timeout_s'])
    campaign = Campaign(config, entry, store, driver)
    final = campaign.run(limit_pairs)
    try:
        commit = subprocess.check_output(['git','rev-parse','HEAD'],cwd=HERE,text=True,stderr=subprocess.DEVNULL).strip()
    except (OSError,subprocess.CalledProcessError):
        commit = None
    manifest = {'run_id':run_id, 'git_commit':commit, 'config':config, 'entry_config':entry,
        'environment':driver.environment, 'limit_pairs':limit_pairs,
        'workflow_implementation':implementation(Path(__file__),Path(core.__file__)),
        'profile_task':campaign.profile.dependency(), 'catalog_task':campaign.catalog.dependency(),
        'recompute':recompute, 'retry_failures':retry_failures}
    staging = results/('.building-'+run_id)
    summary = publish(staging,campaign,final,manifest,previous)
    destination = results/run_id
    os.replace(staging,destination)
    if limit_pairs is None:
        write_json(results/'latest.json', {'run_id':run_id,'scientific_result_sha256':summary['scientific_result_sha256']})
    print('[published] '+str(destination), flush=True)
    print(json.dumps(summary,indent=2,sort_keys=True), flush=True)
    return summary,destination


def main(argv=None):
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--config',type=Path,default=HERE/'mission-campaign.yaml')
    p.add_argument('--cache-dir',type=Path,default=HERE/'.campaign-cache')
    p.add_argument('--results-dir',type=Path,default=HERE/'mission-v5m4-results')
    p.add_argument('--workers',type=int)
    p.add_argument('--limit-pairs',type=int)
    p.add_argument('--retry-failures',action='store_true')
    p.add_argument('--recompute',action='store_true')
    a=p.parse_args(argv)
    try:
        c=yaml.safe_load(a.config.read_text(encoding='utf-8'))
        if a.workers is not None: c['execution']['workers']=a.workers
        if a.limit_pairs is not None and a.limit_pairs<1: raise ValueError('limit-pairs must be positive')
        core.validate(c)
        ep=Path(c['mission']['entry_config'])
        if not ep.is_absolute(): ep=a.config.resolve().parent/ep
        entry=yaml.safe_load(ep.read_text(encoding='utf-8'))
        cache,results=a.cache_dir.resolve(),a.results_dir.resolve()
        if cache==results or cache in results.parents or results in cache.parents:
            raise ValueError('cache and result roots must not contain each other')
        with exclusive_lock(cache),exclusive_lock(results):
            summary,_=execute(c,entry,cache,results,limit_pairs=a.limit_pairs,
                             retry_failures=a.retry_failures,recompute=a.recompute)
        return 1 if summary['failed_task_count'] else 0
    except (ValueError,KeyError,TypeError,OSError,RuntimeError,ImportError) as exc:
        p.error(str(exc))


if __name__=='__main__':
    raise SystemExit(main())
