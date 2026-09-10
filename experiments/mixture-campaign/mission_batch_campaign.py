#!/usr/bin/env python3
"""V5m-4.1: full inventory, bounded-memory batches and cumulative versioned studies.

This orchestrates the unchanged V5m-4 solver. It does not introduce reactive or
porous-flow physics, fabricate reference data, or promote a production coolant.
"""
from __future__ import annotations

import argparse
from collections import Counter
from concurrent.futures import ThreadPoolExecutor, as_completed
from copy import deepcopy
from datetime import datetime, timezone
import json
import os
from pathlib import Path
import subprocess
import uuid

import yaml
import mission_core as core
import mission_batch_analysis as analysis
import mission_batch_archive as archives
from campaign_store import Store, digest, implementation, read_json, write_json, exclusive_lock
from campaign_publication import compressed, table
from campaign_store import encoded
from mission_campaign import Campaign, require

HERE = Path(__file__).resolve().parent
DEFAULT_CONTINUATION = {'batch_pairs': 50, 'max_batches': None, 'evidence_pairs': 2,
    'history_roots': ['mission-v5m4-results', 'campaign-v5m3-results'], 'study_modules': []}


def validate(config):
    core.validate(config)
    if config['catalog']['maximum_pairs'] is not None:
        raise ValueError('V5m-4.1 requires catalog.maximum_pairs: null; use --max-batches for a work budget')
    policy = config.setdefault('continuation', {})
    for key, value in DEFAULT_CONTINUATION.items():
        policy.setdefault(key, deepcopy(value))
    for key in ('batch_pairs', 'evidence_pairs'):
        value = policy[key]
        if isinstance(value, bool) or not isinstance(value, int) or value < (1 if key == 'batch_pairs' else 0):
            raise ValueError('invalid continuation.'+key)
    value = policy['max_batches']
    if value is not None and (isinstance(value, bool) or not isinstance(value, int) or value < 1):
        raise ValueError('max_batches must be null or a positive integer')
    for key in ('history_roots', 'study_modules'):
        if not isinstance(policy[key], list) or any(not isinstance(v, str) for v in policy[key]):
            raise ValueError('continuation.'+key+' must be a list of strings')


class BatchedCampaign:
    def __init__(self, config, entry, cache, results, *, driver_factory=None,
                 retry_failures=False, recompute=False, history_base=HERE):
        self.config, self.entry = deepcopy(config), deepcopy(entry)
        validate(self.config)
        self.cache, self.results = Path(cache), Path(results)
        self.retry, self.recompute = retry_failures, recompute
        self.history_base = Path(history_base)
        self.invocation = datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%SZ')+'-'+uuid.uuid4().hex[:8]
        self.work = self.cache/'work'/('v41-'+self.invocation)
        self.work.mkdir(parents=True, exist_ok=True)
        if driver_factory is None:
            from mission_driver import Driver
            driver_factory = Driver
        self.driver = driver_factory(self.work/'workers', self.config['execution']['timeout_s'])
        self.store = self.make_store('shared')
        self.archive = archives.Archive(self.results, self.cache)
        self.shared = Campaign(self.config, self.entry, self.store, self.driver)
        self.studies = analysis.load_studies(self.config['continuation']['study_modules'])
        self.study_code = analysis.study_signature(self.studies)
        self.protocol_code = implementation(Path(__file__), Path(analysis.__file__), Path(archives.__file__))
        self.registry, self.history, self.history_sources = {}, {}, []
        self.known_pairs = {}
        self.stats = Counter()
        self.last_summary = None

    def make_store(self, label):
        return Store(self.cache, self.work/(label+'.events.jsonl'),
                     retry_failures=self.retry, recompute=self.recompute)

    def prepare(self):
        catalog = self.shared.prepare(None)
        if not self.shared.pairs:
            raise ValueError('no eligible catalog pairs; no successful full sweep can be claimed')
        # Preserve the standard water score and inherited profile validity flags.
        if any(require(a).get('mass_kg') is None for a in self.shared.water.values()):
            raise ValueError('HEOS mission reference incomplete')
        self.pairs = {p['pair_id']: p for p in self.shared.pairs}
        if len(self.pairs) != len(self.shared.pairs):
            raise ValueError('duplicate pair IDs in discovered inventory')
        self.fractions = core.compositions(self.config['composition'])
        previous = archives.read_latest(self.results)
        if previous:
            self.history = previous['registry']['history']
            self.known_pairs = previous['registry'].get('known_pairs', {})
        roots = [Path(v) if Path(v).is_absolute() else self.history_base/v
                 for v in self.config['continuation']['history_roots']]
        incoming, self.history_sources = archives.legacy_history(roots)
        self.history = archives.merge_history(self.history, incoming)
        # Existing screened pairs persist independently of the old random catalog budget.
        known = {**self.known_pairs, **self.history}
        self.out_of_scope = [h for p, h in sorted(known.items()) if p not in self.pairs]
        self.catalog_ref = self.archive.put(catalog)
        self.profile_ref = self.archive.put(require(self.shared.profile))
        self.shared_graph_ref = self.archive.put(self.store.graph())
        self.previous = previous
        self.study_settings = {k: v for k, v in self.config.items()
                               if k not in ('catalog', 'execution', 'publication', 'continuation')}
        self.context = {'study_settings': self.study_settings, 'profile': self.shared.profile.dependency(),
            'fluid': self.config['fluid'], 'thermodynamics': self.config['thermodynamics'],
            'thermo_code': self.driver.thermo_code, 'environment': self.driver.environment,
            'analysis': self.shared.analysis_code, 'suite': self.study_code,
            'protocol': self.protocol_code, 'screening': self.config['screening']}
        # Selection / execution budgets do not alter thermodynamic cache signatures.
        self.context_id = digest(self.context)
        if previous and previous['registry']['context_id'] != self.context_id:
            retained = {}
            for row in previous['scientific'].get('finalists', []):
                if not row.get('qualification', '').startswith('controlled_fixed_load'):
                    continue
                pid = row['pair_id']
                item = retained.setdefault(pid, {'pair_id': pid, 'cas_number': row['cas_number'],
                    'hypothesis': True, 'fractions': [],
                    'sources': ['v41-context:'+previous['registry']['context_id']]})
                item['fractions'].append(row['mass_fraction'])
            self.history = archives.merge_history(self.history, retained)
        self.stats.update('shared_'+e['action'] for e in self.store.events)
        for p in self.ordered_pairs():
            fractions = self.coarse_fractions(p)
            key = self.stage_key(p, fractions, 'coarse')
            saved = self.archive.saved(key, retry_failures=self.retry, recompute=self.recompute)
            if saved:
                self.registry[p['pair_id']] = {'coarse': self.short_pointer(saved)}
                self.stats['coarse_dossier_cache_hits'] += 1
        print(f'[continuation] inventory={len(self.pairs)} current_coarse={len(self.registry)} '
              f'pending={len(self.pairs)-len(self.registry)} historical_out_of_scope={len(self.out_of_scope)}', flush=True)

    def ordered_pairs(self):
        seed = self.config['catalog']['selection_seed']
        return sorted(self.pairs.values(), key=lambda p: (
            0 if self.history.get(p['pair_id'], {}).get('hypothesis') else
            1 if p['pair_id'] in self.history or p['pair_id'] in self.known_pairs else 2,
            digest([seed, p['inchi_key']]), p['pair_id']))

    def coarse_fractions(self, pair):
        lo, hi = self.config['composition']['minimum'], self.config['composition']['maximum']
        old = self.history.get(pair['pair_id'], {}).get('fractions', [])
        return sorted(set(self.fractions) | {w for w in old if lo <= w <= hi and w > 0})

    def stage_key(self, pair, fractions, stage):
        return digest({'context': self.context, 'pair': pair, 'fractions': fractions, 'stage': stage,
                       'nodes': self.shared.nodes[stage].dependency(),
                       'water': self.shared.water[stage].dependency()})

    @staticmethod
    def short_pointer(saved):
        return {k: saved[k] for k in ('key', 'object', 'brief')}

    def evaluate(self, pair, fractions, stage):
        key = self.stage_key(pair, fractions, stage)
        saved = self.archive.saved(key, retry_failures=self.retry, recompute=self.recompute)
        if saved:
            return self.short_pointer(saved), Counter({'dossier_cache_hits': 1})
        store = self.make_store(stage+'-'+key)
        current = Campaign(self.config, self.entry, store, self.driver)
        current.profile, current.nodes, current.water = self.shared.profile, self.shared.nodes, self.shared.water
        # Exactly zero is a control, not an extra discovered formulation.
        records, raw = current.evaluate_pair(pair, sorted(set(fractions) | {0.0}), stage)
        records, controls = analysis.annotate(records, raw, require(current.water[stage])['mass_kg'])
        data = {'pair': pair, 'records': records, 'raw': raw, 'controls': controls, 'stage': stage}
        dependencies = {a.key: a for a in store.used.values()}
        dependencies['nodes'] = current.nodes[stage]
        dependencies['profile'] = current.profile
        bundle = store.run('mission.v41.observations', {'stage': stage}, self.protocol_code,
                           lambda: data, dependencies=dependencies)
        require(bundle)
        studies = {}
        for name, fn in sorted(self.studies.items()):
            def invoke(fn=fn):
                value = fn(deepcopy(data), deepcopy(self.study_settings))
                if not isinstance(value, dict):
                    raise TypeError('study must return a JSON object')
                return value
            result = store.run('mission.v41.study.'+name, {'config': self.study_settings},
                self.study_code[name], invoke, dependencies={'observations': bundle})
            studies[name] = ({**result.data, 'execution_status': 'complete'} if result.outcome['status'] == 'complete'
                             else {'execution_status': 'failed', 'error': result.outcome['error']})
        for required in ('coverage', 'failure_diagnosis', 'zero_limit'):
            if studies[required]['execution_status'] != 'complete':
                raise RuntimeError('required core study failed: '+required)
        brief = analysis.compact(data, studies, self.config)
        brief['failed_numerical_task_keys'] = sorted(a.key for a in store.used.values()
            if a.outcome['status'] == 'failed' and not a.spec['kind'].startswith('mission.v41.study.'))
        brief['failed_control_task_count'] = sum(x['states']['status'] != 'complete' for x in raw if x['mass_fraction'] == 0)
        dossier = {'schema': 'mixture-mission-dossier-v41', 'stage_key': key,
                   'data': data, 'studies': studies, 'graph': store.graph()}
        ref = self.archive.put(dossier)
        tasks = [a.dependency() for a in sorted(store.used.values(), key=lambda a: a.key)]
        self.archive.checkpoint(key, ref, brief, tasks)
        counts = Counter(e['action'] for e in store.events)
        counts['model_states_executed'] = sum(e['kind'] == 'mission.model.states' and e['action'] == 'executed' for e in store.events)
        counts['model_states_cache_hit'] = sum(e['kind'] == 'mission.model.states' and e['action'] == 'cache_hit' for e in store.events)
        return {'key': key, 'object': ref, 'brief': brief}, counts

    def perform(self, design, stage):
        # At most batch_pairs futures exist. Each worker owns one short-lived Store.
        with ThreadPoolExecutor(max_workers=self.config['execution']['workers']) as pool:
            futures = {pool.submit(self.evaluate, self.pairs[p], ws, stage): p for p, ws in design.items()}
            for future in as_completed(futures):
                pair_id = futures[future]
                pointer, stats = future.result()
                self.registry.setdefault(pair_id, {})[stage] = pointer
                self.stats.update(stats)
                brief = pointer['brief']
                print(f'[pair {stage}] {brief["cas_number"]} '
                      f'complete={brief["complete_formulation_count"]}/{brief["formulation_count"]} '
                      f'zero_controls={all(x["status"] == "complete_model_estimate" for x in brief["controls"])}', flush=True)

    def details(self):
        briefs = [v['coarse']['brief'] for v in self.registry.values()]
        chosen = analysis.detail_selection(briefs, self.history, self.config)
        design = {}
        for p, ws in chosen.items():
            grid = self.coarse_fractions(self.pairs[p])
            # Every selection is an evaluated coarse point, including historic compositions.
            design[p] = core.refinement_fractions(grid, ws)
        batch_size = self.config['continuation']['batch_pairs']
        items = sorted(design.items())
        for start in range(0, len(items), batch_size):
            self.perform(dict(items[start:start+batch_size]), 'refined')
        candidates = sorted((self.registry[p]['refined']['brief']['top_records'][0]
                            for p in chosen if self.registry[p]['refined']['brief']['top_records']), key=analysis.merit)
        exact_design = {r['pair_id']: {r['mass_fraction']} for r in candidates[:self.config['screening']['exact_candidates']]}
        lo, hi = self.config['composition']['minimum'], self.config['composition']['maximum']
        for p, h in self.history.items():
            if p in chosen and h['hypothesis']:
                exact_design.setdefault(p, set()).update(w for w in h['fractions'] if lo <= w <= hi and w > 0)
        exact_design = {p: sorted(ws) for p, ws in exact_design.items() if ws}
        items = sorted(exact_design.items())
        for start in range(0, len(items), batch_size):
            self.perform(dict(items[start:start+batch_size]), 'exact')
        finalists = []
        for p in sorted(exact_design):
            detail = self.archive.get(self.registry[p]['refined']['object'])['data']['records']
            prior = {r['mass_fraction']: r for r in detail}
            for row in self.archive.get(self.registry[p]['exact']['object'])['data']['records']:
                finalists.append(analysis.qualify_exact(row, prior.get(row['mass_fraction']), self.config))
        return chosen, exact_design, sorted(finalists, key=analysis.merit)

    def publish(self, chosen, exact_design, finalists, *, interrupted=False):
        """A snapshot contains small cumulative summaries and shared dossier references."""
        batch_id = datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%SZ')+'-'+uuid.uuid4().hex[:8]
        run_id = 'run-'+batch_id
        staging = self.results/('.building-'+run_id)
        staging.mkdir(parents=True)
        baseline = [v['coarse']['brief'] for _, v in sorted(self.registry.items())]
        counts = Counter()
        for b in baseline:
            for r in b['failure_categories']:
                counts[r['category']] += r['failed_model_state_or_task_count']
        summaries = []
        for b in baseline:
            best = b['top_records'][0] if b['top_records'] else {}
            summaries.append({k: b[k] for k in ('pair_id','cas_number','name','formulation_count',
                'complete_formulation_count','unresolved_formulation_count','controlled_formulation_count','failed_task_count')} |
                {'best_mass_fraction': best.get('mass_fraction'), 'best_same_model_ratio': best.get('worst_same_model_mass_ratio'),
                 'best_heos_ratio': best.get('mass_ratio_vs_heos_water'),
                 'refined_selected': b['pair_id'] in chosen, 'exact_selected': b['pair_id'] in exact_design,
                 'historical_pair': b['pair_id'] in self.history,
                 'study_failures': b['study_failures']})
        active_entries = {p: {stage: value for stage, value in stages.items()
                           if stage == 'coarse' or stage == 'refined' and p in chosen or stage == 'exact' and p in exact_design}
                          for p, stages in sorted(self.registry.items())}
        for p in self.registry:
            self.known_pairs[p] = {k: self.pairs[p][k] for k in ('pair_id','cas_number','name')}
        registry = {'schema': 'mixture-mission-registry-v41', 'context_id': self.context_id,
            'history': self.history, 'history_sources': self.history_sources, 'known_pairs': self.known_pairs,
            'entries': active_entries, 'out_of_current_scope': self.out_of_scope}
        current = {p: {stage: value['object'] for stage, value in stages.items()}
                   for p, stages in active_entries.items()}
        active_briefs = [v['brief'] for stages in active_entries.values() for v in stages.values()]
        incomplete_studies = sum(len(b['study_failures']) for b in active_briefs)
        numerical_failures = {key for b in active_briefs for key in b.get('failed_numerical_task_keys', [])}
        all_catalog = len(self.registry) == len(self.pairs)
        summary = {'study_version': 'v5m-4.1', 'context_id': self.context_id,
            'catalog_identity_count': self.shared.catalog.data['scanned_identity_count'],
            'eligible_pair_count': len(self.pairs), 'coarse_evaluated_pair_count': len(self.registry),
            'remaining_pair_count': len(self.pairs)-len(self.registry),
            'complete_formulation_count': sum(b['complete_formulation_count'] for b in baseline),
            'unresolved_formulation_count': sum(b['unresolved_formulation_count'] for b in baseline),
            'total_formulation_count': sum(b['formulation_count'] for b in baseline),
            'zero_controls_count': sum(len(b['controls']) for b in baseline),
            'failure_category_counts': dict(sorted(counts.items())),
            'failed_numerical_task_count': len(numerical_failures),
            'failed_study_count': incomplete_studies,
            'exact_replay_count': len(finalists),
            'controlled_gain_hypothesis_count': sum(r['qualification'].startswith('controlled_fixed_load') for r in finalists),
            'water_reference_mass_kg': require(self.shared.water['exact'])['mass_kg'],
            'mission_profile': require(self.shared.profile)['summary'],
            'known_pair_count': len(self.known_pairs), 'historical_pair_count': len(self.history), 'historical_out_of_scope_count': len(self.out_of_scope),
            'catalog_screening_complete': all_catalog and not interrupted,
            'run_status': 'interrupted' if interrupted else 'catalog_screened' if all_catalog else 'batch_budget_pause',
            'suite_execution_complete': all_catalog and not incomplete_studies and not numerical_failures and not interrupted,
            'rankable_promotions': 0, 'experimental_winners': 0,
            'physical_validation_complete': False, 'continuous_domain_exhaustively_searched': False,
            'interpretation': 'Cumulative fixed-load frozen-species predictions. Incomplete states remain unresolved; zero limits are controls, not validation.'}
        scientific = {'schema': 'mixture-mission-results-v41', 'summary': summary,
            'catalog': self.catalog_ref, 'profile': self.profile_ref,
            'current_dossiers': current, 'finalists': finalists,
            'registry_sha256': digest(registry)}
        digest_result = digest(scientific)
        # Publication details do not participate in model cache keys.
        try:
            commit = subprocess.check_output(['git','rev-parse','HEAD'], cwd=HERE, text=True, stderr=subprocess.DEVNULL).strip()
        except (OSError, subprocess.CalledProcessError):
            commit = None
        manifest = {'run_id': run_id, 'git_commit': commit, 'config': self.config,
            'context': self.context, 'catalog_task': self.shared.catalog.dependency(),
            'profile_task': self.shared.profile.dependency(), 'shared_graph': self.shared_graph_ref,
            'previous_run_id': (self.previous or {}).get('pointer', {}).get('run_id'),
            'invocation_counts': dict(self.stats), 'scientific_result_sha256': digest_result}
        index = []
        for name, value in (('scientific-results', scientific), ('registry', registry), ('manifest', manifest)):
            index.append(compressed(staging, name+'.json.gz', encoded(value)))
        write_json(staging/'summary.json', {**summary, 'run_id': run_id,
            'scientific_result_sha256': digest_result, 'invocation_counts': dict(self.stats)})
        index.append(compressed(staging, 'summary.json.gz', (staging/'summary.json').read_bytes()))
        errors = [{'pair_id': b['pair_id'], 'cas_number': b['cas_number'], **d}
                  for b in baseline for d in b['failure_categories']]
        controls = [{'pair_id': b['pair_id'], 'cas_number': b['cas_number'], **d}
                    for b in baseline for d in b['controls']]
        old = {r['candidate_id']: r for r in (self.previous or {}).get('scientific', {}).get('finalists', [])}
        new = {r['candidate_id']: r for r in finalists}
        changes = [{'candidate_id': k, 'before_qualification': old.get(k, {}).get('qualification'),
                    'after_qualification': new.get(k, {}).get('qualification'),
                    'before_ratio': old.get(k, {}).get('mass_ratio_vs_heos_water'),
                    'after_ratio': new.get(k, {}).get('mass_ratio_vs_heos_water')}
                   for k in sorted(old.keys() | new.keys()) if old.get(k) != new.get(k)]
        for name, rows in (('pair-summary', summaries), ('finalists', finalists),
                           ('failure-diagnosis', errors), ('zero-controls', controls),
                           ('history-out-of-scope', self.out_of_scope), ('changes', changes)):
            index.extend(table(staging, name, rows, self.config['publication']['shard_bytes']))
        write_json(staging/'publication-index.json', {'files': index,
            'object_location': '../objects/<content_sha256>.json.gz relative to a run snapshot; see registry entries',
            'gzip_decoding_by_connector_assumed': False})
        os.replace(staging, self.results/run_id)
        write_json(self.results/'latest.json', {'run_id': run_id, 'scientific_result_sha256': digest_result})
        self.previous = {'pointer': {'run_id': run_id}, 'scientific': scientific, 'registry': registry}
        self.last_summary = {**summary, 'run_id': run_id, 'invocation_counts': dict(self.stats)}
        print('[published] '+str(self.results/run_id), flush=True)
        print(json.dumps(self.last_summary, sort_keys=True), flush=True)
        return self.last_summary

    def run(self):
        self.prepare()
        pending = [p for p in self.ordered_pairs() if p['pair_id'] not in self.registry]
        size = self.config['continuation']['batch_pairs']
        maximum = self.config['continuation']['max_batches']
        batches = 0
        if not pending:
            chosen, exact, finalists = self.details()
            return self.publish(chosen, exact, finalists)
        try:
            for start in range(0, len(pending), size):
                if maximum is not None and batches >= maximum:
                    break
                pairs = pending[start:start+size]
                print(f'[batch] {batches+1} pairs={len(pairs)} pending_before={len(self.pairs)-len(self.registry)}', flush=True)
                self.perform({p['pair_id']: self.coarse_fractions(p) for p in pairs}, 'coarse')
                chosen, exact, finalists = self.details()
                self.publish(chosen, exact, finalists)
                batches += 1
        except KeyboardInterrupt:
            # Every completed pair already has a verified object and atomic checkpoint.
            # An interrupted partial batch is not mislabeled as a finished campaign.
            self.publish({}, {}, [], interrupted=True)
            raise
        return self.last_summary


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--config', type=Path, default=HERE/'mission-campaign-v41.yaml')
    parser.add_argument('--cache-dir', type=Path, default=HERE/'.campaign-cache')
    parser.add_argument('--results-dir', type=Path, default=HERE/'mission-v5m41-results')
    parser.add_argument('--workers', type=int)
    parser.add_argument('--batch-pairs', type=int)
    parser.add_argument('--max-batches', type=int)
    parser.add_argument('--limit-pairs', type=int, help='Compatibility work-budget alias: one batch of this many pending pairs')
    parser.add_argument('--status', action='store_true')
    parser.add_argument('--retry-failures', action='store_true')
    parser.add_argument('--recompute', action='store_true')
    args = parser.parse_args(argv)
    try:
        if args.status:
            snapshot = archives.read_latest(args.results_dir.resolve())
            print(json.dumps(snapshot['scientific']['summary'] if snapshot else {'status': 'not_started'}, indent=2))
            return 0
        config = yaml.safe_load(args.config.read_text(encoding='utf-8'))
        validate(config)
        if args.workers is not None:
            config['execution']['workers'] = args.workers
        if args.batch_pairs is not None:
            config['continuation']['batch_pairs'] = args.batch_pairs
        if args.max_batches is not None:
            config['continuation']['max_batches'] = args.max_batches
        if args.limit_pairs is not None:
            config['continuation']['batch_pairs'] = args.limit_pairs
            config['continuation']['max_batches'] = 1
        validate(config)
        ep = Path(config['mission']['entry_config'])
        if not ep.is_absolute():
            ep = args.config.resolve().parent/ep
        entry = yaml.safe_load(ep.read_text(encoding='utf-8'))
        cache, results = args.cache_dir.resolve(), args.results_dir.resolve()
        if cache == results or cache in results.parents or results in cache.parents:
            raise ValueError('cache and results must not contain each other')
        with exclusive_lock(cache), exclusive_lock(results):
            engine = BatchedCampaign(config, entry, cache, results, retry_failures=args.retry_failures,
                                     recompute=args.recompute, history_base=args.config.resolve().parent)
            summary = engine.run()
        return 1 if summary['failed_numerical_task_count'] or summary['failed_study_count'] else 0
    except KeyboardInterrupt:
        print('Interrupted; completed pairs are checkpointed. Run the same command to continue.', flush=True)
        return 130
    except (ValueError, KeyError, TypeError, OSError, RuntimeError, ImportError) as exc:
        parser.error(str(exc))


if __name__ == '__main__':
    raise SystemExit(main())
