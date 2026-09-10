"""V5m-4.1 orchestration tests. Synthetic fixtures are not thermodynamic evidence."""
from __future__ import annotations

from collections import Counter
from contextlib import redirect_stdout
from copy import deepcopy
import gzip
import importlib.util
import io
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import threading
import unittest
from unittest.mock import patch

ROOT = Path(__file__).resolve().parents[3]
MODULE = ROOT/'experiments/mixture-campaign'
sys.path.insert(0, str(MODULE))
import yaml
import mission_batch_analysis as analysis
import mission_batch_archive as archive
from campaign_store import Store, digest, encoded, read_json, write_json, file_hash
from mission_batch_campaign import BatchedCampaign, validate


def config():
    c = yaml.safe_load((MODULE/'mission-campaign-v41.yaml').read_text())
    c['execution'].update(workers=2, batch_compositions=3)
    c['composition'].update(minimum=.01, transition=.1, maximum=.6, log_points=2, linear_points=3)
    c['screening'].update(coarse_nodes=2, refined_nodes=3, detail_pairs=1, per_pair=1, exact_candidates=1)
    c['continuation'].update(batch_pairs=2, max_batches=None, evidence_pairs=0, history_roots=[])
    return c


def pairs(n=4):
    return [{'pair_id': f'water-P{i:02}', 'cas_number': f'{i+11}-11-1', 'name': f'Fixture {i}',
             'inchi_key': f'KEY{i:02}', 'smiles': 'CCO',
             'models': {'unifac_dortmund': {'fixture': i}}, 'parameter_sha256': digest(i)} for i in range(n)]


class FakeDriver:
    environment = {'library': 'synthetic-v1'}
    thermo_code = {'thermo': 'synthetic-v1'}
    profile_code = {'profile': 'synthetic-v1'}
    catalog_code = {'catalog': 'synthetic-v1'}
    calls = []
    guard = threading.Lock()
    inventory = pairs()
    fail_pair = None
    bad_control = False

    def __init__(self, work, timeout):
        self.work = work

    def call(self, operation, **kw):
        with self.guard:
            self.calls.append((operation, kw.get('pair', {}).get('pair_id'), list(kw.get('fractions', []))))
        if operation == 'catalog':
            p = deepcopy(self.inventory)
            return {'candidates': p, 'eligible_inventory': p, 'rejections': [], 'rejection_counts': {},
                    'eligible_pair_count_before_budget': len(p), 'scanned_identity_count': 12,
                    'budget_limited': False, 'model_coverage_counts': {'unifac_dortmund': len(p)}}
        if operation == 'profile':
            rows = [{'temperature_k': 900., 'pressure_pa': float(101325+1000*i), 'heat_j': 2e6,
                     'heat_flux_w_m2': 1e5, 'time_s': float(i), 'dt_s': 1., 'area_m2': 20.} for i in range(4)]
            return {'demands': rows, 'summary': {'profile_complete': True, 'status': 'terminal_velocity',
                    'required_heat_j': 8e6, 'demand_rows': 4, 'fluid_independent_profile': True}}
        if operation == 'water':
            return {'states': [{**n, 'status': 'ok', 'delta_h_j_kg': 4e6, 'storage_density_kg_m3': 1000.}
                               for n in kw['nodes']]}
        if operation == 'freeze':
            return {'metadata': {'pair': kw['pair']['pair_id'], 'model': kw['model']}}
        if operation == 'states':
            if kw['pair']['pair_id'] == self.fail_pair:
                raise RuntimeError('worker_failed: synthetic failure')
            i = int(kw['pair']['pair_id'][-2:])
            results = []
            for w in kw['fractions']:
                states = []
                for n in kw['nodes']:
                    r = {'temperature_k': n['temperature_k'], 'pressure_pa': n['pressure_pa'],
                         'status': 'ok', 'storage_density_kg_m3': 1000.,
                         'delta_h_j_kg': 4.01e6*(1+(i+1)*.02*w), 'vapor_mole_fraction': 1.}
                    if (w > .3 and i == 1) or (w == 0 and self.bad_control):
                        r.update(status='unsupported_or_blocked_state', reason='no_bounded_flash_or_vapor_bridge_for_requested_state')
                    states.append(r)
                results.append({'mass_fraction': w, 'states': states})
            return {'results': results}
        raise AssertionError(operation)


def run(root, c=None, driver=FakeDriver, **kw):
    with redirect_stdout(io.StringIO()):
        runner = BatchedCampaign(c or config(), {}, root/'cache', root/'results',
                                 driver_factory=driver, history_base=root, **kw)
        result = runner.run()
    return result, runner


class PolicyTests(unittest.TestCase):
    def test_complete_catalog_required(self):
        c=config(); c['catalog']['maximum_pairs']=500
        with self.assertRaises(ValueError): validate(c)
    def test_invalid_budgets(self):
        for value in (-1,0,True,1.5):
            c=config(); c['continuation']['batch_pairs']=value
            with self.subTest(value=value), self.assertRaises(ValueError): validate(c)
    def test_null_invocation_budget_allowed(self):
        validate(config())
    def test_failure_routes_preserve_unknown(self):
        self.assertEqual(analysis.classify('strange'), 'unclassified_requires_diagnosis')
    def test_known_failure_categories(self):
        for reason,category in [('storage_not_single_liquid','storage_phase_model_rejection'),
            ('water_caloric_reference_gate','reference_blocked'),
            ('outside_pure_correlation_range','property_domain_or_parameters_missing'),
            ('no_bounded_flash_or_vapor_bridge_for_requested_state','gas_bridge_or_flash_unsupported'),
            ('worker_timeout','numerical_task_failed'),('injection_pressure_limit','delivery_limit')]:
            self.assertEqual(analysis.classify(reason),category)
    def test_plugin_cannot_replace_core_studies(self):
        with patch.object(analysis.importlib,'import_module',return_value=type('Plugin',(),{'STUDIES':{'coverage':lambda a,b:{}}})):
            with self.assertRaises(ValueError): analysis.load_studies(['fake'])
    def test_unknown_evidence_does_not_fabricate_values(self):
        data={'records':[],'raw':[],'controls':[]}
        for r in analysis.evidence_needs(data,config())['requests']:
            self.assertIsNone(r['reference_value']); self.assertFalse(r['satisfied'])
    def test_no_control_is_not_a_pass(self):
        row={'mass_fraction':.1,'model':'m','score':{'status':'complete_model_estimate','mass_kg':1.}}
        records=[{'mass_fraction':.1,'status':'complete_model_estimate','available_model_count':1}]
        result,_=analysis.annotate(records,[row],1.)
        self.assertIsNone(result[0]['worst_same_model_mass_ratio'])
    def test_bias_removed_only_in_separate_column(self):
        raw=[{'mass_fraction':w,'model':'m','score':{'status':'complete_model_estimate','mass_kg':9.}} for w in (0,.1)]
        original={'mass_fraction':.1,'status':'complete_model_estimate','available_model_count':1,'mass_ratio_vs_heos_water':.9}
        rows,controls=analysis.annotate([original],raw,10.)
        self.assertEqual(rows[0]['mass_ratio_vs_heos_water'],.9)
        self.assertEqual(rows[0]['worst_same_model_mass_ratio'],1.)
        self.assertNotIn('worst_same_model_mass_ratio',original)
    def test_exact_resolution_required_for_hypothesis(self):
        r={'status':'complete_model_estimate','zero_control_status':'complete','mass_ratio_vs_heos_water':.8,
           'worst_same_model_mass_ratio':.8,'worst_model_mass_kg':8.}
        self.assertFalse(analysis.qualify_exact(r,{'worst_model_mass_kg':10.},config())['within_quadrature_tolerance'])
        self.assertEqual(analysis.qualify_exact(r,{'worst_model_mass_kg':10.},config())['qualification'],'numerical_resolution_unresolved')
    def test_both_references_required_for_hypothesis(self):
        r={'status':'complete_model_estimate','zero_control_status':'complete','mass_ratio_vs_heos_water':.8,
           'worst_same_model_mass_ratio':1.,'worst_model_mass_kg':8.}
        self.assertEqual(analysis.qualify_exact(r,r,config())['qualification'],'no_controlled_mission_gain')


class ArchiveTests(unittest.TestCase):
    def test_object_roundtrip(self):
        with tempfile.TemporaryDirectory() as t:
            a=archive.Archive(Path(t)/'results',Path(t)/'cache'); ref=a.put({'value':1})
            self.assertEqual(a.get(ref),{'value':1}); self.assertEqual(a.put({'value':1}),ref)
    def test_corrupt_gzip_rejected(self):
        with tempfile.TemporaryDirectory() as t:
            a=archive.Archive(Path(t)/'r',Path(t)/'c'); ref=a.put({'value':1})
            a.path(ref).write_bytes(b'corrupt')
            with self.assertRaises(ValueError): a.get(ref)
    def test_path_traversal_rejected(self):
        with tempfile.TemporaryDirectory() as t:
            a=archive.Archive(Path(t),Path(t))
            with self.assertRaises(ValueError): a.path({'content_sha256':'../evil'})
    def test_changed_upstream_reference_invalidates_dossier(self):
        with tempfile.TemporaryDirectory() as t:
            root=Path(t); a=archive.Archive(root/'r',root/'c'); store=Store(root/'c',root/'log')
            one=store.run('test',{}, {},lambda:1)
            key=digest('pointer'); a.checkpoint(key,a.put({'value':1}),{'failed_task_count':0},[one.dependency()])
            self.assertIsNotNone(a.saved(key))
            store.put(one.spec,2)
            self.assertIsNone(a.saved(key))
    def test_corrupt_pointer_rejected(self):
        with tempfile.TemporaryDirectory() as t:
            root=Path(t);a=archive.Archive(root/'r',root/'c');key=digest('p')
            a.checkpoint(key,a.put({'x':1}),{},[])
            p=root/'c/mission-v41-pointers'/(key+'.json');v=read_json(p);v['data']['brief']['x']=2;write_json(p,v)
            with self.assertRaises(ValueError):a.saved(key)
    def test_history_merge_retains_old_hypothesis(self):
        a={'x':{'pair_id':'x','cas_number':'1','hypothesis':True,'fractions':[.3],'sources':['a']}}
        b={'x':{'pair_id':'x','cas_number':'1','hypothesis':False,'fractions':[.01],'sources':['b']}}
        result=archive.merge_history(a,b)
        self.assertTrue(result['x']['hypothesis']);self.assertEqual(result['x']['fractions'],[.01,.3])
        self.assertEqual(a['x']['fractions'],[.3])
    def test_conflicting_history_rejected(self):
        a={'x':{'cas_number':'1','hypothesis':False,'fractions':[],'sources':[]}}
        b={'x':{'cas_number':'2','hypothesis':False,'fractions':[],'sources':[]}}
        with self.assertRaises(ValueError):archive.merge_history(a,b)


class CampaignTests(unittest.TestCase):
    def setUp(self):
        FakeDriver.calls=[];FakeDriver.inventory=pairs();FakeDriver.fail_pair=None;FakeDriver.bad_control=False
    def test_full_inventory_and_controls(self):
        with tempfile.TemporaryDirectory() as t:
            summary,runner=run(Path(t))
            self.assertEqual(summary['coarse_evaluated_pair_count'],4)
            self.assertEqual(summary['total_formulation_count'],16)
            self.assertEqual(summary['zero_controls_count'],4)
            self.assertTrue(summary['catalog_screening_complete'])
            self.assertTrue(summary['failure_category_counts'])
            self.assertTrue(all('zero_limit' in v['coarse']['brief']['required_studies'] for v in runner.registry.values()))
    def test_partial_then_continue(self):
        with tempfile.TemporaryDirectory() as t:
            root=Path(t); c=config();c['continuation']['max_batches']=1
            a,_=run(root,c);self.assertEqual(a['remaining_pair_count'],2)
            FakeDriver.calls=[]
            b,_=run(root,c);self.assertTrue(b['catalog_screening_complete'])
            self.assertEqual(b['coarse_evaluated_pair_count'],4)
            self.assertEqual(b['invocation_counts']['coarse_dossier_cache_hits'],2)
    def test_second_completed_run_does_not_call_solver(self):
        with tempfile.TemporaryDirectory() as t:
            root=Path(t);run(root);FakeDriver.calls=[];s,_=run(root)
            self.assertEqual([x for x in FakeDriver.calls if x[0] in ('states','freeze','profile','water')],[])
            self.assertEqual(s['invocation_counts']['coarse_dossier_cache_hits'],4)
    def test_incremental_matches_full(self):
        with tempfile.TemporaryDirectory() as t:
            root=Path(t);run(root/'full')
            c=config();c['continuation']['max_batches']=1;run(root/'inc',c);run(root/'inc',c)
            a=archive.read_latest(root/'full/results')['scientific']
            b=archive.read_latest(root/'inc/results')['scientific']
            self.assertEqual(a,b)
    def test_recompute_matches_science(self):
        with tempfile.TemporaryDirectory() as t:
            root=Path(t);run(root);a=archive.read_latest(root/'results')['scientific']
            run(root,recompute=True);b=archive.read_latest(root/'results')['scientific']
            self.assertEqual(a,b)
    def test_new_pair_does_not_recompute_old_states(self):
        with tempfile.TemporaryDirectory() as t:
            root=Path(t);run(root);FakeDriver.calls=[]
            class More(FakeDriver): inventory=pairs(5);catalog_code={'catalog':'more'}
            summary,_=run(root,driver=More)
            self.assertEqual(summary['eligible_pair_count'],5)
            self.assertEqual(summary['invocation_counts']['coarse_dossier_cache_hits'],4)
    def test_thermo_change_recomputes_states_not_profile(self):
        with tempfile.TemporaryDirectory() as t:
            root=Path(t);run(root);FakeDriver.calls=[]
            class Changed(FakeDriver):thermo_code={'thermo':'new'}
            run(root,driver=Changed)
            self.assertTrue(any(x[0]=='states' for x in FakeDriver.calls))
            self.assertFalse(any(x[0]=='profile' for x in FakeDriver.calls))
    def test_added_study_reuses_states(self):
        with tempfile.TemporaryDirectory() as t:
            root=Path(t);run(root);FakeDriver.calls=[]
            plugin=root/'fixture_plugin.py'
            plugin.write_text('def extra(data, config):\n    return {"fixture_only": True}\nSTUDIES={"extra": extra}\n')
            sys.path.insert(0,str(root))
            try:
                c=config();c['continuation']['study_modules']=['fixture_plugin'];summary,runner=run(root,c)
                self.assertTrue(all('extra' in v['coarse']['brief']['required_studies'] for v in runner.registry.values()))
                self.assertFalse(any(x[0]=='states' for x in FakeDriver.calls))
            finally:
                sys.path.remove(str(root));sys.modules.pop('fixture_plugin',None)
    def test_control_failure_does_not_promote(self):
        with tempfile.TemporaryDirectory() as t:
            FakeDriver.bad_control=True
            s,r=run(Path(t))
            self.assertEqual(s['controlled_gain_hypothesis_count'],0)
            for entry in r.registry.values():
                self.assertTrue(all(x['worst_same_model_mass_ratio'] is None for x in entry['coarse']['brief']['top_records']))
    def test_failure_is_recorded_not_negative_result(self):
        with tempfile.TemporaryDirectory() as t:
            FakeDriver.fail_pair='water-P00';s,r=run(Path(t))
            self.assertGreater(s['failed_numerical_task_count'],0)
            b=r.registry['water-P00']['coarse']['brief'];self.assertEqual(b['complete_formulation_count'],0)
            self.assertEqual(b['unresolved_formulation_count'],4)
    def test_retry_task_failures(self):
        with tempfile.TemporaryDirectory() as t:
            root=Path(t);FakeDriver.fail_pair='water-P00';run(root)
            FakeDriver.fail_pair=None;FakeDriver.calls=[]
            s,_=run(root,retry_failures=True)
            self.assertEqual(s['failed_numerical_task_count'],0)
            self.assertTrue(any(x[0]=='states' and x[1]=='water-P00' for x in FakeDriver.calls))
    def test_interruption_after_completed_pair_preserves_checkpoint(self):
        with tempfile.TemporaryDirectory() as t:
            root=Path(t);old=BatchedCampaign.publish
            def stop(*a,**kw):raise KeyboardInterrupt
            with patch.object(BatchedCampaign,'publish',stop),self.assertRaises(KeyboardInterrupt):run(root)
            FakeDriver.calls=[];s,_=run(root)
            self.assertGreaterEqual(s['invocation_counts']['coarse_dossier_cache_hits'],2)
    def test_published_data_compressed_and_index_valid(self):
        with tempfile.TemporaryDirectory() as t:
            root=Path(t);run(root)
            latest=read_json(root/'results/latest.json');folder=root/'results'/latest['run_id']
            self.assertFalse(list(folder.glob('*.csv')))
            self.assertTrue(list(folder.glob('*.csv.gz')))
            for f in read_json(folder/'publication-index.json')['files']:
                self.assertEqual(file_hash(folder/f['file']),f['sha256'])
            self.assertTrue((folder/'registry.json.gz').exists())
            self.assertTrue(list((root/'results/objects').glob('*.json.gz')))
    def test_old_snapshots_not_overwritten(self):
        with tempfile.TemporaryDirectory() as t:
            root=Path(t);run(root);name=read_json(root/'results/latest.json')['run_id'];path=root/'results'/name/'summary.json';raw=path.read_bytes()
            run(root);self.assertEqual(path.read_bytes(),raw)
    def test_old_v4_states_are_reused_by_coarse_pass(self):
        from mission_campaign import execute as old_execute
        with tempfile.TemporaryDirectory() as t:
            root=Path(t);c=config();legacy=deepcopy(c);legacy['catalog']['maximum_pairs']=500
            with redirect_stdout(io.StringIO()):
                old_execute(legacy,{},root/'cache',root/'old',driver_factory=FakeDriver)
                engine=BatchedCampaign(c,{},root/'cache',root/'results',driver_factory=FakeDriver,history_base=root)
                engine.prepare();FakeDriver.calls=[]
                for pair in engine.ordered_pairs():engine.evaluate(pair,engine.coarse_fractions(pair),'coarse')
            requested=[w for op,_,ws in FakeDriver.calls if op=='states' for w in ws]
            self.assertTrue(requested)
            self.assertTrue(all(w==0 for w in requested),requested)
    def test_historical_hypothesis_gets_exact_replay(self):
        from campaign_publication import compressed
        with tempfile.TemporaryDirectory() as t:
            root=Path(t);folder=root/'legacy/run-old';folder.mkdir(parents=True)
            pair=pairs()[0]
            science={'pairs':[{'pair_id':pair['pair_id'],'cas_number':pair['cas_number'],
                'qualification':'local_model_hypothesis_requires_validation','best_point':[.3,350.,101325.]}]}
            write_json(folder/'summary.json',{'scientific_result_sha256':digest(science)})
            compressed(folder,'scientific-results.json.gz',encoded(science))
            item=compressed(folder,'catalog.json.gz',encoded({'candidates':[pair]}))
            write_json(folder/'publication-index.json',{'files':[item]})
            c=config();c['continuation']['history_roots']=['legacy']
            _,engine=run(root,c)
            self.assertTrue(engine.history[pair['pair_id']]['hypothesis'])
            data=engine.archive.get(engine.registry[pair['pair_id']]['exact']['object'])
            self.assertTrue(any(r['mass_fraction']==.3 for r in data['data']['records']))
    def test_history_hash_mismatch_is_not_silently_ignored(self):
        from campaign_publication import compressed
        with tempfile.TemporaryDirectory() as t:
            root=Path(t);folder=root/'legacy/run-old';folder.mkdir(parents=True)
            write_json(folder/'summary.json',{'scientific_result_sha256':digest({'different':1})})
            compressed(folder,'scientific-results.json.gz',encoded({'pairs':[]}))
            c=config();c['continuation']['history_roots']=['legacy']
            with self.assertRaises(ValueError):run(root,c)
    def test_removed_pair_remains_in_historical_registry(self):
        with tempfile.TemporaryDirectory() as t:
            root=Path(t);run(root)
            class Smaller(FakeDriver):inventory=pairs(3);catalog_code={'catalog':'smaller'}
            s,engine=run(root,driver=Smaller)
            self.assertEqual(s['historical_out_of_scope_count'],1)
            self.assertEqual(engine.out_of_scope[0]['pair_id'],'water-P03')
    def test_plain_text_previews_are_bounded(self):
        with tempfile.TemporaryDirectory() as t:
            root=Path(t);run(root)
            for path in (root/'results').rglob('preview-*.md'):
                self.assertLessEqual(path.stat().st_size,12288)
                self.assertIn('limited preview',path.read_text())
    def test_failed_plugin_is_reported_incomplete(self):
        with tempfile.TemporaryDirectory() as t:
            root=Path(t);plugin=root/'broken_plugin.py'
            plugin.write_text('def fail(data, config):\n    raise ValueError("fixture")\nSTUDIES={"broken": fail}\n')
            sys.path.insert(0,str(root))
            try:
                c=config();c['continuation']['study_modules']=['broken_plugin'];s,r=run(root,c)
                self.assertGreater(s['failed_study_count'],0)
                self.assertFalse(s['suite_execution_complete'])
                self.assertTrue(all('broken' in x['coarse']['brief']['study_failures'] for x in r.registry.values()))
            finally:
                sys.path.remove(str(root));sys.modules.pop('broken_plugin',None)
    def test_catalog_completion_is_not_physical_validation(self):
        with tempfile.TemporaryDirectory() as t:
            summary,_=run(Path(t))
            self.assertTrue(summary['catalog_screening_complete'])
            self.assertFalse(summary['physical_validation_complete'])
            self.assertEqual(summary['rankable_promotions'],0)
    def test_status_without_numerical_dependencies(self):
        with tempfile.TemporaryDirectory() as t:
            result=subprocess.run([sys.executable,str(MODULE/'run_campaign.py'),'--status','--results-dir',t],capture_output=True,text=True)
            self.assertEqual(result.returncode,0,result.stderr);self.assertIn('not_started',result.stdout)
    def test_help_modes(self):
        for mode in ('mission','mission-v4','properties'):
            result=subprocess.run([sys.executable,str(MODULE/'run_campaign.py'),'--mode',mode,'--help'],capture_output=True,text=True)
            self.assertEqual(result.returncode,0,result.stderr)
    def test_gitignore_tracks_objects_not_work(self):
        with tempfile.TemporaryDirectory() as t:
            root=Path(t);subprocess.run(['git','init','-q',t],check=True)
            (root/'.gitignore').write_bytes((MODULE/'.gitignore').read_bytes())
            for path,expected in [('mission-v5m41-results/objects/'+'a'*64+'.json.gz',False),
                 ('mission-v5m41-results/run-x/registry.json.gz',False),('mission-v5m41-results/run-x/report-x-0001.csv.gz',False),
                 ('mission-v5m41-results/run-x/report-x-0001.csv',True),('.campaign-cache/mission-v41-pointers/a.json',True),
                 ('mission-v5m41-results/objects/secret.log',True),('mission-v5m41-results/.building-run-x/summary.json',True)]:
                r=subprocess.run(['git','check-ignore','-q',path],cwd=root)
                self.assertEqual(r.returncode==0,expected,path)


HAVE_THERMO = all(importlib.util.find_spec(x) is not None for x in ('thermo','chemicals','CoolProp','rdkit'))
class NumericalIntegrationTests(unittest.TestCase):
    @unittest.skipUnless(HAVE_THERMO,'requires the pinned project thermodynamic environment')
    def test_same_binary_ethanol_zero_limit_and_dilution(self):
        from campaign_backend import legacy_import
        from chemicals.identifiers import get_pubchem_db
        from thermo.interaction_parameters import IPDB
        from mission_thermo import MissionModel
        legacy=legacy_import();c=config()
        item=get_pubchem_db().search_CAS('64-17-5',autoload=True)
        pair=legacy.molecular_metadata(item,'64-17-5',c['catalog'])
        pair['models']={'chemsep_nrtl':legacy.nrtl_parameters(IPDB,'64-17-5')}
        pair['parameter_sha256']=digest(pair['models'])
        obj=MissionModel(pair,'chemsep_nrtl',c)
        nodes=[{'temperature_k':900.,'pressure_pa':101325.}]
        zero=obj.evaluate(0.,nodes)[0];dilute=obj.evaluate(1e-8,nodes)[0]
        self.assertEqual(zero['status'],'ok',zero)
        self.assertEqual(dilute['status'],'ok',dilute)
        self.assertAlmostEqual(zero['delta_h_j_kg']/dilute['delta_h_j_kg'],1.,places=5)


if __name__=='__main__': unittest.main()
