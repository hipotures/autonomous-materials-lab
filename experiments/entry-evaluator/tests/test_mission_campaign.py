"""Unit, incremental, publication and optional real-library mission tests."""
from __future__ import annotations
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
import unittest
import yaml

ROOT = Path(__file__).resolve().parents[3]
MODULE = ROOT/'experiments'/'mixture-campaign'
sys.path.insert(0,str(MODULE))
import mission_core as core
import mission_campaign as workflow
from campaign_store import digest, read_json
from mission_profile import entry_configuration


def config():
    c=yaml.safe_load((MODULE/'mission-campaign.yaml').read_text())
    c['execution']['workers']=1
    c['composition'].update(log_points=2,linear_points=3)
    c['screening'].update(coarse_nodes=2,refined_nodes=3,detail_pairs=2,per_pair=1,exact_candidates=2)
    c['fluid'].update(porous_delta_p_pa=0.,injection_pressure_margin=1.,storage_pressure_pa=200000.)
    return c


def profile():
    rows=[{'time_s':float(i), 'dt_s':1., 'temperature_k':500.+i*10,
           'pressure_pa':100000.+i*1000, 'heat_j':1000.+100*i,
           'heat_flux_w_m2':1000.+100*i, 'area_m2':1.} for i in range(5)]
    rows.append({**rows[2],'heat_j':0.,'area_m2':0.,'heat_flux_w_m2':2500.})
    return {'demands':rows,'history':[],'summary':{'profile_complete':True,'status':'terminal_altitude',
        'required_heat_j':sum(r['heat_j'] for r in rows),'coolant_mass_kg':None}}


def pair(i=1):
    return {'pair_id':f'water-{i}', 'cas_number':f'10-{i}-0','name':f'candidate {i}',
            'inchi_key':f'key-{i}','smiles':'CO', 'models':{'unifac_dortmund':{'value':i}},
            'parameter_sha256':str(i)}


class FakeDriver:
    environment={'fake':True}
    thermo_code={'version':'1'}
    catalog_code={'version':'1'}
    profile_code={'version':'1'}
    calls=[]
    def __init__(self,work,timeout): pass
    def call(self,operation,**v):
        self.calls.append(operation)
        if operation=='catalog':
            rows=[pair(1),pair(2)]
            return {'candidates':rows,'eligible_inventory':rows,'rejections':[],
                'scanned_identity_count':100,'eligible_pair_count_before_budget':2,
                'selected_pair_count':2,'budget_limited':False,
                'model_coverage_counts':{'unifac_dortmund':2},'rejection_counts':{}}
        if operation=='profile': return profile()
        if operation=='freeze': return {'metadata':{'model':v['model'],'pair':v['pair']['pair_id']}}
        if operation=='water':
            return {'states':[{'temperature_k':n['temperature_k'],'pressure_pa':n['pressure_pa'],
                             'status':'ok','delta_h_j_kg':1000.,'storage_density_kg_m3':1000.} for n in v['nodes']]}
        if operation=='states':
            return {'results':[{'mass_fraction':w,'states':[
                {'temperature_k':n['temperature_k'],'pressure_pa':n['pressure_pa'], 'status':'ok',
                 'delta_h_j_kg':1000.*(1+w),'storage_density_kg_m3':1000.} for n in v['nodes']]} for w in v['fractions']]}
        raise ValueError(operation)


class CoreTests(unittest.TestCase):
    def test_default_design_has_40_compositions(self):
        c=yaml.safe_load((MODULE/'mission-campaign.yaml').read_text())
        core.validate(c)
        xs=core.compositions(c['composition'])
        self.assertEqual(len(xs),40); self.assertAlmostEqual(xs[0],.0001); self.assertAlmostEqual(xs[-1],.6)
    def test_invalid_schema(self):
        c=config();c['schema']='bad'
        with self.assertRaises(ValueError):core.validate(c)
    def test_invalid_pump(self):
        c=config();c['fluid']['pump_efficiency']=1.1
        with self.assertRaises(ValueError):core.validate(c)
    def test_missing_budget_is_not_infinite(self):
        c=config();c['screening']['exact_candidates']=0
        with self.assertRaises(ValueError):core.validate(c)
    def test_nan_rejected(self):
        with self.assertRaises(ValueError):core.finite(float('nan'),'x')
    def test_boolean_not_number(self):
        with self.assertRaises(ValueError):core.finite(True,'x')
    def test_equal_composition_bounds_rejected(self):
        c=config();c['composition']['minimum']=.01
        with self.assertRaises(ValueError):core.validate(c)
    def test_catalog_sampling_is_order_invariant(self):
        rows=[pair(i) for i in range(20)]
        self.assertEqual(core.select_catalog(rows,4,9),core.select_catalog(rows[::-1],4,9))
    def test_catalog_all_without_budget(self):
        self.assertEqual(len(core.select_catalog([pair(i) for i in range(7)],None,1)),7)
    def test_quadrature_conserves_heat(self):
        p=profile()
        for n in (1,2,3,99):
            q=core.quadrature(p,n)
            self.assertAlmostEqual(sum(r['heat_j'] for r in q['nodes']),p['summary']['required_heat_j'])
            self.assertTrue({core.coord(r) for r in q['nodes']} <= {core.coord(r) for r in p['demands']})
    def test_exact_preserves_zero_area_peak(self):
        n=core.exact_nodes(profile())
        self.assertEqual(n[2]['peak_flux_w_m2'],2500.)
    def test_empty_profile_rejected(self):
        with self.assertRaises(ValueError):core.quadrature({'demands':[]},2)
    def states(self,nodes):
        return [{'temperature_k':n['temperature_k'],'pressure_pa':n['pressure_pa'],
                 'status':'ok','delta_h_j_kg':1000.,'storage_density_kg_m3':1000.} for n in nodes]
    def test_mass_integral(self):
        p=profile();n=core.exact_nodes(p)
        r=core.score(n,self.states(n),config()['fluid'],profile=p)
        self.assertAlmostEqual(r['mass_kg'],6.);self.assertAlmostEqual(r['peak_mass_flow_kg_s'],1.4)
    def test_missing_state_blocks_score(self):
        n=core.exact_nodes(profile())
        with self.assertRaises(ValueError):core.score(n,self.states(n)[:-1],config()['fluid'])
    def test_duplicate_state_blocks_score(self):
        n=core.exact_nodes(profile());s=self.states(n);s[-1]=s[0]
        with self.assertRaises(ValueError):core.score(n,s,config()['fluid'])
    def test_unsupported_heat_not_renormalized(self):
        n=core.exact_nodes(profile());s=self.states(n);s[0]['status']='unsupported'
        r=core.score(n,s,config()['fluid'])
        self.assertIsNone(r['mass_kg']);self.assertLess(r['heat_coverage_fraction'],1)
    def test_pressure_limit_blocks(self):
        n=core.exact_nodes(profile());f=config()['fluid'];f['max_injection_pressure_pa']=100.
        self.assertIsNone(core.score(n,self.states(n),f)['mass_kg'])
    def test_zero_area_stagnation_can_fail(self):
        n=core.exact_nodes(profile());f=config()['fluid'];f['maximum_mass_flux_kg_m2_s']=2.
        self.assertIsNone(core.score(n,self.states(n),f,profile=profile())['mass_kg'])
    def test_pump_work_reduces_wall_capacity(self):
        n=core.exact_nodes(profile());f=config()['fluid'];f['storage_pressure_pa']=100000.;f['porous_delta_p_pa']=100000.
        r=core.score(n,self.states(n),f)
        self.assertGreater(r['mass_kg'],6.);self.assertGreater(r['pump_work_j'],0)
    def test_single_model_is_not_rejected_for_missing_nrtl(self):
        r=core.combine(pair(),.1,{'unifac_dortmund':{'status':'complete_model_estimate','mass_kg':5.,'heat_coverage_fraction':1.}},6.)
        self.assertEqual(r['status'],'complete_model_estimate');self.assertEqual(r['evidence_tier'],'single_model_exploratory')
    def test_failed_second_model_is_not_ignored(self):
        p=pair();p['models']['chemsep_nrtl']={}
        r=core.combine(p,.1,{'unifac_dortmund':{'status':'complete_model_estimate','mass_kg':5.}},6.)
        self.assertIsNone(r['mass_ratio_vs_heos_water'])
    def test_refinement_generic_not_named(self):
        self.assertEqual(core.refinement_fractions([.1,.2,.3],[.2]),[.15,.2,.25])
    def test_no_unvalidated_promotions(self):
        n=core.exact_nodes(profile());r=core.score(n,self.states(n),config()['fluid'])
        self.assertFalse(r['rankable']);self.assertFalse(r['experimental_winner'])
    def test_high_pressure_is_rejected_before_any_flash(self):
        from mission_thermo import MissionModel
        obj=object.__new__(MissionModel)
        obj.config=config()
        with self.assertRaisesRegex(ValueError,'outside_declared_mission_thermo_domain'):
            obj.state(.1,900.,1e9)
    def test_unknown_entry_override_rejected(self):
        c=config()
        with self.assertRaises(ValueError):entry_configuration({},c)


class WorkflowTests(unittest.TestCase):
    def run_campaign(self,root,c=None,driver=FakeDriver,**kw):
        with redirect_stdout(io.StringIO()):
            return workflow.execute(c or config(),{},root/'cache',root/'results',driver_factory=driver,**kw)
    def test_end_to_end(self):
        with tempfile.TemporaryDirectory() as t:
            summary,path=self.run_campaign(Path(t))
            self.assertEqual(summary['screened_pair_count'],2)
            self.assertEqual(summary['exact_replay_candidates'],2)
            self.assertEqual(summary['failed_task_count'],0)
            self.assertEqual(summary['experimental_winners'],0)
            self.assertTrue((path/'mission-profile.json.gz').is_file())
    def test_second_run_uses_cached_numerical_work(self):
        with tempfile.TemporaryDirectory() as t:
            root=Path(t);self.run_campaign(root);FakeDriver.calls=[]
            summary,_=self.run_campaign(root)
            self.assertNotIn('states',FakeDriver.calls);self.assertNotIn('profile',FakeDriver.calls)
            self.assertGreater(summary['cached_task_count'],0)
    def test_full_recompute_matches_scientific_results(self):
        with tempfile.TemporaryDirectory() as t:
            root=Path(t);a,_=self.run_campaign(root);b,_=self.run_campaign(root,recompute=True)
            self.assertEqual(a['scientific_result_sha256'],b['scientific_result_sha256'])
    def test_thermo_change_invalidates_states_not_profile(self):
        class Changed(FakeDriver):thermo_code={'version':'2'}
        with tempfile.TemporaryDirectory() as t:
            root=Path(t);self.run_campaign(root);FakeDriver.calls=[]
            self.run_campaign(root,driver=Changed)
            self.assertIn('states',FakeDriver.calls);self.assertNotIn('profile',FakeDriver.calls)
    def test_published_shards_are_bounded_and_gzip_roundtrips(self):
        with tempfile.TemporaryDirectory() as t:
            _,path=self.run_campaign(Path(t))
            for p in path.glob('report-*.csv'): self.assertLessEqual(p.stat().st_size,32768)
            self.assertEqual(json.loads(gzip.decompress((path/'summary.json.gz').read_bytes())),read_json(path/'summary.json'))
    def test_limit_does_not_replace_latest(self):
        with tempfile.TemporaryDirectory() as t:
            root=Path(t);self.run_campaign(root);old=(root/'results/latest.json').read_bytes()
            self.run_campaign(root,limit_pairs=1)
            self.assertEqual(old,(root/'results/latest.json').read_bytes())
    def test_snapshots_not_overwritten(self):
        with tempfile.TemporaryDirectory() as t:
            root=Path(t);_,a=self.run_campaign(root);before=(a/'summary.json').read_bytes()
            _,b=self.run_campaign(root)
            self.assertNotEqual(a,b);self.assertEqual(before,(a/'summary.json').read_bytes())
    def test_cli_help(self):
        r=subprocess.run([sys.executable,str(MODULE/'run_campaign.py'),'--help'],capture_output=True,text=True)
        self.assertEqual(r.returncode,0,r.stderr);self.assertIn('--limit-pairs',r.stdout)
    def test_legacy_cli_remains(self):
        r=subprocess.run([sys.executable,str(MODULE/'run_campaign.py'),'--mode','properties','--help'],capture_output=True,text=True)
        self.assertEqual(r.returncode,0,r.stderr)
    def test_gitignore_allows_reports_but_not_workers(self):
        with tempfile.TemporaryDirectory() as t:
            r=Path(t);subprocess.run(['git','init','-q',str(r)],check=True)
            (r/'.gitignore').write_bytes((MODULE/'.gitignore').read_bytes())
            for path,ignored in [('mission-v5m4-results/run-x/summary.json',False),('mission-v5m4-results/run-x/pair-0001.json.gz',False),
                                  ('mission-v5m4-results/run-x/worker.log',True),('.campaign-cache/work/x/request.json',True)]:
                result=subprocess.run(['git','check-ignore','-q',path],cwd=r)
                self.assertEqual(result.returncode==0,ignored,path)


HAVE_THERMO=all(importlib.util.find_spec(x) is not None for x in ('thermo','chemicals','CoolProp','rdkit'))
class IntegrationTests(unittest.TestCase):
    @unittest.skipUnless(HAVE_THERMO,'thermo/chemicals/CoolProp not installed')
    def test_ethanol_liquid_storage_and_frozen_gas_900k(self):
        from campaign_backend import legacy_import
        from chemicals.identifiers import get_pubchem_db
        from thermo.interaction_parameters import IPDB
        from mission_thermo import MissionModel
        legacy=legacy_import();c=config();c['fluid']['storage_pressure_pa']=101325.
        item=get_pubchem_db().search_CAS('64-17-5',autoload=True)
        p=legacy.molecular_metadata(item,'64-17-5',c['catalog'])
        p['models']={'chemsep_nrtl':legacy.nrtl_parameters(IPDB,'64-17-5')}
        p['parameter_sha256']=digest(p['models'])
        obj=MissionModel(p,'chemsep_nrtl',c)
        r=obj.evaluate(.1,[{'temperature_k':900.,'pressure_pa':101325.}])[0]
        self.assertEqual(r['status'],'ok',r)
        self.assertEqual(r['gas_branch'],'frozen_species_continuation')
        self.assertFalse(r['chemistry_validated'])
    @unittest.skipUnless(HAVE_THERMO,'thermo/chemicals/CoolProp not installed')
    def test_heos_reference_supports_gas_mission_endpoint(self):
        from mission_thermo import heos_states
        r=heos_states([{'temperature_k':900.,'pressure_pa':101325.}],config()['fluid'])[0]
        self.assertEqual(r['status'],'ok');self.assertGreater(r['delta_h_j_kg'],1e6)
    @unittest.skipUnless(HAVE_THERMO and importlib.util.find_spec('pymsis') is not None,'numerical entry dependencies unavailable')
    def test_full_demand_heat_matches_existing_evaluator(self):
        from mission_profile import build_profile
        c=config();c['fluid']['storage_pressure_pa']=101325.
        base=yaml.safe_load((ROOT/'experiments/entry-evaluator/config.yaml').read_text())
        c['mission']['overrides']['numerics.dt_s']=.2
        result=build_profile(base,c)
        self.assertTrue(result['summary']['profile_complete'])
        sys.path.insert(0,str(ROOT/'experiments/entry-evaluator'))
        from evaluator import evaluate
        cfg=entry_configuration(base,c);cfg['coolant']['coolprop_name']='Water'
        old=evaluate(cfg)
        self.assertAlmostEqual(result['summary']['required_heat_j']/1e6,old['vehicle_coolant_heat_mj'],places=6)


if __name__=='__main__':unittest.main()
