"""Isolated numerical batches for catalog, entry demands and thermodynamic states."""
from __future__ import annotations
from importlib import import_module
from importlib.metadata import version
from pathlib import Path
import os
import subprocess
import sys
import uuid
from campaign_backend import runtime, LEGACY
from campaign_store import implementation, file_hash, digest, read_json, write_json

HERE = Path(__file__).resolve().parent
ENTRY = HERE.parent / 'entry-evaluator'


class Driver:
    def __init__(self, work, timeout):
        self.work, self.timeout = Path(work), timeout
        self.environment = runtime()
        root = Path(import_module('pymsis').__file__).resolve().parent
        self.environment['packages']['pymsis'] = version('pymsis')
        self.environment['library_tree_sha256']['pymsis'] = digest({str(p.relative_to(root)): file_hash(p)
            for p in sorted(root.rglob('*')) if p.is_file() and '__pycache__' not in p.parts})
        self.thermo_code = implementation(HERE/'mission_thermo.py', HERE/'mission_core.py', Path(__file__),
            LEGACY/'v5m2_backend.py', LEGACY/'v5m2_core.py')
        self.catalog_code = implementation(HERE/'mission_catalog.py', HERE/'mission_core.py',
            LEGACY/'v5m2_backend.py', LEGACY/'v5m2_core.py', Path(__file__))
        self.profile_code = implementation(HERE/'mission_profile.py', Path(__file__),
            *[ENTRY/p for p in ('atmosphere.py', 'trajectory.py', 'heating.py', 'surface.py', 'wall.py',
                               'coolant.py', 'chemistry.py', 'property_provider.py')])

    def call(self, operation, **values):
        folder = self.work/uuid.uuid4().hex
        folder.mkdir(parents=True)
        inp, out = folder/'request.json', folder/'response.json'
        write_json(inp, {'operation': operation, **values})
        env = {**os.environ, 'OMP_NUM_THREADS': '1', 'OPENBLAS_NUM_THREADS': '1', 'MKL_NUM_THREADS': '1'}
        with (folder/'worker.log').open('w', encoding='utf-8') as log:
            try:
                proc = subprocess.run([sys.executable, str(Path(__file__).resolve()), str(inp), str(out)],
                    stdout=log, stderr=subprocess.STDOUT, timeout=self.timeout, env=env, check=False)
            except subprocess.TimeoutExpired as exc:
                raise RuntimeError(f'worker_timeout: {folder / "worker.log"}') from exc
        if proc.returncode or not out.is_file():
            raise RuntimeError(f'worker_failed: {folder / "worker.log"}')
        return read_json(out)


def worker(r):
    operation = r['operation']
    if operation == 'catalog':
        from mission_catalog import discover
        return discover(r['config'])
    if operation == 'profile':
        from mission_profile import build_profile
        return build_profile(r['entry'], r['config'])
    if operation == 'water':
        from mission_thermo import heos_states
        return {'states': heos_states(r['nodes'], r['fluid'])}
    from mission_thermo import MissionModel
    obj = MissionModel(r['pair'], r['model'], r['config'])
    if operation == 'freeze':
        return {'metadata': obj.metadata}
    if operation != 'states':
        raise ValueError('unknown mission worker operation')
    if obj.metadata != r['frozen']['metadata']:
        raise ValueError('selected correlations changed after freeze')
    return {'results': [{'mass_fraction': w, 'states': obj.evaluate(w, r['nodes'])} for w in r['fractions']]}


if __name__ == '__main__':
    if len(sys.argv) != 3:
        raise SystemExit('internal mission worker requires input and output files')
    write_json(Path(sys.argv[2]), worker(read_json(Path(sys.argv[1]))))
