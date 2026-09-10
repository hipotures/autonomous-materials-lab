"""Immutable compressed dossiers, resumable stage pointers and verified history.

Only small summaries remain resident across pairs. Scientific cache objects from
V5m-4 are reused by its unchanged Store; publication is not a numerical cache.
"""
from __future__ import annotations

import gzip
import hashlib
import re
from pathlib import Path

from campaign_store import atomic_bytes, digest, encoded, file_hash, read_json, write_json

SHA = re.compile(r'[0-9a-f]{64}')


def checked_sha(value):
    if not isinstance(value, str) or not SHA.fullmatch(value):
        raise ValueError('invalid content hash')
    return value


class Archive:
    def __init__(self, results: Path, cache: Path):
        self.results, self.cache = Path(results), Path(cache)

    def path(self, ref):
        return self.results/'objects'/(checked_sha(ref['content_sha256'])+'.json.gz')

    def put(self, data):
        raw = encoded(data)
        sha = hashlib.sha256(raw).hexdigest()
        target = self.results/'objects'/(sha+'.json.gz')
        if target.exists():
            if read_json(target) != data:
                raise ValueError('immutable dossier collision or corruption')
        else:
            packed = gzip.compress(raw, compresslevel=6, mtime=0)
            if gzip.decompress(packed) != raw:
                raise ValueError('dossier gzip verification failed')
            atomic_bytes(target, packed)
        return {'content_sha256': sha, 'gzip_sha256': file_hash(target),
                'bytes': target.stat().st_size, 'uncompressed_bytes': len(raw)}

    def verify(self, ref):
        target = self.path(ref)
        if not target.is_file() or file_hash(target) != checked_sha(ref['gzip_sha256']):
            raise ValueError('missing or corrupt dossier: '+str(target))
        return target

    def get(self, ref):
        value = read_json(self.verify(ref))
        if digest(value) != ref['content_sha256']:
            raise ValueError('dossier uncompressed hash mismatch')
        return value

    def checkpoint(self, key, ref, brief, tasks):
        body = {'key': checked_sha(key), 'object': ref, 'brief': brief, 'tasks': tasks}
        write_json(self.cache/'mission-v41-pointers'/(key+'.json'),
                   {'data': body, 'sha256': digest(body)})

    def saved(self, key, *, retry_failures=False, recompute=False):
        if recompute:
            return None
        p = self.cache/'mission-v41-pointers'/(checked_sha(key)+'.json')
        if not p.is_file():
            return None
        envelope = read_json(p)
        value = envelope['data']
        if envelope['sha256'] != digest(value) or value['key'] != key:
            raise ValueError('corrupt pair checkpoint')
        if retry_failures and (value['brief'].get('failed_task_count', 0) or value['brief'].get('failed_control_task_count', 0)):
            return None
        if retry_failures and value['brief'].get('study_failures'):
            return None
        # A result for the same numerical key may have changed after an explicit retry.
        for dependency in value['tasks']:
            rp = self.cache/'refs'/(checked_sha(dependency['task_key'])+'.json')
            if not rp.is_file() or read_json(rp).get('result_hash') != dependency['result_hash']:
                return None
        self.verify(value['object'])
        return value


def read_latest(results):
    pointer = Path(results)/'latest.json'
    if not pointer.exists():
        return None
    p = read_json(pointer)
    run = p['run_id']
    if not isinstance(run, str) or not re.fullmatch(r'run-[A-Za-z0-9-]+', run):
        raise ValueError('unsafe historical snapshot name')
    data = read_json(Path(results)/run/'scientific-results.json.gz')
    if digest(data) != checked_sha(p['scientific_result_sha256']):
        raise ValueError('historical snapshot hash mismatch')
    registry = read_json(Path(results)/run/'registry.json.gz')
    if digest(registry) != data['registry_sha256']:
        raise ValueError('historical registry hash mismatch')
    return {'pointer': p, 'scientific': data, 'registry': registry}


def legacy_history(roots):
    """Read identities/compositions only; never import old scores as current evidence."""
    found, sources = {}, []

    def add(row, source, hypothesis=False, fraction=None):
        pid = row.get('pair_id')
        if not pid or not row.get('cas_number'):
            return
        item = found.setdefault(pid, {'pair_id': pid, 'cas_number': row['cas_number'],
                'hypothesis': False, 'fractions': [], 'sources': []})
        if item['cas_number'] != row['cas_number']:
            raise ValueError('conflicting CAS for historical pair identity')
        item['hypothesis'] |= bool(hypothesis)
        if source not in item['sources']:
            item['sources'].append(source)
        if fraction is not None:
            if isinstance(fraction, bool) or not isinstance(fraction, (int, float)) or not 0 <= fraction <= 1:
                raise ValueError('invalid historical composition')
            if fraction not in item['fractions']:
                item['fractions'].append(fraction)

    for root in roots:
        root = Path(root)
        if not root.exists():
            continue
        for run in sorted(root.glob('run-*')):
            if not run.is_dir() or not (run/'summary.json').is_file():
                continue
            source = str(run)
            summary = read_json(run/'summary.json')
            scientific = read_json(run/'scientific-results.json.gz')
            if digest(scientific) != checked_sha(summary['scientific_result_sha256']):
                raise ValueError('legacy scientific hash mismatch: '+source)
            # Both the old property suite and mission V5m-4 write the selected catalog.
            index = read_json(run/'publication-index.json')
            entry = next((r for r in index['files'] if r['file'] == 'catalog.json.gz'), None)
            if entry is None or file_hash(run/'catalog.json.gz') != entry['sha256']:
                raise ValueError('legacy catalog integrity failure: '+source)
            catalog = read_json(run/'catalog.json.gz')
            for row in catalog.get('candidates', []):
                add(row, source)
            for row in scientific.get('pairs', []):
                q = str(row.get('qualification', ''))
                w = (row.get('best_point') or [None])[0]
                add(row, source, 'hypothesis' in q and not q.startswith('no_'), w)
            for row in scientific.get('finalists', []):
                q = str(row.get('qualification', ''))
                add(row, source, 'hypothesis' in q and not q.startswith('no_'), row.get('mass_fraction'))
            sources.append({'path': source, 'scientific_result_sha256': digest(scientific)})
    for item in found.values():
        item['fractions'].sort(); item['sources'].sort()
    return found, sources


def merge_history(old, new):
    result = {k: {**v, 'fractions': list(v['fractions']), 'sources': list(v['sources'])} for k, v in old.items()}
    for key, value in new.items():
        if key not in result:
            result[key] = {**value, 'fractions': list(value['fractions']), 'sources': list(value['sources'])}
        else:
            before = result[key]
            if before['cas_number'] != value['cas_number']:
                raise ValueError('historical identity conflict')
            before['hypothesis'] |= value['hypothesis']
            before['fractions'] = sorted(set(before['fractions']) | set(value['fractions']))
            before['sources'] = sorted(set(before['sources']) | set(value['sources']))
    return result
