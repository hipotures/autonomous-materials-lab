"""Execution-scheduling tests. No thermodynamic or physical evidence is generated."""
from __future__ import annotations

from pathlib import Path
import sys
import unittest

ROOT = Path(__file__).resolve().parents[3]
MODULE = ROOT / 'experiments/mixture-campaign'
sys.path.insert(0, str(MODULE))

import runtime_tuning as tuning


class RuntimeTuningTests(unittest.TestCase):
    def test_recommendation_reserves_cpu_and_caps_concurrency(self):
        self.assertEqual(tuning.recommend_workers(
            cpus=64, memory_bytes=128*tuning.GIB, max_workers=32,
            worker_memory_mib=768, environ={}), 32)
        self.assertEqual(tuning.recommend_workers(
            cpus=8, memory_bytes=32*tuning.GIB, max_workers=32,
            worker_memory_mib=768, environ={}), 6)

    def test_memory_can_limit_workers(self):
        self.assertEqual(tuning.recommend_workers(
            cpus=64, memory_bytes=8*tuning.GIB, max_workers=64,
            worker_memory_mib=1024, environ={}), 7)

    def test_default_batch_removes_old_fifty_pair_barrier(self):
        self.assertEqual(tuning.recommend_batch_pairs(4), 64)
        self.assertEqual(tuning.recommend_batch_pairs(32), 512)

    def test_defaults_are_injected(self):
        argv, info = tuning.tune_mission_argv([], cpus=64,
            memory_bytes=128*tuning.GIB, environ={})
        self.assertEqual(argv, ['--workers', '32', '--batch-pairs', '512'])
        self.assertEqual(info['injected'], {'workers': 32, 'batch_pairs': 512})

    def test_explicit_cli_is_never_overwritten(self):
        argv, info = tuning.tune_mission_argv(
            ['--workers', '7', '--batch-pairs=123'], cpus=64,
            memory_bytes=128*tuning.GIB, environ={})
        self.assertEqual(argv, ['--workers', '7', '--batch-pairs=123'])
        self.assertEqual(info['injected'], {})

    def test_limit_pairs_keeps_compatibility_budget(self):
        argv, info = tuning.tune_mission_argv(
            ['--limit-pairs', '12'], cpus=16, memory_bytes=32*tuning.GIB, environ={})
        self.assertIn('--workers', argv)
        self.assertNotIn('--batch-pairs', argv)
        self.assertNotIn('batch_pairs', info['injected'])

    def test_environment_overrides_are_explicit(self):
        argv, info = tuning.tune_mission_argv([], cpus=64,
            memory_bytes=128*tuning.GIB,
            environ={'AML_WORKERS': '20', 'AML_BATCH_PAIRS': '700'})
        self.assertEqual(argv, ['--workers', '20', '--batch-pairs', '700'])
        self.assertEqual(info['injected'], {'workers': 20, 'batch_pairs': 700})

    def test_invalid_environment_override_is_rejected(self):
        with self.assertRaises(ValueError):
            tuning.tune_mission_argv([], cpus=8, memory_bytes=8*tuning.GIB,
                                    environ={'AML_WORKERS': 'zero'})


if __name__ == '__main__':
    unittest.main()
