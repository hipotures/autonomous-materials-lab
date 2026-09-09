from __future__ import annotations

import copy
import hashlib
from pathlib import Path
import sys
from types import SimpleNamespace
import unittest

import yaml

ENTRY = Path(__file__).resolve().parents[1]
ROOT = ENTRY.parent.parent
V5E = ROOT / "experiments" / "active-learning-v5e"
if str(V5E) not in sys.path:
    sys.path.insert(0, str(V5E))

import gap_acquisition as g
import run_gap_acquisition as runner

CONFIG = yaml.safe_load((V5E / "config-v5e3.yaml").read_text())


class GapAccountingTests(unittest.TestCase):
    def test_current_v5b3_counts_produce_expected_margin_deficits(self):
        summary = {
            "property_results": {
                g.TC: {"calibration_group_count": 6, "evaluation_group_count": 3},
                g.PC: {"calibration_group_count": 6, "evaluation_group_count": 3},
                g.PS: {"calibration_group_count": 8, "evaluation_group_count": 3},
                g.HV: {"calibration_group_count": 6, "evaluation_group_count": 2},
            }
        }
        deficits = g.property_deficits(summary, CONFIG)
        self.assertEqual((deficits[g.TC]["calibration"], deficits[g.TC]["evaluation"]), (6, 2))
        self.assertEqual((deficits[g.PS]["calibration"], deficits[g.PS]["evaluation"]), (4, 2))
        self.assertEqual((deficits[g.HV]["calibration"], deficits[g.HV]["evaluation"]), (6, 3))

    def test_assignment_exactly_matches_v5b3_hash_contract(self):
        stats = {"split_seed": "v5b3-test", "evaluation_fraction": 0.30}
        group = "LFQSCWFLJHTTHZ"
        digest = hashlib.sha256((stats["split_seed"] + ":" + group).encode()).hexdigest()
        expected = "evaluation" if int(digest[:16], 16) / 2**64 < .30 else "calibration"
        self.assertEqual(g.assignment(group, stats), expected)

    def test_strict_anchors_ignore_diagnostics(self):
        audited = {"targets": [{"observations": [
            {"property": g.TC, "smiles": "CCO", "benchmark_eligible": True},
            {"property": g.TC, "smiles": "CCC", "benchmark_eligible": False},
        ]}]}
        self.assertEqual(g.strict_property_anchors(audited, [g.TC])[g.TC], ["CCO"])


class CandidateTests(unittest.TestCase):
    def test_candidate_pool_excludes_prior_acquisition_and_nonexploratory(self):
        payload = {"passes": [
            {"candidate_id": "a", "smiles": "CCO", "lane": "exploratory_domain_expansion",
             "property_screen": {"property_priority_score": .5}, "provenance": [{"family": "oxygenated"}]},
            {"candidate_id": "b", "smiles": "CCC", "lane": "rankable",
             "property_screen": {"property_priority_score": .8}, "provenance": [{"family": "alkane"}]},
            {"candidate_id": "c", "smiles": "CCCC", "lane": "exploratory_domain_expansion",
             "property_screen": {"property_priority_score": .1}, "provenance": [{"family": "alkane"}]},
        ]}
        self.assertEqual(g.candidate_pool(payload, excluded_smiles={"CCO"}, config=CONFIG), [])

    def test_preprobe_feature_uses_candidate_split_not_global_average(self):
        identity = g.canonical_identity("CCO")
        row = {"candidate_id": "a", "smiles": "CCO", "_identity": identity,
               "_property_priority_score": .5, "_family": "oxygenated"}
        stats = {"split_seed": "seed", "evaluation_fraction": .3}
        split = g.assignment(identity["molecule_group"], stats)
        deficits = {prop: {"calibration": 0, "evaluation": 0,
                           "target_calibration": 12, "target_evaluation": 5}
                    for prop in CONFIG["priority_properties"]}
        deficits[g.TC][split] = 2
        features = g.prepare_preprobe_features([row], deficits=deficits,
                                               anchors={p: [] for p in deficits},
                                               statistics=stats, config=CONFIG)
        self.assertGreater(features[0].preprobe_score, 0.0)
        self.assertEqual(features[0].split, split)


class FakeBackend:
    def __init__(self, key):
        self.key = key
        self.constant_methods = {
            g.TC: lambda cas: ["IUPAC"],
            g.PC: lambda cas: ["IUPAC"],
            "normal_boiling_temperature_k": lambda cas: ["CRC_ORG"],
        }
    def identity(self, cas):
        return {"cas": cas, "inchi_key": self.key, "smiles": "CCO", "formula": "C2H6O", "mw": 46.069}
    def constant(self, cas, prop, method):
        return {g.TC: 514.0, g.PC: 6.14e6, "normal_boiling_temperature_k": 351.44}[prop]
    def object(self, cas, prop, packet):
        return SimpleNamespace(all_methods={"ANTOINE_POLING"} if prop == g.PS else {"CRC_HVAP_TB"})
    def native_point(self, cas, prop, method, packet):
        return 351.44, 38600.0
    def series_value(self, cas, prop, method, temperature_k, packet):
        return 10000.0, [280.0, 400.0]


def fake_allowed(prop, method):
    return method in {"IUPAC", "CRC_ORG", "ANTOINE_POLING", "CRC_HVAP_TB"}


def fake_structural_identity(smiles):
    return g.canonical_identity(smiles)


class AvailabilityTests(unittest.TestCase):
    def test_exact_identity_and_local_cas_expose_strict_property_inventory(self):
        identity = g.canonical_identity("CCO")
        pubchem = {"inchi_key": identity["inchi_key"], "molecular_formula": identity["molecular_formula"],
                   "molecular_weight_g_mol": identity["molecular_weight_g_mol"], "cas_numbers": ["64-17-5"]}
        result = g.empirical_availability(candidate_identity=identity, pubchem=pubchem,
                                          backend=FakeBackend(identity["inchi_key"]), allowed_fn=fake_allowed,
                                          structural_identity_fn=fake_structural_identity, config=CONFIG)
        self.assertTrue(result["identity_verified"])
        self.assertEqual(result["selected_cas_number"], "64-17-5")
        self.assertEqual(set(result["available_properties"]), {g.TC, g.PC, g.PS, g.HV})

    def test_local_cas_structure_mismatch_cannot_be_selected(self):
        identity = g.canonical_identity("CCO")
        pubchem = {"inchi_key": identity["inchi_key"], "molecular_formula": identity["molecular_formula"],
                   "molecular_weight_g_mol": identity["molecular_weight_g_mol"], "cas_numbers": ["64-17-5"]}
        result = g.empirical_availability(candidate_identity=identity, pubchem=pubchem,
                                          backend=FakeBackend("AAAAAAAAAAAAAA-UHFFFAOYSA-N"), allowed_fn=fake_allowed,
                                          structural_identity_fn=fake_structural_identity, config=CONFIG)
        self.assertIsNone(result["selected_cas_number"])
        self.assertEqual(result["available_properties"], [])


class SelectionTests(unittest.TestCase):
    def setUp(self):
        self.config = copy.deepcopy(CONFIG)
        self.config["acquisition"]["target_count"] = 4
        self.config["acquisition"]["maximum_family_fraction"] = .5
        self.deficits = {prop: {"current_calibration": 0, "current_evaluation": 0,
                                "target_calibration": 2, "target_evaluation": 1,
                                "calibration": 2, "evaluation": 1}
                         for prop in self.config["priority_properties"]}
        self.anchors = {prop: ["CCO"] for prop in self.config["priority_properties"]}

    def row(self, candidate_id, smiles, family, split, props, merit):
        ident = g.canonical_identity(smiles)
        return {"candidate_id": candidate_id, "smiles": ident["canonical_smiles"], "family": family,
                "split": split, "molecule_group": ident["molecule_group"], "property_percentile": merit,
                "identity_verified": True, "selected_cas_number": "64-17-5",
                "available_properties": props, "property_methods": {p: ["IUPAC"] for p in props},
                "inchi_key": ident["inchi_key"], "canonical_smiles": ident["canonical_smiles"]}

    def test_multi_gap_candidate_is_preferred(self):
        rows = [
            self.row("multi", "CCCC", "a", "calibration", [g.TC, g.PC, g.PS, g.HV], .5),
            self.row("single", "CCCCC", "b", "calibration", [g.TC], 1.0),
        ]
        selected, _ = g.select_gap_targets(rows, deficits=self.deficits, anchors=self.anchors, config=self.config)
        self.assertEqual(selected[0]["candidate_id"], "multi")

    def test_selection_reduces_only_candidate_split_deficits(self):
        rows = [self.row("multi", "CCCC", "a", "evaluation", [g.TC, g.PC], .5)]
        _, remaining = g.select_gap_targets(rows, deficits=self.deficits, anchors=self.anchors, config=self.config)
        self.assertEqual(remaining[g.TC]["evaluation"], 0)
        self.assertEqual(remaining[g.TC]["calibration"], 2)

    def test_family_cap_is_respected(self):
        rows = [
            self.row("a1", "CCCC", "same", "calibration", [g.TC, g.PC, g.PS, g.HV], .9),
            self.row("a2", "CCCCC", "same", "calibration", [g.TC, g.PC, g.PS, g.HV], .8),
            self.row("a3", "CCCCCC", "same", "evaluation", [g.TC, g.PC, g.PS, g.HV], .7),
            self.row("b1", "CCCO", "other", "evaluation", [g.TC, g.PC, g.PS, g.HV], .6),
        ]
        selected, _ = g.select_gap_targets(rows, deficits=self.deficits, anchors=self.anchors, config=self.config)
        counts = {}
        for row in selected:
            counts[row["family"]] = counts.get(row["family"], 0) + 1
        self.assertLessEqual(counts.get("same", 0), 2)

    def test_projection_is_conditional_not_claimed_validation(self):
        remaining = {prop: {"calibration": 1, "evaluation": 0} for prop in self.deficits}
        projection = g.gap_projection(self.deficits, remaining)
        self.assertIn("projected_after_if_all_selected_sources_audit_clean", projection[g.TC])
        self.assertEqual(projection[g.TC]["remaining_gap_after_selection"]["calibration"], 1)


class RunnerContractTests(unittest.TestCase):
    def test_default_config_is_valid(self):
        runner._validate_config(CONFIG)

    def test_invalid_weight_sum_is_rejected(self):
        config = copy.deepcopy(CONFIG)
        config["acquisition"]["weights"]["gap_closure"] = 0.51
        with self.assertRaises(ValueError):
            runner._validate_config(config)

    def test_checkpoint_detects_corruption(self):
        import tempfile
        with tempfile.TemporaryDirectory() as folder:
            path = Path(folder) / "check.json"
            runner._checkpoint(path, "sig", {"value": 1})
            saved = runner._load_json(path)
            saved["data"]["value"] = 2
            runner._write_json(path, saved)
            with self.assertRaises(ValueError):
                runner._checkpoint(path, "sig")


if __name__ == "__main__":
    unittest.main()
