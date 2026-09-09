from __future__ import annotations

import copy
from pathlib import Path
import sys
import unittest

import yaml

ENTRY = Path(__file__).resolve().parents[1]
ROOT = ENTRY.parent.parent
V5E = ROOT / "experiments" / "active-learning-v5e"
if str(V5E) not in sys.path:
    sys.path.insert(0, str(V5E))

from empirical_reference import (  # noqa: E402
    build_empirical_packet,
    build_v5b3_property_handoff,
    classify_packet,
    ordered_allowed_methods,
    select_constant,
    select_temperature_series,
    selected_method_is_forbidden,
    selector_matches,
    source_family,
)
from run_reference_resolution_v5e21 import (  # noqa: E402
    _load_targets,
    _validate_config,
)


class FakeTemperatureProperty:
    def __init__(self):
        self.all_methods = {
            "JOBACK",
            "CRCSTD",
            "WEBBOOK_SHOMATE",
        }

    def test_method_validity(self, temperature_k, method):
        if method == "CRCSTD":
            return 250.0 <= temperature_k <= 350.0
        if method == "WEBBOOK_SHOMATE":
            return 250.0 <= temperature_k <= 500.0
        return True

    def calculate(self, temperature_k, method):
        if method == "JOBACK":
            return 999.0
        if method == "CRCSTD":
            return 100.0
        if method == "WEBBOOK_SHOMATE":
            return 120.0 + 0.01 * temperature_k
        raise ValueError(method)


class V5e21PolicyTests(unittest.TestCase):
    def test_selector_matching_is_exact_unless_wildcarded(self):
        self.assertTrue(selector_matches("EOS", "EOS"))
        self.assertFalse(selector_matches("HEOS_FIT", "EOS"))
        self.assertTrue(selector_matches("RACKETTFIT", "RACKETT*"))
        self.assertFalse(selector_matches("HEOS_FIT", "EOS*"))

    def test_ordered_allowlist_does_not_admit_unlisted_methods(self):
        methods = ordered_allowed_methods(
            ["JOBACK", "WEBBOOK_SHOMATE", "CRCSTD"],
            ["WEBBOOK_SHOMATE*", "CRCSTD"],
        )
        self.assertEqual(methods, ["WEBBOOK_SHOMATE", "CRCSTD"])
        self.assertNotIn("JOBACK", methods)

    def test_forbidden_selector_does_not_false_positive_heos(self):
        selectors = ["JOBACK", "EOS", "RACKETT*"]
        self.assertFalse(
            selected_method_is_forbidden("HEOS_FIT", selectors)
        )
        self.assertTrue(
            selected_method_is_forbidden("EOS", selectors)
        )
        self.assertTrue(
            selected_method_is_forbidden("RACKETTFIT", selectors)
        )

    def test_source_family_preserves_reference_provenance(self):
        self.assertEqual(source_family("WEBBOOK"), "NIST")
        self.assertEqual(source_family("ANTOINE_WEBBOOK"), "NIST")
        self.assertEqual(source_family("IUPAC"), "IUPAC")
        self.assertEqual(source_family("DIPPR_PERRY_8E"), "DIPPR_PERRY")


class V5e21SelectionTests(unittest.TestCase):
    def test_constant_selector_skips_predicted_method(self):
        available = lambda cas: ["JOBACK", "WEBBOOK"]  # noqa: E731

        def getter(cas, method=None):
            return {
                "JOBACK": 600.0,
                "WEBBOOK": 500.0,
            }[method]

        record, audit = select_constant(
            property_name="critical_temperature_k",
            unit="K",
            available_methods=available,
            getter=getter,
            cas_number="64-17-5",
            allowed_methods=["WEBBOOK"],
        )
        self.assertEqual(record["method"], "WEBBOOK")
        joback = next(row for row in audit if row["method"] == "JOBACK")
        self.assertFalse(joback["allowed"])
        self.assertFalse(joback["selected"])

    def test_temperature_series_prefers_grid_coverage_then_priority(self):
        obj = FakeTemperatureProperty()
        record, audit = select_temperature_series(
            property_name="liquid_cp_j_kg_k",
            unit="J/kg/K",
            obj=obj,
            temperatures_k=[293.15, 350.0, 400.0, 500.0],
            allowed_methods=["CRCSTD", "WEBBOOK_SHOMATE*"],
            state_basis="liquid",
        )
        self.assertEqual(record["method"], "WEBBOOK_SHOMATE")
        self.assertEqual(len(record["points"]), 4)
        joback = next(row for row in audit if row["method"] == "JOBACK")
        self.assertFalse(joback["allowed"])


class V5e21GateTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.config = yaml.safe_load(
            (V5E / "config-v5e21.yaml").read_text(encoding="utf-8")
        )

    def _packet(self):
        return {
            "constants": {
                "critical_temperature_k": {
                    "value": 500.0,
                    "method": "IUPAC",
                    "source_family": "IUPAC",
                },
                "critical_pressure_pa": {
                    "value": 5.0e6,
                    "method": "IUPAC",
                    "source_family": "IUPAC",
                },
                "normal_boiling_temperature_k": {
                    "value": 350.0,
                    "method": "WEBBOOK",
                    "source_family": "NIST",
                },
            },
            "series": {
                "vapor_pressure_pa": {
                    "method": "ANTOINE_WEBBOOK",
                    "source_family": "NIST",
                    "points": [{"temperature_k": 293.15, "value": 50000.0}],
                },
                "liquid_density_kg_m3": {
                    "method": "DIPPR_PERRY_8E",
                    "source_family": "DIPPR_PERRY",
                    "points": [{"temperature_k": 293.15, "value": 800.0}],
                },
                "liquid_cp_j_kg_k": None,
                "latent_heat_vaporization_j_kg": None,
            },
        }

    def test_medium_packet_passes_calibration_gate(self):
        result = classify_packet(self._packet(), self.config)
        self.assertTrue(result["ready_for_v5b3_property_calibration"])
        self.assertIn(result["reference_quality"], {"medium", "high"})
        self.assertEqual(result["dynamic_property_count"], 2)
        self.assertGreaterEqual(result["source_family_count"], 2)

    def test_missing_pc_blocks_packet(self):
        packet = self._packet()
        packet["constants"]["critical_pressure_pa"] = None
        result = classify_packet(packet, self.config)
        self.assertFalse(result["ready_for_v5b3_property_calibration"])
        self.assertIn(
            "critical_pressure_pa",
            result["missing_required_properties"],
        )

    def test_forbidden_selected_method_blocks_packet(self):
        packet = self._packet()
        packet["series"]["liquid_density_kg_m3"]["method"] = "RACKETTFIT"
        result = classify_packet(packet, self.config)
        self.assertFalse(result["ready_for_v5b3_property_calibration"])
        self.assertEqual(result["forbidden_selected_methods"], ["RACKETTFIT"])


class V5e21RunnerTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.config = yaml.safe_load(
            (V5E / "config-v5e21.yaml").read_text(encoding="utf-8")
        )

    def test_default_config_is_valid(self):
        _validate_config(self.config)

    def test_weaker_high_quality_gate_is_rejected(self):
        config = copy.deepcopy(self.config)
        config["calibration_gate"]["high_quality"][
            "minimum_dynamic_properties"
        ] = 1
        with self.assertRaises(ValueError):
            _validate_config(config)

    def test_target_loader_is_deterministic(self):
        payload = {
            "targets": [
                {"candidate_id": "b", "acquisition_rank": 2},
                {"candidate_id": "a", "acquisition_rank": 1},
            ]
        }
        rows = _load_targets(payload)
        self.assertEqual(
            [row["candidate_id"] for row in rows],
            ["a", "b"],
        )

    def test_handoff_includes_only_property_ready_anchors(self):
        config = copy.deepcopy(self.config)
        config["handoff"]["minimum_ready_targets"] = 1
        rows = [
            {
                "candidate_id": "ready",
                "smiles": "CCO",
                "canonical_smiles": "CCO",
                "inchi_key": "KEY",
                "selected_cas_number": "64-17-5",
                "family": "oxygenated",
                "reference_quality": "medium",
                "ready_for_v5b3_property_calibration": True,
                "entry_reference_ready": False,
                "classification": {
                    "source_families": ["IUPAC", "NIST"],
                    "present_properties": [
                        "critical_temperature_k",
                        "critical_pressure_pa",
                        "vapor_pressure_pa",
                        "liquid_density_kg_m3",
                    ],
                },
            },
            {
                "candidate_id": "partial",
                "smiles": "NCO",
                "canonical_smiles": "NCO",
                "inchi_key": "OTHER",
                "selected_cas_number": None,
                "family": "mixed",
                "reference_quality": "partial",
                "ready_for_v5b3_property_calibration": False,
                "entry_reference_ready": False,
                "classification": {
                    "source_families": [],
                    "present_properties": [],
                },
            },
        ]
        handoff = build_v5b3_property_handoff(rows, config)
        self.assertTrue(handoff["executable"])
        self.assertEqual(len(handoff["anchors"]), 1)
        self.assertEqual(handoff["anchors"][0]["id"], "ready")
        self.assertFalse(
            handoff["anchors"][0]["entry_reference_available"]
        )


class V5e21InstalledDataTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.config = yaml.safe_load(
            (V5E / "config-v5e21.yaml").read_text(encoding="utf-8")
        )

    def test_ethanol_packet_uses_only_allowed_nonpredictive_methods(self):
        try:
            import chemicals  # noqa: F401
            import thermo  # noqa: F401
        except ImportError:
            self.skipTest("chemicals/thermo not installed in this environment")

        packet = build_empirical_packet(
            cas_number="64-17-5",
            molecular_weight_g_mol=46.06844,
            config=self.config,
        )
        classification = classify_packet(packet, self.config)
        self.assertIn(
            "critical_temperature_k",
            classification["present_properties"],
        )
        self.assertIn(
            "critical_pressure_pa",
            classification["present_properties"],
        )
        self.assertEqual(classification["forbidden_selected_methods"], [])
        selected = [
            record["method"]
            for record in list(packet["constants"].values())
            + list(packet["series"].values())
            if record
        ]
        self.assertNotIn("JOBACK", selected)


if __name__ == "__main__":
    unittest.main()
