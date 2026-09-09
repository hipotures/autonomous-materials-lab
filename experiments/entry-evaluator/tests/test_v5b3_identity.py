"""Regression coverage for missing CAS metadata in the V5b-3 audit."""
from __future__ import annotations

import importlib.util
from pathlib import Path
import sys
from types import SimpleNamespace
import unittest
from unittest.mock import Mock

V5B = Path(__file__).resolve().parents[2] / "property-predictor-v5b"
if str(V5B) not in sys.path:
    sys.path.insert(0, str(V5B))
import v5b3_audit as audit

CAS = "64-17-5"
KEY = "LFQSCWFLJHTTHZ-UHFFFAOYSA-N"
AUDIT_CONFIG = {
    "identity_mw_relative_tolerance": 0.002,
    "reproduction_relative_tolerance": 0.0001,
    "constant_conflict_relative_tolerance": 0.10,
}


def target(candidate_id="metadata-test"):
    return {
        "candidate_id": candidate_id,
        "smiles": "CCO",
        "inchi_key": KEY,
        "selected_cas_number": CAS,
        "packet": {
            "cas_number": CAS,
            "molecular_weight_g_mol": 46.069,
            "constants": {
                audit.TC: {"value": 514.0, "unit": "K", "method": "IUPAC"},
            },
            "series": {},
        },
    }


def metadata_record():
    return SimpleNamespace(
        CASs=CAS, InChI_key=KEY, smiles="CCO", formula="C2H6O", MW=46.069,
    )


def backend_with_result(result):
    # Exercise the production identity adapter without importing reference tables.
    backend = audit.LocalReferenceBackend.__new__(audit.LocalReferenceBackend)
    backend.metadata = SimpleNamespace(search_CAS=Mock(return_value=result))
    backend.constant_getters = {audit.TC: lambda cas, method: 514.0}
    backend.constant_methods = {audit.TC: lambda cas: ["IUPAC"]}
    return backend


class LocalMetadataTests(unittest.TestCase):
    def test_false_lookup_is_reported_as_missing_metadata(self):
        backend = backend_with_result(False)
        with self.assertRaisesRegex(ValueError, "^cas_not_in_local_metadata$"):
            backend.identity(CAS)
        backend.metadata.search_CAS.assert_called_once_with(CAS, autoload=True)

    def test_none_lookup_remains_reported_as_missing_metadata(self):
        with self.assertRaisesRegex(ValueError, "^cas_not_in_local_metadata$"):
            backend_with_result(None).identity(CAS)

    def test_found_metadata_is_preserved(self):
        result = backend_with_result(metadata_record()).identity(CAS)
        self.assertEqual(result, {
            "cas": CAS, "inchi_key": KEY, "smiles": "CCO",
            "formula": "C2H6O", "mw": 46.069,
        })

    def test_unexpected_library_error_is_not_disguised_as_missing_data(self):
        backend = backend_with_result(None)
        backend.metadata.search_CAS.side_effect = AttributeError("unexpected API failure")
        with self.assertRaisesRegex(AttributeError, "unexpected API failure"):
            backend.identity(CAS)


@unittest.skipUnless(importlib.util.find_spec("rdkit"), "RDKit is not installed")
class MissingMetadataAuditTests(unittest.TestCase):
    def test_false_lookup_quarantines_target_without_accepting_properties(self):
        result = audit.audit_target(target(), AUDIT_CONFIG, backend_with_result(False))
        self.assertFalse(result["identity_verified"])
        self.assertEqual(result["identity_reason"], "cas_not_in_local_metadata")
        self.assertEqual(result["observations"], [])
        self.assertFalse(result["rankable"])

    def test_none_lookup_quarantines_target_without_accepting_properties(self):
        result = audit.audit_target(target(), AUDIT_CONFIG, backend_with_result(None))
        self.assertFalse(result["identity_verified"])
        self.assertEqual(result["identity_reason"], "cas_not_in_local_metadata")
        self.assertEqual(result["observations"], [])

    def test_next_target_is_audited_after_missing_metadata(self):
        backend = backend_with_result(None)
        backend.metadata.search_CAS.side_effect = [False, metadata_record()]
        results = [audit.audit_target(target(name), AUDIT_CONFIG, backend)
                   for name in ("missing-target", "found-target")]
        self.assertFalse(results[0]["identity_verified"])
        self.assertTrue(results[1]["identity_verified"])
        self.assertEqual(results[0]["observations"], [])
        self.assertEqual(len(results[1]["observations"]), 1)
        self.assertTrue(results[1]["observations"][0]["benchmark_eligible"])

    def test_cas_identity_mismatch_is_still_rejected(self):
        metadata = metadata_record()
        metadata.CASs = "67-64-1"
        result = audit.audit_target(target(), AUDIT_CONFIG, backend_with_result(metadata))
        self.assertFalse(result["identity_verified"])
        self.assertEqual(result["identity_reason"], "cas_identity_mismatch")
        self.assertEqual(result["observations"], [])


@unittest.skipUnless(importlib.util.find_spec("chemicals"), "chemicals is not installed")
class InstalledMetadataTests(unittest.TestCase):
    def test_real_database_missing_lookup_is_handled(self):
        from chemicals.identifiers import get_pubchem_db
        backend = audit.LocalReferenceBackend.__new__(audit.LocalReferenceBackend)
        backend.metadata = get_pubchem_db()
        # A synthetic lookup key, not a claim that a real compound has this CAS.
        cas = "9999999-99-5"
        result = backend.metadata.search_CAS(cas, autoload=True)
        self.assertTrue(result is False or result is None)
        with self.assertRaisesRegex(ValueError, "^cas_not_in_local_metadata$"):
            backend.identity(cas)
        self.assertEqual(backend.identity(CAS)["inchi_key"], KEY)


if __name__ == "__main__":
    unittest.main()
