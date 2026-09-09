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

from reference_resolver import (  # noqa: E402
    CoolPropReference,
    HttpClient,
    SourceUnavailable,
    classify_resolution,
    extract_cas_numbers,
    identity_consistency,
    match_coolprop_fluid,
    resolve_coolprop_reference,
    resolve_nist_reference,
    resolve_pubchem_identity,
    valid_cas,
)
from run_reference_resolution import (  # noqa: E402
    _targets_from_resolution_template,
    _v5b3_handoff,
    _validate_config,
)


class FakeHttp(HttpClient):
    def __init__(
        self,
        *,
        json_by_url=None,
        text_by_url=None,
    ):
        self.json_by_url = (
            json_by_url or {}
        )
        self.text_by_url = (
            text_by_url or {}
        )

    def get_json(self, url: str):
        for needle, payload in (
            self.json_by_url.items()
        ):
            if needle in url:
                return copy.deepcopy(
                    payload
                )
        raise SourceUnavailable(
            f"no fake JSON for {url}"
        )

    def get_text(self, url: str):
        for needle, payload in (
            self.text_by_url.items()
        ):
            if needle in url:
                return payload
        raise SourceUnavailable(
            f"no fake text for {url}"
        )


class FakeCoolProp:
    def __init__(self):
        self.fluids = {
            "TestFluid": {
                "CAS": "64-17-5",
                "INCHI_Key": (
                    "LFQSCWFLJHTTHZ-"
                    "UHFFFAOYSA-N"
                ),
            },
            "OtherFluid": {
                "CAS": "67-64-1",
                "INCHI_Key": (
                    "CSCPPACGZOOCGX-"
                    "UHFFFAOYSA-N"
                ),
            },
        }

    def get_global_param_string(
        self,
        key,
    ):
        if key == "FluidsList":
            return ",".join(
                self.fluids
            )
        if key == "version":
            return "8.0.0-test"
        if key == "gitrevision":
            return "test-revision"
        raise ValueError(key)

    def get_fluid_param_string(
        self,
        fluid,
        key,
    ):
        return self.fluids[fluid][key]

    def PropsSI(
        self,
        output,
        *args,
    ):
        if len(args) == 1:
            values = {
                "Tcrit": 514.0,
                "Pcrit": 6.14e6,
                "acentric": 0.64,
            }
            return values[output]

        key1, value1, key2, value2, fluid = args
        self._check(fluid)
        pair = {
            key1: value1,
            key2: value2,
        }
        if (
            "Q" in pair
            and "T" in pair
            and output == "P"
        ):
            return 8.0e4
        if (
            "Q" in pair
            and "P" in pair
        ):
            if output == "T":
                return 351.4
            if output == "Hmass":
                return (
                    1.0e5
                    if int(pair["Q"]) == 0
                    else 9.0e5
                )
        if (
            "T" in pair
            and "P" in pair
        ):
            temperature = float(
                pair["T"]
            )
            if output == "Dmass":
                return 789.0
            if output == "Cpmass":
                return 2400.0
            if output == "Hmass":
                return (
                    temperature
                    * 2400.0
                )
        raise ValueError(
            (output, args)
        )

    def PhaseSI(
        self,
        key1,
        value1,
        key2,
        value2,
        fluid,
    ):
        self._check(fluid)
        return "liquid"

    def _check(self, fluid):
        if fluid not in self.fluids:
            raise ValueError(fluid)


class V5e2IdentityTests(
    unittest.TestCase
):
    def test_cas_validation_and_extraction(
        self,
    ):
        self.assertTrue(
            valid_cas("64-17-5")
        )
        self.assertTrue(
            valid_cas("67-64-1")
        )
        self.assertFalse(
            valid_cas("64-17-6")
        )
        self.assertEqual(
            extract_cas_numbers(
                [
                    "ethanol",
                    "64-17-5",
                    "64-17-6",
                ]
            ),
            ("64-17-5",),
        )

    def test_pubchem_exact_inchikey_resolution(
        self,
    ):
        key = (
            "LFQSCWFLJHTTHZ-"
            "UHFFFAOYSA-N"
        )
        http = FakeHttp(
            json_by_url={
                "/property/": {
                    "PropertyTable": {
                        "Properties": [
                            {
                                "CID": 702,
                                "Title": (
                                    "Ethanol"
                                ),
                                "IUPACName": (
                                    "ethanol"
                                ),
                                "CanonicalSMILES": (
                                    "CCO"
                                ),
                                "IsomericSMILES": (
                                    "CCO"
                                ),
                                "InChIKey": key,
                                "MolecularFormula": (
                                    "C2H6O"
                                ),
                                "MolecularWeight": (
                                    "46.07"
                                ),
                            }
                        ]
                    }
                },
                "/synonyms/": {
                    "InformationList": {
                        "Information": [
                            {
                                "Synonym": [
                                    "Ethanol",
                                    "64-17-5",
                                ]
                            }
                        ]
                    }
                },
            }
        )
        identity = (
            resolve_pubchem_identity(
                key,
                http,
            )
        )
        self.assertIsNotNone(
            identity
        )
        self.assertEqual(
            identity.cid,
            702,
        )
        self.assertTrue(
            identity.exact_inchi_key_match
        )
        self.assertEqual(
            identity.cas_numbers,
            ("64-17-5",),
        )

        consistency = (
            identity_consistency(
                expected_inchi_key=key,
                expected_formula="C2H6O",
                expected_molecular_weight_g_mol=46.069,
                pubchem=identity,
                max_molecular_weight_relative_error=0.002,
            )
        )
        self.assertTrue(
            consistency[
                "identity_verified"
            ]
        )

    def test_identity_mismatch_cannot_be_verified(
        self,
    ):
        consistency = (
            identity_consistency(
                expected_inchi_key=(
                    "AAAAAAAAAAAAAA-"
                    "BBBBBBBBBB-C"
                ),
                expected_formula="C2H6O",
                expected_molecular_weight_g_mol=46.069,
                pubchem=None,
                max_molecular_weight_relative_error=0.002,
            )
        )
        self.assertFalse(
            consistency[
                "identity_verified"
            ]
        )


class V5e2CoolPropTests(
    unittest.TestCase
):
    def setUp(self):
        self.cp = FakeCoolProp()

    def test_identifier_match_prefers_exact_inchikey(
        self,
    ):
        match = match_coolprop_fluid(
            inchi_key=(
                "LFQSCWFLJHTTHZ-"
                "UHFFFAOYSA-N"
            ),
            cas_numbers=[
                "64-17-5"
            ],
            cp=self.cp,
        )
        self.assertEqual(
            match,
            (
                "TestFluid",
                "inchi_key_exact",
            ),
        )

    def test_coolprop_reference_computes_liquid_storage_and_enthalpy_window(
        self,
    ):
        reference = (
            resolve_coolprop_reference(
                inchi_key=(
                    "LFQSCWFLJHTTHZ-"
                    "UHFFFAOYSA-N"
                ),
                cas_numbers=[
                    "64-17-5"
                ],
                cp=self.cp,
                storage_temperature_k=293.15,
                default_storage_pressure_pa=50000.0,
                max_storage_pressure_pa=5.0e6,
                saturation_pressure_margin=1.10,
                max_exit_temperature_k=500.0,
            )
        )
        self.assertIsNotNone(
            reference
        )
        self.assertAlmostEqual(
            reference.storage_pressure_pa,
            88000.0,
        )
        self.assertTrue(
            reference.liquid_storage_feasible
        )
        self.assertTrue(
            reference.required_property_complete
        )
        self.assertGreater(
            reference.enthalpy_window_j_kg,
            0.0,
        )
        self.assertAlmostEqual(
            reference.latent_heat_vaporization_j_kg,
            8.0e5,
        )


class V5e2NistTests(
    unittest.TestCase
):
    def test_nist_phase_and_fluid_tables_are_normalized(
        self,
    ):
        compound_html = """
        <html>
        <head>
        <title>Ethanol - NIST Chemistry WebBook</title>
        </head>
        <body>
        CAS Registry Number: 64-17-5
        <table>
          <tr><th>Tboil</th><td>351.44</td><td>K</td></tr>
          <tr><th>Tc</th><td>514.0</td><td>K</td></tr>
          <tr><th>Pc</th><td>61.4</td><td>bar</td></tr>
          <tr><th>Enthalpy of vaporization</th><td>38.6</td><td>kJ/mol</td></tr>
        </table>
        </body>
        </html>
        """
        fluid_html = """
        <html><body><table>
          <tr>
            <th>Temperature (K)</th>
            <th>Pressure (MPa)</th>
            <th>Density (kg/m3)</th>
            <th>Enthalpy (kJ/kg)</th>
            <th>Cp</th>
            <th>Phase</th>
          </tr>
          <tr>
            <td>293.15</td>
            <td>0.101325</td>
            <td>789.2</td>
            <td>105.0</td>
            <td>2.44</td>
            <td>liquid</td>
          </tr>
        </table></body></html>
        """
        http = FakeHttp(
            text_by_url={
                "cbook.cgi": (
                    compound_html
                ),
                "fluid.cgi": fluid_html,
            }
        )
        ref = resolve_nist_reference(
            cas_number="64-17-5",
            http=http,
            storage_temperature_k=293.15,
            storage_pressure_pa=101325.0,
        )
        self.assertIsNotNone(ref)
        self.assertTrue(
            ref.compound_page_verified
        )
        self.assertTrue(
            ref.fluid_table_verified
        )
        self.assertAlmostEqual(
            ref.critical_pressure_pa,
            6.14e6,
        )
        self.assertAlmostEqual(
            ref.storage_density_kg_m3,
            789.2,
        )
        self.assertAlmostEqual(
            ref.storage_cp_j_kg_k,
            2440.0,
        )


class V5e2ClassificationTests(
    unittest.TestCase
):
    def test_only_complete_independent_reference_is_auto_usable(
        self,
    ):
        complete = CoolPropReference(
            fluid_name="TestFluid",
            cas_number="64-17-5",
            inchi_key="KEY",
            match_method="cas_exact",
            version="8",
            git_revision="rev",
            storage_temperature_k=293.15,
            storage_pressure_pa=101325.0,
            saturation_pressure_pa=50000.0,
            storage_phase="liquid",
            max_exit_temperature_k=500.0,
            critical_temperature_k=500.0,
            critical_pressure_pa=5e6,
            acentric_factor=0.5,
            normal_boiling_temperature_k=350.0,
            storage_density_kg_m3=800.0,
            storage_cp_j_kg_k=2000.0,
            latent_heat_vaporization_j_kg=5e5,
            storage_enthalpy_j_kg=1e5,
            exit_enthalpy_j_kg=8e5,
            enthalpy_window_j_kg=7e5,
            liquid_storage_feasible=True,
            required_property_complete=True,
            errors=(),
        )
        status, quality, usable = (
            classify_resolution(
                identity_verified=True,
                coolprop=complete,
                nist=None,
            )
        )
        self.assertEqual(
            status,
            "resolved",
        )
        self.assertEqual(
            quality,
            "high",
        )
        self.assertTrue(usable)

        status, _, usable = (
            classify_resolution(
                identity_verified=False,
                coolprop=complete,
                nist=None,
            )
        )
        self.assertEqual(
            status,
            "unresolved_identity",
        )
        self.assertFalse(usable)


class V5e2RunnerTests(
    unittest.TestCase
):
    @classmethod
    def setUpClass(cls):
        cls.config = yaml.safe_load(
            (
                V5E
                / "config-v5e2.yaml"
            ).read_text(
                encoding="utf-8"
            )
        )

    def test_default_config_is_valid(
        self,
    ):
        _validate_config(
            self.config
        )

    def test_invalid_reference_pressure_margin_is_rejected(
        self,
    ):
        config = copy.deepcopy(
            self.config
        )
        config["properties"][
            "saturation_pressure_margin"
        ] = 1.0
        with self.assertRaises(
            ValueError
        ):
            _validate_config(config)

    def test_target_template_requires_identity_fields(
        self,
    ):
        payload = {
            "targets": [
                {
                    "candidate_id": "a",
                    "smiles": "CCO",
                    "canonical_smiles": (
                        "CCO"
                    ),
                    "inchi_key": "KEY",
                    "molecular_formula": (
                        "C2H6O"
                    ),
                    "molecular_weight_g_mol": (
                        46.0
                    ),
                    "family": (
                        "oxygenated"
                    ),
                    "acquisition_rank": (
                        1
                    ),
                }
            ]
        }
        targets = (
            _targets_from_resolution_template(
                payload
            )
        )
        self.assertEqual(
            targets[0][
                "candidate_id"
            ],
            "a",
        )

    def test_v5b3_handoff_excludes_nonusable_targets(
        self,
    ):
        config = copy.deepcopy(
            self.config
        )
        config["handoff"][
            "minimum_ready_targets"
        ] = 1
        rows = [
            {
                "candidate_id": "ready",
                "smiles": "CCO",
                "canonical_smiles": (
                    "CCO"
                ),
                "family": (
                    "oxygenated"
                ),
                "reference_quality": (
                    "high"
                ),
                "usable_for_v5b_calibration": (
                    True
                ),
                "pubchem": {
                    "cid": 702,
                    "inchi_key": (
                        "KEY"
                    ),
                },
                "coolprop": {
                    "fluid_name": (
                        "Ethanol"
                    ),
                    "storage_temperature_k": (
                        293.15
                    ),
                    "storage_pressure_pa": (
                        101325.0
                    ),
                    "max_exit_temperature_k": (
                        500.0
                    ),
                    "cas_number": (
                        "64-17-5"
                    ),
                    "match_method": (
                        "cas_exact"
                    ),
                },
            },
            {
                "candidate_id": "review",
                "smiles": "NCO",
                "canonical_smiles": (
                    "NCO"
                ),
                "family": "mixed",
                "reference_quality": (
                    "low"
                ),
                "usable_for_v5b_calibration": (
                    False
                ),
                "pubchem": {
                    "cid": 0,
                    "inchi_key": (
                        "OTHER"
                    ),
                },
                "coolprop": None,
            },
        ]
        handoff = _v5b3_handoff(
            rows,
            config,
        )
        self.assertTrue(
            handoff["executable"]
        )
        self.assertEqual(
            len(
                handoff["holdouts"]
            ),
            1,
        )
        self.assertEqual(
            handoff["holdouts"][0][
                "id"
            ],
            "ready",
        )


if __name__ == "__main__":
    unittest.main()
