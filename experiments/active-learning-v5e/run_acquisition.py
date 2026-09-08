#!/usr/bin/env python3
"""V5e-1: select active-learning targets for calibration-domain expansion."""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import platform
import subprocess
import sys
from collections import Counter
from importlib.metadata import version
from pathlib import Path
from typing import Any

import yaml

HERE = Path(__file__).resolve().parent
EXPERIMENTS = HERE.parent
V5B = EXPERIMENTS / "property-predictor-v5b"
if str(V5B) not in sys.path:
    sys.path.insert(0, str(V5B))

from applicability import canonical_smiles  # noqa: E402
from acquisition import (  # noqa: E402
    hypothetical_domain_expansion,
    prepare_features,
    select_targets,
)


def _validate_config(config: dict[str, Any]) -> None:
    acquisition = config["acquisition"]
    target_count = int(acquisition["target_count"])
    minimum_per_family = int(acquisition["minimum_per_family"])
    maximum_family_fraction = float(
        acquisition["maximum_family_fraction"]
    )
    if target_count <= 0:
        raise ValueError("acquisition.target_count must be positive")
    if minimum_per_family < 0:
        raise ValueError(
            "acquisition.minimum_per_family must be non-negative"
        )
    if not 0.0 < maximum_family_fraction <= 1.0:
        raise ValueError(
            "acquisition.maximum_family_fraction must be in (0, 1]"
        )

    weights = {
        key: float(value)
        for key, value in acquisition["weights"].items()
    }
    if any(value < 0.0 for value in weights.values()):
        raise ValueError("acquisition weights must be non-negative")
    if abs(sum(weights.values()) - 1.0) > 1.0e-9:
        raise ValueError("acquisition weights must sum to 1")

    bands = acquisition["novelty_bands"]
    bridge = float(bands["bridge_min_similarity"])
    near = float(bands["near_domain_min_similarity"])
    if not 0.0 <= bridge <= near <= 1.0:
        raise ValueError(
            "novelty-band thresholds must satisfy 0 <= bridge <= near <= 1"
        )

    uncertainty = config["uncertainty"]
    edge = float(uncertainty["edge_similarity"])
    in_domain = float(uncertainty["in_domain_similarity"])
    if not 0.0 <= edge <= in_domain <= 1.0:
        raise ValueError(
            "uncertainty thresholds must satisfy 0 <= edge <= in_domain <= 1"
        )


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _write_json(path: Path, value: Any) -> None:
    with path.open("w", encoding="utf-8") as handle:
        json.dump(value, handle, indent=2, sort_keys=True, allow_nan=False)


def _write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    if not rows:
        path.write_text("", encoding="utf-8")
        return
    fields = sorted({key for row in rows for key in row})
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def _resolve_input(config_path: Path, value: str) -> Path:
    path = Path(value)
    if path.is_absolute():
        return path
    return (config_path.parent / path).resolve()


def _calibration_smiles(summary: dict[str, Any]) -> list[str]:
    result: set[str] = set()
    for row in summary.get("candidates", []):
        smiles = row.get("smiles")
        if not smiles:
            continue
        result.add(canonical_smiles(str(smiles)))
    if not result:
        raise ValueError("V5b summary contains no calibration SMILES")
    return sorted(result)


def _identity_metadata(smiles: str) -> dict[str, Any]:
    from rdkit import Chem
    from rdkit.Chem import Descriptors, rdMolDescriptors

    molecule = Chem.MolFromSmiles(smiles)
    if molecule is None:
        raise ValueError(f"invalid SMILES: {smiles}")
    canonical = Chem.MolToSmiles(molecule, canonical=True)
    try:
        inchi_key = Chem.MolToInchiKey(molecule)
    except Exception:
        inchi_key = None
    return {
        "canonical_smiles": canonical,
        "inchi_key": inchi_key,
        "molecular_formula": rdMolDescriptors.CalcMolFormula(
            molecule
        ),
        "molecular_weight_g_mol": float(
            Descriptors.MolWt(molecule)
        ),
    }


def _entry_by_smiles(payload: dict[str, Any]) -> dict[str, dict[str, Any]]:
    result: dict[str, dict[str, Any]] = {}
    for key in ("ranked", "exploratory"):
        for row in payload.get(key, []):
            smiles = row.get("smiles")
            if not smiles:
                continue
            result[canonical_smiles(str(smiles))] = row
    return result


def _candidate_pool(
    payload: dict[str, Any],
    config: dict[str, Any],
) -> list[dict[str, Any]]:
    acquisition = config["acquisition"]
    require_exploratory = bool(
        acquisition.get("require_exploratory_lane", True)
    )
    minimum_property = float(
        acquisition.get("minimum_property_priority_score", 0.0)
    )

    output: dict[str, dict[str, Any]] = {}
    for row in payload.get("passes", []):
        if require_exploratory and row.get("lane") != "exploratory_domain_expansion":
            continue
        metrics = row.get("property_screen") or {}
        property_score = metrics.get("property_priority_score")
        if not isinstance(property_score, (int, float)):
            continue
        if float(property_score) < minimum_property:
            continue

        smiles = canonical_smiles(str(row["smiles"]))
        normalized = {
            **row,
            "smiles": smiles,
        }
        output.setdefault(smiles, normalized)

    return [
        output[key]
        for key in sorted(output)
    ]


def _family_summary(rows: list[dict[str, Any]]) -> dict[str, int]:
    counts = Counter(str(row["family"]) for row in rows)
    return dict(sorted(counts.items()))


def _band_summary(rows: list[dict[str, Any]]) -> dict[str, int]:
    counts = Counter(str(row["novelty_band"]) for row in rows)
    return dict(sorted(counts.items()))


def _role_summary(rows: list[dict[str, Any]]) -> dict[str, int]:
    counts = Counter(str(row["acquisition_role"]) for row in rows)
    return dict(sorted(counts.items()))


def _pool_csv(features) -> list[dict[str, Any]]:
    return [
        {
            "candidate_id": row.candidate_id,
            "smiles": row.smiles,
            "family": row.family,
            "property_priority_score": row.property_priority_score,
            "property_percentile": row.property_percentile,
            "nearest_calibration_similarity": (
                row.nearest_calibration_similarity
            ),
            "calibration_novelty": row.calibration_novelty,
            "entry_predicted_ratio_vs_water": (
                row.entry_predicted_ratio_vs_water
            ),
            "entry_percentile": row.entry_percentile,
            "acquisition_base_score": row.acquisition_base_score,
            "domain_status": row.domain_status,
            "lane_reason": row.lane_reason,
        }
        for row in features
    ]


def _target_csv(selected: list[dict[str, Any]]) -> list[dict[str, Any]]:
    fields = (
        "acquisition_rank",
        "candidate_id",
        "smiles",
        "canonical_smiles",
        "inchi_key",
        "molecular_formula",
        "molecular_weight_g_mol",
        "family",
        "property_priority_score",
        "property_percentile",
        "nearest_calibration_similarity",
        "calibration_novelty",
        "entry_predicted_ratio_vs_water",
        "entry_percentile",
        "acquisition_base_score",
        "minimum_distance_to_earlier_selected",
        "acquisition_score_at_selection",
        "novelty_band",
        "acquisition_role",
        "domain_status",
        "lane_reason",
        "reference_status",
    )
    return [
        {
            key: row.get(key)
            for key in fields
        }
        for row in selected
    ]


def _reference_resolution_template(
    selected: list[dict[str, Any]],
    config: dict[str, Any],
) -> dict[str, Any]:
    resolution = config["reference_resolution"]
    return {
        "version": "v5e-1",
        "contract": {
            "preferred_sources": list(
                resolution["preferred_sources"]
            ),
            "require_independent_reference_before_v5b_ingest": bool(
                resolution[
                    "require_independent_reference_before_v5b_ingest"
                ]
            ),
            "instruction": (
                "Fill reference_status/reference_source/"
                "reference_identifier only after the molecular identity and "
                "trusted property source have been independently verified. "
                "Do not use the FeOS prediction being calibrated as reference."
            ),
        },
        "targets": [
            {
                "acquisition_rank": row["acquisition_rank"],
                "candidate_id": row["candidate_id"],
                "smiles": row["smiles"],
                "canonical_smiles": row["canonical_smiles"],
                "inchi_key": row["inchi_key"],
                "molecular_formula": row["molecular_formula"],
                "molecular_weight_g_mol": row[
                    "molecular_weight_g_mol"
                ],
                "family": row["family"],
                "reference_status": str(
                    resolution.get(
                        "default_status",
                        "unresolved",
                    )
                ),
                "reference_source": None,
                "reference_identifier": None,
                "reference_fluid_name": None,
                "identity_verified": False,
                "property_reference_verified": False,
                "notes": None,
            }
            for row in selected
        ],
    }


def _v5b_expansion_template(
    selected: list[dict[str, Any]],
) -> dict[str, Any]:
    return {
        "version": "v5e-1-v5b-expansion-template",
        "executable": False,
        "important": (
            "This becomes executable only after reference-resolution.yaml "
            "contains independently verified reference identities and sources."
        ),
        "holdouts": [
            {
                "id": row["candidate_id"],
                "smiles": row["smiles"],
                "inchi_key": row["inchi_key"],
                "molecular_formula": row["molecular_formula"],
                "challenge_family": row["family"],
                "acquisition_rank": row["acquisition_rank"],
                "storage_temperature_k": 293.15,
                "storage_pressure_pa": 101325.0,
                "max_exit_temperature_k": 500.0,
                "reference_backend": None,
                "reference_identifier": None,
                "reference_status": "unresolved",
            }
            for row in selected
        ],
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--config",
        type=Path,
        default=HERE / "config.yaml",
    )
    parser.add_argument(
        "--v5d-full-prescreen",
        type=Path,
        default=None,
    )
    parser.add_argument(
        "--v5d-entry-results",
        type=Path,
        default=None,
    )
    parser.add_argument(
        "--v5b-summary",
        type=Path,
        default=None,
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("active-learning-v5e1-results"),
    )
    args = parser.parse_args()

    if not args.config.is_file():
        parser.error(f"missing config: {args.config}")

    config = yaml.safe_load(
        args.config.read_text(encoding="utf-8")
    )
    try:
        _validate_config(config)
    except (KeyError, TypeError, ValueError) as exc:
        parser.error(f"invalid V5e config: {exc}")
    configured_inputs = config["inputs"]

    full_prescreen_path = (
        args.v5d_full_prescreen
        if args.v5d_full_prescreen is not None
        else _resolve_input(
            args.config,
            str(configured_inputs["v5d_full_prescreen"]),
        )
    )
    entry_results_path = (
        args.v5d_entry_results
        if args.v5d_entry_results is not None
        else _resolve_input(
            args.config,
            str(configured_inputs["v5d_entry_results"]),
        )
    )
    v5b_summary_path = (
        args.v5b_summary
        if args.v5b_summary is not None
        else _resolve_input(
            args.config,
            str(configured_inputs["v5b_summary"]),
        )
    )

    for path in (
        full_prescreen_path,
        entry_results_path,
        v5b_summary_path,
    ):
        if not path.is_file():
            parser.error(f"missing required input: {path}")

    if args.output_dir.exists() and any(args.output_dir.iterdir()):
        parser.error(
            "output directory must be empty; use a new directory"
        )
    args.output_dir.mkdir(parents=True, exist_ok=True)

    full_prescreen = json.loads(
        full_prescreen_path.read_text(encoding="utf-8")
    )
    entry_results = json.loads(
        entry_results_path.read_text(encoding="utf-8")
    )
    v5b_summary = json.loads(
        v5b_summary_path.read_text(encoding="utf-8")
    )

    pool = _candidate_pool(full_prescreen, config)
    if not pool:
        parser.error("no exploratory V5d candidates satisfy acquisition gate")

    calibration_smiles = _calibration_smiles(v5b_summary)
    entry_map = _entry_by_smiles(entry_results)
    features = prepare_features(
        pool,
        calibration_smiles,
        entry_map,
        config,
    )
    if not features:
        parser.error("no acquisition features could be constructed")

    selected = select_targets(features, config)
    if not selected:
        parser.error("active-learning selector returned no targets")
    selected = [
        {
            **row,
            **_identity_metadata(str(row["smiles"])),
        }
        for row in selected
    ]

    expansion = hypothetical_domain_expansion(
        features,
        calibration_smiles,
        selected,
        config,
    )

    pool_rows = _pool_csv(features)
    target_rows = _target_csv(selected)
    _write_csv(
        args.output_dir / "acquisition-pool.csv",
        pool_rows,
    )
    _write_csv(
        args.output_dir / "calibration-targets.csv",
        target_rows,
    )

    _write_json(
        args.output_dir / "acquisition-set.json",
        {
            "version": "v5e-1",
            "targets": selected,
        },
    )
    _write_json(
        args.output_dir / "domain-expansion-whatif.json",
        expansion,
    )

    reference_template = _reference_resolution_template(
        selected,
        config,
    )
    with (
        args.output_dir / "reference-resolution.yaml"
    ).open("w", encoding="utf-8") as handle:
        yaml.safe_dump(
            reference_template,
            handle,
            sort_keys=False,
            allow_unicode=True,
        )

    expansion_template = _v5b_expansion_template(selected)
    with (
        args.output_dir / "v5b-expansion-template.yaml"
    ).open("w", encoding="utf-8") as handle:
        yaml.safe_dump(
            expansion_template,
            handle,
            sort_keys=False,
            allow_unicode=True,
        )

    selected_with_entry = sum(
        row.get("entry_predicted_ratio_vs_water") is not None
        for row in selected
    )
    summary = {
        "study_version": "v5e-1",
        "study_complete": True,
        "candidate_pool_count": len(features),
        "target_count_requested": int(
            config["acquisition"]["target_count"]
        ),
        "target_count_selected": len(selected),
        "target_family_counts": _family_summary(selected),
        "target_novelty_band_counts": _band_summary(selected),
        "target_role_counts": _role_summary(selected),
        "selected_with_v5d_entry_evaluation_count": (
            selected_with_entry
        ),
        "domain_expansion_whatif": expansion,
        "top_targets": selected[:10],
        "handoff_status": (
            "reference_resolution_required_before_v5b_ingest"
        ),
        "interpretation": (
            "V5e-1 selects structurally diverse calibration targets. "
            "Selection alone does not validate their properties and does not "
            "expand the calibrated uncertainty domain. Trusted independent "
            "reference data must be resolved before V5b ingestion."
        ),
    }
    _write_json(
        args.output_dir / "summary.json",
        summary,
    )

    try:
        git_commit = subprocess.check_output(
            ["git", "rev-parse", "HEAD"],
            cwd=HERE,
            text=True,
        ).strip()
        git_dirty = bool(
            subprocess.check_output(
                ["git", "status", "--porcelain"],
                cwd=HERE,
                text=True,
            ).strip()
        )
    except (OSError, subprocess.CalledProcessError):
        git_commit = None
        git_dirty = None

    manifest = {
        "study_version": "v5e-1",
        "git_commit": git_commit,
        "git_dirty": git_dirty,
        "python": platform.python_version(),
        "packages": {
            package: version(package)
            for package in ("rdkit", "PyYAML")
        },
        "inputs": {
            "config_sha256": _sha256(args.config),
            "v5d_full_prescreen": str(full_prescreen_path),
            "v5d_full_prescreen_sha256": _sha256(
                full_prescreen_path
            ),
            "v5d_entry_results": str(entry_results_path),
            "v5d_entry_results_sha256": _sha256(
                entry_results_path
            ),
            "v5b_summary": str(v5b_summary_path),
            "v5b_summary_sha256": _sha256(v5b_summary_path),
        },
        "acquisition_contract": {
            "candidate_source": (
                "V5d full-prescreen exploratory pass set"
            ),
            "primary_terms": [
                "thermodynamic property merit percentile",
                "novelty relative to V5b calibration",
                "greedy diversity relative to already selected targets",
                "small bonus for evaluated nominal entry merit",
            ],
            "family_balance": {
                "minimum_per_family": int(
                    config["acquisition"]["minimum_per_family"]
                ),
                "maximum_family_fraction": float(
                    config["acquisition"][
                        "maximum_family_fraction"
                    ]
                ),
            },
            "reference_data": (
                "unresolved at acquisition time; independent trusted source "
                "required before calibration ingest"
            ),
        },
    }
    _write_json(
        args.output_dir / "manifest.json",
        manifest,
    )

    print(
        json.dumps(
            {
                "study_complete": True,
                "candidate_pool_count": len(features),
                "target_count_selected": len(selected),
                "target_family_counts": summary[
                    "target_family_counts"
                ],
                "target_novelty_band_counts": summary[
                    "target_novelty_band_counts"
                ],
                "target_role_counts": summary[
                    "target_role_counts"
                ],
                "selected_with_v5d_entry_evaluation_count": (
                    selected_with_entry
                ),
                "domain_expansion_whatif": expansion,
                "top_targets": [
                    {
                        "rank": row["acquisition_rank"],
                        "candidate_id": row["candidate_id"],
                        "smiles": row["smiles"],
                        "family": row["family"],
                        "property_priority_score": row[
                            "property_priority_score"
                        ],
                        "nearest_calibration_similarity": row[
                            "nearest_calibration_similarity"
                        ],
                        "entry_predicted_ratio_vs_water": row[
                            "entry_predicted_ratio_vs_water"
                        ],
                        "acquisition_role": row[
                            "acquisition_role"
                        ],
                    }
                    for row in selected[:10]
                ],
            },
            indent=2,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
