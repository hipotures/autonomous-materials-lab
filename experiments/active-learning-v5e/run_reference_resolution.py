#!/usr/bin/env python3
"""V5e-2: automatically resolve trusted reference identities and properties."""
from __future__ import annotations

import argparse
import csv
import hashlib
from importlib.metadata import PackageNotFoundError, version
import json
import platform
from pathlib import Path
import subprocess
from typing import Any

import yaml

from reference_resolver import CachedHttpClient, resolve_target


HERE = Path(__file__).resolve().parent


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _write_json(path: Path, value: Any) -> None:
    with path.open("w", encoding="utf-8") as handle:
        json.dump(
            value,
            handle,
            indent=2,
            sort_keys=True,
            allow_nan=False,
        )
        handle.write("\n")


def _write_csv(
    path: Path,
    rows: list[dict[str, Any]],
) -> None:
    if not rows:
        path.write_text("", encoding="utf-8")
        return
    fields = sorted(
        {key for row in rows for key in row}
    )
    with path.open(
        "w",
        newline="",
        encoding="utf-8",
    ) as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=fields,
        )
        writer.writeheader()
        writer.writerows(rows)


def _resolve_input(
    config_path: Path,
    value: str,
) -> Path:
    path = Path(value)
    if path.is_absolute():
        return path
    return (config_path.parent / path).resolve()


def _validate_config(config: dict[str, Any]) -> None:
    if str(config.get("version")) != "v5e-2":
        raise ValueError("config.version must be v5e-2")

    network = config["network"]
    if float(network["timeout_s"]) <= 0.0:
        raise ValueError(
            "network.timeout_s must be positive"
        )
    if int(network["retries"]) < 0:
        raise ValueError(
            "network.retries must be non-negative"
        )
    if float(network["minimum_interval_s"]) < 0.0:
        raise ValueError(
            "network.minimum_interval_s must be non-negative"
        )

    identity = config["identity"]
    mw_error = float(
        identity[
            "max_molecular_weight_relative_error"
        ]
    )
    if not 0.0 <= mw_error <= 0.1:
        raise ValueError(
            "identity.max_molecular_weight_relative_error "
            "must be in [0, 0.1]"
        )

    properties = config["properties"]
    storage_t = float(
        properties["storage_temperature_k"]
    )
    storage_p = float(
        properties["default_storage_pressure_pa"]
    )
    max_p = float(
        properties["max_storage_pressure_pa"]
    )
    margin = float(
        properties["saturation_pressure_margin"]
    )
    exit_t = float(
        properties["max_exit_temperature_k"]
    )
    if storage_t <= 0.0 or exit_t <= storage_t:
        raise ValueError(
            "invalid storage/exit temperatures"
        )
    if storage_p <= 0.0 or max_p < storage_p:
        raise ValueError(
            "invalid storage pressure limits"
        )
    if margin <= 1.0:
        raise ValueError(
            "saturation_pressure_margin must be > 1"
        )

    quality = config["quality"]
    disagreement = float(
        quality[
            "maximum_cross_source_relative_difference"
        ]
    )
    if not 0.0 <= disagreement <= 1.0:
        raise ValueError(
            "quality.maximum_cross_source_relative_difference "
            "must be in [0, 1]"
        )


def _targets_from_resolution_template(
    payload: dict[str, Any],
) -> list[dict[str, Any]]:
    targets = payload.get("targets")
    if not isinstance(targets, list) or not targets:
        raise ValueError(
            "V5e-1 reference-resolution input "
            "contains no targets"
        )
    required = {
        "candidate_id",
        "smiles",
        "canonical_smiles",
        "inchi_key",
        "molecular_formula",
        "molecular_weight_g_mol",
        "family",
    }
    output: list[dict[str, Any]] = []
    seen: set[str] = set()
    for index, row in enumerate(
        targets,
        start=1,
    ):
        if not isinstance(row, dict):
            raise ValueError(
                f"target {index} is not a mapping"
            )
        missing = sorted(
            required - set(row)
        )
        if missing:
            raise ValueError(
                f"target "
                f"{row.get('candidate_id', index)} "
                f"missing: {missing}"
            )
        candidate_id = str(
            row["candidate_id"]
        )
        if candidate_id in seen:
            raise ValueError(
                f"duplicate candidate_id: "
                f"{candidate_id}"
            )
        seen.add(candidate_id)
        output.append(dict(row))
    output.sort(
        key=lambda row: (
            int(
                row.get("acquisition_rank")
                or 10**9
            ),
            str(row["candidate_id"]),
        )
    )
    return output


def _flatten_reference_row(
    row: dict[str, Any],
) -> dict[str, Any]:
    pubchem = row.get("pubchem") or {}
    coolprop = row.get("coolprop") or {}
    nist = row.get("nist") or {}
    identity = row.get("identity") or {}
    return {
        "acquisition_rank": row.get(
            "acquisition_rank"
        ),
        "candidate_id": row.get(
            "candidate_id"
        ),
        "smiles": row.get("smiles"),
        "family": row.get("family"),
        "preferred_name": row.get(
            "preferred_name"
        ),
        "resolution_status": row.get(
            "resolution_status"
        ),
        "reference_quality": row.get(
            "reference_quality"
        ),
        "usable_for_v5b_calibration": row.get(
            "usable_for_v5b_calibration"
        ),
        "identity_verified": identity.get(
            "identity_verified"
        ),
        "pubchem_cid": pubchem.get("cid"),
        "inchi_key": (
            pubchem.get("inchi_key")
            or row.get("expected_inchi_key")
        ),
        "cas_numbers": ";".join(
            pubchem.get("cas_numbers") or []
        ),
        "coolprop_fluid_name": coolprop.get(
            "fluid_name"
        ),
        "coolprop_match_method": coolprop.get(
            "match_method"
        ),
        "storage_temperature_k": coolprop.get(
            "storage_temperature_k"
        ),
        "storage_pressure_pa": coolprop.get(
            "storage_pressure_pa"
        ),
        "storage_phase": coolprop.get(
            "storage_phase"
        ),
        "storage_density_kg_m3": coolprop.get(
            "storage_density_kg_m3"
        ),
        "storage_cp_j_kg_k": coolprop.get(
            "storage_cp_j_kg_k"
        ),
        "enthalpy_window_j_kg": coolprop.get(
            "enthalpy_window_j_kg"
        ),
        "normal_boiling_temperature_k": (
            coolprop.get(
                "normal_boiling_temperature_k"
            )
        ),
        "critical_temperature_k": coolprop.get(
            "critical_temperature_k"
        ),
        "critical_pressure_pa": coolprop.get(
            "critical_pressure_pa"
        ),
        "latent_heat_vaporization_j_kg": (
            coolprop.get(
                "latent_heat_vaporization_j_kg"
            )
        ),
        "nist_compound_page_verified": (
            nist.get(
                "compound_page_verified"
            )
        ),
        "nist_fluid_table_verified": (
            nist.get(
                "fluid_table_verified"
            )
        ),
        "conflict_fields": ";".join(
            row.get("conflict_fields") or []
        ),
    }


def _unresolved_reason(
    row: dict[str, Any],
) -> str:
    reasons: list[str] = []
    identity = row.get("identity") or {}
    reasons.extend(
        identity.get("reasons") or []
    )
    coolprop = row.get("coolprop") or {}
    reasons.extend(
        coolprop.get("errors") or []
    )
    nist = row.get("nist") or {}
    reasons.extend(
        nist.get("errors") or []
    )
    reasons.extend(
        f"cross_source_conflict:{field}"
        for field in (
            row.get("conflict_fields")
            or []
        )
    )
    if not reasons:
        status = row.get(
            "resolution_status"
        )
        if status:
            reasons.append(str(status))
    return ";".join(
        dict.fromkeys(
            str(item)
            for item in reasons
        )
    )


def _audit_row(
    row: dict[str, Any],
) -> dict[str, Any]:
    flattened = _flatten_reference_row(row)
    flattened[
        "unresolved_or_review_reason"
    ] = _unresolved_reason(row)
    pubchem = row.get("pubchem") or {}
    nist = row.get("nist") or {}
    flattened.update(
        {
            "pubchem_source_url": (
                pubchem.get("source_url")
            ),
            "nist_compound_url": (
                nist.get("compound_url")
            ),
            "nist_fluid_url": (
                nist.get("fluid_url")
            ),
            "predictor_reused_as_reference": (
                row.get(
                    "reference_independence"
                )
                or {}
            ).get(
                "predictor_reused_as_reference"
            ),
        }
    )
    return flattened


def _v5b3_handoff(
    rows: list[dict[str, Any]],
    config: dict[str, Any],
) -> dict[str, Any]:
    ready = [
        row
        for row in rows
        if row.get(
            "usable_for_v5b_calibration"
        )
    ]
    holdouts: list[dict[str, Any]] = []
    for row in ready:
        cp = row["coolprop"]
        holdouts.append(
            {
                "id": str(
                    row["candidate_id"]
                ),
                "smiles": str(
                    row["canonical_smiles"]
                    or row["smiles"]
                ),
                "reference_coolprop_name": str(
                    cp["fluid_name"]
                ),
                "storage_temperature_k": float(
                    cp["storage_temperature_k"]
                ),
                "storage_pressure_pa": float(
                    cp["storage_pressure_pa"]
                ),
                "max_exit_temperature_k": float(
                    cp["max_exit_temperature_k"]
                ),
                "challenge_family": str(
                    row.get("family")
                    or "v5e2"
                ),
                "expected_prediction_outcome": (
                    "success"
                ),
                "reference_identity": {
                    "pubchem_cid": (
                        row.get("pubchem")
                        or {}
                    ).get("cid"),
                    "inchi_key": (
                        row.get("pubchem")
                        or {}
                    ).get("inchi_key"),
                    "cas_number": cp.get(
                        "cas_number"
                    ),
                    "coolprop_match_method": (
                        cp.get("match_method")
                    ),
                },
                "v5e2_reference_quality": (
                    row.get(
                        "reference_quality"
                    )
                ),
            }
        )

    minimum = int(
        config["handoff"][
            "minimum_ready_targets"
        ]
    )
    return {
        "version": (
            "v5e-2-v5b3-calibration-input"
        ),
        "executable": (
            len(holdouts) >= minimum
        ),
        "minimum_ready_targets": minimum,
        "ready_target_count": len(
            holdouts
        ),
        "important": (
            "Only targets with exact identity "
            "resolution and an independent "
            "CoolProp reference are included. "
            "FeOS/Joback predictions are never "
            "used as reference data. V5b-3 must "
            "still rerun blind validation and "
            "uncertainty calibration before these "
            "molecules become rankable."
        ),
        "holdouts": holdouts,
    }


def _package_version(
    name: str,
) -> str | None:
    try:
        return version(name)
    except PackageNotFoundError:
        return None


def main() -> int:
    parser = argparse.ArgumentParser(
        description=__doc__
    )
    parser.add_argument(
        "--config",
        type=Path,
        default=(
            HERE / "config-v5e2.yaml"
        ),
    )
    parser.add_argument(
        "--v5e1-reference-resolution",
        type=Path,
        default=None,
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path(
            "reference-v5e2-results"
        ),
    )
    parser.add_argument(
        "--limit",
        type=int,
        default=None,
        help=(
            "Resolve only the first N targets; "
            "intended for smoke tests."
        ),
    )
    args = parser.parse_args()

    if not args.config.is_file():
        parser.error(
            f"missing config: {args.config}"
        )
    config = yaml.safe_load(
        args.config.read_text(
            encoding="utf-8"
        )
    )
    try:
        _validate_config(config)
    except (
        KeyError,
        TypeError,
        ValueError,
    ) as exc:
        parser.error(
            f"invalid V5e-2 config: {exc}"
        )

    configured_input = _resolve_input(
        args.config,
        str(
            config["inputs"][
                "v5e1_reference_resolution"
            ]
        ),
    )
    input_path = (
        args.v5e1_reference_resolution
        if (
            args.v5e1_reference_resolution
            is not None
        )
        else configured_input
    )
    if not input_path.is_file():
        parser.error(
            "missing V5e-1 "
            "reference-resolution input: "
            f"{input_path}"
        )

    if (
        args.limit is not None
        and args.limit <= 0
    ):
        parser.error(
            "--limit must be positive"
        )
    if (
        args.output_dir.exists()
        and any(
            args.output_dir.iterdir()
        )
    ):
        parser.error(
            "output directory must be empty; "
            "use a new directory"
        )
    args.output_dir.mkdir(
        parents=True,
        exist_ok=True,
    )

    try:
        import CoolProp.CoolProp as CP
    except ImportError:
        parser.error(
            "CoolProp is required; install "
            "experiments/entry-evaluator/"
            "requirements.txt"
        )

    input_payload = yaml.safe_load(
        input_path.read_text(
            encoding="utf-8"
        )
    )
    try:
        targets = (
            _targets_from_resolution_template(
                input_payload
            )
        )
    except ValueError as exc:
        parser.error(str(exc))
    if args.limit is not None:
        targets = targets[: args.limit]

    network = config["network"]
    http = CachedHttpClient(
        cache_dir=(
            args.output_dir
            / "source-cache"
        ),
        user_agent=str(
            network["user_agent"]
        ),
        timeout_s=float(
            network["timeout_s"]
        ),
        retries=int(
            network["retries"]
        ),
        minimum_interval_s=float(
            network[
                "minimum_interval_s"
            ]
        ),
    )

    resolved: list[
        dict[str, Any]
    ] = []
    for index, target in enumerate(
        targets,
        start=1,
    ):
        candidate_id = str(
            target["candidate_id"]
        )
        print(
            f"[{index}/{len(targets)}] "
            f"resolving {candidate_id}",
            flush=True,
        )
        row = resolve_target(
            target,
            http=http,
            cp=CP,
            config=config,
        )
        resolved.append(row)

    ready = [
        row
        for row in resolved
        if row[
            "usable_for_v5b_calibration"
        ]
    ]
    review = [
        row
        for row in resolved
        if not row[
            "usable_for_v5b_calibration"
        ]
    ]

    _write_csv(
        args.output_dir
        / "resolved-targets.csv",
        [
            _flatten_reference_row(row)
            for row in ready
        ],
    )
    _write_csv(
        args.output_dir
        / "unresolved-targets.csv",
        [
            {
                **_flatten_reference_row(
                    row
                ),
                "reason": (
                    _unresolved_reason(
                        row
                    )
                ),
            }
            for row in review
        ],
    )
    _write_csv(
        args.output_dir
        / "reference-audit.csv",
        [
            _audit_row(row)
            for row in resolved
        ],
    )

    _write_json(
        args.output_dir
        / "reference-properties.json",
        {
            "version": "v5e-2",
            "targets": resolved,
        },
    )

    provenance = {
        "version": "v5e-2",
        "source_contract": {
            "identity": {
                "source": (
                    "PubChem PUG REST"
                ),
                "method": (
                    "exact InChIKey lookup "
                    "plus formula/MW "
                    "consistency"
                ),
            },
            "primary_thermophysical_reference": {
                "source": "CoolProp",
                "method": (
                    "exact InChIKey or CAS "
                    "match only"
                ),
            },
            "secondary_reference": {
                "source": (
                    "NIST Chemistry "
                    "WebBook SRD 69"
                ),
                "method": (
                    "CAS-based compound "
                    "and fluid-system lookup"
                ),
            },
            "forbidden_reference": (
                "FeOS GC-PC-SAFT/Joback "
                "predictor under calibration"
            ),
        },
        "source_cache": "source-cache/",
        "target_sources": [
            {
                "candidate_id": (
                    row["candidate_id"]
                ),
                "pubchem_url": (
                    row.get("pubchem")
                    or {}
                ).get("source_url"),
                "nist_compound_url": (
                    row.get("nist")
                    or {}
                ).get("compound_url"),
                "nist_fluid_url": (
                    row.get("nist")
                    or {}
                ).get("fluid_url"),
                "coolprop_fluid_name": (
                    row.get("coolprop")
                    or {}
                ).get("fluid_name"),
                "resolution_status": (
                    row[
                        "resolution_status"
                    ]
                ),
                "reference_quality": (
                    row[
                        "reference_quality"
                    ]
                ),
            }
            for row in resolved
        ],
    }
    _write_json(
        args.output_dir
        / "provenance.json",
        provenance,
    )

    handoff = _v5b3_handoff(
        resolved,
        config,
    )
    for filename in (
        "calibration-expansion.yaml",
        "v5b3-calibration-input.yaml",
    ):
        with (
            args.output_dir
            / filename
        ).open(
            "w",
            encoding="utf-8",
        ) as handle:
            yaml.safe_dump(
                handoff,
                handle,
                sort_keys=False,
                allow_unicode=True,
            )

    status_counts: dict[
        str,
        int,
    ] = {}
    quality_counts: dict[
        str,
        int,
    ] = {}
    for row in resolved:
        status = str(
            row["resolution_status"]
        )
        quality = str(
            row["reference_quality"]
        )
        status_counts[status] = (
            status_counts.get(
                status,
                0,
            )
            + 1
        )
        quality_counts[quality] = (
            quality_counts.get(
                quality,
                0,
            )
            + 1
        )

    summary = {
        "study_version": "v5e-2",
        "study_complete": True,
        "target_count_attempted": len(
            resolved
        ),
        "ready_for_v5b3_count": len(
            ready
        ),
        "review_or_unresolved_count": (
            len(review)
        ),
        "resolution_status_counts": dict(
            sorted(
                status_counts.items()
            )
        ),
        "reference_quality_counts": dict(
            sorted(
                quality_counts.items()
            )
        ),
        "v5b3_handoff_executable": bool(
            handoff["executable"]
        ),
        "minimum_ready_targets_for_handoff": int(
            config["handoff"][
                "minimum_ready_targets"
            ]
        ),
        "manual_work_required": (
            "Only unresolved/review rows "
            "require human follow-up; ready "
            "rows are emitted directly in "
            "v5b3-calibration-input.yaml."
        ),
        "interpretation": (
            "V5e-2 resolves identity and "
            "independent reference evidence. "
            "A ready target is not yet rankable: "
            "V5b-3 must rerun blind prediction, "
            "error calibration, applicability-"
            "domain analysis and uncertainty "
            "validation using the expanded "
            "reference set."
        ),
    }
    _write_json(
        args.output_dir
        / "summary.json",
        summary,
    )

    try:
        git_commit = (
            subprocess.check_output(
                [
                    "git",
                    "rev-parse",
                    "HEAD",
                ],
                cwd=HERE,
                text=True,
            ).strip()
        )
        git_dirty = bool(
            subprocess.check_output(
                [
                    "git",
                    "status",
                    "--porcelain",
                ],
                cwd=HERE,
                text=True,
            ).strip()
        )
    except (
        OSError,
        subprocess.CalledProcessError,
    ):
        git_commit = None
        git_dirty = None

    manifest = {
        "study_version": "v5e-2",
        "git_commit": git_commit,
        "git_dirty": git_dirty,
        "python": (
            platform.python_version()
        ),
        "packages": {
            "CoolProp": (
                _package_version(
                    "CoolProp"
                )
            ),
            "PyYAML": (
                _package_version(
                    "PyYAML"
                )
            ),
            "rdkit": (
                _package_version(
                    "rdkit"
                )
            ),
        },
        "inputs": {
            "config": str(args.config),
            "config_sha256": (
                _sha256(args.config)
            ),
            "v5e1_reference_resolution": (
                str(input_path)
            ),
            "v5e1_reference_resolution_sha256": (
                _sha256(input_path)
            ),
        },
        "network_contract": {
            "user_agent": str(
                network["user_agent"]
            ),
            "timeout_s": float(
                network["timeout_s"]
            ),
            "retries": int(
                network["retries"]
            ),
            "minimum_interval_s": float(
                network[
                    "minimum_interval_s"
                ]
            ),
            "responses_cached": True,
        },
        "calibration_gate": {
            "identity_source": "PubChem",
            "accepted_property_backend": (
                "CoolProp"
            ),
            "corroborating_source": (
                "NIST Chemistry WebBook"
            ),
            "cross_source_conflict_threshold": float(
                config["quality"][
                    "maximum_cross_source_relative_difference"
                ]
            ),
            "predictor_reused_as_reference": (
                False
            ),
        },
    }
    _write_json(
        args.output_dir
        / "manifest.json",
        manifest,
    )

    print(
        json.dumps(
            summary,
            indent=2,
            sort_keys=True,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
