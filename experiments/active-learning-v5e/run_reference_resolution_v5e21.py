#!/usr/bin/env python3
"""V5e-2.1: resolve multi-source empirical property packets for V5b-3."""
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

from empirical_reference import (
    build_v5b3_property_handoff,
    resolve_empirical_target,
)


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


def _package_version(name: str) -> str | None:
    try:
        return version(name)
    except PackageNotFoundError:
        return None


def _validate_config(config: dict[str, Any]) -> None:
    if str(config.get("version")) != "v5e-2.1":
        raise ValueError("config.version must be v5e-2.1")

    temperatures = [
        float(value)
        for value in config["properties"]["temperatures_k"]
    ]
    if not temperatures or any(value <= 0.0 for value in temperatures):
        raise ValueError("properties.temperatures_k must contain positive values")
    if temperatures != sorted(set(temperatures)):
        raise ValueError(
            "properties.temperatures_k must be unique and sorted ascending"
        )

    policy = config["source_policy"]
    required_policy_keys = (
        "critical_temperature_methods",
        "critical_pressure_methods",
        "normal_boiling_temperature_methods",
        "vapor_pressure_methods",
        "liquid_volume_methods",
        "liquid_heat_capacity_methods",
        "enthalpy_vaporization_methods",
    )
    for key in required_policy_keys:
        values = policy.get(key)
        if not isinstance(values, list) or not values:
            raise ValueError(f"source_policy.{key} must be a non-empty list")

    gate = config["calibration_gate"]
    required = list(gate["required_properties"])
    dynamic = list(gate["dynamic_properties"])
    if not required:
        raise ValueError("calibration_gate.required_properties cannot be empty")
    if not dynamic:
        raise ValueError("calibration_gate.dynamic_properties cannot be empty")
    if int(gate["minimum_dynamic_properties"]) < 1:
        raise ValueError(
            "calibration_gate.minimum_dynamic_properties must be positive"
        )
    if int(gate["minimum_total_properties"]) < len(required):
        raise ValueError(
            "calibration_gate.minimum_total_properties must cover required properties"
        )
    if int(gate["minimum_source_families"]) < 1:
        raise ValueError(
            "calibration_gate.minimum_source_families must be positive"
        )
    forbidden = gate["forbidden_method_selectors"]
    if not isinstance(forbidden, list) or not forbidden:
        raise ValueError(
            "calibration_gate.forbidden_method_selectors must be non-empty"
        )
    high = gate["high_quality"]
    if int(high["minimum_dynamic_properties"]) < int(
        gate["minimum_dynamic_properties"]
    ):
        raise ValueError(
            "high_quality minimum_dynamic_properties cannot be weaker than gate"
        )
    if int(high["minimum_total_properties"]) < int(
        gate["minimum_total_properties"]
    ):
        raise ValueError(
            "high_quality minimum_total_properties cannot be weaker than gate"
        )
    if int(high["minimum_source_families"]) < int(
        gate["minimum_source_families"]
    ):
        raise ValueError(
            "high_quality minimum_source_families cannot be weaker than gate"
        )

    if int(config["handoff"]["minimum_ready_targets"]) < 1:
        raise ValueError("handoff.minimum_ready_targets must be positive")


def _load_targets(payload: dict[str, Any]) -> list[dict[str, Any]]:
    targets = payload.get("targets")
    if not isinstance(targets, list) or not targets:
        raise ValueError("V5e-2 reference-properties input contains no targets")

    output: list[dict[str, Any]] = []
    seen: set[str] = set()
    for index, row in enumerate(targets, start=1):
        if not isinstance(row, dict):
            raise ValueError(f"target {index} is not a mapping")
        candidate_id = row.get("candidate_id")
        if not candidate_id:
            raise ValueError(f"target {index} has no candidate_id")
        candidate_id = str(candidate_id)
        if candidate_id in seen:
            raise ValueError(f"duplicate candidate_id: {candidate_id}")
        seen.add(candidate_id)
        output.append(dict(row))

    output.sort(
        key=lambda row: (
            int(row.get("acquisition_rank") or 10**9),
            str(row["candidate_id"]),
        )
    )
    return output


def _compact_row(row: dict[str, Any]) -> dict[str, Any]:
    classification = row.get("classification") or {}
    packet = row.get("packet") or {}
    constants = packet.get("constants") or {}
    series = packet.get("series") or {}

    def method_of(name: str) -> str | None:
        record = constants.get(name) or series.get(name)
        return record.get("method") if record else None

    def point_count(name: str) -> int:
        record = series.get(name) or {}
        return len(record.get("points") or [])

    return {
        "acquisition_rank": row.get("acquisition_rank"),
        "candidate_id": row.get("candidate_id"),
        "smiles": row.get("smiles"),
        "family": row.get("family"),
        "preferred_name": row.get("preferred_name"),
        "inchi_key": row.get("inchi_key"),
        "selected_cas_number": row.get("selected_cas_number"),
        "resolution_status": row.get("resolution_status"),
        "reference_quality": row.get("reference_quality"),
        "ready_for_v5b3_property_calibration": row.get(
            "ready_for_v5b3_property_calibration"
        ),
        "entry_reference_ready": row.get("entry_reference_ready"),
        "dynamic_property_count": classification.get("dynamic_property_count"),
        "total_property_count": classification.get("total_property_count"),
        "source_family_count": classification.get("source_family_count"),
        "source_families": ";".join(
            classification.get("source_families") or []
        ),
        "present_properties": ";".join(
            classification.get("present_properties") or []
        ),
        "reasons": ";".join(classification.get("reasons") or []),
        "tc_method": method_of("critical_temperature_k"),
        "pc_method": method_of("critical_pressure_pa"),
        "tb_method": method_of("normal_boiling_temperature_k"),
        "vapor_pressure_method": method_of("vapor_pressure_pa"),
        "density_method": method_of("liquid_density_kg_m3"),
        "cp_method": method_of("liquid_cp_j_kg_k"),
        "hvap_method": method_of("latent_heat_vaporization_j_kg"),
        "vapor_pressure_grid_points": point_count("vapor_pressure_pa"),
        "density_grid_points": point_count("liquid_density_kg_m3"),
        "cp_grid_points": point_count("liquid_cp_j_kg_k"),
        "hvap_grid_points": point_count("latent_heat_vaporization_j_kg"),
    }


def _source_audit_rows(row: dict[str, Any]) -> list[dict[str, Any]]:
    packet = row.get("packet") or {}
    result: list[dict[str, Any]] = []
    for audit in packet.get("source_audit") or []:
        result.append({
            "candidate_id": row.get("candidate_id"),
            "cas_number": row.get("selected_cas_number"),
            **audit,
        })
    return result


def _property_coverage(rows: list[dict[str, Any]]) -> dict[str, int]:
    counts: dict[str, int] = {}
    for row in rows:
        for name in (
            (row.get("classification") or {}).get("present_properties") or []
        ):
            counts[str(name)] = counts.get(str(name), 0) + 1
    return dict(sorted(counts.items()))


def _source_family_counts(rows: list[dict[str, Any]]) -> dict[str, int]:
    counts: dict[str, int] = {}
    for row in rows:
        families = (
            (row.get("classification") or {}).get("source_families") or []
        )
        for family in families:
            counts[str(family)] = counts.get(str(family), 0) + 1
    return dict(sorted(counts.items()))


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--config",
        type=Path,
        default=HERE / "config-v5e21.yaml",
    )
    parser.add_argument(
        "--v5e2-reference-properties",
        type=Path,
        default=None,
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("reference-v5e21-results"),
    )
    parser.add_argument(
        "--limit",
        type=int,
        default=None,
        help="Resolve only the first N targets; intended for smoke tests.",
    )
    args = parser.parse_args()

    if not args.config.is_file():
        parser.error(f"missing config: {args.config}")
    config = yaml.safe_load(args.config.read_text(encoding="utf-8"))
    try:
        _validate_config(config)
    except (KeyError, TypeError, ValueError) as exc:
        parser.error(f"invalid V5e-2.1 config: {exc}")

    input_path = (
        args.v5e2_reference_properties
        if args.v5e2_reference_properties is not None
        else _resolve_input(
            args.config,
            str(config["inputs"]["v5e2_reference_properties"]),
        )
    )
    if not input_path.is_file():
        parser.error(f"missing V5e-2 reference-properties input: {input_path}")

    if args.limit is not None and args.limit <= 0:
        parser.error("--limit must be positive")
    if args.output_dir.exists() and any(args.output_dir.iterdir()):
        parser.error("output directory must be empty; use a new directory")
    args.output_dir.mkdir(parents=True, exist_ok=True)

    try:
        import chemicals  # noqa: F401
        import thermo  # noqa: F401
    except ImportError:
        parser.error(
            "chemicals and thermo are required; install "
            "experiments/entry-evaluator/requirements-v5a.txt"
        )

    payload = json.loads(input_path.read_text(encoding="utf-8"))
    try:
        targets = _load_targets(payload)
    except ValueError as exc:
        parser.error(str(exc))
    if args.limit is not None:
        targets = targets[: args.limit]

    rows: list[dict[str, Any]] = []
    for index, target in enumerate(targets, start=1):
        candidate_id = str(target["candidate_id"])
        print(
            f"[{index}/{len(targets)}] empirical resolve {candidate_id}",
            flush=True,
        )
        row = resolve_empirical_target(target, config)
        rows.append(row)
        classification = row.get("classification") or {}
        print(
            "  "
            f"{row['resolution_status']} "
            f"quality={row['reference_quality']} "
            f"properties={classification.get('total_property_count', 0)} "
            f"dynamic={classification.get('dynamic_property_count', 0)} "
            f"sources={classification.get('source_family_count', 0)}",
            flush=True,
        )

    ready = [
        row
        for row in rows
        if row.get("ready_for_v5b3_property_calibration")
    ]
    review = [
        row
        for row in rows
        if not row.get("ready_for_v5b3_property_calibration")
    ]

    _write_csv(
        args.output_dir / "calibration-ready.csv",
        [_compact_row(row) for row in ready],
    )
    _write_csv(
        args.output_dir / "partial-or-unresolved.csv",
        [_compact_row(row) for row in review],
    )
    _write_csv(
        args.output_dir / "reference-audit.csv",
        [_compact_row(row) for row in rows],
    )
    _write_csv(
        args.output_dir / "source-method-audit.csv",
        [
            audit
            for row in rows
            for audit in _source_audit_rows(row)
        ],
    )

    _write_json(
        args.output_dir / "empirical-reference-properties.json",
        {
            "version": "v5e-2.1",
            "targets": rows,
        },
    )

    handoff = build_v5b3_property_handoff(rows, config)
    with (
        args.output_dir / "v5b3-property-calibration-input.yaml"
    ).open("w", encoding="utf-8") as handle:
        yaml.safe_dump(
            handoff,
            handle,
            sort_keys=False,
            allow_unicode=True,
        )

    quality_counts: dict[str, int] = {}
    status_counts: dict[str, int] = {}
    for row in rows:
        quality = str(row["reference_quality"])
        status = str(row["resolution_status"])
        quality_counts[quality] = quality_counts.get(quality, 0) + 1
        status_counts[status] = status_counts.get(status, 0) + 1

    forbidden_selected = sum(
        len(
            (row.get("classification") or {}).get(
                "forbidden_selected_methods"
            ) or []
        )
        for row in rows
    )
    entry_ready = sum(bool(row.get("entry_reference_ready")) for row in rows)

    summary = {
        "study_version": "v5e-2.1",
        "study_complete": True,
        "target_count_attempted": len(rows),
        "identity_verified_count": sum(
            row["resolution_status"] != "unresolved_identity"
            for row in rows
        ),
        "ready_for_v5b3_property_calibration_count": len(ready),
        "review_or_unresolved_count": len(review),
        "entry_reference_ready_count": entry_ready,
        "reference_quality_counts": dict(sorted(quality_counts.items())),
        "resolution_status_counts": dict(sorted(status_counts.items())),
        "property_coverage_counts": _property_coverage(rows),
        "source_family_candidate_counts": _source_family_counts(rows),
        "forbidden_selected_method_count": forbidden_selected,
        "v5b3_property_handoff_executable": bool(handoff["executable"]),
        "minimum_ready_targets_for_handoff": int(
            config["handoff"]["minimum_ready_targets"]
        ),
        "interpretation": (
            "V5e-2.1 broadens reference acquisition beyond CoolProp using "
            "explicitly whitelisted empirical/review/tabulated/data-fitted "
            "sources exposed by chemicals/thermo. Ready anchors support "
            "property-model calibration and applicability-domain expansion; "
            "they do not become entry-level ground truth unless an independent "
            "entry-capable reference backend also exists."
        ),
    }
    _write_json(args.output_dir / "summary.json", summary)

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
        "study_version": "v5e-2.1",
        "git_commit": git_commit,
        "git_dirty": git_dirty,
        "python": platform.python_version(),
        "packages": {
            name: _package_version(name)
            for name in ("chemicals", "thermo", "fluids", "scipy", "PyYAML")
        },
        "inputs": {
            "config": str(args.config),
            "config_sha256": _sha256(args.config),
            "v5e2_reference_properties": str(input_path),
            "v5e2_reference_properties_sha256": _sha256(input_path),
        },
        "reference_contract": {
            "identity_source": "V5e-2 exact PubChem identity resolution",
            "source_selection": "explicit allowlist only",
            "source_inventory_preserved": True,
            "forbidden_method_selectors": list(
                config["calibration_gate"]["forbidden_method_selectors"]
            ),
            "predictor_under_calibration": "FeOS GC-PC-SAFT + Joback",
            "predictor_reused_as_reference": False,
            "entry_reference_semantics": (
                "property packet is not automatically entry-level ground truth"
            ),
        },
    }
    _write_json(args.output_dir / "manifest.json", manifest)

    print(json.dumps(summary, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
