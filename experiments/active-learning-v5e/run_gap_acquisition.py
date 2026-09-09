#!/usr/bin/env python3
"""V5e-3: acquire a small empirical-reference batch to close V5b-3 property gaps."""
from __future__ import annotations

import argparse
from collections import Counter
import csv
import hashlib
from importlib.metadata import PackageNotFoundError, version
import json
import platform
from pathlib import Path
import subprocess
import sys
from typing import Any

import yaml

HERE = Path(__file__).resolve().parent
EXPERIMENTS = HERE.parent
V5B = EXPERIMENTS / "property-predictor-v5b"
if str(V5B) not in sys.path:
    sys.path.insert(0, str(V5B))
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

from gap_acquisition import (  # noqa: E402
    candidate_pool,
    canonical_identity,
    empirical_availability,
    gap_projection,
    prepare_preprobe_features,
    property_deficits,
    select_gap_targets,
    select_probe_shortlist,
    strict_property_anchors,
)
from reference_resolver import CachedHttpClient, resolve_pubchem_identity  # noqa: E402
from v5b3_audit import LocalReferenceBackend, allowed, structural_identity  # noqa: E402


def _sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def _json_digest(value: Any) -> str:
    return hashlib.sha256(
        json.dumps(value, sort_keys=True, separators=(",", ":"), allow_nan=False).encode()
    ).hexdigest()


def _write_json(path: Path, value: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temp = path.with_suffix(path.suffix + ".tmp")
    temp.write_text(
        json.dumps(value, indent=2, sort_keys=True, allow_nan=False) + "\n",
        encoding="utf-8",
    )
    temp.replace(path)


def _load_json(path: Path) -> Any:
    def reject(text: str):
        raise ValueError("nonfinite_json_constant:" + text)
    return json.loads(path.read_text(encoding="utf-8"), parse_constant=reject)


def _write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    fields = sorted({key for row in rows for key in row}) or ["status"]
    temp = path.with_suffix(path.suffix + ".tmp")
    with temp.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader()
        for row in rows:
            writer.writerow(
                {
                    key: json.dumps(value, sort_keys=True, allow_nan=False)
                    if isinstance(value, (dict, list))
                    else value
                    for key, value in row.items()
                }
            )
    temp.replace(path)


def _resolve(config_path: Path, value: str) -> Path:
    path = Path(value)
    return path if path.is_absolute() else (config_path.parent / path).resolve()


def _validate_config(config: dict[str, Any]) -> None:
    if config.get("version") != "v5e-3":
        raise ValueError("config.version must be v5e-3")
    priorities = list(config["priority_properties"])
    if not priorities or len(priorities) != len(set(priorities)):
        raise ValueError("priority_properties must be a non-empty unique list")
    for prop in priorities:
        if prop not in config["property_weights"]:
            raise ValueError(f"missing property weight: {prop}")
        if float(config["property_weights"][prop]) <= 0.0:
            raise ValueError(f"property weight must be positive: {prop}")
    targets = config["gap_targets"]
    if int(targets["calibration_groups"]) < 9:
        raise ValueError("gap target must preserve at least nine calibration groups")
    if int(targets["evaluation_groups"]) < 3:
        raise ValueError("gap target must preserve at least three evaluation groups")
    weights = {k: float(v) for k, v in config["preprobe_weights"].items()}
    if any(v < 0.0 for v in weights.values()) or abs(sum(weights.values()) - 1.0) > 1e-9:
        raise ValueError("preprobe_weights must be nonnegative and sum to 1")
    weights = {k: float(v) for k, v in config["acquisition"]["weights"].items()}
    if any(v < 0.0 for v in weights.values()) or abs(sum(weights.values()) - 1.0) > 1e-9:
        raise ValueError("acquisition.weights must be nonnegative and sum to 1")
    if int(config["probe"]["candidate_count"]) < int(config["acquisition"]["target_count"]):
        raise ValueError("probe.candidate_count must be >= acquisition.target_count")
    for section in ("probe", "acquisition"):
        fraction = float(config[section]["maximum_family_fraction"])
        if not 0.0 < fraction <= 1.0:
            raise ValueError(f"{section}.maximum_family_fraction must be in (0, 1]")
    if float(config["identity"]["mw_relative_tolerance"]) <= 0.0:
        raise ValueError("identity.mw_relative_tolerance must be positive")
    network = config["network"]
    if float(network["timeout_s"]) <= 0.0 or int(network["retries"]) < 0:
        raise ValueError("invalid network timeout/retries")


def _excluded_smiles(acquisition_set: dict[str, Any]) -> set[str]:
    result: set[str] = set()
    for row in acquisition_set.get("targets", []):
        smiles = row.get("canonical_smiles") or row.get("smiles")
        if smiles:
            result.add(canonical_identity(str(smiles))["canonical_smiles"])
    return result


def _feature_csv(features) -> list[dict[str, Any]]:
    return [row.to_dict() for row in features]


def _resolved_csv(rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    return [
        {
            "probe_rank": row.get("probe_rank"),
            "candidate_id": row.get("candidate_id"),
            "smiles": row.get("smiles"),
            "family": row.get("family"),
            "split": row.get("split"),
            "property_priority_score": row.get("property_priority_score"),
            "property_percentile": row.get("property_percentile"),
            "identity_verified": row.get("identity_verified"),
            "pubchem_cid": (row.get("pubchem") or {}).get("cid"),
            "inchi_key": row.get("inchi_key"),
            "selected_cas_number": row.get("selected_cas_number"),
            "available_properties": ";".join(row.get("available_properties") or []),
            "property_methods": row.get("property_methods") or {},
            "identity_reasons": ";".join(row.get("identity_reasons") or []),
            "preprobe_score": row.get("preprobe_score"),
        }
        for row in rows
    ]


def _target_csv(rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    return [
        {
            "acquisition_rank": row["acquisition_rank"],
            "candidate_id": row["candidate_id"],
            "smiles": row["smiles"],
            "family": row["family"],
            "split": row["split"],
            "property_priority_score": row["property_priority_score"],
            "property_percentile": row["property_percentile"],
            "inchi_key": row["inchi_key"],
            "selected_cas_number": row["selected_cas_number"],
            "available_properties": ";".join(row.get("available_properties") or []),
            "gap_properties_covered_at_selection": ";".join(row.get("gap_properties_covered_at_selection") or []),
            "gap_score_at_selection": row["gap_score_at_selection"],
            "property_structural_gain_at_selection": row["property_structural_gain_at_selection"],
            "minimum_distance_to_earlier_selected": row["minimum_distance_to_earlier_selected"],
        }
        for row in rows
    ]


def _checkpoint(path: Path, signature: str, data: Any | None = None) -> Any | None:
    if data is not None:
        _write_json(path, {"signature": signature, "data_sha256": _json_digest(data), "data": data})
        return data
    if not path.is_file():
        return None
    saved = _load_json(path)
    if saved.get("signature") != signature:
        raise ValueError(f"checkpoint signature mismatch: {path}")
    if saved.get("data_sha256") != _json_digest(saved.get("data")):
        raise ValueError(f"checkpoint corruption: {path}")
    return saved["data"]


def _package_version(name: str) -> str | None:
    try:
        return version(name)
    except PackageNotFoundError:
        return None


def _run_contract(v5b3_summary: dict[str, Any]) -> None:
    if not bool(v5b3_summary.get("study_complete")):
        raise ValueError("V5b-3 study is incomplete")
    if not bool(v5b3_summary.get("benchmark_has_comparable_results")):
        raise ValueError("V5b-3 has no comparable strict references")
    if int(v5b3_summary.get("strict_failed_prediction_count") or 0) != 0:
        raise ValueError("V5b-3 has strict prediction failures; fix numerical failures before acquisition")
    if bool(v5b3_summary.get("predictor_parameters_changed")):
        raise ValueError("V5b-3 predictor parameters unexpectedly changed")
    if int(v5b3_summary.get("rankable_promotions") or 0) != 0:
        raise ValueError("V5b-3 unexpectedly promoted candidates")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, default=HERE / "config-v5e3.yaml")
    parser.add_argument("--v5b3-summary", type=Path, default=None)
    parser.add_argument("--v5b3-audited-reference", type=Path, default=None)
    parser.add_argument("--v5b3-config", type=Path, default=None)
    parser.add_argument("--v5d-full-prescreen", type=Path, default=None)
    parser.add_argument("--v5e1-acquisition-set", type=Path, default=None)
    parser.add_argument("--output-dir", type=Path, default=Path("active-learning-v5e3-results"))
    parser.add_argument("--resume", action="store_true")
    parser.add_argument("--limit-probes", type=int, default=None,
                        help="Smoke-test only: resolve at most N probes before final selection.")
    args = parser.parse_args()

    if not args.config.is_file():
        parser.error(f"missing config: {args.config}")
    config = yaml.safe_load(args.config.read_text(encoding="utf-8"))
    try:
        _validate_config(config)
    except (KeyError, TypeError, ValueError) as exc:
        parser.error(f"invalid V5e-3 config: {exc}")

    configured = config["inputs"]
    paths = {
        "v5b3_summary": args.v5b3_summary or _resolve(args.config, configured["v5b3_summary"]),
        "v5b3_audited_reference": args.v5b3_audited_reference or _resolve(args.config, configured["v5b3_audited_reference"]),
        "v5b3_config": args.v5b3_config or _resolve(args.config, configured["v5b3_config"]),
        "v5d_full_prescreen": args.v5d_full_prescreen or _resolve(args.config, configured["v5d_full_prescreen"]),
        "v5e1_acquisition_set": args.v5e1_acquisition_set or _resolve(args.config, configured["v5e1_acquisition_set"]),
    }
    missing = [str(path) for path in paths.values() if not path.is_file()]
    if missing:
        parser.error("missing required inputs: " + ", ".join(missing))
    if args.limit_probes is not None and args.limit_probes <= 0:
        parser.error("--limit-probes must be positive")
    if args.output_dir.exists() and any(args.output_dir.iterdir()) and not args.resume:
        parser.error("output directory is nonempty; use --resume or a new directory")
    args.output_dir.mkdir(parents=True, exist_ok=True)

    v5b3_summary = _load_json(paths["v5b3_summary"])
    audited_reference = _load_json(paths["v5b3_audited_reference"])
    v5b3_config = yaml.safe_load(paths["v5b3_config"].read_text(encoding="utf-8"))
    full_prescreen = _load_json(paths["v5d_full_prescreen"])
    v5e1_acquisition = _load_json(paths["v5e1_acquisition_set"])
    try:
        _run_contract(v5b3_summary)
    except ValueError as exc:
        parser.error(str(exc))
    if v5b3_config.get("version") != "v5b-3":
        parser.error("expected V5b-3 benchmark configuration")
    statistics = v5b3_config["statistics"]

    input_hashes = {name: _sha256(path) for name, path in paths.items()}
    code_hashes = {
        "run_gap_acquisition.py": _sha256(Path(__file__)),
        "gap_acquisition.py": _sha256(HERE / "gap_acquisition.py"),
    }
    signature = _json_digest({
        "study": "v5e-3",
        "config_sha256": _sha256(args.config),
        "inputs": input_hashes,
        "code": code_hashes,
    })
    progress_path = args.output_dir / "progress.json"
    if args.resume and progress_path.is_file():
        old = _load_json(progress_path)
        if old.get("signature") != signature:
            parser.error("existing V5e-3 output was produced by a different input/config/code signature")

    deficits = property_deficits(v5b3_summary, config)
    anchors = strict_property_anchors(audited_reference, config["priority_properties"])
    excluded = _excluded_smiles(v5e1_acquisition)
    pool = candidate_pool(full_prescreen, excluded_smiles=excluded, config=config)
    if not pool:
        parser.error("no eligible V5d exploratory candidates remain after exclusions")
    features = prepare_preprobe_features(
        pool, deficits=deficits, anchors=anchors, statistics=statistics, config=config)
    shortlist = select_probe_shortlist(features, config)
    if args.limit_probes is not None:
        shortlist = shortlist[: args.limit_probes]
    if not shortlist:
        parser.error("property-gap preprobe produced no candidates")

    _write_csv(args.output_dir / "preprobe-pool.csv", _feature_csv(features))
    _write_csv(args.output_dir / "probe-shortlist.csv", shortlist)
    _write_json(progress_path, {"signature": signature, "phase": "probing", "completed": 0,
                                "total": len(shortlist), "property_deficits": deficits})

    network = config["network"]
    http = CachedHttpClient(
        cache_dir=args.output_dir / "source-cache",
        user_agent=str(network["user_agent"]),
        timeout_s=float(network["timeout_s"]),
        retries=int(network["retries"]),
        minimum_interval_s=float(network["minimum_interval_s"]),
    )
    try:
        backend = LocalReferenceBackend()
    except Exception as exc:
        parser.error(f"could not initialize local empirical reference backend: {exc}")

    resolved: list[dict[str, Any]] = []
    for index, feature in enumerate(shortlist, start=1):
        candidate_id = str(feature["candidate_id"])
        checkpoint_path = args.output_dir / "checkpoints" / f"{candidate_id}.json"
        cached = False
        row = _checkpoint(checkpoint_path, signature)
        if row is None:
            identity = canonical_identity(str(feature["smiles"]))
            pubchem_obj = resolve_pubchem_identity(identity["inchi_key"], http)
            pubchem = pubchem_obj.to_dict() if pubchem_obj is not None else None
            availability = empirical_availability(
                candidate_identity=identity,
                pubchem=pubchem,
                backend=backend,
                allowed_fn=allowed,
                structural_identity_fn=structural_identity,
                config=config,
            )
            row = {
                **feature,
                **identity,
                "pubchem": pubchem,
                **availability,
            }
            _checkpoint(checkpoint_path, signature, row)
        else:
            cached = True
        resolved.append(row)
        available = ",".join(row.get("available_properties") or []) or "-"
        print(
            f"[probe {index}/{len(shortlist)}] {candidate_id} "
            f"identity={row.get('identity_verified', False)} cas={row.get('selected_cas_number') or '-'} "
            f"split={row.get('split')} properties={available} cached={cached}",
            flush=True,
        )
        _write_json(progress_path, {"signature": signature, "phase": "probing", "completed": index,
                                    "total": len(shortlist), "property_deficits": deficits})

    _write_csv(args.output_dir / "resolved-probe.csv", _resolved_csv(resolved))
    _write_json(args.output_dir / "resolved-probe.json", {"version": "v5e-3", "probes": resolved})

    selected, remaining = select_gap_targets(
        resolved, deficits=deficits, anchors=anchors, config=config)
    projection = gap_projection(deficits, remaining)
    _write_csv(args.output_dir / "calibration-targets.csv", _target_csv(selected))
    _write_json(args.output_dir / "property-gap-projection.json", projection)
    _write_json(args.output_dir / "acquisition-set.json", {"version": "v5e-3", "targets": selected})
    _write_json(
        args.output_dir / "v5b3-expansion-targets.json",
        {
            "version": "v5e-3-v5b3-expansion-targets",
            "reference_values_included": False,
            "important": (
                "Targets contain exact identity and local empirical-source availability only. "
                "Reference values must be reconstructed and audited by the next V5b benchmark; "
                "availability does not make a candidate rankable."
            ),
            "targets": [
                {
                    "id": row["candidate_id"],
                    "smiles": row["canonical_smiles"],
                    "inchi_key": row["inchi_key"],
                    "cas_number": row["selected_cas_number"],
                    "family": row["family"],
                    "split": row["split"],
                    "expected_strict_properties": row["available_properties"],
                    "available_reference_methods": row["property_methods"],
                    "acquisition_rank": row["acquisition_rank"],
                }
                for row in selected
            ],
        },
    )

    property_selected_counts = Counter(
        prop
        for row in selected
        for prop in row.get("gap_properties_covered_at_selection") or []
    )
    split_counts = Counter(str(row["split"]) for row in selected)
    identity_count = sum(bool(row.get("identity_verified")) for row in resolved)
    local_cas_count = sum(bool(row.get("selected_cas_number")) for row in resolved)
    smoke_test = args.limit_probes is not None
    summary = {
        "study_version": "v5e-3",
        "study_complete": not smoke_test,
        "run_status": "smoke_test" if smoke_test else "complete",
        "candidate_pool_count": len(features),
        "probe_count": len(shortlist),
        "probe_identity_verified_count": identity_count,
        "probe_local_exact_cas_count": local_cas_count,
        "selected_target_count": len(selected),
        "selected_split_counts": dict(sorted(split_counts.items())),
        "selected_property_gap_contribution_counts": dict(sorted(property_selected_counts.items())),
        "property_gap_projection": projection,
        "all_configured_gaps_projected_closed": all(
            values["remaining_gap_after_selection"][split] == 0
            for values in projection.values()
            for split in ("calibration", "evaluation")
        ),
        "rankable_promotions": 0,
        "reference_values_fabricated": False,
        "interpretation": (
            "V5e-3 is a property-gap acquisition plan. It prioritizes exact-identity molecules "
            "with local empirical source support that are predicted to close V5b-3 calibration/"
            "evaluation deficits. Projected gap closure is conditional on the next strict audit."
        ),
    }
    _write_json(args.output_dir / "summary.json", summary)

    try:
        git_commit = subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=HERE, text=True).strip()
        git_dirty = bool(subprocess.check_output(["git", "status", "--porcelain"], cwd=HERE, text=True).strip())
    except (OSError, subprocess.CalledProcessError):
        git_commit = None
        git_dirty = None
    manifest = {
        "study_version": "v5e-3",
        "signature": signature,
        "git_commit": git_commit,
        "git_dirty": git_dirty,
        "python": platform.python_version(),
        "packages": {name: _package_version(name) for name in ("rdkit", "chemicals", "thermo", "PyYAML")},
        "inputs": {name: {"path": str(paths[name]), "sha256": digest} for name, digest in input_hashes.items()},
        "config_sha256": _sha256(args.config),
        "code_sha256": code_hashes,
        "network": {
            "identity_source": "PubChem exact InChIKey lookup",
            "timeout_s": float(network["timeout_s"]),
            "retries": int(network["retries"]),
            "minimum_interval_s": float(network["minimum_interval_s"]),
            "responses_cached": True,
        },
        "reference_availability_contract": {
            "local_identity": "chemicals exact CAS + full InChIKey + structure match",
            "source_policy": "V5b-3 audited allowlist",
            "reference_values_used_for_selection": False,
            "next_audit_required": True,
        },
    }
    _write_json(args.output_dir / "manifest.json", manifest)
    _write_json(progress_path, {"signature": signature, "phase": summary["run_status"],
                                "selected_target_count": len(selected)})
    print(json.dumps(summary, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
