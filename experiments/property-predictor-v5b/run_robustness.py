#!/usr/bin/env python3
"""V5b-2.1: repeated-split and leave-one-family-out uncertainty robustness."""
from __future__ import annotations

import argparse
import json
import math
import statistics
from pathlib import Path
from typing import Any

import yaml

from applicability import (
    assess_domain,
    calibrate_metric,
    deterministic_split,
    predict_relative_uncertainty,
)
from failure_taxonomy import (
    SUCCESS,
    classify_prediction,
    expected_prediction_outcome,
)
from run_holdout import write_summary

HERE = Path(__file__).resolve().parent
METRIC_KEYS = (
    "entry_error",
    "storage_density_error",
    "delta_h_error",
    "cp_error",
    "saturation_error",
)


def _load_records(summary: dict[str, Any]) -> list[dict[str, Any]]:
    records: list[dict[str, Any]] = []
    for row in summary.get("candidates", []):
        candidate_id = row.get("candidate_id")
        smiles = row.get("smiles")
        if not candidate_id or not smiles:
            continue
        records.append(
            {
                "candidate_id": str(candidate_id),
                "smiles": str(smiles),
                **{
                    key: row.get(key)
                    for key in METRIC_KEYS
                },
            }
        )
    return records


def _coverage(rows: list[dict[str, Any]], metric: str) -> dict[str, Any]:
    eligible = [
        row
        for row in rows
        if row["domain_status"] == "in_domain"
        and isinstance(row.get(f"{metric}_uncertainty"), (int, float))
        and isinstance(row.get(metric), (int, float))
    ]
    covered = sum(
        float(row[metric])
        <= float(row[f"{metric}_uncertainty"])
        for row in eligible
    )
    return {
        "eligible_count": len(eligible),
        "covered_count": covered,
        "coverage": (
            covered / len(eligible)
            if eligible
            else None
        ),
    }


def _evaluate_split(
    records: list[dict[str, Any]],
    cfg: dict[str, Any],
    *,
    seed: str,
) -> dict[str, Any]:
    calibration_ids, evaluation_ids = deterministic_split(
        [row["candidate_id"] for row in records],
        calibration_fraction=float(cfg["calibration_fraction"]),
        seed=seed,
    )
    calibration_id_set = set(calibration_ids)
    evaluation_id_set = set(evaluation_ids)
    calibration_records = [
        row
        for row in records
        if row["candidate_id"] in calibration_id_set
    ]
    calibrations = {
        metric: calibrate_metric(
            calibration_records,
            error_key=metric,
            coverage=float(cfg["coverage_target"]),
            k_neighbors=int(cfg["k_neighbors"]),
            similarity_floor=float(cfg["similarity_floor"]),
        )
        for metric in METRIC_KEYS
    }

    evaluation_rows: list[dict[str, Any]] = []
    for record in records:
        if record["candidate_id"] not in evaluation_id_set:
            continue
        domain = assess_domain(
            record["smiles"],
            calibration_records,
            in_domain_similarity=float(
                cfg["in_domain_similarity"]
            ),
            edge_similarity=float(cfg["edge_similarity"]),
            minimum_neighbors=int(cfg["minimum_neighbors"]),
        )
        result = {
            **record,
            "domain_status": domain.status,
            "nearest_similarity": domain.nearest_similarity,
            "domain_neighbor_count": domain.neighbor_count,
            "uncertainty_supported": domain.uncertainty_valid,
        }
        for metric in METRIC_KEYS:
            result[f"{metric}_uncertainty"] = (
                predict_relative_uncertainty(
                    record["smiles"],
                    calibration_records,
                    calibrations[metric],
                    error_key=metric,
                    k_neighbors=int(cfg["k_neighbors"]),
                    similarity_floor=float(
                        cfg["similarity_floor"]
                    ),
                )
            )
        evaluation_rows.append(result)

    coverage = {
        metric: _coverage(evaluation_rows, metric)
        for metric in METRIC_KEYS
    }
    in_domain_count = sum(
        row["domain_status"] == "in_domain"
        for row in evaluation_rows
    )
    return {
        "seed": seed,
        "calibration_candidate_ids": calibration_ids,
        "evaluation_candidate_ids": evaluation_ids,
        "calibration_count": len(calibration_ids),
        "evaluation_count": len(evaluation_ids),
        "in_domain_evaluation_count": in_domain_count,
        "edge_evaluation_count": sum(
            row["domain_status"] == "edge"
            for row in evaluation_rows
        ),
        "out_of_domain_evaluation_count": sum(
            row["domain_status"] == "out_of_domain"
            for row in evaluation_rows
        ),
        "valid_for_coverage": (
            in_domain_count
            >= int(cfg["minimum_in_domain_evaluation_candidates"])
        ),
        "coverage": coverage,
        "metric_calibrations": calibrations,
    }


def _stats(values: list[float]) -> dict[str, Any]:
    finite = [
        float(value)
        for value in values
        if math.isfinite(float(value))
    ]
    if not finite:
        return {
            "count": 0,
            "min": None,
            "median": None,
            "max": None,
        }
    return {
        "count": len(finite),
        "min": min(finite),
        "median": statistics.median(finite),
        "max": max(finite),
    }


def _aggregate_repeated(
    splits: list[dict[str, Any]],
) -> dict[str, Any]:
    valid = [
        split
        for split in splits
        if split["valid_for_coverage"]
    ]
    metrics: dict[str, Any] = {}
    for metric in METRIC_KEYS:
        coverage_values = [
            split["coverage"][metric]["coverage"]
            for split in valid
            if split["coverage"][metric]["coverage"] is not None
        ]
        eligible_total = sum(
            split["coverage"][metric]["eligible_count"]
            for split in valid
        )
        covered_total = sum(
            split["coverage"][metric]["covered_count"]
            for split in valid
        )
        metrics[metric] = {
            "split_coverage": _stats(coverage_values),
            "pooled_eligible_count": eligible_total,
            "pooled_covered_count": covered_total,
            "pooled_coverage": (
                covered_total / eligible_total
                if eligible_total
                else None
            ),
        }

    return {
        "split_count": len(splits),
        "valid_split_count": len(valid),
        "valid_split_fraction": (
            len(valid) / len(splits)
            if splits
            else 0.0
        ),
        "in_domain_evaluation_count": _stats(
            [
                split["in_domain_evaluation_count"]
                for split in splits
            ]
        ),
        "metrics": metrics,
    }


def _family_challenge(
    records: list[dict[str, Any]],
    family_by_id: dict[str, str],
    cfg: dict[str, Any],
    *,
    minimum_family_size: int,
) -> dict[str, Any]:
    families = sorted(
        {
            family_by_id[row["candidate_id"]]
            for row in records
            if row["candidate_id"] in family_by_id
            and family_by_id[row["candidate_id"]] != "probe"
        }
    )
    result: dict[str, Any] = {}
    for family in families:
        evaluation = [
            row
            for row in records
            if family_by_id.get(row["candidate_id"]) == family
        ]
        if len(evaluation) < minimum_family_size:
            continue
        calibration = [
            row
            for row in records
            if family_by_id.get(row["candidate_id"]) != family
        ]
        calibrations = {
            metric: calibrate_metric(
                calibration,
                error_key=metric,
                coverage=float(cfg["coverage_target"]),
                k_neighbors=int(cfg["k_neighbors"]),
                similarity_floor=float(cfg["similarity_floor"]),
            )
            for metric in METRIC_KEYS
        }
        rows: list[dict[str, Any]] = []
        for record in evaluation:
            domain = assess_domain(
                record["smiles"],
                calibration,
                in_domain_similarity=float(
                    cfg["in_domain_similarity"]
                ),
                edge_similarity=float(cfg["edge_similarity"]),
                minimum_neighbors=int(cfg["minimum_neighbors"]),
            )
            item = {
                **record,
                "domain_status": domain.status,
                "nearest_similarity": domain.nearest_similarity,
                "domain_neighbor_count": domain.neighbor_count,
            }
            for metric in METRIC_KEYS:
                item[f"{metric}_uncertainty"] = (
                    predict_relative_uncertainty(
                        record["smiles"],
                        calibration,
                        calibrations[metric],
                        error_key=metric,
                        k_neighbors=int(cfg["k_neighbors"]),
                        similarity_floor=float(
                            cfg["similarity_floor"]
                        ),
                    )
                )
            rows.append(item)

        result[family] = {
            "evaluation_count": len(rows),
            "calibration_count": len(calibration),
            "in_domain_count": sum(
                row["domain_status"] == "in_domain"
                for row in rows
            ),
            "edge_count": sum(
                row["domain_status"] == "edge"
                for row in rows
            ),
            "out_of_domain_count": sum(
                row["domain_status"] == "out_of_domain"
                for row in rows
            ),
            "coverage": {
                metric: _coverage(rows, metric)
                for metric in METRIC_KEYS
            },
            "candidates": [
                {
                    "candidate_id": row["candidate_id"],
                    "domain_status": row["domain_status"],
                    "nearest_similarity": row[
                        "nearest_similarity"
                    ],
                }
                for row in rows
            ],
        }
    return result


def _taxonomy_check(
    source_summary: dict[str, Any],
    benchmark: dict[str, Any],
) -> dict[str, Any]:
    candidate_by_id = {
        str(candidate["id"]): candidate
        for candidate in benchmark["holdouts"]
    }
    source_probes = (
        source_summary.get("applicability_domain", {})
        .get("domain_probe_results", [])
    )
    source_by_id = {
        str(probe["candidate_id"]): probe
        for probe in source_probes
    }
    expected_probes = [
        candidate
        for candidate in benchmark["holdouts"]
        if expected_prediction_outcome(candidate) != SUCCESS
    ]

    rows: list[dict[str, Any]] = []
    for candidate in expected_probes:
        candidate_id = str(candidate["id"])
        probe = source_by_id.get(candidate_id)
        expected = expected_prediction_outcome(candidate)
        if probe is None:
            rows.append(
                {
                    "candidate_id": candidate_id,
                    "expected_prediction_outcome": expected,
                    "observed_prediction_outcome": "missing_probe_result",
                    "matches": False,
                    "failure_reason": None,
                }
            )
            continue

        observed = probe.get("observed_prediction_outcome")
        if observed is None:
            observed = classify_prediction(
                result=None,
                failure_reason=probe.get("failure_reason"),
            )
        rows.append(
            {
                "candidate_id": candidate_id,
                "expected_prediction_outcome": expected,
                "observed_prediction_outcome": observed,
                "matches": observed == expected,
                "failure_reason": probe.get("failure_reason"),
            }
        )
    return {
        "probe_count": len(rows),
        "matching_count": sum(row["matches"] for row in rows),
        "mismatch_count": sum(not row["matches"] for row in rows),
        "probes": rows,
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--input-dir",
        type=Path,
        default=Path("property-v5b2-results"),
        help=(
            "Existing V5b-2 result directory. No FeOS/CoolProp "
            "trajectory rerun is required."
        ),
    )
    parser.add_argument(
        "--study",
        type=Path,
        default=HERE / "benchmark-v5b2.yaml",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("property-v5b21-results"),
    )
    args = parser.parse_args()

    if args.output_dir.exists() and any(args.output_dir.iterdir()):
        parser.error("output directory must be empty; use a new directory")

    summary_path = args.input_dir / "summary.json"
    if not summary_path.is_file():
        parser.error(f"missing V5b-2 summary: {summary_path}")

    summary = json.loads(summary_path.read_text(encoding="utf-8"))
    study = yaml.safe_load(args.study.read_text(encoding="utf-8"))
    records = _load_records(summary)
    if len(records) < int(
        study["uncertainty"]["minimum_comparable_candidates"]
    ):
        parser.error(
            "not enough comparable V5b-2 candidates for robustness analysis"
        )

    args.output_dir.mkdir(parents=True, exist_ok=True)
    uncertainty_cfg = study["uncertainty"]
    robustness_cfg = study["robustness"]

    split_count = int(robustness_cfg["repeated_split_count"])
    base_seed = str(uncertainty_cfg["split_seed"])
    repeated = [
        _evaluate_split(
            records,
            uncertainty_cfg,
            seed=f"{base_seed}:repeat:{index:03d}",
        )
        for index in range(split_count)
    ]
    repeated_aggregate = _aggregate_repeated(repeated)

    family_by_id = {
        str(candidate["id"]): str(
            candidate.get("challenge_family", "unassigned")
        )
        for candidate in study["holdouts"]
        if expected_prediction_outcome(candidate) == SUCCESS
    }
    family = _family_challenge(
        records,
        family_by_id,
        uncertainty_cfg,
        minimum_family_size=int(
            robustness_cfg["minimum_family_challenge_size"]
        ),
    )
    taxonomy = _taxonomy_check(summary, study)

    entry_summary = repeated_aggregate["metrics"]["entry_error"]
    pooled_entry = entry_summary["pooled_coverage"]
    median_entry = entry_summary["split_coverage"]["median"]
    valid_fraction = repeated_aggregate["valid_split_fraction"]

    family_entry_checks = {
        family_name: data["coverage"]["entry_error"]
        for family_name, data in family.items()
        if data["coverage"]["entry_error"]["eligible_count"] >= 2
    }
    family_entry_coverage_pass = all(
        isinstance(value["coverage"], (int, float))
        and float(value["coverage"])
        >= float(robustness_cfg["minimum_pooled_entry_coverage"])
        for value in family_entry_checks.values()
    )

    v5c_ready = (
        valid_fraction
        >= float(robustness_cfg["minimum_valid_split_fraction"])
        and isinstance(pooled_entry, (int, float))
        and pooled_entry
        >= float(robustness_cfg["minimum_pooled_entry_coverage"])
        and isinstance(median_entry, (int, float))
        and median_entry
        >= float(robustness_cfg["minimum_median_entry_coverage"])
        and family_entry_coverage_pass
        and taxonomy["mismatch_count"] == 0
    )

    report = {
        "study_version": "v5b-2.1",
        "study_complete": len(repeated) == split_count,
        "v5c_ready": v5c_ready,
        "source_result_dir": str(args.input_dir),
        "source_candidate_count": len(records),
        "repeated_split": repeated_aggregate,
        "family_challenge": family,
        "failure_taxonomy": taxonomy,
        "v5c_readiness_gate": {
            "minimum_valid_split_fraction": float(
                robustness_cfg["minimum_valid_split_fraction"]
            ),
            "minimum_pooled_entry_coverage": float(
                robustness_cfg["minimum_pooled_entry_coverage"]
            ),
            "minimum_median_entry_coverage": float(
                robustness_cfg["minimum_median_entry_coverage"]
            ),
            "observed_valid_split_fraction": valid_fraction,
            "observed_pooled_entry_coverage": pooled_entry,
            "observed_median_entry_coverage": median_entry,
            "taxonomy_mismatch_count": taxonomy["mismatch_count"],
            "family_entry_coverage_checks": family_entry_checks,
            "family_entry_coverage_pass": family_entry_coverage_pass,
        },
        "interpretation": (
            "Repeated splits reuse the same finite holdout population. "
            "They test split sensitivity and screening robustness; they do "
            "not create new independent experimental evidence."
        ),
    }
    write_summary(args.output_dir / "summary.json", report)
    write_summary(
        args.output_dir / "repeated_splits.json",
        {"splits": repeated},
    )
    write_summary(
        args.output_dir / "family_challenge.json",
        family,
    )

    print(json.dumps(
        {
            "study_complete": report["study_complete"],
            "v5c_ready": report["v5c_ready"],
            "repeated_split": repeated_aggregate,
            "failure_taxonomy": taxonomy,
            "family_challenge": {
                family_name: {
                    "evaluation_count": data["evaluation_count"],
                    "in_domain_count": data["in_domain_count"],
                    "edge_count": data["edge_count"],
                    "out_of_domain_count": data[
                        "out_of_domain_count"
                    ],
                    "entry_coverage": data["coverage"][
                        "entry_error"
                    ],
                }
                for family_name, data in family.items()
            },
        },
        indent=2,
    ))
    return 0 if report["study_complete"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
