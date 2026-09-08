#!/usr/bin/env python3
"""V5b-2: broad blind holdout, applicability domain and calibrated uncertainty."""
from __future__ import annotations

import argparse
import csv
import json
import math
import statistics
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path
from typing import Any

import yaml

from applicability import (
    assess_domain,
    calibrate_metric,
    deterministic_split,
    domain_to_dict,
    predict_relative_uncertainty,
)
from run_holdout import (
    _compare_candidate,
    _predict_candidate,
    _reference_candidate,
    _sha256,
    _write_csv,
    write_summary,
    ENTRY,
)

SUCCESS_STATUSES = {"terminal_velocity", "terminal_altitude"}
HERE = Path(__file__).resolve().parent


def _prediction_task(args: tuple[dict[str, Any], dict[str, Any], dict[str, Any], dict[str, Any]]):
    candidate, study, base, physical = args
    try:
        return {
            "candidate": candidate,
            "result": _predict_candidate(candidate, study, base, physical),
            "failure_reason": None,
        }
    except Exception as exc:
        return {
            "candidate": candidate,
            "result": None,
            "failure_reason": str(exc),
        }


def _reference_task(args: tuple[dict[str, Any], dict[str, Any], dict[str, Any], dict[str, Any]]):
    candidate, study, base, physical = args
    try:
        return {
            "candidate": candidate,
            "result": _reference_candidate(candidate, study, base, physical),
            "failure_reason": None,
        }
    except Exception as exc:
        return {
            "candidate": candidate,
            "result": None,
            "failure_reason": str(exc),
        }


def _parallel_map(fn, tasks: list[tuple], workers: int) -> list[dict[str, Any]]:
    if workers <= 1:
        return [fn(task) for task in tasks]
    rows: list[dict[str, Any]] = []
    with ProcessPoolExecutor(max_workers=workers) as pool:
        futures = [pool.submit(fn, task) for task in tasks]
        for future in as_completed(futures):
            rows.append(future.result())
    return rows


def _median_abs(values: list[float | None]) -> float | None:
    finite = [
        abs(float(value))
        for value in values
        if isinstance(value, (int, float))
        and math.isfinite(float(value))
    ]
    return statistics.median(finite) if finite else None


def _metric_records(
    comparisons: list[dict[str, Any]],
    property_rows: list[dict[str, Any]],
) -> list[dict[str, Any]]:
    by_candidate: dict[str, list[dict[str, Any]]] = {}
    for row in property_rows:
        by_candidate.setdefault(str(row["candidate_id"]), []).append(row)

    result: list[dict[str, Any]] = []
    for comparison in comparisons:
        candidate_id = str(comparison["candidate_id"])
        rows = by_candidate.get(candidate_id, [])
        saturation_errors = [
            item.get("relative_error")
            for item in comparison.get("saturation", [])
        ]
        result.append(
            {
                "candidate_id": candidate_id,
                "smiles": str(comparison["smiles"]),
                "entry_error": comparison.get(
                    "entry_coolant_mass_abs_relative_error"
                ),
                "storage_density_error": (
                    abs(float(comparison["storage_density_relative_error"]))
                    if isinstance(
                        comparison.get("storage_density_relative_error"),
                        (int, float),
                    )
                    else None
                ),
                "delta_h_error": _median_abs(
                    [row.get("delta_h_relative_error") for row in rows]
                ),
                "cp_error": _median_abs(
                    [row.get("cp_relative_error") for row in rows]
                ),
                "saturation_error": _median_abs(saturation_errors),
                "entry_comparable": bool(comparison.get("entry_comparable")),
                "predicted_coolant_kg": comparison.get("predicted_coolant_kg"),
                "reference_coolant_kg": comparison.get("reference_coolant_kg"),
                "predicted_storage_density_kg_m3": comparison.get(
                    "predicted_storage_density_kg_m3"
                ),
                "reference_storage_density_kg_m3": comparison.get(
                    "reference_storage_density_kg_m3"
                ),
            }
        )
    return result


def _interval(value: float | None, relative_half_width: float | None):
    if (
        not isinstance(value, (int, float))
        or not isinstance(relative_half_width, (int, float))
        or not math.isfinite(float(value))
        or not math.isfinite(float(relative_half_width))
    ):
        return None
    value = float(value)
    width = abs(value) * max(0.0, float(relative_half_width))
    return {
        "lower": max(0.0, value - width),
        "upper": value + width,
        "half_width": width,
        "relative_half_width": float(relative_half_width),
    }


def _observed_coverage(
    candidate_rows: list[dict[str, Any]],
    *,
    split_ids: set[str],
    error_key: str,
    uncertainty_key: str,
) -> dict[str, Any]:
    eligible = [
        row
        for row in candidate_rows
        if row["candidate_id"] in split_ids
        and row["domain"]["status"] == "in_domain"
        and row["domain"]["uncertainty_valid"]
        and isinstance(row.get(error_key), (int, float))
        and isinstance(row.get(uncertainty_key), (int, float))
    ]
    if not eligible:
        return {
            "eligible_count": 0,
            "covered_count": 0,
            "coverage": None,
        }
    covered = sum(
        float(row[error_key]) <= float(row[uncertainty_key])
        for row in eligible
    )
    return {
        "eligible_count": len(eligible),
        "covered_count": covered,
        "coverage": covered / len(eligible),
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--study",
        type=Path,
        default=HERE / "benchmark-v5b2.yaml",
    )
    parser.add_argument(
        "--base-config",
        type=Path,
        default=ENTRY / "config.yaml",
    )
    parser.add_argument(
        "--physical-study",
        type=Path,
        default=ENTRY / "physical-sensitivity.yaml",
    )
    parser.add_argument("--workers", type=int, default=16)
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("property-v5b2-results"),
    )
    args = parser.parse_args()

    if args.workers <= 0:
        parser.error("--workers must be positive")
    if args.output_dir.exists() and any(args.output_dir.iterdir()):
        parser.error("output directory must be empty; use a new directory")
    args.output_dir.mkdir(parents=True, exist_ok=True)

    study = yaml.safe_load(args.study.read_text(encoding="utf-8"))
    base = yaml.safe_load(args.base_config.read_text(encoding="utf-8"))
    physical = yaml.safe_load(args.physical_study.read_text(encoding="utf-8"))
    candidates = list(study["holdouts"])
    uncertainty_cfg = study["uncertainty"]

    prediction_tasks = [
        (candidate, study, base, physical)
        for candidate in candidates
    ]
    prediction_results = _parallel_map(
        _prediction_task,
        prediction_tasks,
        args.workers,
    )
    prediction_results.sort(
        key=lambda row: str(row["candidate"]["id"])
    )

    predictions: dict[str, dict[str, Any]] = {}
    prediction_failures: list[dict[str, Any]] = []
    domain_probe_results: list[dict[str, Any]] = []
    for row in prediction_results:
        candidate = row["candidate"]
        candidate_id = str(candidate["id"])
        expected_support = bool(candidate.get("expected_model_support", True))
        if row["result"] is not None:
            predictions[candidate_id] = row["result"]
            print(
                f"predicted {candidate_id}: "
                f"{row['result']['entry'].get('status')} "
                f"coolant={row['result']['entry'].get('coolant_used_kg')}",
                flush=True,
            )
            if not expected_support:
                domain_probe_results.append(
                    {
                        "candidate_id": candidate_id,
                        "expected_model_support": False,
                        "observed_model_support": True,
                        "failure_reason": None,
                    }
                )
        else:
            failure = {
                "candidate_id": candidate_id,
                "smiles": candidate["smiles"],
                "expected_model_support": expected_support,
                "failure_reason": row["failure_reason"],
            }
            prediction_failures.append(failure)
            print(
                f"prediction failed {candidate_id}: {row['failure_reason']}",
                flush=True,
            )
            if not expected_support:
                domain_probe_results.append(
                    {
                        "candidate_id": candidate_id,
                        "expected_model_support": False,
                        "observed_model_support": False,
                        "failure_reason": row["failure_reason"],
                    }
                )

    # Blind protocol boundary: reference backend has not been opened yet.
    prediction_artifact = {
        "study_version": "v5b-2",
        "reference_properties_used_during_prediction": False,
        "predictions": list(predictions.values()),
        "failures": prediction_failures,
        "domain_probe_results": domain_probe_results,
    }
    write_summary(
        args.output_dir / "predictions_before_reference.json",
        prediction_artifact,
    )

    candidate_by_id = {
        str(candidate["id"]): candidate
        for candidate in candidates
    }
    reference_candidates = [
        candidate_by_id[candidate_id]
        for candidate_id in predictions
    ]
    reference_tasks = [
        (candidate, study, base, physical)
        for candidate in reference_candidates
    ]
    reference_results = _parallel_map(
        _reference_task,
        reference_tasks,
        args.workers,
    )
    references: dict[str, dict[str, Any]] = {}
    reference_failures: list[dict[str, Any]] = []
    for row in reference_results:
        candidate_id = str(row["candidate"]["id"])
        if row["result"] is not None:
            references[candidate_id] = row["result"]
        else:
            reference_failures.append(
                {
                    "candidate_id": candidate_id,
                    "failure_reason": row["failure_reason"],
                }
            )

    property_rows: list[dict[str, Any]] = []
    comparisons: list[dict[str, Any]] = []
    for candidate_id in sorted(predictions):
        reference = references.get(candidate_id)
        if reference is None:
            continue
        rows, comparison = _compare_candidate(
            predictions[candidate_id],
            reference,
        )
        property_rows.extend(rows)
        comparisons.append(comparison)

    _write_csv(
        args.output_dir / "property_comparison.csv",
        property_rows,
    )
    _write_csv(
        args.output_dir / "entry_comparison.csv",
        comparisons,
    )

    metric_records = _metric_records(comparisons, property_rows)
    for record in metric_records:
        candidate = candidate_by_id[str(record["candidate_id"])]
        record["expected_model_support"] = bool(
            candidate.get("expected_model_support", True)
        )

    eligible_metric_records = [
        row
        for row in metric_records
        if row["entry_comparable"]
        and row["expected_model_support"]
    ]

    calibration_ids, evaluation_ids = deterministic_split(
        [row["candidate_id"] for row in eligible_metric_records],
        calibration_fraction=float(
            uncertainty_cfg["calibration_fraction"]
        ),
        seed=str(uncertainty_cfg["split_seed"]),
    )
    calibration_id_set = set(calibration_ids)
    evaluation_id_set = set(evaluation_ids)
    calibration_records = [
        row
        for row in eligible_metric_records
        if row["candidate_id"] in calibration_id_set
    ]

    metric_keys = [
        "entry_error",
        "storage_density_error",
        "delta_h_error",
        "cp_error",
        "saturation_error",
    ]
    metric_calibrations = {
        key: calibrate_metric(
            calibration_records,
            error_key=key,
            coverage=float(uncertainty_cfg["coverage_target"]),
            k_neighbors=int(uncertainty_cfg["k_neighbors"]),
            similarity_floor=float(uncertainty_cfg["similarity_floor"]),
        )
        for key in metric_keys
    }

    candidate_uncertainty: list[dict[str, Any]] = []
    for record in eligible_metric_records:
        candidate_id = str(record["candidate_id"])
        exclude_id = (
            candidate_id
            if candidate_id in calibration_id_set
            else None
        )
        domain = assess_domain(
            str(record["smiles"]),
            calibration_records,
            in_domain_similarity=float(
                uncertainty_cfg["in_domain_similarity"]
            ),
            edge_similarity=float(
                uncertainty_cfg["edge_similarity"]
            ),
            minimum_neighbors=int(
                uncertainty_cfg["minimum_neighbors"]
            ),
            exclude_id=exclude_id,
        )
        row = {
            **record,
            "split": (
                "calibration"
                if candidate_id in calibration_id_set
                else "evaluation"
            ),
            "domain": domain_to_dict(domain),
        }
        for key in metric_keys:
            uncertainty = predict_relative_uncertainty(
                str(record["smiles"]),
                calibration_records,
                metric_calibrations[key],
                error_key=key,
                k_neighbors=int(uncertainty_cfg["k_neighbors"]),
                similarity_floor=float(
                    uncertainty_cfg["similarity_floor"]
                ),
                exclude_id=exclude_id,
            )
            # A numerical bound outside the calibrated domain is retained as a
            # diagnostic, but is explicitly not certified.
            row[f"{key}_relative_uncertainty"] = uncertainty
            row[f"{key}_uncertainty_certified"] = bool(
                domain.uncertainty_valid
                and domain.status == "in_domain"
                and uncertainty is not None
            )

        row["predicted_coolant_interval_kg"] = _interval(
            record.get("predicted_coolant_kg"),
            row.get("entry_error_relative_uncertainty"),
        )
        row["predicted_storage_density_interval_kg_m3"] = _interval(
            record.get("predicted_storage_density_kg_m3"),
            row.get("storage_density_error_relative_uncertainty"),
        )
        candidate_uncertainty.append(row)

    coverage = {
        key: _observed_coverage(
            candidate_uncertainty,
            split_ids=evaluation_id_set,
            error_key=key,
            uncertainty_key=f"{key}_relative_uncertainty",
        )
        for key in metric_keys
    }

    unexpected_prediction_failures = [
        row
        for row in prediction_failures
        if row["expected_model_support"]
    ]
    comparable_count = len(eligible_metric_records)
    minimum_comparable = int(
        uncertainty_cfg["minimum_comparable_candidates"]
    )
    entry_calibration_ready = (
        metric_calibrations["entry_error"]["conformal_factor"]
        is not None
    )
    in_domain_evaluation_count = sum(
        row["split"] == "evaluation"
        and row["domain"]["status"] == "in_domain"
        and row["domain"]["uncertainty_valid"]
        for row in candidate_uncertainty
    )
    minimum_in_domain_evaluation = int(
        uncertainty_cfg["minimum_in_domain_evaluation_candidates"]
    )
    study_complete = (
        comparable_count >= minimum_comparable
        and len(evaluation_ids) >= 3
        and in_domain_evaluation_count >= minimum_in_domain_evaluation
        and entry_calibration_ready
    )

    calibration_model = {
        "study_version": "v5b-2",
        "method": (
            "Morgan radius-2 2048-bit Tanimoto neighborhood + "
            "split-conformal factor on local leave-one-out error scale"
        ),
        "settings": uncertainty_cfg,
        "calibration_candidate_ids": calibration_ids,
        "evaluation_candidate_ids": evaluation_ids,
        "calibration_records": calibration_records,
        "metric_calibrations": metric_calibrations,
        "usage_contract": {
            "in_domain": (
                "relative uncertainty may be used as calibrated screening bound"
            ),
            "edge": (
                "bound is diagnostic only; escalate validation before ranking"
            ),
            "out_of_domain": (
                "do not use uncertainty as calibrated; reject or escalate model"
            ),
        },
    }
    write_summary(
        args.output_dir / "uncertainty_model.json",
        calibration_model,
    )

    uncertainty_rows_flat = []
    for row in candidate_uncertainty:
        uncertainty_rows_flat.append(
            {
                "candidate_id": row["candidate_id"],
                "smiles": row["smiles"],
                "split": row["split"],
                "domain_status": row["domain"]["status"],
                "nearest_similarity": row["domain"]["nearest_similarity"],
                "domain_neighbor_count": row["domain"]["neighbor_count"],
                "uncertainty_valid": row["domain"]["uncertainty_valid"],
                "entry_error": row["entry_error"],
                "entry_relative_uncertainty": row[
                    "entry_error_relative_uncertainty"
                ],
                "storage_density_error": row["storage_density_error"],
                "storage_density_relative_uncertainty": row[
                    "storage_density_error_relative_uncertainty"
                ],
                "delta_h_error": row["delta_h_error"],
                "delta_h_relative_uncertainty": row[
                    "delta_h_error_relative_uncertainty"
                ],
                "cp_error": row["cp_error"],
                "cp_relative_uncertainty": row[
                    "cp_error_relative_uncertainty"
                ],
                "saturation_error": row["saturation_error"],
                "saturation_relative_uncertainty": row[
                    "saturation_error_relative_uncertainty"
                ],
            }
        )
    _write_csv(
        args.output_dir / "candidate_uncertainty.csv",
        uncertainty_rows_flat,
    )

    report = {
        "study_complete": study_complete,
        "candidate_count": len(candidates),
        "prediction_success_count": len(predictions),
        "prediction_failure_count": len(prediction_failures),
        "unexpected_prediction_failure_count": len(
            unexpected_prediction_failures
        ),
        "reference_failure_count": len(reference_failures),
        "entry_comparable_count": comparable_count,
        "minimum_comparable_candidates": minimum_comparable,
        "in_domain_evaluation_count": in_domain_evaluation_count,
        "minimum_in_domain_evaluation_candidates": (
            minimum_in_domain_evaluation
        ),
        "reference_revealed_after_prediction": True,
        "structure_only_prediction": True,
        "calibration_split": {
            "calibration_candidate_ids": calibration_ids,
            "evaluation_candidate_ids": evaluation_ids,
        },
        "applicability_domain": {
            "in_domain_count": sum(
                row["domain"]["status"] == "in_domain"
                for row in candidate_uncertainty
            ),
            "edge_count": sum(
                row["domain"]["status"] == "edge"
                for row in candidate_uncertainty
            ),
            "out_of_domain_count": sum(
                row["domain"]["status"] == "out_of_domain"
                for row in candidate_uncertainty
            ),
            "domain_probe_results": domain_probe_results,
        },
        "uncertainty": {
            "coverage_target": float(
                uncertainty_cfg["coverage_target"]
            ),
            "metric_calibrations": metric_calibrations,
            "evaluation_in_domain_observed_coverage": coverage,
            "status": (
                "split-conformal factor calibrated on local structural "
                "error scale; certified only for in-domain candidates"
            ),
        },
        "candidates": candidate_uncertainty,
        "unexpected_prediction_failures": unexpected_prediction_failures,
        "reference_failures": reference_failures,
    }
    write_summary(
        args.output_dir / "summary.json",
        report,
    )

    manifest = {
        "study_version": "v5b-2",
        "study_sha256": _sha256(args.study),
        "v5b1_prediction_contract_preserved": True,
        "reference_properties_used_during_prediction": False,
        "applicability_descriptor": (
            "RDKit Morgan fingerprint radius=2 fpSize=2048"
        ),
        "similarity": "Tanimoto",
        "uncertainty_method": calibration_model["method"],
        "important_limitations": [
            "Calibration set is small and chemically nonuniform.",
            "Observed coverage is descriptive, not proof of exchangeability.",
            "Applicability thresholds are screening thresholds, not physical laws.",
            "Transport, decomposition chemistry, mixtures and system mass remain outside V5b-2.",
        ],
    }
    write_summary(
        args.output_dir / "manifest.json",
        manifest,
    )

    print(json.dumps(
        {
            "study_complete": report["study_complete"],
            "candidate_count": report["candidate_count"],
            "prediction_success_count": report["prediction_success_count"],
            "unexpected_prediction_failure_count": report[
                "unexpected_prediction_failure_count"
            ],
            "entry_comparable_count": report["entry_comparable_count"],
            "applicability_domain": report["applicability_domain"],
            "evaluation_in_domain_observed_coverage": coverage,
        },
        indent=2,
    ))
    return 0 if study_complete else 1


if __name__ == "__main__":
    raise SystemExit(main())
