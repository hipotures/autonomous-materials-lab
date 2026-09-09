"""Molecule-grouped benchmark metrics; diagnostic envelopes, not deployment gates."""
from __future__ import annotations

from collections import defaultdict
import hashlib
import math
from statistics import median
from typing import Any


def assignment(group: str, config: dict[str, Any]) -> str:
    digest = hashlib.sha256((config["split_seed"] + ":" + group).encode()).hexdigest()
    return "evaluation" if int(digest[:16], 16) / 2**64 < config["evaluation_fraction"] else "calibration"


def comparison_rows(observations: list[dict[str, Any]], predictions: list[dict[str, Any]],
                    config: dict[str, Any]) -> list[dict[str, Any]]:
    by_id = {}
    for p in predictions:
        if p["observation_id"] in by_id:
            raise ValueError("duplicate_prediction_observation_id")
        by_id[p["observation_id"]] = p
    result = []
    for obs in observations:
        pred = by_id.get(obs["observation_id"], {})
        value = pred.get("predicted_value")
        ok = (not isinstance(value, bool) and isinstance(value, (float, int))
              and math.isfinite(value) and value > 0)
        relative = value / obs["reference_value"] - 1 if ok else None
        result.append({**obs, "split": assignment(obs["molecule_group"], config),
                       "predicted_value": value if ok else None,
                       "status": "ok" if ok else pred.get("status", "prediction_missing"),
                       "prediction_error": pred.get("error"),
                       "signed_relative_error": relative,
                       "absolute_relative_error": abs(relative) if relative is not None else None})
    return result


def molecule_metrics(rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    grouped = defaultdict(list)
    for row in rows:
        grouped[(row["molecule_group"], row["property"], row["benchmark_eligible"])].append(row)
    results = []
    for (group, prop, strict), points in sorted(grouped.items()):
        errors = [x["absolute_relative_error"] for x in points if x["absolute_relative_error"] is not None]
        results.append({"molecule_group": group, "property": prop, "benchmark_eligible": strict,
                        "candidate_ids": sorted({x["candidate_id"] for x in points}),
                        "smiles": points[0]["smiles"], "family": points[0]["family"], "split": points[0]["split"],
                        "point_count": len(points), "successful_point_count": len(errors),
                        "failed_point_count": len(points) - len(errors),
                        "median_absolute_relative_error": median(errors) if errors else None,
                        "maximum_absolute_relative_error": max(errors) if errors else None,
                        "all_points_succeeded": len(errors) == len(points)})
    return results


def summarize(rows: list[dict[str, Any]], config: dict[str, Any]) -> dict[str, Any]:
    molecules = molecule_metrics(rows)
    output = {}
    for prop in sorted({r["property"] for r in molecules}):
        strict = [r for r in molecules if r["property"] == prop and r["benchmark_eligible"]]
        calibration = [r for r in strict if r["split"] == "calibration"]
        evaluation = [r for r in strict if r["split"] == "evaluation"]
        calibration_ok = [r for r in calibration if r["all_points_succeeded"]]
        # Max error per molecule prevents treating temperatures as independent data.
        errors = sorted(r["maximum_absolute_relative_error"] for r in calibration_ok)
        k = math.ceil((len(errors) + 1) * config["coverage_target"])
        enough = (len(errors) >= config["minimum_calibration_groups"] and 1 <= k <= len(errors)
                  and len(calibration_ok) == len(calibration))
        bound = errors[k - 1] if enough else None
        covered = sum(r["all_points_succeeded"] and r["maximum_absolute_relative_error"] <= bound
                      for r in evaluation) if bound is not None else None
        enough_eval = len(evaluation) >= config["minimum_evaluation_groups"]
        output[prop] = {
            "strict_molecule_group_count": len(strict),
            "diagnostic_only_group_count": sum(r["property"] == prop and not r["benchmark_eligible"] for r in molecules),
            "calibration_group_count": len(calibration), "evaluation_group_count": len(evaluation),
            "calibration_complete_prediction_group_count": len(calibration_ok),
            "failed_point_count": sum(r["failed_point_count"] for r in strict),
            "evaluation_median_molecule_error": median(
                [r["median_absolute_relative_error"] for r in evaluation if r["median_absolute_relative_error"] is not None]
            ) if any(r["median_absolute_relative_error"] is not None for r in evaluation) else None,
            "maximum_error_envelope": bound, "finite_sample_order_statistic": k,
            "evaluation_coverage_including_failures": covered / len(evaluation)
                if covered is not None and evaluation else None,
            "envelope_status": ("diagnostic_evaluated" if enough_eval else "insufficient_evaluation_groups")
                if bound is not None else "insufficient_calibration_or_prediction_failures",
            "coverage_guarantee_claimed": False,
            "selection_bias_note": "Actively selected molecules are not an exchangeable deployment sample.",
            "production_uncertainty_updated": False,
        }
    return {"properties": output, "molecule_metrics": molecules, "rankable_promotions": 0,
            "predictor_parameters_changed": False, "entry_uncertainty_updated": False}


def applicability(rows: list[dict[str, Any]], config: dict[str, Any]) -> list[dict[str, Any]]:
    from rdkit import Chem, DataStructs
    from rdkit.Chem import rdFingerprintGenerator
    generator = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
    cache = {}
    def fp(smiles):
        if smiles not in cache:
            mol = Chem.MolFromSmiles(smiles)
            if mol is None:
                raise ValueError("invalid_benchmark_smiles")
            cache[smiles] = generator.GetFingerprint(mol)
        return cache[smiles]
    groups = molecule_metrics(rows)
    output = []
    for row in groups:
        if row["split"] != "evaluation" or not row["benchmark_eligible"]:
            continue
        anchors = [r for r in groups if r["split"] == "calibration" and r["benchmark_eligible"]
                   and r["property"] == row["property"] and r["all_points_succeeded"]
                   and r["molecule_group"] != row["molecule_group"]]
        similarities = [(float(DataStructs.TanimotoSimilarity(fp(row["smiles"]), fp(a["smiles"]))),
                         a["molecule_group"]) for a in anchors]
        similarities.sort(reverse=True)
        nearest = similarities[0][0] if similarities else None
        neighbors = sum(s >= config["neighbor_similarity"] for s, _ in similarities)
        output.append({"molecule_group": row["molecule_group"], "property": row["property"],
                       "nearest_calibration_similarity": nearest, "neighbor_count": neighbors,
                       "structural_support_only": bool(nearest is not None and nearest >= config["in_domain_similarity"]
                                                       and neighbors >= config["minimum_neighbors"]),
                       "nearest_groups": similarities[:5], "rankable": False,
                       "calibrated_domain_expanded": False})
    return output
