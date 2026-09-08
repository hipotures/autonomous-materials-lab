"""V5b-2 structural applicability domain and calibrated local error model."""
from __future__ import annotations

from dataclasses import dataclass, asdict
from functools import lru_cache
import hashlib
import math
from typing import Any, Iterable


@dataclass(frozen=True)
class DomainAssessment:
    nearest_similarity: float | None
    neighbor_count: int
    status: str
    uncertainty_valid: bool
    reason: str


@lru_cache(maxsize=4096)
def canonical_smiles(smiles: str) -> str:
    from rdkit import Chem

    molecule = Chem.MolFromSmiles(str(smiles))
    if molecule is None:
        raise ValueError(f"invalid SMILES: {smiles}")
    return Chem.MolToSmiles(molecule, canonical=True)


@lru_cache(maxsize=4096)
def fingerprint(smiles: str):
    from rdkit import Chem
    from rdkit.Chem import rdFingerprintGenerator

    molecule = Chem.MolFromSmiles(canonical_smiles(smiles))
    generator = rdFingerprintGenerator.GetMorganGenerator(
        radius=2,
        fpSize=2048,
    )
    return generator.GetFingerprint(molecule)


def tanimoto(smiles_a: str, smiles_b: str) -> float:
    from rdkit import DataStructs

    return float(
        DataStructs.TanimotoSimilarity(
            fingerprint(smiles_a),
            fingerprint(smiles_b),
        )
    )


def deterministic_split(
    candidate_ids: Iterable[str],
    *,
    calibration_fraction: float,
    seed: str,
) -> tuple[list[str], list[str]]:
    ids = sorted({str(candidate_id) for candidate_id in candidate_ids})
    if not ids:
        return [], []
    if not 0.0 < calibration_fraction < 1.0:
        raise ValueError("calibration_fraction must be in (0, 1)")

    ranked = sorted(
        ids,
        key=lambda candidate_id: hashlib.sha256(
            f"{seed}:{candidate_id}".encode("utf-8")
        ).hexdigest(),
    )
    calibration_count = max(
        1,
        min(
            len(ranked) - 1,
            int(math.ceil(len(ranked) * calibration_fraction)),
        ),
    )
    return ranked[:calibration_count], ranked[calibration_count:]


def similarity_neighbors(
    target_smiles: str,
    records: Iterable[dict[str, Any]],
    *,
    exclude_id: str | None = None,
) -> list[dict[str, Any]]:
    result: list[dict[str, Any]] = []
    for record in records:
        candidate_id = str(record["candidate_id"])
        if exclude_id is not None and candidate_id == exclude_id:
            continue
        similarity = tanimoto(target_smiles, str(record["smiles"]))
        result.append(
            {
                **record,
                "similarity": similarity,
            }
        )
    result.sort(
        key=lambda row: (
            -float(row["similarity"]),
            str(row["candidate_id"]),
        )
    )
    return result


def assess_domain(
    target_smiles: str,
    calibration_records: Iterable[dict[str, Any]],
    *,
    in_domain_similarity: float,
    edge_similarity: float,
    minimum_neighbors: int,
    exclude_id: str | None = None,
) -> DomainAssessment:
    neighbors = similarity_neighbors(
        target_smiles,
        calibration_records,
        exclude_id=exclude_id,
    )
    if not neighbors:
        return DomainAssessment(
            nearest_similarity=None,
            neighbor_count=0,
            status="out_of_domain",
            uncertainty_valid=False,
            reason="no calibrated structural neighbors",
        )

    nearest = float(neighbors[0]["similarity"])
    useful_count = sum(
        float(row["similarity"]) >= edge_similarity
        for row in neighbors
    )

    if nearest >= in_domain_similarity and useful_count >= minimum_neighbors:
        return DomainAssessment(
            nearest_similarity=nearest,
            neighbor_count=useful_count,
            status="in_domain",
            uncertainty_valid=True,
            reason="sufficient calibrated structural neighborhood",
        )
    if nearest >= edge_similarity:
        return DomainAssessment(
            nearest_similarity=nearest,
            neighbor_count=useful_count,
            status="edge",
            uncertainty_valid=False,
            reason=(
                "some structural support exists, but the calibrated "
                "neighborhood is too weak for certified uncertainty"
            ),
        )
    return DomainAssessment(
        nearest_similarity=nearest,
        neighbor_count=useful_count,
        status="out_of_domain",
        uncertainty_valid=False,
        reason="nearest calibrated structure is below similarity threshold",
    )


def local_error_scale(
    target_smiles: str,
    calibration_records: Iterable[dict[str, Any]],
    *,
    error_key: str,
    k_neighbors: int,
    similarity_floor: float,
    exclude_id: str | None = None,
) -> float | None:
    neighbors = [
        row
        for row in similarity_neighbors(
            target_smiles,
            calibration_records,
            exclude_id=exclude_id,
        )
        if isinstance(row.get(error_key), (int, float))
        and math.isfinite(float(row[error_key]))
    ]
    if not neighbors:
        return None

    selected = neighbors[: max(1, int(k_neighbors))]
    weighted_error = 0.0
    total_weight = 0.0
    for row in selected:
        similarity = max(
            float(row["similarity"]),
            float(similarity_floor),
        )
        weight = similarity * similarity
        weighted_error += weight * abs(float(row[error_key]))
        total_weight += weight

    if total_weight <= 0.0:
        return None
    return weighted_error / total_weight


def conformal_quantile(
    scores: Iterable[float],
    *,
    coverage: float,
) -> float | None:
    values = sorted(
        float(value)
        for value in scores
        if math.isfinite(float(value)) and float(value) >= 0.0
    )
    if not values:
        return None
    if not 0.0 < coverage < 1.0:
        raise ValueError("coverage must be in (0, 1)")

    # Finite-sample split-conformal order statistic:
    # ceil((n + 1) * coverage), clipped to the available sample.
    rank = int(math.ceil((len(values) + 1) * coverage))
    rank = max(1, min(rank, len(values)))
    return values[rank - 1]


def calibrate_metric(
    calibration_records: list[dict[str, Any]],
    *,
    error_key: str,
    coverage: float,
    k_neighbors: int,
    similarity_floor: float,
) -> dict[str, Any]:
    nonconformity: list[float] = []
    usable_ids: list[str] = []

    for record in calibration_records:
        error = record.get(error_key)
        if not isinstance(error, (int, float)):
            continue
        error = abs(float(error))
        if not math.isfinite(error):
            continue
        scale = local_error_scale(
            str(record["smiles"]),
            calibration_records,
            error_key=error_key,
            k_neighbors=k_neighbors,
            similarity_floor=similarity_floor,
            exclude_id=str(record["candidate_id"]),
        )
        if scale is None:
            continue
        scale = max(scale, 1.0e-12)
        nonconformity.append(error / scale)
        usable_ids.append(str(record["candidate_id"]))

    factor = conformal_quantile(
        nonconformity,
        coverage=coverage,
    )
    return {
        "error_key": error_key,
        "coverage_target": coverage,
        "calibration_count": len(nonconformity),
        "calibration_candidate_ids": usable_ids,
        "conformal_factor": factor,
        "nonconformity_scores": nonconformity,
    }


def predict_relative_uncertainty(
    target_smiles: str,
    calibration_records: list[dict[str, Any]],
    calibration: dict[str, Any],
    *,
    error_key: str,
    k_neighbors: int,
    similarity_floor: float,
    exclude_id: str | None = None,
) -> float | None:
    factor = calibration.get("conformal_factor")
    if not isinstance(factor, (int, float)) or not math.isfinite(float(factor)):
        return None
    scale = local_error_scale(
        target_smiles,
        calibration_records,
        error_key=error_key,
        k_neighbors=k_neighbors,
        similarity_floor=similarity_floor,
        exclude_id=exclude_id,
    )
    if scale is None:
        return None
    return float(factor) * scale


def domain_to_dict(value: DomainAssessment) -> dict[str, Any]:
    return asdict(value)
