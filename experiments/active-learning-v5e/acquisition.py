"""V5e-1 active-learning acquisition for calibration-domain expansion."""
from __future__ import annotations

from dataclasses import dataclass, asdict
from functools import lru_cache
import math
from typing import Any, Iterable


@dataclass(frozen=True)
class AcquisitionFeature:
    candidate_id: str
    smiles: str
    family: str
    property_priority_score: float
    property_percentile: float
    nearest_calibration_similarity: float
    calibration_novelty: float
    entry_predicted_ratio_vs_water: float | None
    entry_percentile: float | None
    acquisition_base_score: float
    domain_status: str | None
    lane_reason: str | None

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)


@lru_cache(maxsize=262144)
def _fingerprint(smiles: str):
    from rdkit import Chem
    from rdkit.Chem import rdFingerprintGenerator

    molecule = Chem.MolFromSmiles(smiles)
    if molecule is None:
        raise ValueError(f"invalid SMILES: {smiles}")
    generator = rdFingerprintGenerator.GetMorganGenerator(
        radius=2,
        fpSize=2048,
    )
    return generator.GetFingerprint(molecule)


def tanimoto(smiles_a: str, smiles_b: str) -> float:
    from rdkit import DataStructs

    return float(
        DataStructs.TanimotoSimilarity(
            _fingerprint(smiles_a),
            _fingerprint(smiles_b),
        )
    )


def nearest_similarity(
    smiles: str,
    references: Iterable[str],
) -> float:
    values = [
        tanimoto(smiles, reference)
        for reference in references
    ]
    return max(values) if values else 0.0


def _percentiles(
    values: dict[str, float],
    *,
    higher_is_better: bool,
) -> dict[str, float]:
    if not values:
        return {}
    ordered = sorted(
        values.items(),
        key=lambda item: (
            item[1] if higher_is_better else -item[1],
            item[0],
        ),
    )
    if len(ordered) == 1:
        return {ordered[0][0]: 1.0}

    result: dict[str, float] = {}
    for index, (key, _) in enumerate(ordered):
        result[key] = index / (len(ordered) - 1)
    return result


def _finite_positive(value: Any) -> float | None:
    if not isinstance(value, (int, float)):
        return None
    value = float(value)
    if not math.isfinite(value) or value <= 0.0:
        return None
    return value


def prepare_features(
    candidates: list[dict[str, Any]],
    calibration_smiles: list[str],
    entry_by_smiles: dict[str, dict[str, Any]],
    config: dict[str, Any],
) -> list[AcquisitionFeature]:
    weights = config["acquisition"]["weights"]
    property_values: dict[str, float] = {}
    entry_values: dict[str, float] = {}
    raw: list[dict[str, Any]] = []

    seen: set[str] = set()
    for row in candidates:
        candidate_id = str(row["candidate_id"])
        smiles = str(row["smiles"])
        if smiles in seen:
            continue
        seen.add(smiles)

        metrics = row.get("property_screen") or {}
        property_score = _finite_positive(
            metrics.get("property_priority_score")
        )
        if property_score is None:
            continue

        family = str(
            (row.get("provenance") or [{}])[0].get(
                "family",
                "unknown",
            )
        )
        entry = entry_by_smiles.get(smiles)
        entry_ratio = (
            _finite_positive(entry.get("predicted_ratio_vs_water"))
            if entry is not None
            else None
        )

        property_values[candidate_id] = property_score
        if entry_ratio is not None:
            # Lower coolant ratio is better. Percentile helper expects a value
            # that can be ranked; use the ratio directly with lower-is-better.
            entry_values[candidate_id] = entry_ratio

        raw.append(
            {
                "candidate_id": candidate_id,
                "smiles": smiles,
                "family": family,
                "property_priority_score": property_score,
                "entry_predicted_ratio_vs_water": entry_ratio,
                "domain_status": (
                    (row.get("domain") or {}).get("status")
                ),
                "lane_reason": row.get("lane_reason"),
            }
        )

    property_pct = _percentiles(
        property_values,
        higher_is_better=True,
    )
    # _percentiles(lower-is-better) returns 0 for the worst and 1 for best.
    entry_pct = _percentiles(
        entry_values,
        higher_is_better=False,
    )

    features: list[AcquisitionFeature] = []
    for row in raw:
        candidate_id = row["candidate_id"]
        similarity = nearest_similarity(
            row["smiles"],
            calibration_smiles,
        )
        novelty = 1.0 - similarity
        entry_percentile = entry_pct.get(candidate_id)

        score = (
            float(weights["property_merit"])
            * property_pct[candidate_id]
            + float(weights["calibration_novelty"])
            * novelty
        )
        if entry_percentile is not None:
            score += (
                float(weights["entry_merit"])
                * entry_percentile
            )

        features.append(
            AcquisitionFeature(
                candidate_id=candidate_id,
                smiles=row["smiles"],
                family=row["family"],
                property_priority_score=row[
                    "property_priority_score"
                ],
                property_percentile=property_pct[candidate_id],
                nearest_calibration_similarity=similarity,
                calibration_novelty=novelty,
                entry_predicted_ratio_vs_water=row[
                    "entry_predicted_ratio_vs_water"
                ],
                entry_percentile=entry_percentile,
                acquisition_base_score=score,
                domain_status=row["domain_status"],
                lane_reason=row["lane_reason"],
            )
        )

    features.sort(
        key=lambda row: (
            -row.acquisition_base_score,
            row.candidate_id,
        )
    )
    return features


def _min_selected_distance(
    smiles: str,
    selected: list[AcquisitionFeature],
) -> float:
    if not selected:
        return 1.0
    return min(
        1.0 - tanimoto(smiles, row.smiles)
        for row in selected
    )


def novelty_band(
    nearest_similarity_value: float,
    config: dict[str, Any],
) -> str:
    bands = config["acquisition"]["novelty_bands"]
    if nearest_similarity_value >= float(bands["bridge_min_similarity"]):
        if nearest_similarity_value < float(
            bands["near_domain_min_similarity"]
        ):
            return "bridge"
        return "near_domain"
    return "frontier"


def _family_limits(
    features: list[AcquisitionFeature],
    config: dict[str, Any],
) -> tuple[dict[str, int], int]:
    acquisition = config["acquisition"]
    target_count = int(acquisition["target_count"])
    minimum = int(acquisition["minimum_per_family"])
    max_fraction = float(acquisition["maximum_family_fraction"])
    max_per_family = max(
        minimum,
        int(math.ceil(target_count * max_fraction)),
    )

    counts: dict[str, int] = {}
    for feature in features:
        counts[feature.family] = counts.get(feature.family, 0) + 1

    floors = {
        family: min(minimum, count)
        for family, count in counts.items()
    }
    return floors, max_per_family


def select_targets(
    features: list[AcquisitionFeature],
    config: dict[str, Any],
) -> list[dict[str, Any]]:
    acquisition = config["acquisition"]
    target_count = min(
        int(acquisition["target_count"]),
        len(features),
    )
    if target_count <= 0:
        return []

    diversity_weight = float(
        acquisition["weights"]["selected_diversity"]
    )
    floors, max_per_family = _family_limits(
        features,
        config,
    )
    if sum(floors.values()) > target_count:
        raise ValueError(
            "minimum_per_family cannot be satisfied within target_count"
        )
    available_by_family: dict[str, int] = {}
    for feature in features:
        available_by_family[feature.family] = (
            available_by_family.get(feature.family, 0) + 1
        )
    maximum_selectable = sum(
        min(count, max_per_family)
        for count in available_by_family.values()
    )
    if maximum_selectable < target_count:
        raise ValueError(
            "maximum_family_fraction prevents selecting target_count"
        )

    selected: list[AcquisitionFeature] = []
    selected_ids: set[str] = set()
    family_counts: dict[str, int] = {}

    def choose_one(
        pool: list[AcquisitionFeature],
        *,
        family: str | None = None,
    ) -> AcquisitionFeature | None:
        best = None
        best_key = None
        for feature in pool:
            if feature.candidate_id in selected_ids:
                continue
            if family is not None and feature.family != family:
                continue
            if family_counts.get(feature.family, 0) >= max_per_family:
                continue

            diversity = _min_selected_distance(
                feature.smiles,
                selected,
            )
            dynamic_score = (
                feature.acquisition_base_score
                + diversity_weight * diversity
            )
            key = (
                dynamic_score,
                diversity,
                feature.property_percentile,
                feature.calibration_novelty,
                feature.candidate_id,
            )
            if best is None or key > best_key:
                best = feature
                best_key = key
        return best

    # Phase 1: guarantee family coverage where the pool permits it.
    for family in sorted(floors):
        for _ in range(floors[family]):
            feature = choose_one(features, family=family)
            if feature is None or len(selected) >= target_count:
                break
            selected.append(feature)
            selected_ids.add(feature.candidate_id)
            family_counts[family] = (
                family_counts.get(family, 0) + 1
            )

    # Phase 2: globally optimize merit + novelty + incremental diversity while
    # keeping a cap on any single family.
    while len(selected) < target_count:
        feature = choose_one(features)
        if feature is None:
            break
        selected.append(feature)
        selected_ids.add(feature.candidate_id)
        family_counts[feature.family] = (
            family_counts.get(feature.family, 0) + 1
        )

    output: list[dict[str, Any]] = []
    previous: list[AcquisitionFeature] = []
    for index, feature in enumerate(selected, start=1):
        min_distance = _min_selected_distance(
            feature.smiles,
            previous,
        )
        output.append(
            {
                **feature.to_dict(),
                "acquisition_rank": index,
                "minimum_distance_to_earlier_selected": min_distance,
                "acquisition_score_at_selection": (
                    feature.acquisition_base_score
                    + diversity_weight * min_distance
                ),
                "novelty_band": novelty_band(
                    feature.nearest_calibration_similarity,
                    config,
                ),
                "acquisition_role": (
                    "exploit"
                    if feature.property_percentile >= 0.75
                    else (
                        "frontier"
                        if feature.calibration_novelty >= 0.75
                        else "bridge"
                    )
                ),
                "reference_status": "unresolved",
            }
        )
        previous.append(feature)

    return output


def hypothetical_domain_expansion(
    pool: list[AcquisitionFeature],
    calibration_smiles: list[str],
    selected: list[dict[str, Any]],
    config: dict[str, Any],
) -> dict[str, Any]:
    thresholds = {
        "in_domain_similarity": float(
            config["uncertainty"]["in_domain_similarity"]
        ),
        "edge_similarity": float(
            config["uncertainty"]["edge_similarity"]
        ),
    }
    selected_smiles = [
        str(row["smiles"])
        for row in selected
    ]

    before: list[float] = []
    after: list[float] = []
    before_unselected: list[float] = []
    after_unselected: list[float] = []
    selected_set = set(selected_smiles)
    for feature in pool:
        before_value = nearest_similarity(
            feature.smiles,
            calibration_smiles,
        )
        after_value = nearest_similarity(
            feature.smiles,
            calibration_smiles + selected_smiles,
        )
        before.append(before_value)
        after.append(after_value)
        if feature.smiles not in selected_set:
            before_unselected.append(before_value)
            after_unselected.append(after_value)

    def summarize(values: list[float]) -> dict[str, Any]:
        if not values:
            return {
                "count": 0,
                "median_nearest_similarity": None,
                "fraction_ge_edge": None,
                "fraction_ge_in_domain_similarity": None,
            }
        ordered = sorted(values)
        middle = len(ordered) // 2
        median = (
            ordered[middle]
            if len(ordered) % 2 == 1
            else 0.5 * (
                ordered[middle - 1]
                + ordered[middle]
            )
        )
        return {
            "count": len(values),
            "median_nearest_similarity": median,
            "fraction_ge_edge": (
                sum(
                    value >= thresholds["edge_similarity"]
                    for value in values
                )
                / len(values)
            ),
            "fraction_ge_in_domain_similarity": (
                sum(
                    value >= thresholds["in_domain_similarity"]
                    for value in values
                )
                / len(values)
            ),
        }

    return {
        "important_limitation": (
            "This is a structural what-if calculation only. Selected targets "
            "do not expand the calibrated uncertainty domain until trusted "
            "reference data are acquired and V5b is rerun."
        ),
        "thresholds": thresholds,
        "before": summarize(before),
        "after_if_all_selected_are_validated": summarize(after),
        "before_excluding_selected_targets": summarize(
            before_unselected
        ),
        "after_excluding_selected_targets_if_validated": summarize(
            after_unselected
        ),
    }
