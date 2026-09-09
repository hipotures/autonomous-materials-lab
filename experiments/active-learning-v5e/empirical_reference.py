"""V5e-2.1 multi-source empirical thermophysical reference resolution.

This stage extends V5e-2 beyond CoolProp without weakening the blind-reference
contract. Only explicitly whitelisted experimental, reviewed, tabulated, or
data-fitted methods from chemicals/thermo may be selected. Structure-only,
group-contribution, CSP, corresponding-state, or generic EOS estimation methods
remain visible in the audit but are never calibration evidence.
"""
from __future__ import annotations

import math
from typing import Any, Callable, Iterable


QUALITY_RANK = {
    "insufficient": 0,
    "partial": 1,
    "medium": 2,
    "high": 3,
}


def _finite_positive(value: Any) -> float | None:
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        return None
    value = float(value)
    if not math.isfinite(value) or value <= 0.0:
        return None
    return value


def selector_matches(method: str, selector: str) -> bool:
    method = str(method)
    selector = str(selector)
    if selector.endswith("*"):
        return method.startswith(selector[:-1])
    return method == selector


def ordered_allowed_methods(
    available: Iterable[str],
    selectors: Iterable[str],
) -> list[str]:
    available_set = {str(method) for method in available}
    result: list[str] = []
    seen: set[str] = set()
    for selector in selectors:
        matches = sorted(
            method
            for method in available_set
            if selector_matches(method, str(selector))
        )
        for method in matches:
            if method not in seen:
                seen.add(method)
                result.append(method)
    return result


def source_family(method: str) -> str:
    value = str(method).upper()
    if value.startswith("DIRECT_NIST") or "WEBBOOK" in value:
        return "NIST"
    if value.startswith("HEOS") or value == "COOLPROP":
        return "HEOS_REFPROP_COOLPROP"
    if value.startswith("IUPAC"):
        return "IUPAC"
    if value.startswith("MATTHEWS"):
        return "MATTHEWS"
    if value == "PD":
        return "PD_COMPILATION"
    if value.startswith("CRC"):
        return "CRC"
    if value.startswith("COMMON_CHEMISTRY"):
        return "CAS_COMMON_CHEMISTRY"
    if value.startswith("ZABRANSKY"):
        return "ZABRANSKY"
    if value.startswith("VDI"):
        return "VDI"
    if value.startswith("DIPPR"):
        return "DIPPR_PERRY"
    if value.startswith("POLING") or value.startswith("WAGNER_POLING"):
        return "POLING"
    if value.startswith("WAGNER_MCGARRY"):
        return "MCGARRY"
    if value.startswith("LANDOLT"):
        return "LANDOLT"
    if value.startswith("ALCOCK"):
        return "ALCOCK"
    return value


def evidence_class(method: str) -> str:
    value = str(method).upper()
    if value.startswith("HEOS") or value == "COOLPROP":
        return "reviewed_eos_or_fit"
    if (
        value.startswith("IUPAC")
        or value.startswith("MATTHEWS")
        or value.startswith("DIRECT_NIST")
    ):
        return "experimental_review_or_compilation"
    if (
        value.startswith("CRC")
        or value == "PD"
        or "WEBBOOK" in value
        or value.startswith("COMMON_CHEMISTRY")
    ):
        return "processed_or_compiled_experimental"
    if (
        value.startswith("ZABRANSKY")
        or value.startswith("VDI")
        or value.startswith("DIPPR")
        or value.startswith("WAGNER")
        or value.startswith("LANDOLT")
        or value.startswith("ALCOCK")
    ):
        return "correlation_fit_to_reference_data"
    if value.startswith("POLING"):
        return "published_reference_data"
    return "approved_reference_method"


def selected_method_is_forbidden(
    method: str,
    forbidden_selectors: Iterable[str],
) -> bool:
    upper = str(method).upper()
    return any(
        selector_matches(upper, str(selector).upper())
        for selector in forbidden_selectors
    )


def select_constant(
    *,
    property_name: str,
    unit: str,
    available_methods: Callable[[str], Iterable[str]],
    getter: Callable[..., Any],
    cas_number: str,
    allowed_methods: Iterable[str],
) -> tuple[dict[str, Any] | None, list[dict[str, Any]]]:
    try:
        available = [str(item) for item in available_methods(cas_number)]
    except Exception as exc:
        return None, [{
            "property": property_name,
            "method": None,
            "available": False,
            "allowed": False,
            "selected": False,
            "reason": f"method_inventory_failed:{type(exc).__name__}",
        }]

    ordered = ordered_allowed_methods(available, allowed_methods)
    selected: dict[str, Any] | None = None
    selected_method: str | None = None
    failures: dict[str, str] = {}
    for method in ordered:
        try:
            value = _finite_positive(
                getter(cas_number, method=method)
            )
        except Exception as exc:
            failures[method] = f"calculation_failed:{type(exc).__name__}"
            continue
        if value is None:
            failures[method] = "nonpositive_or_nonfinite"
            continue
        selected_method = method
        selected = {
            "property": property_name,
            "value": value,
            "unit": unit,
            "method": method,
            "source_family": source_family(method),
            "evidence_class": evidence_class(method),
        }
        break

    audit: list[dict[str, Any]] = []
    allowed_set = set(ordered)
    for method in sorted(set(available)):
        is_allowed = method in allowed_set
        audit.append({
            "property": property_name,
            "method": method,
            "available": True,
            "allowed": is_allowed,
            "selected": method == selected_method,
            "reason": (
                "selected"
                if method == selected_method
                else failures.get(method)
                if is_allowed
                else "method_not_in_empirical_allowlist"
            ),
        })
    return selected, audit


def _method_valid_at(obj: Any, temperature_k: float, method: str) -> bool:
    try:
        return bool(obj.test_method_validity(float(temperature_k), method))
    except Exception:
        return False


def _transform_value(
    raw: Any,
    transform: Callable[[float], float],
) -> float | None:
    value = _finite_positive(raw)
    if value is None:
        return None
    try:
        converted = _finite_positive(transform(value))
    except Exception:
        return None
    return converted


def select_temperature_series(
    *,
    property_name: str,
    unit: str,
    obj: Any,
    temperatures_k: Iterable[float],
    allowed_methods: Iterable[str],
    transform: Callable[[float], float] = lambda value: value,
    state_basis: str,
) -> tuple[dict[str, Any] | None, list[dict[str, Any]]]:
    available = sorted(
        str(item)
        for item in getattr(obj, "all_methods", set())
    )
    ordered = ordered_allowed_methods(available, allowed_methods)
    candidates: list[tuple[int, int, str, list[dict[str, float]]]] = []
    failures: dict[str, str] = {}

    for priority, method in enumerate(ordered):
        points: list[dict[str, float]] = []
        for temperature_k in temperatures_k:
            temperature_k = float(temperature_k)
            if not _method_valid_at(obj, temperature_k, method):
                continue
            try:
                raw = obj.calculate(temperature_k, method)
            except Exception:
                continue
            value = _transform_value(raw, transform)
            if value is None:
                continue
            points.append({
                "temperature_k": temperature_k,
                "value": value,
            })
        if points:
            candidates.append(
                (len(points), -priority, method, points)
            )
        else:
            failures[method] = "no_valid_grid_points"

    selected: dict[str, Any] | None = None
    selected_method: str | None = None
    if candidates:
        _, _, selected_method, points = max(
            candidates,
            key=lambda item: (item[0], item[1], item[2]),
        )
        selected = {
            "property": property_name,
            "unit": unit,
            "method": selected_method,
            "source_family": source_family(selected_method),
            "evidence_class": evidence_class(selected_method),
            "state_basis": state_basis,
            "points": points,
            "temperature_min_k": min(
                point["temperature_k"] for point in points
            ),
            "temperature_max_k": max(
                point["temperature_k"] for point in points
            ),
        }

    allowed_set = set(ordered)
    audit = [
        {
            "property": property_name,
            "method": method,
            "available": True,
            "allowed": method in allowed_set,
            "selected": method == selected_method,
            "reason": (
                "selected"
                if method == selected_method
                else failures.get(method)
                if method in allowed_set
                else "method_not_in_empirical_allowlist"
            ),
        }
        for method in available
    ]
    return selected, audit


def _direct_nist_constant(
    direct_nist: dict[str, Any] | None,
    *,
    cas_number: str,
    key: str,
    property_name: str,
    unit: str,
) -> dict[str, Any] | None:
    if not direct_nist:
        return None
    if str(direct_nist.get("cas_number") or "") != cas_number:
        return None
    if not bool(direct_nist.get("compound_page_verified")):
        return None
    value = _finite_positive(direct_nist.get(key))
    if value is None:
        return None
    method = "DIRECT_NIST_V5E2"
    return {
        "property": property_name,
        "value": value,
        "unit": unit,
        "method": method,
        "source_family": source_family(method),
        "evidence_class": evidence_class(method),
    }


def _safe_construct(factory: Callable[[], Any]) -> tuple[Any | None, str | None]:
    try:
        return factory(), None
    except Exception as exc:
        return None, f"{type(exc).__name__}:{exc}"


def build_empirical_packet(
    *,
    cas_number: str,
    molecular_weight_g_mol: float,
    config: dict[str, Any],
    direct_nist: dict[str, Any] | None = None,
) -> dict[str, Any]:
    from chemicals import critical, phase_change
    from thermo import (
        EnthalpyVaporization,
        HeatCapacityLiquid,
        VaporPressure,
        VolumeLiquid,
    )

    policy = config["source_policy"]
    grid = [
        float(value)
        for value in config["properties"]["temperatures_k"]
    ]
    mw = _finite_positive(molecular_weight_g_mol)
    if mw is None:
        raise ValueError("molecular_weight_g_mol must be positive")

    audit: list[dict[str, Any]] = []
    constants: dict[str, dict[str, Any] | None] = {}

    tc, rows = select_constant(
        property_name="critical_temperature_k",
        unit="K",
        available_methods=critical.Tc_methods,
        getter=critical.Tc,
        cas_number=cas_number,
        allowed_methods=policy["critical_temperature_methods"],
    )
    audit.extend(rows)
    if tc is None:
        tc = _direct_nist_constant(
            direct_nist,
            cas_number=cas_number,
            key="critical_temperature_k",
            property_name="critical_temperature_k",
            unit="K",
        )
    constants["critical_temperature_k"] = tc

    pc, rows = select_constant(
        property_name="critical_pressure_pa",
        unit="Pa",
        available_methods=critical.Pc_methods,
        getter=critical.Pc,
        cas_number=cas_number,
        allowed_methods=policy["critical_pressure_methods"],
    )
    audit.extend(rows)
    if pc is None:
        pc = _direct_nist_constant(
            direct_nist,
            cas_number=cas_number,
            key="critical_pressure_pa",
            property_name="critical_pressure_pa",
            unit="Pa",
        )
    constants["critical_pressure_pa"] = pc

    tb, rows = select_constant(
        property_name="normal_boiling_temperature_k",
        unit="K",
        available_methods=phase_change.Tb_methods,
        getter=phase_change.Tb,
        cas_number=cas_number,
        allowed_methods=policy["normal_boiling_temperature_methods"],
    )
    audit.extend(rows)
    if tb is None:
        tb = _direct_nist_constant(
            direct_nist,
            cas_number=cas_number,
            key="normal_boiling_temperature_k",
            property_name="normal_boiling_temperature_k",
            unit="K",
        )
    constants["normal_boiling_temperature_k"] = tb

    tc_value = tc["value"] if tc else None
    pc_value = pc["value"] if pc else None
    tb_value = tb["value"] if tb else None

    series: dict[str, dict[str, Any] | None] = {}
    construction_errors: dict[str, str] = {}

    vp_obj, error = _safe_construct(
        lambda: VaporPressure(
            Tb=tb_value,
            Tc=tc_value,
            Pc=pc_value,
            CASRN=cas_number,
            extrapolation=None,
        )
    )
    if error:
        construction_errors["vapor_pressure_pa"] = error
    if vp_obj is not None:
        record, rows = select_temperature_series(
            property_name="vapor_pressure_pa",
            unit="Pa",
            obj=vp_obj,
            temperatures_k=grid,
            allowed_methods=policy["vapor_pressure_methods"],
            state_basis="saturation_pressure",
        )
        series["vapor_pressure_pa"] = record
        audit.extend(rows)
    else:
        series["vapor_pressure_pa"] = None

    volume_obj, error = _safe_construct(
        lambda: VolumeLiquid(
            MW=mw,
            Tb=tb_value,
            Tc=tc_value,
            Pc=pc_value,
            CASRN=cas_number,
            extrapolation=None,
        )
    )
    if error:
        construction_errors["liquid_density_kg_m3"] = error
    if volume_obj is not None:
        record, rows = select_temperature_series(
            property_name="liquid_density_kg_m3",
            unit="kg/m^3",
            obj=volume_obj,
            temperatures_k=grid,
            allowed_methods=policy["liquid_volume_methods"],
            transform=lambda vm: (mw / 1000.0) / vm,
            state_basis="low_pressure_or_saturation_liquid",
        )
        series["liquid_density_kg_m3"] = record
        audit.extend(rows)
    else:
        series["liquid_density_kg_m3"] = None

    cp_obj, error = _safe_construct(
        lambda: HeatCapacityLiquid(
            CASRN=cas_number,
            MW=mw,
            Tc=tc_value,
            extrapolation=None,
        )
    )
    if error:
        construction_errors["liquid_cp_j_kg_k"] = error
    if cp_obj is not None:
        record, rows = select_temperature_series(
            property_name="liquid_cp_j_kg_k",
            unit="J/kg/K",
            obj=cp_obj,
            temperatures_k=grid,
            allowed_methods=policy["liquid_heat_capacity_methods"],
            transform=lambda cp_molar: cp_molar * 1000.0 / mw,
            state_basis="low_pressure_or_saturation_liquid",
        )
        series["liquid_cp_j_kg_k"] = record
        audit.extend(rows)
    else:
        series["liquid_cp_j_kg_k"] = None

    hvap_obj, error = _safe_construct(
        lambda: EnthalpyVaporization(
            CASRN=cas_number,
            Tb=tb_value,
            Tc=tc_value,
            Pc=pc_value,
            Psat=vp_obj,
            extrapolation=None,
        )
    )
    if error:
        construction_errors["latent_heat_vaporization_j_kg"] = error
    if hvap_obj is not None:
        record, rows = select_temperature_series(
            property_name="latent_heat_vaporization_j_kg",
            unit="J/kg",
            obj=hvap_obj,
            temperatures_k=grid,
            allowed_methods=policy["enthalpy_vaporization_methods"],
            transform=lambda hvap_molar: hvap_molar * 1000.0 / mw,
            state_basis="saturation_curve",
        )
        series["latent_heat_vaporization_j_kg"] = record
        audit.extend(rows)
    else:
        series["latent_heat_vaporization_j_kg"] = None

    return {
        "cas_number": cas_number,
        "molecular_weight_g_mol": mw,
        "constants": constants,
        "series": series,
        "source_audit": audit,
        "construction_errors": construction_errors,
    }


def classify_packet(
    packet: dict[str, Any],
    config: dict[str, Any],
) -> dict[str, Any]:
    gate = config["calibration_gate"]
    constants = packet.get("constants") or {}
    series = packet.get("series") or {}

    present_constants = sorted(
        key for key, value in constants.items() if value
    )
    present_series = sorted(
        key
        for key, value in series.items()
        if value and value.get("points")
    )
    present_properties = present_constants + present_series

    required = [
        str(value)
        for value in gate["required_properties"]
    ]
    dynamic = [
        str(value)
        for value in gate["dynamic_properties"]
    ]
    missing_required = [
        name
        for name in required
        if name not in present_properties
    ]
    dynamic_count = sum(
        name in present_properties
        for name in dynamic
    )

    selected_records = [
        value for value in constants.values() if value
    ] + [
        value
        for value in series.values()
        if value and value.get("points")
    ]
    families = sorted({
        str(record["source_family"])
        for record in selected_records
        if record.get("source_family")
    })
    forbidden = sorted({
        str(record["method"])
        for record in selected_records
        if selected_method_is_forbidden(
            str(record["method"]),
            gate["forbidden_method_selectors"],
        )
    })

    total_count = len(present_properties)
    ready = (
        not missing_required
        and dynamic_count >= int(gate["minimum_dynamic_properties"])
        and total_count >= int(gate["minimum_total_properties"])
        and len(families) >= int(gate["minimum_source_families"])
        and not forbidden
    )

    if ready:
        high = gate["high_quality"]
        if (
            dynamic_count >= int(high["minimum_dynamic_properties"])
            and total_count >= int(high["minimum_total_properties"])
            and len(families) >= int(high["minimum_source_families"])
        ):
            quality = "high"
        else:
            quality = "medium"
    elif total_count > 0:
        quality = "partial"
    else:
        quality = "insufficient"

    reasons: list[str] = []
    if missing_required:
        reasons.append(
            "missing_required:" + ",".join(missing_required)
        )
    if dynamic_count < int(gate["minimum_dynamic_properties"]):
        reasons.append(
            f"dynamic_coverage:{dynamic_count}/"
            f"{int(gate['minimum_dynamic_properties'])}"
        )
    if total_count < int(gate["minimum_total_properties"]):
        reasons.append(
            f"total_coverage:{total_count}/"
            f"{int(gate['minimum_total_properties'])}"
        )
    if len(families) < int(gate["minimum_source_families"]):
        reasons.append(
            f"source_families:{len(families)}/"
            f"{int(gate['minimum_source_families'])}"
        )
    if forbidden:
        reasons.append(
            "forbidden_selected_method:" + ",".join(forbidden)
        )

    return {
        "ready_for_v5b3_property_calibration": ready,
        "reference_quality": quality,
        "present_properties": present_properties,
        "missing_required_properties": missing_required,
        "dynamic_property_count": dynamic_count,
        "total_property_count": total_count,
        "source_families": families,
        "source_family_count": len(families),
        "forbidden_selected_methods": forbidden,
        "reasons": reasons,
    }


def _attempt_score(
    classification: dict[str, Any],
) -> tuple[int, int, int, int, int]:
    return (
        int(bool(classification["ready_for_v5b3_property_calibration"])),
        QUALITY_RANK.get(str(classification["reference_quality"]), -1),
        int(classification["dynamic_property_count"]),
        int(classification["total_property_count"]),
        int(classification["source_family_count"]),
    )


def resolve_empirical_target(
    target: dict[str, Any],
    config: dict[str, Any],
) -> dict[str, Any]:
    identity = target.get("identity") or {}
    pubchem = target.get("pubchem") or {}
    base = {
        "acquisition_rank": target.get("acquisition_rank"),
        "candidate_id": target.get("candidate_id"),
        "smiles": target.get("smiles"),
        "canonical_smiles": target.get("canonical_smiles"),
        "family": target.get("family"),
        "preferred_name": target.get("preferred_name"),
        "inchi_key": (
            pubchem.get("inchi_key")
            or target.get("expected_inchi_key")
        ),
        "pubchem_cid": pubchem.get("cid"),
        "v5e2_resolution_status": target.get("resolution_status"),
        "v5e2_reference_quality": target.get("reference_quality"),
        "entry_reference_ready": bool(
            target.get("usable_for_v5b_calibration")
        ),
    }
    if not bool(identity.get("identity_verified")):
        return {
            **base,
            "resolution_status": "unresolved_identity",
            "reference_quality": "insufficient",
            "ready_for_v5b3_property_calibration": False,
            "selected_cas_number": None,
            "classification": {
                "ready_for_v5b3_property_calibration": False,
                "reference_quality": "insufficient",
                "present_properties": [],
                "missing_required_properties": list(
                    config["calibration_gate"]["required_properties"]
                ),
                "dynamic_property_count": 0,
                "total_property_count": 0,
                "source_families": [],
                "source_family_count": 0,
                "forbidden_selected_methods": [],
                "reasons": ["identity_not_verified"],
            },
            "packet": None,
            "cas_attempts": [],
        }

    cas_numbers = [
        str(value)
        for value in pubchem.get("cas_numbers") or []
        if value
    ]
    nist = target.get("nist") or {}
    nist_cas = nist.get("cas_number")
    if nist_cas and str(nist_cas) not in cas_numbers:
        cas_numbers.append(str(nist_cas))
    cas_numbers = sorted(set(cas_numbers))

    mw = _finite_positive(
        pubchem.get("molecular_weight_g_mol")
        or target.get("expected_molecular_weight_g_mol")
    )
    if not cas_numbers or mw is None:
        return {
            **base,
            "resolution_status": (
                "resolved_identity_no_cas"
                if not cas_numbers
                else "resolved_identity_no_molecular_weight"
            ),
            "reference_quality": "insufficient",
            "ready_for_v5b3_property_calibration": False,
            "selected_cas_number": None,
            "classification": {
                "ready_for_v5b3_property_calibration": False,
                "reference_quality": "insufficient",
                "present_properties": [],
                "missing_required_properties": list(
                    config["calibration_gate"]["required_properties"]
                ),
                "dynamic_property_count": 0,
                "total_property_count": 0,
                "source_families": [],
                "source_family_count": 0,
                "forbidden_selected_methods": [],
                "reasons": [
                    "no_valid_cas" if not cas_numbers else "no_molecular_weight"
                ],
            },
            "packet": None,
            "cas_attempts": [],
        }

    attempts: list[dict[str, Any]] = []
    for cas_number in cas_numbers:
        try:
            packet = build_empirical_packet(
                cas_number=cas_number,
                molecular_weight_g_mol=mw,
                config=config,
                direct_nist=nist,
            )
            classification = classify_packet(packet, config)
            attempts.append({
                "cas_number": cas_number,
                "packet": packet,
                "classification": classification,
                "error": None,
            })
        except Exception as exc:
            attempts.append({
                "cas_number": cas_number,
                "packet": None,
                "classification": {
                    "ready_for_v5b3_property_calibration": False,
                    "reference_quality": "insufficient",
                    "present_properties": [],
                    "missing_required_properties": list(
                        config["calibration_gate"]["required_properties"]
                    ),
                    "dynamic_property_count": 0,
                    "total_property_count": 0,
                    "source_families": [],
                    "source_family_count": 0,
                    "forbidden_selected_methods": [],
                    "reasons": [f"resolver_exception:{type(exc).__name__}"],
                },
                "error": f"{type(exc).__name__}:{exc}",
            })

    chosen = max(
        attempts,
        key=lambda row: (
            _attempt_score(row["classification"]),
            row["cas_number"],
        ),
    )
    classification = chosen["classification"]
    ready = bool(
        classification["ready_for_v5b3_property_calibration"]
    )
    quality = str(classification["reference_quality"])
    if ready:
        status = "empirical_reference_ready"
    elif classification["total_property_count"] > 0:
        status = "empirical_reference_partial"
    else:
        status = "empirical_reference_insufficient"

    return {
        **base,
        "resolution_status": status,
        "reference_quality": quality,
        "ready_for_v5b3_property_calibration": ready,
        "selected_cas_number": chosen["cas_number"],
        "classification": classification,
        "packet": chosen["packet"],
        "cas_attempts": [
            {
                "cas_number": row["cas_number"],
                "classification": row["classification"],
                "error": row["error"],
            }
            for row in attempts
        ],
    }


def build_v5b3_property_handoff(
    rows: list[dict[str, Any]],
    config: dict[str, Any],
) -> dict[str, Any]:
    ready = [
        row
        for row in rows
        if row.get("ready_for_v5b3_property_calibration")
    ]
    minimum = int(config["handoff"]["minimum_ready_targets"])
    anchors = [
        {
            "id": str(row["candidate_id"]),
            "smiles": str(
                row.get("canonical_smiles")
                or row.get("smiles")
            ),
            "inchi_key": row.get("inchi_key"),
            "cas_number": row.get("selected_cas_number"),
            "challenge_family": row.get("family"),
            "reference_quality": row.get("reference_quality"),
            "source_families": (
                row.get("classification") or {}
            ).get("source_families", []),
            "present_properties": (
                row.get("classification") or {}
            ).get("present_properties", []),
            "reference_packet_file": "empirical-reference-properties.json",
            "reference_packet_key": str(row["candidate_id"]),
            "reference_mode": "multi_source_empirical_property_packet",
            "entry_reference_available": bool(
                row.get("entry_reference_ready")
            ),
        }
        for row in ready
    ]
    return {
        "version": "v5e-2.1-v5b3-property-calibration-input",
        "executable": len(anchors) >= minimum,
        "minimum_ready_targets": minimum,
        "ready_target_count": len(anchors),
        "reference_mode": "multi_source_empirical_property_packets",
        "entry_level_reference_is_required_for_anchor": False,
        "important": (
            "These anchors are calibration-grade property packets, not complete "
            "entry-trajectory ground truth. V5b-3 must use them to expand and "
            "validate the property/applicability model while retaining the "
            "original CoolProp-backed entry holdouts for entry-level error "
            "calibration. Predicted/GC/CSP sources are excluded by policy."
        ),
        "anchors": anchors,
    }
