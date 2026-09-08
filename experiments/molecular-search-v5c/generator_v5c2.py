"""Expanded deterministic molecular generator for V5c-2."""
from __future__ import annotations

from generator import GeneratedCandidate, generate_candidates as generate_v5c1


def _linear_aldehyde(carbon_count: int) -> str:
    if carbon_count < 2:
        raise ValueError("aldehyde requires at least 2 carbons")
    return "C" * (carbon_count - 1) + "C=O"


def _terminal_alkyne(carbon_count: int) -> str:
    if carbon_count < 3:
        raise ValueError("terminal alkyne requires at least 3 carbons")
    return "C" * (carbon_count - 2) + "C#C"


def _ester(acyl_carbons: int, alkoxy_carbons: int) -> str:
    if acyl_carbons < 1 or alkoxy_carbons < 1:
        raise ValueError("ester fragments require at least one carbon")
    acyl = (
        "C(=O)O"
        if acyl_carbons == 1
        else "C" * (acyl_carbons - 1) + "C(=O)O"
    )
    return acyl + "C" * alkoxy_carbons


def _branched_primary_alcohol(
    total_carbons: int,
    branch_position: int,
) -> str:
    main_chain = total_carbons - 1
    if main_chain < 3:
        raise ValueError(
            "branched primary alcohol requires at least 4 carbons"
        )
    if not 2 <= branch_position <= main_chain - 1:
        raise ValueError("branch must leave terminal CH2OH intact")
    return (
        "C" * branch_position
        + "(C)"
        + "C" * (main_chain - branch_position)
        + "O"
    )


def _primary_amine(carbon_count: int) -> str:
    if carbon_count < 2:
        raise ValueError(
            "pinned NH2 SMARTS excludes methylamine; use C2+"
        )
    return "C" * carbon_count + "N"


def generate_candidates(config: dict) -> list[GeneratedCandidate]:
    """Return the V5c-1 space plus new parameter-supported motifs."""
    candidates = list(generate_v5c1(config))
    limits = config["generation"]

    aldehydes = limits["aldehydes"]
    for carbon_count in range(
        int(aldehydes["min_c"]),
        int(aldehydes["max_c"]) + 1,
    ):
        candidates.append(
            GeneratedCandidate(
                smiles=_linear_aldehyde(carbon_count),
                family="aldehydes",
                operation="linear_aldehyde_chain_length",
                parameters={"carbon_count": carbon_count},
            )
        )

    alkynes = limits["terminal_alkynes"]
    for carbon_count in range(
        int(alkynes["min_c"]),
        int(alkynes["max_c"]) + 1,
    ):
        candidates.append(
            GeneratedCandidate(
                smiles=_terminal_alkyne(carbon_count),
                family="alkynes",
                operation="terminal_alkyne_chain_length",
                parameters={"carbon_count": carbon_count},
            )
        )

    esters = limits["esters"]
    min_total = int(esters["min_total_c"])
    max_total = int(esters["max_total_c"])
    max_acyl = int(esters.get("max_acyl_c", max_total - 1))
    for total_carbons in range(min_total, max_total + 1):
        for acyl_carbons in range(
            1,
            min(max_acyl, total_carbons - 1) + 1,
        ):
            alkoxy_carbons = total_carbons - acyl_carbons
            candidates.append(
                GeneratedCandidate(
                    smiles=_ester(
                        acyl_carbons,
                        alkoxy_carbons,
                    ),
                    family="esters",
                    operation="ester_acyl_alkoxy_partition",
                    parameters={
                        "total_carbons": total_carbons,
                        "acyl_carbons": acyl_carbons,
                        "alkoxy_carbons": alkoxy_carbons,
                    },
                )
            )

    alcohols = limits["branched_primary_alcohols"]
    for total_carbons in range(
        int(alcohols["min_total_c"]),
        int(alcohols["max_total_c"]) + 1,
    ):
        main_chain = total_carbons - 1
        for branch_position in range(2, main_chain):
            candidates.append(
                GeneratedCandidate(
                    smiles=_branched_primary_alcohol(
                        total_carbons,
                        branch_position,
                    ),
                    family="oxygenated",
                    operation="branched_primary_alcohol",
                    parameters={
                        "total_carbons": total_carbons,
                        "branch_position": branch_position,
                    },
                )
            )

    amines = limits["primary_amines"]
    for carbon_count in range(
        int(amines["min_c"]),
        int(amines["max_c"]) + 1,
    ):
        candidates.append(
            GeneratedCandidate(
                smiles=_primary_amine(carbon_count),
                family="amines",
                operation="primary_amine_chain_length",
                parameters={"carbon_count": carbon_count},
            )
        )

    return candidates
