"""Deterministic V5c-1 molecular candidate generator.

This is intentionally a constrained template generator, not an unrestricted
molecular-design model. Every generated structure has explicit provenance.
"""
from __future__ import annotations

from dataclasses import dataclass, asdict
@dataclass(frozen=True)
class GeneratedCandidate:
    smiles: str
    family: str
    operation: str
    parameters: dict[str, int | str]

    def to_dict(self) -> dict:
        return asdict(self)


def _linear_alkane(carbon_count: int) -> str:
    return "C" * carbon_count


def _methyl_branched_alkane(total_carbons: int, branch_position: int) -> str:
    main_chain = total_carbons - 1
    if main_chain < 3:
        raise ValueError("branched alkane requires at least 4 total carbons")
    if not 2 <= branch_position <= main_chain - 1:
        raise ValueError("branch position must be internal")
    return (
        "C" * branch_position
        + "(C)"
        + "C" * (main_chain - branch_position)
    )


def _linear_alkene(carbon_count: int, bond_position: int) -> str:
    if carbon_count < 3:
        raise ValueError("alkene requires at least 3 carbons")
    if not 1 <= bond_position < carbon_count:
        raise ValueError("invalid alkene bond position")
    left = "C" * bond_position
    right = "C" * (carbon_count - bond_position)
    return f"{left}={right}"


def _primary_alcohol(carbon_count: int) -> str:
    if carbon_count < 1:
        raise ValueError("alcohol requires at least one carbon")
    return "C" * carbon_count + "O"


def _linear_ketone(carbon_count: int, carbonyl_position: int) -> str:
    if carbon_count < 3:
        raise ValueError("ketone requires at least 3 carbons")
    if not 2 <= carbonyl_position <= carbon_count - 1:
        raise ValueError("carbonyl must be internal")
    left = "C" * (carbonyl_position - 1)
    right = "C" * (carbon_count - carbonyl_position)
    return f"{left}C(=O){right}"


def _linear_ether(left_carbons: int, right_carbons: int) -> str:
    if left_carbons < 1 or right_carbons < 1:
        raise ValueError("ether sides must contain carbon")
    return "C" * left_carbons + "O" + "C" * right_carbons


def generate_candidates(config: dict) -> list[GeneratedCandidate]:
    """Generate deterministic template candidates from a V5c config."""
    limits = config["generation"]
    candidates: list[GeneratedCandidate] = []

    for carbon_count in range(
        int(limits["alkanes"]["min_c"]),
        int(limits["alkanes"]["max_c"]) + 1,
    ):
        candidates.append(
            GeneratedCandidate(
                smiles=_linear_alkane(carbon_count),
                family="alkanes",
                operation="linear_chain_length",
                parameters={"carbon_count": carbon_count},
            )
        )

    for total_carbons in range(
        int(limits["branched_alkanes"]["min_total_c"]),
        int(limits["branched_alkanes"]["max_total_c"]) + 1,
    ):
        main_chain = total_carbons - 1
        for position in range(2, (main_chain + 1) // 2 + 1):
            candidates.append(
                GeneratedCandidate(
                    smiles=_methyl_branched_alkane(
                        total_carbons,
                        position,
                    ),
                    family="alkanes",
                    operation="single_methyl_branch",
                    parameters={
                        "total_carbons": total_carbons,
                        "branch_position": position,
                    },
                )
            )

    for carbon_count in range(
        int(limits["alkenes"]["min_c"]),
        int(limits["alkenes"]["max_c"]) + 1,
    ):
        position = 1
        candidates.append(
            GeneratedCandidate(
                smiles=_linear_alkene(carbon_count, position),
                family="alkenes",
                operation="terminal_alkene_chain_length",
                parameters={
                    "carbon_count": carbon_count,
                    "bond_position": position,
                },
            )
        )

    for carbon_count in range(
        int(limits["primary_alcohols"]["min_c"]),
        int(limits["primary_alcohols"]["max_c"]) + 1,
    ):
        candidates.append(
            GeneratedCandidate(
                smiles=_primary_alcohol(carbon_count),
                family="oxygenated",
                operation="primary_alcohol_chain_length",
                parameters={"carbon_count": carbon_count},
            )
        )

    for carbon_count in range(
        int(limits["ketones"]["min_c"]),
        int(limits["ketones"]["max_c"]) + 1,
    ):
        for position in range(
            2,
            (carbon_count + 1) // 2 + 1,
        ):
            candidates.append(
                GeneratedCandidate(
                    smiles=_linear_ketone(carbon_count, position),
                    family="oxygenated",
                    operation="ketone_chain_and_position",
                    parameters={
                        "carbon_count": carbon_count,
                        "carbonyl_position": position,
                    },
                )
            )

    max_ether_c = int(limits["ethers"]["max_total_c"])
    min_ether_c = int(limits["ethers"]["min_total_c"])
    for total_c in range(min_ether_c, max_ether_c + 1):
        for left_c in range(1, total_c // 2 + 1):
            right_c = total_c - left_c
            candidates.append(
                GeneratedCandidate(
                    smiles=_linear_ether(left_c, right_c),
                    family="oxygenated",
                    operation="ether_chain_partition",
                    parameters={
                        "left_carbons": left_c,
                        "right_carbons": right_c,
                    },
                )
            )

    for chain_c in range(
        int(limits["alkyl_benzenes"]["min_sidechain_c"]),
        int(limits["alkyl_benzenes"]["max_sidechain_c"]) + 1,
    ):
        candidates.append(
            GeneratedCandidate(
                smiles="C" * chain_c + "c1ccccc1",
                family="aromatics",
                operation="benzene_linear_alkyl_substitution",
                parameters={"sidechain_carbons": chain_c},
            )
        )

    for ring_size in range(
        int(limits["cycloalkanes"]["min_ring"]),
        int(limits["cycloalkanes"]["max_ring"]) + 1,
    ):
        ring = "C1" + "C" * (ring_size - 1) + "1"
        candidates.append(
            GeneratedCandidate(
                smiles=ring,
                family="cyclic_hydrocarbons",
                operation="cycloalkane_ring_size",
                parameters={"ring_size": ring_size},
            )
        )
        if bool(limits["cycloalkanes"].get("methyl_substitution", True)):
            candidates.append(
                GeneratedCandidate(
                    smiles="C" + ring,
                    family="cyclic_hydrocarbons",
                    operation="methyl_cycloalkane",
                    parameters={"ring_size": ring_size},
                )
            )

    return candidates
