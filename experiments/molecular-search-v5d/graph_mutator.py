"""Deterministic graph mutations for V5d-1 high-throughput search."""
from __future__ import annotations

from dataclasses import dataclass, asdict
import hashlib
import random
from typing import Any


@dataclass(frozen=True)
class Mutation:
    smiles: str
    parent_smiles: str
    operator: str
    generation: int

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)


def _canonical_mol(mol):
    from rdkit import Chem

    try:
        Chem.SanitizeMol(mol)
    except Exception:
        return None
    smiles = Chem.MolToSmiles(mol, canonical=True)
    if "." in smiles:
        return None
    check = Chem.MolFromSmiles(smiles)
    if check is None:
        return None
    return check


def _canonical_smiles(mol) -> str | None:
    from rdkit import Chem

    clean = _canonical_mol(mol)
    if clean is None:
        return None
    return Chem.MolToSmiles(clean, canonical=True)


def _seed_int(seed: str, parent: str, generation: int, index: int) -> int:
    digest = hashlib.sha256(
        f"{seed}|{generation}|{index}|{parent}".encode("utf-8")
    ).digest()
    return int.from_bytes(digest[:8], "big", signed=False)


def _available_nonaromatic_bonds(mol) -> list[int]:
    return [
        bond.GetIdx()
        for bond in mol.GetBonds()
        if not bond.GetIsAromatic()
    ]


def _mutate_add_terminal(mol, rng: random.Random, atomic_numbers: list[int]):
    from rdkit import Chem

    if mol.GetNumAtoms() == 0:
        return None
    rw = Chem.RWMol(mol)
    parent_idx = rng.randrange(rw.GetNumAtoms())
    new_idx = rw.AddAtom(Chem.Atom(rng.choice(atomic_numbers)))
    rw.AddBond(parent_idx, new_idx, Chem.BondType.SINGLE)
    return rw.GetMol()


def _mutate_delete_terminal(mol, rng: random.Random):
    from rdkit import Chem

    candidates = [
        atom.GetIdx()
        for atom in mol.GetAtoms()
        if atom.GetDegree() == 1 and not atom.IsInRing()
    ]
    if not candidates or mol.GetNumAtoms() <= 2:
        return None
    rw = Chem.RWMol(mol)
    rw.RemoveAtom(rng.choice(candidates))
    return rw.GetMol()


def _mutate_change_atom(mol, rng: random.Random, atomic_numbers: list[int]):
    from rdkit import Chem

    candidates = [
        atom.GetIdx()
        for atom in mol.GetAtoms()
        if not atom.GetIsAromatic()
    ]
    if not candidates:
        return None
    rw = Chem.RWMol(mol)
    idx = rng.choice(candidates)
    current = rw.GetAtomWithIdx(idx).GetAtomicNum()
    alternatives = [z for z in atomic_numbers if z != current]
    if not alternatives:
        return None
    rw.GetAtomWithIdx(idx).SetAtomicNum(rng.choice(alternatives))
    return rw.GetMol()


def _mutate_change_bond(mol, rng: random.Random):
    from rdkit import Chem

    bond_indices = _available_nonaromatic_bonds(mol)
    if not bond_indices:
        return None
    rw = Chem.RWMol(mol)
    bond = rw.GetBondWithIdx(rng.choice(bond_indices))
    current = bond.GetBondType()
    options = [
        Chem.BondType.SINGLE,
        Chem.BondType.DOUBLE,
        Chem.BondType.TRIPLE,
    ]
    options = [option for option in options if option != current]
    bond.SetBondType(rng.choice(options))
    return rw.GetMol()


def _mutate_insert_atom(mol, rng: random.Random, atomic_numbers: list[int]):
    from rdkit import Chem

    bond_indices = _available_nonaromatic_bonds(mol)
    if not bond_indices:
        return None
    source = mol.GetBondWithIdx(rng.choice(bond_indices))
    begin = source.GetBeginAtomIdx()
    end = source.GetEndAtomIdx()

    rw = Chem.RWMol(mol)
    rw.RemoveBond(begin, end)
    new_idx = rw.AddAtom(Chem.Atom(rng.choice(atomic_numbers)))
    rw.AddBond(begin, new_idx, Chem.BondType.SINGLE)
    rw.AddBond(new_idx, end, Chem.BondType.SINGLE)
    return rw.GetMol()


def _mutate_close_ring(mol, rng: random.Random):
    from rdkit import Chem

    if mol.GetNumAtoms() < 5:
        return None
    candidates: list[tuple[int, int]] = []
    for left in range(mol.GetNumAtoms()):
        for right in range(left + 1, mol.GetNumAtoms()):
            if mol.GetBondBetweenAtoms(left, right) is not None:
                continue
            try:
                path = Chem.GetShortestPath(mol, left, right)
            except Exception:
                continue
            ring_size = len(path)
            if ring_size in {5, 6}:
                candidates.append((left, right))
    if not candidates:
        return None

    left, right = rng.choice(candidates)
    rw = Chem.RWMol(mol)
    rw.AddBond(left, right, Chem.BondType.SINGLE)
    return rw.GetMol()


def _mutate_remove_ring_bond(mol, rng: random.Random):
    from rdkit import Chem

    candidates = [
        bond.GetIdx()
        for bond in mol.GetBonds()
        if bond.IsInRing() and not bond.GetIsAromatic()
    ]
    if not candidates:
        return None
    rw = Chem.RWMol(mol)
    bond = rw.GetBondWithIdx(rng.choice(candidates))
    rw.RemoveBond(bond.GetBeginAtomIdx(), bond.GetEndAtomIdx())
    return rw.GetMol()


_OPERATORS = {
    "add_terminal": _mutate_add_terminal,
    "delete_terminal": _mutate_delete_terminal,
    "change_atom": _mutate_change_atom,
    "change_bond": _mutate_change_bond,
    "insert_atom": _mutate_insert_atom,
    "close_ring_5_6": _mutate_close_ring,
    "remove_ring_bond": _mutate_remove_ring_bond,
}


def mutate_once(
    parent_smiles: str,
    *,
    operator: str,
    seed: int,
    allowed_atomic_numbers: list[int],
) -> str | None:
    from rdkit import Chem

    mol = Chem.MolFromSmiles(parent_smiles)
    if mol is None:
        return None
    function = _OPERATORS[operator]
    rng = random.Random(int(seed))

    if operator in {"add_terminal", "change_atom", "insert_atom"}:
        candidate = function(mol, rng, allowed_atomic_numbers)
    else:
        candidate = function(mol, rng)
    if candidate is None:
        return None
    return _canonical_smiles(candidate)


def generate_mutations(
    parent_smiles: str,
    *,
    generation: int,
    count: int,
    seed: str,
    allowed_atomic_numbers: list[int],
    operators: list[str],
) -> list[Mutation]:
    if count <= 0:
        return []
    unknown = sorted(set(operators) - set(_OPERATORS))
    if unknown:
        raise ValueError(f"unknown mutation operators: {unknown}")

    output: list[Mutation] = []
    for index in range(count):
        local_seed = _seed_int(seed, parent_smiles, generation, index)
        rng = random.Random(local_seed)
        operator = rng.choice(operators)
        smiles = mutate_once(
            parent_smiles,
            operator=operator,
            seed=local_seed,
            allowed_atomic_numbers=allowed_atomic_numbers,
        )
        if smiles is None or smiles == parent_smiles:
            continue
        output.append(
            Mutation(
                smiles=smiles,
                parent_smiles=parent_smiles,
                operator=operator,
                generation=generation,
            )
        )
    return output


def classify_family(smiles: str) -> str:
    """Assign a conservative family label used for search-policy gating."""
    from rdkit import Chem

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"invalid SMILES: {smiles}")

    patterns = {
        "amine": Chem.MolFromSmarts("[NX3H2][#6]"),
        "ester": Chem.MolFromSmarts("[CX3](=O)[OX2][#6]"),
        "aldehyde": Chem.MolFromSmarts("[CX3H1](=O)[#6]"),
        "ketone": Chem.MolFromSmarts("[CX3](=O)([#6])[#6]"),
        "alcohol": Chem.MolFromSmarts("[OX2H][#6]"),
        "ether": Chem.MolFromSmarts("[OD2]([#6])[#6]"),
        "alkyne": Chem.MolFromSmarts("[#6]#[#6]"),
        "alkene": Chem.MolFromSmarts("[#6]=[#6]"),
    }
    hits = {
        name
        for name, pattern in patterns.items()
        if pattern is not None and mol.HasSubstructMatch(pattern)
    }
    # The generic ether pattern also sees the single-bond ester oxygen.
    # Ester is the more specific functional-group assignment.
    if "ester" in hits:
        hits.discard("ether")

    functional_hits = hits & {
        "amine",
        "ester",
        "aldehyde",
        "ketone",
        "alcohol",
        "ether",
        "alkyne",
    }
    has_aromatic = any(atom.GetIsAromatic() for atom in mol.GetAtoms())
    has_ring = mol.GetRingInfo().NumRings() > 0
    if len(functional_hits) > 1:
        return "mixed_functional"
    if functional_hits and (has_aromatic or has_ring or "alkene" in hits):
        return "mixed_functional"
    if "alkene" in hits and (has_aromatic or has_ring):
        return "mixed_functional"

    if "amine" in hits:
        return "amines"
    if "ester" in hits:
        return "esters"
    if "aldehyde" in hits:
        return "aldehydes"
    if "alkyne" in hits:
        return "alkynes"
    if functional_hits & {"ketone", "alcohol", "ether"}:
        return "oxygenated"
    if has_aromatic:
        return "aromatics"
    if has_ring:
        return "cyclic_hydrocarbons"
    if "alkene" in hits:
        return "alkenes"
    if all(atom.GetAtomicNum() == 6 for atom in mol.GetAtoms()):
        return "alkanes"
    return "mixed_functional"


def structural_bucket(smiles: str) -> tuple[str, int, int, int]:
    from rdkit import Chem

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"invalid SMILES: {smiles}")
    heavy = mol.GetNumHeavyAtoms()
    hetero = sum(
        atom.GetAtomicNum() not in {1, 6}
        for atom in mol.GetAtoms()
    )
    rings = mol.GetRingInfo().NumRings()
    return classify_family(smiles), heavy, hetero, rings
