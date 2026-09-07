# MVP Plan

## 1. Objective

Build the smallest system that proves the closed-loop architecture works.

The MVP should **not** begin by searching arbitrary unknown chemistry. It should first reproduce a controlled benchmark where a good answer is already known.

This makes it possible to measure whether the AI-guided policy adds value.

## 2. MVP research question

Suggested first task:

> Within a deliberately bounded family of inorganic crystals, can the system recover low-energy structures while using fewer DFT calculations than naive screening?

The exact chemistry can be selected later.

Important properties of the benchmark:

- small unit cells;
- inexpensive elements;
- known reference structures;
- manageable Quantum ESPRESSO calculations;
- compatible with an available ML potential;
- enough alternative structures to make selection non-trivial.

## 3. Phase 0 — repository design

Deliverables:

- architecture;
- candidate lifecycle;
- AI-controller contract;
- compute strategy;
- experiment schema;
- benchmark definition.

Status: current phase.

## 4. Phase 1 — deterministic pipeline

No AI controller yet.

Implement:

```text
structure input
    |
    v
normalize
    |
    v
deduplicate
    |
    v
ML relax
    |
    v
rank
    |
    v
DFT selected candidates
    |
    v
database
```

Goal: prove all tools are callable and results are reproducible.

## 5. Phase 2 — baseline search policies

Implement several simple selectors:

- random;
- top-N predicted energy;
- uncertainty sampling;
- diversity sampling;
- fixed mixture policy.

These become the baselines the AI must beat.

## 6. Phase 3 — AI controller

Add the reasoning model above the same tool interface.

The AI controller receives only structured summaries and can choose:

- search subspace;
- generation batch;
- promotion set;
- retry actions;
- budget allocation.

Every action is validated by the deterministic orchestrator.

## 7. Phase 4 — active learning

After enough DFT data is collected:

```text
DFT observations
      |
      v
surrogate update
      |
      v
recalculate uncertainty
      |
      v
new acquisition
      |
      v
next DFT batch
```

Measure whether the updated surrogate reduces expensive reference calculations.

## 8. Phase 5 — genuinely novel candidates

Only after benchmark success:

- expand composition space;
- introduce stronger novelty pressure;
- compare against known-material databases;
- perform convex-hull analysis;
- run phonon validation;
- escalate the best candidates to HPC.

## 9. Minimal software components

A pragmatic first implementation could contain:

```text
src/
  orchestrator/
  candidates/
  screening/
  dft/
  controller/
  storage/

configs/
  objectives/
  dft-presets/
  models/

data/
  .gitkeep

tests/
```

Likely Python ecosystem:

- ASE;
- pymatgen;
- PyTorch;
- one ML interatomic potential;
- Quantum ESPRESSO;
- SQLite initially;
- Pydantic for typed schemas.

AiiDA can be introduced when workflow/provenance complexity justifies it rather than on day one.

## 10. Minimal data model

First entities:

```text
Objective
Candidate
Structure
Prediction
Calculation
Observation
Decision
Hypothesis
BudgetLedger
```

This is sufficient to reconstruct why any DFT job was launched.

## 11. First acceptance test

The MVP passes when:

1. at least 100 candidate structures are ingested/generated;
2. all are normalized and deduplicated;
3. all valid candidates receive an ML score;
4. at least two selection policies choose DFT subsets;
5. Quantum ESPRESSO calculations complete for selected candidates;
6. results are parsed into the database;
7. the full lineage from objective to final observation is queryable;
8. rerunning the experiment from the same configuration is reproducible.

## 12. AI acceptance test

The AI layer should not be considered useful merely because it produces plausible scientific prose.

It passes only if, across repeated benchmark runs, it demonstrates measurable benefit such as:

- better candidate found for the same DFT budget;
- equivalent candidate found with fewer DFT jobs;
- fewer redundant calculations;
- better recovery from failed calculations;
- broader useful exploration without excessive compute cost.

## 13. First practical milestone on 2 × RTX 4090

A reasonable workstation demonstration is:

```text
1000 candidate structures
        |
        v
ML relax / score on 2 GPUs
        |
        v
cluster + uncertainty analysis
        |
        v
select 20-50 candidates
        |
        v
small DFT validation batch
        |
        v
compare selection policies
```

The exact numbers should be adjusted to atom count and DFT cost.

## 14. What comes after the MVP

Once the benchmark is trustworthy:

- integrate a crystal generator;
- add active-learning surrogate updates;
- add external materials databases;
- add phase-diagram construction;
- add phonon workflows;
- add HPC submission;
- add an AI hypothesis ledger;
- add finite-temperature validation;
- eventually test inverse-design objectives.

The system should grow by adding new tools beneath the same controller interface, not by making the controller itself increasingly unconstrained.
