# MVP Plan

## 1. Objective

Build the smallest system that proves the closed-loop architecture works.

The MVP should **not** begin by searching arbitrary unknown chemistry. It should first reproduce a controlled benchmark where a good answer is already known.

This makes it possible to measure whether the AI-guided policy adds value.

## 2. MVP research question

The project now has a preferred first **system benchmark**:

> Can a transparent reduced-order model determine the minimum water mass required for an atmospheric-entry thermal-protection concept using transpiration cooling and, later, directed vapor microjets?

This benchmark deliberately avoids novel-material generation at first. Water properties are known well enough that the project can validate the coupled evaluator before adding chemical discovery.

The initial benchmark progression is:

~~~text
W0: passive reference
W1: water transpiration
W2: water directed microjets
W3: optimized water system
~~~

The primary result is M_water_min: the minimum water mass satisfying the selected thermal and trajectory constraints.

After that benchmark is reproducible, the discovery loop can compare:

- known pure liquids;
- known mixtures;
- optimized mixtures;
- novel molecular candidates;
- joint fluid + nozzle + control designs.

A bounded inorganic-crystal benchmark remains a useful secondary test for the generic materials-discovery engine, especially for DFT/ML active-learning workflows.

## 3. Phase 0 — repository design

Deliverables:

- architecture;
- candidate lifecycle;
- AI-controller contract;
- compute strategy;
- experiment schema;
- benchmark definition.

Status: current phase.

## 4. Phase 1 — deterministic water benchmark

No AI controller and no novel-fluid generation yet.

Implement:

~~~text
mission definition
    |
    v
water property model
    |
    v
reduced thermal model
    |
    v
transpiration / nozzle model
    |
    v
vehicle force and trajectory model
    |
    v
minimum-water-mass search
    |
    v
provenance database
~~~

Goal: establish M_water_min and show that the complete calculation is reproducible.

The first implementation may use prescribed heat-load and drag histories before moving to fully coupled aerothermodynamics.

## 5. Phase 2 — known-fluid comparison

Add a small library of known liquids and mixtures using trusted thermodynamic data.

For every fluid, optimize the same system variables under the same mission assumptions and compare:

~~~text
mass_ratio = M_candidate_min / M_water_min
~~~

This prevents the project from attributing gains from nozzle geometry or control policy to fluid chemistry.

## 6. Phase 3 — optimization and AI controller

Add the reasoning model above the same tool interface.

The AI controller receives only structured summaries and can choose:

- fluid or mixture subspace;
- geometry/control search region;
- high-fidelity promotion set;
- retry actions;
- budget allocation.

Every action is validated by the deterministic orchestrator.

## 7. Phase 4 — active learning

After enough high-fidelity CFD, MD, chemistry, or DFT data is collected:

```text
high-fidelity observations
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
next high-fidelity batch
```

Measure whether the updated surrogate reduces expensive reference calculations while preserving system-level accuracy.

## 8. Phase 5 — genuinely novel candidates

Only after the water and known-fluid benchmarks are trustworthy:

- generate novel molecular candidates or mixtures;
- predict missing thermodynamic/transport properties;
- validate selected candidates with MD / ab-initio MD / chemistry calculations;
- introduce CFD for the most promising designs;
- escalate expensive coupled simulations to HPC.

The crystal-discovery track can independently add composition generation, convex-hull analysis, phonons and DFT validation.

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
