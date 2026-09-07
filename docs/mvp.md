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

A pragmatic first implementation for the water benchmark could contain:

~~~text
src/
  orchestrator/
  mission/
  fluids/
  thermal/
  nozzle/
  trajectory/
  optimization/
  controller/
  storage/

configs/
  missions/
  fluids/
  solver-presets/

benchmarks/
  water/

data/
  .gitkeep

tests/
~~~

Likely Python ecosystem:

- NumPy / SciPy;
- CoolProp or another validated property backend if selected after review;
- Cantera where reacting-gas chemistry becomes relevant;
- Pydantic for typed schemas;
- SQLite initially;
- plotting/reporting tools for benchmark inspection.

CFD, MD, DFT and AiiDA should be introduced only when the reduced-order benchmark is stable enough to justify higher-fidelity calculations.

## 10. Minimal data model

First entities:

~~~text
Objective
Mission
Fluid
FluidModel
Geometry
ControlPolicy
Simulation
Observation
Decision
Hypothesis
BudgetLedger
~~~

This is sufficient to reconstruct why a given design was evaluated and how its mass score was obtained.

The generic platform may later add domain-specific entities such as Structure and DFTCalculation for the crystal-discovery track.

## 11. First acceptance test

The water-benchmark MVP passes when:

1. a fixed mission definition can be loaded from configuration;
2. water properties are supplied by a versioned, testable property model;
3. the reduced thermal model computes a stable wall-temperature history;
4. the transpiration-only case W1 runs end to end;
5. the directed-microjet case W2 runs end to end;
6. an optimizer can determine a reproducible estimate of M_water_min;
7. all assumptions, solver settings and outputs are stored with provenance;
8. rerunning the benchmark from the same configuration reproduces the result within a defined numerical tolerance.

## 12. AI acceptance test

The AI layer should not be considered useful merely because it produces plausible scientific prose.

For the liquid-TPS track, it passes only if it demonstrates measurable benefit over simpler optimization/search baselines, for example:

- lower feasible system mass for the same simulation budget;
- equivalent design found with fewer high-fidelity simulations;
- better allocation between fluid, geometry and control exploration;
- fewer redundant CFD/MD evaluations;
- better recovery from failed or numerically unstable simulations.

The same principle applies to the crystal-discovery track with DFT budget replacing CFD/MD budget.

## 13. First practical milestone on 2 × RTX 4090

The first practical milestone does not need both GPUs heavily.

It should be:

~~~text
fixed entry mission
      |
      v
W0 passive reference
      |
      v
W1 water transpiration
      |
      v
W2 water microjets
      |
      v
optimize water mass / flow policy
      |
      v
produce reproducible M_water_min
~~~

The GPUs become more important in later phases when the project adds:

- surrogate models;
- molecular dynamics;
- learned fluid/property models;
- CFD acceleration where supported;
- high-throughput candidate evaluation.

This sequencing avoids using GPU compute merely because it is available.

## 14. What comes after the MVP

Once the water benchmark is trustworthy:

- compare known pure liquids;
- compare known mixtures;
- introduce system-mass accounting for tanks/manifolds/nozzles;
- add higher-fidelity aerothermal models;
- add CFD and reacting-gas chemistry;
- add molecular and ab-initio validation for selected fluids;
- add surrogate models and active learning;
- add HPC submission;
- add an AI hypothesis ledger;
- search novel mixtures and molecular candidates.

The crystal-discovery track can be added in parallel using the same controller, provenance and budget infrastructure.

The system should grow by adding new domain tools beneath the same controller interface, not by making the controller itself increasingly unconstrained.
