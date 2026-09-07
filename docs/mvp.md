# MVP Plan

## 1. Objective

Build the smallest system that can **compare candidate working fluids fairly**.

The MVP is not a water-cooling project. Water is used only to validate the evaluator and provide a reference score.

The primary MVP research question is:

> Can the system rank multiple fluids by the working-fluid mass required to satisfy the same thermal and mission constraints?

## 2. MVP stages

### Phase 0 — evaluator contract

Define:

- mission inputs;
- hardware/comparison assumptions;
- thermodynamic boundaries;
- heat-transfer convention;
- flow/nozzle convention;
- validity limits;
- numerical verification gates;
- provenance schema.

Use the [generic fluid evaluator contract](benchmarks/fluid-evaluator-contract.md).

### Phase 1 — evaluator verification with water

Use water because its properties are well characterized.

Goals:

- verify thermophysical-property access;
- verify mass and energy conservation;
- verify phase-change handling;
- verify numerical convergence;
- produce a reference mass score.

This phase produces:

~~~text
M_water_reference
~~~

It does **not** define the project goal as optimizing a water-cooled heat shield.

### Phase 2 — known-fluid screening

Immediately evaluate multiple known candidate fluids using the same geometry, mission and control assumptions.

~~~text
water
fluid A
fluid B
fluid C
...
    |
    v
same evaluator
    |
    v
required mass
~~~

Primary comparison:

~~~text
mass_ratio =
    M_candidate_required /
    M_water_reference
~~~

The important milestone is not merely reproducing water. It is demonstrating that the same evaluator can rank non-water fluids consistently.

### Phase 3 — mixtures

Add known or deliberately constructed mixtures.

Requirements:

- explicit composition basis;
- validated mixture property model;
- no naive averaging of pure-fluid phase properties;
- uncertainty recorded.

Search mixture ratios with conventional numerical optimization before adding an LLM.

### Phase 4 — property prediction and molecular simulation

For candidates without reliable tabulated properties, add:

- molecular dynamics;
- ab-initio MD where justified;
- learned property models;
- uncertainty estimation;
- chemical/decomposition analysis.

The goal is to fill missing inputs for the same evaluator, not to replace it.

### Phase 5 — AI-guided candidate selection

The AI controller may now:

- choose chemical families;
- propose mixture regions;
- choose candidates for expensive MD / ab-initio / CFD;
- identify uncertain properties that dominate the score;
- balance exploration and exploitation.

The AI must be compared against simpler acquisition/search strategies.

### Phase 6 — coupled jet / aerodynamic optimization

Once the fluid-ranking path is stable, add increasingly realistic models of:

- directed vapor discharge;
- nozzle reaction force;
- injection-dependent drag;
- injection-dependent heat flux;
- coupled trajectory dynamics.

Only then evaluate the full benefit of braking and cooling together.

### Phase 7 — joint system optimization

Finally allow simultaneous optimization of:

- fluid chemistry;
- mixture ratio;
- mass-flow schedule;
- pressure schedule;
- nozzle geometry;
- nozzle distribution;
- passive/active TPS split.

At this stage use total system mass rather than fluid mass alone.

## 3. Fair comparison hierarchy

To avoid confusing fluid quality with hardware quality:

~~~text
Experiment A
fixed hardware
fixed control law
compare fluids

Experiment B
fixed hardware bounds
optimize control per fluid

Experiment C
co-optimize fluid + hardware + control
include hardware mass
~~~

A fluid should not receive credit for a better nozzle geometry that was never offered to the water reference.

## 4. Minimal software components

A pragmatic first implementation:

~~~text
src/
  orchestrator/
  mission/
  fluids/
  properties/
  thermal/
  evaluator/
  optimization/
  storage/

configs/
  missions/
  fluids/
  evaluator/

benchmarks/
  references/
    water/

tests/
  conservation/
  properties/
  evaluator/
~~~

Initial Python ecosystem:

- Python 3.12;
- NumPy / SciPy;
- validated thermophysical-property backend;
- Pydantic or equivalent typed schemas;
- SQLite + artifact files;
- plotting/reporting tools.

No DFT engine, crystal generator, GPU stack or LLM is required to compare the first set of known fluids.

## 5. Minimal data model

~~~text
Objective
Mission
FluidCandidate
FluidPropertyModel
HardwareConfiguration
ControlPolicy
SimulationRequest
SimulationAttempt
Observation
Comparison
Decision
BudgetLedger
~~~

The crystal-discovery track may later add Structure, DFTCalculation and related entities.

## 6. First acceptance test

The first useful MVP passes when:

1. a fixed mission/comparison case loads with explicit units and constraints;
2. water property calculations pass independent reference checks;
3. mass and energy balances pass numerical verification;
4. the evaluator returns a reproducible water reference result;
5. at least **three non-water fluids** run through the same evaluator;
6. the system produces a ranked comparison with constraint margins;
7. failures and unsupported states cannot appear as feasible solutions;
8. all model versions, assumptions and outputs are stored with provenance.

This deliberately prevents “water benchmark completed” from being mistaken for the project milestone.

## 7. First scientific output

A first meaningful result should look like:

~~~text
Fixed mission / hardware / control assumptions

Water       required mass: 1.000 reference
Fluid A     mass ratio:    0.91
Fluid B     mass ratio:    1.14
Fluid C     mass ratio:    0.78
~~~

with:

- thermal constraint margins;
- property-model validity;
- uncertainty;
- model fidelity;
- explicit statement of whether braking effects are included.

## 8. First GPU milestone

The 2 × RTX 4090 become useful when the project begins to evaluate candidates whose properties are not already available.

Possible workloads:

- molecular dynamics;
- ML property prediction;
- uncertainty ensembles;
- learned interatomic potentials;
- candidate-generation models;
- later CFD acceleration where supported.

The project should not invent a GPU workload merely to use the hardware.

## 9. AI acceptance test

The AI layer passes only if, at equal physical-evaluation budget, it measurably improves discovery.

Examples:

- finds a lower-mass feasible fluid sooner;
- identifies better mixture regions;
- reduces expensive property calculations;
- avoids redundant high-fidelity evaluations;
- explores high-uncertainty chemical regions without collapsing onto one family.

Plausible prose is not a success metric.

## 10. Parallel crystal track

Crystal discovery remains a valid second domain:

~~~text
candidate structures
    |
    v
ML screening
    |
    v
DFT
    |
    v
phase / phonon validation
~~~

It shares orchestration and AI-selection infrastructure with the fluid track, but it is not a prerequisite for the current working-fluid problem.
