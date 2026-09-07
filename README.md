# Autonomous Materials Lab

An experimental architecture for AI-guided inverse design of materials and material-enabled systems.

The project is intentionally starting as a **design repository**, not as a software implementation. Its purpose is to build a reproducible search-and-validation system in which AI chooses what to evaluate next while numerical physics remains the source of quantitative evidence.

## Current primary use case

The first concrete research track is **atmospheric-entry working-fluid discovery**.

The goal is **not** to design a water-cooled heat shield.

The goal is to discover or formulate a liquid / liquid mixture that minimizes the working-fluid mass required to:

- keep the vehicle within thermal limits;
- provide any useful additional braking through directed vapor discharge;
- remain physically storable, deliverable and compatible with the system.

Water is used only as a **reference and evaluator-validation fluid** because its properties are well characterized.

Read [Project Goal](docs/project-goal.md) first.

## Core liquid-discovery loop

~~~text
mission + architecture
        |
        v
candidate fluid / mixture
        |
        v
known or predicted properties
        |
        v
fluid evaluator
        |
        +---- thermal contribution
        |
        +---- vapor / nozzle contribution
        |
        +---- aerodynamic interaction
        |
        v
required working-fluid mass
        |
        v
compare with water reference
        |
        v
AI / optimizer chooses next candidate
        |
        +--------------------> repeat
~~~

The same evaluator should accept water, known liquids, mixtures and later hypothetical molecular candidates without changing the project objective.

## Water's role

Water is a calibration point, not the desired end product.

~~~text
water
  -> verify property model
  -> verify conservation
  -> obtain reference score

candidate X
  -> same evaluator
  -> compare required mass with water
~~~

A scientifically meaningful result is of the form:

> Under the stated mission, hardware assumptions and model fidelity, candidate X requires Y% less working-fluid mass than the water reference while satisfying the same constraints.

## Other research tracks

The architecture remains general enough to support other domains, including crystal discovery.

~~~text
Autonomous Materials Lab
        |
        +-- fluid / mixture inverse design
        |     -> thermodynamics
        |     -> MD / chemistry
        |     -> aerothermal system evaluator
        |
        +-- crystal discovery
              -> structure generation
              -> ML screening
              -> DFT
              -> phonons / stability
~~~

These tracks share:

- orchestration;
- provenance;
- budget accounting;
- AI experiment selection;
- uncertainty handling;
- high-fidelity escalation.

They do not share the same physical evaluator.

## Initial target hardware

The first prototype is designed around:

- 2 × NVIDIA RTX 4090;
- multicore x86_64 CPU;
- local NVMe storage;
- optional university/HPC access.

The first reduced fluid evaluator may run almost entirely on CPU. The GPUs become important for high-throughput ML, molecular simulation, surrogate models and later expensive candidate screening.

## Proposed software stack

Domain-dependent candidates include:

### Fluid / TPS track

- NumPy / SciPy — reduced-order models and optimization;
- validated thermophysical-property backend;
- Cantera — thermodynamics / reacting-gas chemistry where applicable;
- LAMMPS — molecular dynamics;
- CP2K — selected ab-initio MD / molecular reference calculations;
- SU2 or another open CFD solver — later aerothermal interaction modelling;
- custom trajectory / thermal / injection evaluator.

### Crystal track

- ASE;
- pymatgen;
- MatterGen or another generator;
- MACE / CHGNet;
- Quantum ESPRESSO;
- CP2K;
- Phonopy.

### Shared

- Python 3.12 production baseline;
- SQLite/PostgreSQL + artifact storage;
- Pydantic or equivalent typed schemas;
- LLM / reasoning model as scientific controller;
- optional AiiDA when remote-workflow complexity justifies it.

## Design principles

1. **The search target must be explicit.**
2. **Water is a benchmark, not the design objective.**
3. **Physics outranks language-model preference.**
4. **Every rejected candidate has an explicit reason.**
5. **A numerical failure is not a physical rejection.**
6. **Expensive calculations are treated as budgeted resources.**
7. **Every result is reproducible from stored inputs, versions and settings.**
8. **The controller must explain why a calculation was selected.**
9. **Candidate comparison must use the same evaluator and comparable system assumptions.**

## Documentation

- [Project goal](docs/project-goal.md)
- [Architecture](docs/architecture.md)
- [Liquid TPS / working-fluid use case](docs/use-cases/liquid-transpiration-tps.md)
- [Generic fluid evaluator contract](docs/benchmarks/fluid-evaluator-contract.md)
- [Water reference benchmark](docs/benchmarks/water-reference.md)
- [MVP plan](docs/mvp.md)
- [Discovery loop](docs/discovery-loop.md)
- [AI controller](docs/ai-controller.md)
- [Candidate lifecycle](docs/candidate-lifecycle.md)
- [Compute strategy](docs/compute-strategy.md)
- [Development and compute environment](docs/environment.md)
- [Source development and optimization](docs/development.md)
- [Preflight checklist](docs/preflight.md)

## First executable experiment

Start here:

- [Ambient-liquid heat-sink screening](experiments/fluid-heat-sink/README.md)
- [Cryogenic / storage-state enthalpy screening](experiments/fluid-enthalpy-window/README.md)
- [Reactive fluid + air equilibrium check](experiments/reactive-fluid-air/README.md)
- [Hydrogen / air ignition-delay sweep](experiments/hydrogen-ignition-delay/README.md)
- [Low-fidelity end-to-end Earth entry evaluator](experiments/entry-evaluator/README.md)

The first experiments isolate thermodynamics and chemical kinetics. The entry evaluator is the first system-level calculation: it propagates a configurable Earth-entry trajectory, evaluates atmosphere, stagnation-point convective and radiative heating, wall thermal response, coolant enthalpy and required coolant mass flow. Independent cases can be distributed across CPU processes. It remains intentionally below CFD / DSMC fidelity and does not yet model detailed boundary-layer blowing, porous flow or nozzle geometry.

## Status

**Phase 1 — reduced-order end-to-end entry evaluator and working-fluid comparison.**

The repository now contains local property / chemistry screens and a trajectory-level evaluator. The staged entry work includes [V1 numerical verification](experiments/entry-evaluator/V1.md), [V2 physical-sensitivity screening](experiments/entry-evaluator/V2.md), [V3 cross-model Earth aeroheating comparison](experiments/entry-evaluator/V3.md), and [V4 ignition-delay constrained hydrogen cooling](experiments/entry-evaluator/V4.md). V4 is the final planned calibration stage before V5 starts explicit liquid / mixture candidate modeling. None of these reduced-order stages establishes flight-level physical validity.
