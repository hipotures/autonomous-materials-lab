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

The first experiment compares fluids that are liquid at a common ambient state. The second gives every candidate its own initial `T0/P0`, allowing cryogenic liquids such as methane, oxygen, nitrogen and hydrogen to enter the comparison. The third uses Cantera equilibrium to estimate the chemical oxidation potential of hot H2 / CH4 / NH3 after mixing with air. The fourth adds finite-rate H2 / air chemistry and sweeps homogeneous ignition delay versus temperature, pressure and equivalence ratio. None of these experiments includes CFD, trajectory or nozzle geometry.

## Status

**Phase 0 — scientific objective definition, evaluator contract, workstation preparation, and first thermodynamic screening experiment.**

The next implementation should validate a generic fluid evaluator with a known reference fluid, then immediately begin comparing multiple fluids. Water verification is a prerequisite for confidence in the evaluator, not the scientific endpoint.
