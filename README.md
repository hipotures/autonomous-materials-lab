# Autonomous Materials Lab

An experimental architecture for AI-guided discovery of hypothetical materials.

The project is intentionally starting as a **design repository**, not as a software implementation. The first goal is to define a rigorous closed-loop workflow in which an AI scientific controller decides where to spend computation, while physical models remain the source of numerical truth.

## First milestone

The immediate target is **W1-R: a CPU-only water/wall thermal benchmark with prescribed heating**, using a validated property backend, NumPy/SciPy and lightweight provenance. No crystal generator, DFT engine, GPU or LLM is required to begin. Microjets and coupled trajectories are later, separately validated extensions.

Read the [design review](docs/design-review.md) and [physical/numerical contract](docs/benchmarks/model-contract.md) before implementing the evaluator. No water-mass optimum has yet been computed.

## Core idea — later crystal discovery track

```text
Research objective
      |
      v
AI scientific controller
      |
      v
Chemical-space selection
      |
      v
Candidate structure generation
      |
      v
Deterministic filters
      |
      v
ML potential screening
      |
      v
Candidate ranking + uncertainty
      |
      +-------------------+
      |                   |
      v                   v
   exploit             explore
      |                   |
      +---------+---------+
                |
                v
          DFT validation
                |
                v
          results database
                |
                v
       model / policy update
                |
                +----> next iteration
```

The AI controller does **not** replace density-functional theory, atomistic simulation, or uncertainty-aware surrogate models. Its job is to choose the next useful computation, interpret failures, allocate a finite compute budget, and decide when to broaden or narrow the search.

## Initial target hardware

The first prototype is designed around a workstation with:

- 2 × NVIDIA RTX 4090
- conventional multicore CPU
- local SSD storage
- optional access to an external university/HPC cluster for expensive validation

The local GPUs are best treated primarily as accelerators for structure generation, machine-learning interatomic potentials, surrogate models, and molecular dynamics. High-accuracy DFT can run locally for small jobs and be escalated to HPC when justified.

## Proposed stack

These are candidates, not hard dependencies. For the water MVP use NumPy/SciPy, a validated water property backend such as CoolProp, typed schemas and SQLite plus artifact files. The following stack primarily serves the later atomistic/crystal track:

- **ASE** — common atomistic workflow interface
- **pymatgen** — structures, phase diagrams, materials analysis
- **MatterGen or another crystal generator** — candidate structure generation
- **MACE / CHGNet / related ML potential** — fast relaxation and screening
- **Quantum ESPRESSO** — open-source DFT reference backend
- **CP2K** — molecular dynamics and larger-system DFT workflows
- **Phonopy** — phonon workflows
- **AiiDA or a lightweight custom orchestrator** — provenance and workflow execution
- **SQLite/PostgreSQL + object storage** — experiment metadata and large artifacts
- **LLM / reasoning model** — scientific controller and experiment planner

## Design principles

1. **Physics outranks language-model preference.**
2. **Every rejected candidate has an explicit reason.**
3. **AI recommendations are soft unless backed by deterministic constraints.**
4. **Expensive calculations are treated as a budgeted resource.**
5. **Exploration must be preserved; unfamiliar chemistry must not be discarded merely because it is unfamiliar.**
6. **Every result is reproducible from stored inputs, code versions, model versions, and calculation settings.**
7. **The controller should be able to explain why a particular computation was selected.**

## Documentation

- [Design review and software assessment](docs/design-review.md)
- [Water physical and numerical contract](docs/benchmarks/model-contract.md)
- [Architecture](docs/architecture.md)
- [Discovery loop](docs/discovery-loop.md)
- [AI controller](docs/ai-controller.md)
- [Candidate lifecycle](docs/candidate-lifecycle.md)
- [Compute strategy](docs/compute-strategy.md)
- [Development and compute environment](docs/environment.md)
- [Source development and optimization](docs/development.md)
- [Preflight checklist](docs/preflight.md)
- [Liquid transpiration / microjet TPS use case](docs/use-cases/liquid-transpiration-tps.md)
- [Water reference benchmark](docs/benchmarks/water-reference.md)
- [MVP plan](docs/mvp.md)

## Status

**Phase 0 — architecture definition, workstation preparation, and first system benchmark definition.**

No claim is made yet that the proposed pipeline can discover a synthesizable material. The first concrete system benchmark is a water-based liquid transpiration / directed-microjet thermal-protection model. The immediate objective is to reproduce known physics with a transparent baseline before searching novel fluids, mixtures, crystals, or other material classes.
