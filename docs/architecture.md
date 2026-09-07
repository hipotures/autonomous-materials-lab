# Architecture

## 1. Purpose

The system is a closed-loop computational materials-discovery platform. It should search hypothetical chemical and structural space, estimate which candidates deserve more expensive computation, validate selected candidates with higher-fidelity physics, and use the resulting evidence to choose the next search direction.

The architecture deliberately separates:

- **decision intelligence** — what should be tried next;
- **surrogate physics** — what can be estimated cheaply;
- **reference physics** — what should be trusted numerically;
- **provenance** — how every result can be reproduced.

## 2. Functional layers

### Layer A — Research objective

A user defines a target in machine-readable form.

Example:

```yaml
objective:
  family: ternary_nitrides
  maximize:
    - bulk_modulus
  constraints:
    density_g_cm3:
      max: 4.0
    elements_forbidden:
      - Hg
      - Cd
      - Pb
    energy_above_hull_meV_atom:
      target_max: 50
  compute_budget:
    local_gpu_hours: 200
    local_cpu_hours: 300
    hpc_jobs: 50
```

The objective should distinguish:

- hard constraints;
- soft preferences;
- target properties;
- uncertainty tolerance;
- total compute budget.

### Layer B — AI scientific controller

The controller is a reasoning model with access to structured tools. It does not directly calculate energies or forces.

Responsibilities:

- choose chemical subspaces to explore;
- allocate budget between exploration and exploitation;
- request candidate generation;
- decide which candidates advance to higher-fidelity stages;
- interpret failed calculations;
- propose retries with changed numerical settings;
- identify families that deserve deeper investigation;
- stop unproductive branches;
- produce a machine-readable rationale for every decision.

The controller must not silently convert a heuristic judgment into a physical fact.

### Layer C — Candidate generation

Possible generators include:

- composition generators;
- prototype substitution;
- symmetry-constrained random structures;
- evolutionary crystal-structure prediction;
- diffusion/generative crystal models;
- perturbations of promising structures.

A generated candidate consists minimally of:

```text
candidate_id
composition
lattice
atomic species
fractional coordinates
generation_method
generation_seed
parent_candidate_ids
generator_version
```

### Layer D — Deterministic filters

Cheap rules should remove candidates that violate explicit requirements before any expensive calculation.

Examples:

- forbidden elements;
- impossible composition syntax;
- duplicate structures;
- minimum interatomic-distance violations;
- cell-size limits;
- charge-balance heuristics when explicitly enabled;
- density bounds;
- user-defined chemistry exclusions.

These filters are auditable and may produce a hard rejection.

### Layer E — ML screening

Machine-learning models estimate the expensive quantities needed to prioritize candidates.

Typical tasks:

- geometry relaxation;
- energy prediction;
- forces and stress;
- rough mechanical properties;
- uncertainty estimation;
- short molecular-dynamics stability checks.

The ML stage should be treated as a **surrogate**, never as final evidence of thermodynamic stability.

### Layer F — Acquisition and ranking

Candidates receive a priority based on several signals rather than a single predicted score.

A generic acquisition score can contain:

```text
predicted utility
+ uncertainty bonus
+ novelty bonus
+ diversity bonus
- compute cost
- known-risk penalties
```

A formal implementation might resemble:

```text
score(x) =
    w1 * expected_utility(x)
  + w2 * epistemic_uncertainty(x)
  + w3 * novelty(x)
  + w4 * diversity_contribution(x)
  - w5 * estimated_compute_cost(x)
```

The AI controller may adjust the policy, but the numerical inputs and final decision must be recorded.

### Layer G — Reference calculations

Selected candidates move to higher-fidelity calculations.

Initial reference backend:

- Quantum ESPRESSO for DFT.

Potential later backends:

- CP2K;
- ABINIT;
- Elk/exciting for all-electron cross-checks;
- commercial codes where licenses are available.

Reference stages may include:

1. geometry optimization;
2. converged total-energy calculation;
3. competing magnetic states where relevant;
4. formation energy;
5. convex-hull analysis against reference phases;
6. phonon stability;
7. elastic response;
8. finite-temperature molecular dynamics;
9. higher-level electronic-structure calculations if justified.

### Layer H — Knowledge and provenance store

All computational evidence is written to a persistent store.

The store should distinguish:

- immutable raw calculation artifacts;
- parsed numerical observables;
- model predictions;
- controller decisions;
- human annotations;
- lineage between candidates.

Every observation should include:

```text
code
code_version
input_parameters
pseudopotential identifiers
model checkpoint
hardware
wall time
timestamp
parent calculation
status
failure category
```

## 3. Trust hierarchy

The system should maintain an explicit hierarchy of evidence.

```text
heuristic / LLM suggestion
        <
ML surrogate prediction
        <
converged DFT result
        <
cross-code / higher-level validation
        <
experimental result
```

A lower level can decide what to inspect next. It should not overwrite a contradictory higher-level result.

## 4. Control loop

```text
objective
   |
   v
controller
   |
   v
generate candidates
   |
   v
hard filters
   |
   v
ML relaxation / prediction
   |
   v
rank by utility + uncertainty + novelty
   |
   v
select reference calculations
   |
   v
DFT / phonons / MD
   |
   v
store observations
   |
   v
update surrogate and search policy
   |
   +---------------------> controller
```

## 5. Failure is data

A failed DFT or relaxation job should not be represented only as `failed=true`.

Failure classes should include, for example:

- SCF non-convergence;
- geometry divergence;
- pathological cell collapse;
- unstable magnetic state;
- pseudopotential incompatibility;
- insufficient cutoff;
- memory exhaustion;
- wall-time exhaustion;
- scheduler failure.

The controller can then decide whether the failure is scientific, numerical, or infrastructural.

## 6. Safety against search collapse

The controller must preserve some exploration budget.

A simple initial rule:

- 60% — exploit promising families;
- 25% — explore high-uncertainty candidates;
- 15% — deliberately sample chemically or structurally novel regions.

These fractions should become adaptive later.

## 7. Non-goals for the first version

The first version does not attempt to:

- predict arbitrary material properties from language alone;
- prove synthesizability;
- automate laboratory synthesis;
- replace expert review;
- search all possible chemistry;
- produce publication-grade melting temperatures immediately.

The first milestone is a small, reproducible closed loop.
