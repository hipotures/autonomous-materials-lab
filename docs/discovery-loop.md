# Discovery Loop

## 1. Principle

The discovery loop should answer one question repeatedly:

> Given everything observed so far and a finite compute budget, which calculation should be performed next?

The output of one iteration becomes the evidence for the next.

## 2. Iteration state

Each iteration has a state object:

```json
{
  "iteration": 12,
  "objective_id": "ternary-nitride-001",
  "budget_remaining": {
    "gpu_hours": 71.5,
    "cpu_hours": 122.0,
    "hpc_jobs": 31
  },
  "candidate_pool_size": 4820,
  "surrogate_version": "model-007",
  "controller_policy_version": "policy-003"
}
```

## 3. Step-by-step loop

### Step 1 — Select a search region

The controller selects one or more regions such as:

- a chemical system;
- a stoichiometry family;
- a prototype family;
- a pressure range;
- a magnetic hypothesis;
- a structural motif.

The decision includes an explicit reason.

Example:

```text
Explore B-Si-N compositions because:
1. previous candidates show low predicted formation energies,
2. uncertainty remains high for B-rich structures,
3. the family satisfies the density constraint,
4. only 8% of the current budget has been spent in this region.
```

### Step 2 — Generate a batch

The generator creates a bounded batch, for example 500 structures.

Generation should be reproducible from:

- generator name;
- checkpoint/version;
- random seed;
- conditioning parameters.

### Step 3 — Normalize and deduplicate

Structures are standardized and duplicates removed.

Suggested tools:

- pymatgen structure matching;
- symmetry normalization;
- composition canonicalization.

### Step 4 — Apply hard filters

Candidates violating explicit constraints are rejected.

A hard rejection always stores a reason code.

Example:

```text
REJECT_FORBIDDEN_ELEMENT
REJECT_DENSITY_LIMIT
REJECT_ATOMIC_OVERLAP
REJECT_CELL_TOO_LARGE
REJECT_DUPLICATE
```

### Step 5 — ML relaxation

Remaining candidates are relaxed with an ML interatomic potential.

Outputs:

- relaxed geometry;
- predicted energy;
- forces;
- stress;
- model uncertainty if available;
- convergence status;
- compute cost.

### Step 6 — Fast descriptors

Calculate cheap descriptors useful for search policy:

- composition features;
- density;
- symmetry;
- coordination environment;
- structural fingerprints;
- distance to known candidates;
- distance to training distribution.

### Step 7 — Score candidates

The acquisition layer computes multiple components.

Example:

```json
{
  "candidate_id": "cand-019283",
  "utility": 0.78,
  "uncertainty": 0.64,
  "novelty": 0.81,
  "diversity": 0.55,
  "estimated_dft_cost": 0.31,
  "aggregate_priority": 0.73
}
```

The aggregate score is useful, but individual components must be preserved.

### Step 8 — Controller selection

The AI controller chooses a bounded set for DFT.

It may select candidates for different reasons:

```text
candidate A -> best predicted objective
candidate B -> highest epistemic uncertainty
candidate C -> representative of a new cluster
candidate D -> near known phase boundary
candidate E -> retry after interpretable numerical failure
```

This is preferable to simply taking the top N by one score.

### Step 9 — DFT relaxation and validation

The first DFT stage should be deliberately modest.

Initial sequence:

1. coarse but safe relaxation;
2. tighter relaxation;
3. static converged energy;
4. basic magnetic alternatives when relevant.

Only candidates surviving this stage should proceed to more expensive calculations.

### Step 10 — Stability analysis

For surviving candidates:

- compute formation energies;
- compare against relevant competing phases;
- estimate energy above the convex hull.

Important distinction:

```text
low total energy != thermodynamic stability
```

The comparison set is therefore part of the scientific result.

### Step 11 — Dynamic stability

Promising candidates can move to phonon calculations.

Strong imaginary modes should trigger one of several actions:

- reject;
- distort along the unstable mode and relax again;
- increase the supercell;
- inspect whether the instability is numerical.

The controller should not automatically equate every small imaginary mode with a physically unstable material.

### Step 12 — Feed evidence back

New DFT observations are added to the data store.

Possible responses:

- fine-tune or retrain a surrogate;
- recalibrate uncertainty;
- change sampling weights;
- open a new chemical branch;
- terminate an exhausted branch.

## 4. Exploration versus exploitation

A practical batch can intentionally mix policies.

Example for a 40-job DFT batch:

```text
20 jobs: best predicted utility
10 jobs: highest uncertainty
 6 jobs: structural/chemical diversity
 4 jobs: controller-selected scientific hypotheses
```

This prevents the system from repeatedly sampling the same local optimum.

## 5. Retry policy

The system should separate:

### Numerical retry

Examples:

- increase SCF iterations;
- change mixing;
- alter smearing;
- increase cutoff;
- change k-point density;
- restart from a saved charge density.

### Scientific branch

Examples:

- test alternative magnetic order;
- distort along an unstable phonon mode;
- test a nearby stoichiometry;
- test pressure dependence.

A numerical failure can trigger a bounded retry policy. A scientific branch consumes research budget and should require an explicit controller decision.

## 6. Stopping conditions

A search branch may stop when any of the following is true:

- budget exhausted;
- no improvement after N iterations;
- uncertainty becomes uniformly low;
- all nearby candidates fail stability criteria;
- search has converged to already-known structures;
- user-defined objective reached.

## 7. Minimum viable loop

The first implementation should be much smaller than the final architecture:

```text
known composition family
        |
        v
generate / perturb structures
        |
        v
ML relax
        |
        v
rank
        |
        v
DFT top + uncertain subset
        |
        v
store results
        |
        v
repeat
```

The MVP should first demonstrate that it can rediscover a known answer. Only after this benchmark should it be allowed to spend meaningful compute on genuinely novel candidates.
