# Candidate Lifecycle

## 1. Candidate as a persistent entity

A candidate is not just a CIF file. It is an entity with lineage, structures at different fidelity levels, predictions, reference calculations, and state transitions.

Suggested identity:

```text
candidate_id = immutable project-local UUID
```

A candidate may have several associated structures:

- generated structure;
- ML-relaxed structure;
- DFT-relaxed structure;
- distorted structure after a phonon instability;
- pressure-dependent structure.

## 2. Suggested states

```text
GENERATED
  |
  v
NORMALIZED
  |
  v
FILTERED
  |
  +----> HARD_REJECTED
  |
  v
ML_QUEUED
  |
  v
ML_RELAXED
  |
  +----> ML_FAILED
  |
  v
RANKED
  |
  +----> SOFT_REJECTED
  |
  +----> HELD
  |
  v
DFT_QUEUED
  |
  v
DFT_RELAXED
  |
  +----> DFT_FAILED
  |
  v
DFT_VALIDATED
  |
  v
STABILITY_ANALYSIS
  |
  +----> UNSTABLE
  |
  v
PROMISING
  |
  v
ADVANCED_VALIDATION
```

States should describe process state, not scientific certainty.

## 3. Hard versus soft rejection

### Hard rejection

Candidate cannot proceed under the current research objective.

Examples:

- forbidden element;
- atomic overlap below a defined threshold;
- cell exceeds explicit atom-count limit;
- duplicate of an existing candidate within configured tolerances;
- violates a strict density bound after defined relaxation stage.

### Soft rejection

Candidate is currently low priority but remains in the database.

Examples:

- mediocre surrogate score;
- redundant structural cluster;
- low novelty;
- weak objective fit;
- budget insufficient this iteration.

A soft-rejected candidate may be revived after new evidence appears.

## 4. Provenance

Every state transition records:

```json
{
  "candidate_id": "cand-019283",
  "from": "ML_RELAXED",
  "to": "DFT_QUEUED",
  "actor": "controller",
  "reason_code": "HIGH_UNCERTAINTY_AND_HIGH_UTILITY",
  "decision_id": "decision-881",
  "timestamp": "..."
}
```

No candidate should disappear without an auditable transition.

## 5. Structure lineage

When a new structure is derived from another candidate, the relation is explicit.

Examples:

```text
generated_from
substituted_from
mutated_from
distorted_from
pressure_relaxed_from
phonon_mode_distortion_of
```

This allows the system to discover families rather than isolated points.

## 6. Scientific evidence model

A candidate can accumulate claims such as:

```text
predicted_density
predicted_formation_energy
DFT_total_energy
DFT_formation_energy
energy_above_hull
phonon_stability
bulk_modulus
band_gap
melting_estimate
```

Each claim stores:

- value;
- units;
- method;
- fidelity level;
- uncertainty if applicable;
- calculation ID;
- input structure ID.

There should be no single mutable field called simply `energy`.

## 7. Known versus novel

Novelty is not binary.

Possible dimensions:

- novel composition;
- novel stoichiometry;
- novel structure for a known composition;
- novel local coordination motif;
- far from known database structures;
- novel predicted property combination.

The system should store the basis for a novelty score so that the controller can distinguish novelty from mere model uncertainty.

## 8. Promotion criteria

Promotion from ML screening to DFT can use configurable gates.

Example:

```text
Promote if any of the following:
- top 5% expected utility
- top 5% uncertainty
- representative of a novel cluster
- selected by explicit hypothesis test

And:
- no hard constraint violated
- DFT cost estimate below configured maximum
```

## 9. Advanced validation

A candidate should be called "promising" only after a defined evidence bundle is complete.

Example minimal bundle:

- converged DFT relaxation;
- converged static energy;
- formation energy computed against reference states;
- convex-hull analysis;
- no severe structural pathology;
- basic dynamic-stability assessment.

For different research goals, this bundle can be changed.

## 10. Final status

Even the highest internal status should avoid language such as "discovered material" until experimental evidence exists.

Suggested labels:

```text
computational_candidate
computationally_promising
DFT_stable_candidate
dynamically_stable_candidate
experimentally_reported
experimentally_validated
```
