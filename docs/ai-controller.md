# AI Scientific Controller

## 1. Role

The AI controller is the policy layer of the system.

It is responsible for deciding **what to compute next**, not for inventing numerical material properties. Energies, forces, stresses, phonons, elastic constants, and thermodynamic quantities must come from dedicated physical or surrogate models.

The controller operates over structured evidence and emits structured actions.

## 2. Inputs

At each decision point the controller may receive:

- research objective;
- remaining compute budget;
- candidate summaries;
- previous DFT results;
- surrogate predictions and uncertainty;
- structural and chemical clusters;
- failed-job diagnostics;
- branch history;
- relevant literature or database-derived facts when explicitly provided;
- available computational tools.

The controller should not be given enormous raw output files when a parsed representation is sufficient.

## 3. Outputs

Every decision should be machine-readable.

Example:

```json
{
  "action": "submit_dft_batch",
  "candidate_ids": [
    "cand-1021",
    "cand-1582",
    "cand-2044"
  ],
  "selection_strategy": {
    "exploit": 1,
    "uncertainty": 1,
    "novelty": 1
  },
  "rationale": [
    "cand-1021 has the best predicted objective value",
    "cand-1582 is high-value but outside the surrogate training domain",
    "cand-2044 represents a structurally novel cluster"
  ],
  "maximum_budget": {
    "gpu_hours": 0,
    "cpu_hours": 18,
    "hpc_jobs": 3
  }
}
```

The execution layer validates the action before running it.

## 4. Allowed decision classes

### Search-direction decisions

- expand a chemical family;
- narrow a stoichiometry range;
- alter generator conditioning;
- increase diversity;
- increase uncertainty sampling;
- explore a structurally distinct family.

### Candidate-routing decisions

- promote to DFT;
- retain for later;
- soft reject;
- request more surrogate evaluation;
- request a specific physical test.

### Diagnostic decisions

- retry a failed SCF calculation;
- test an alternative magnetic initialization;
- increase convergence parameters;
- distinguish likely numerical failure from scientific instability.

### Budget decisions

- spend more DFT budget in a promising region;
- stop an exhausted branch;
- reserve compute for high-uncertainty candidates;
- escalate selected calculations to HPC.

## 5. Decisions the controller must not make alone

The controller must not produce a hard rejection solely because:

- a composition looks unusual;
- a structure is absent from known databases;
- the chemistry is unfamiliar;
- the candidate contradicts typical textbook intuition;
- a language model assigns a low subjective plausibility.

Hard rejection should require either:

- violation of an explicit user constraint;
- deterministic structural invalidity;
- a defined numerical/stability criterion;
- a higher-fidelity physical result.

## 6. Decision policy

A useful initial policy is a mixture of four modes.

### Exploit

Select candidates with the strongest expected target properties.

### Uncertainty sampling

Select candidates where the surrogate is least certain.

### Diversity sampling

Select representatives far from already-tested structures or compositions.

### Scientific hypothesis sampling

Allow the controller to reserve a small fraction of the batch for interpretable hypotheses.

Example batch policy:

```text
50% exploit
25% uncertainty
15% diversity
10% scientific hypotheses
```

The percentages should be configuration, not hard-coded behavior.

## 7. Controller memory

The controller should not rely on conversational memory as the scientific record.

Persistent state belongs in the project database.

Useful records include:

```text
hypothesis
supporting evidence
contradicting evidence
status
candidate IDs
calculation IDs
created_at
updated_at
```

Example hypothesis:

```text
H-004:
B-rich Si-B-N structures may form a low-density high-modulus family.

Support:
- 8 ML-relaxed candidates below threshold
- 2 DFT candidates with favorable formation energies

Contradiction:
- 1 candidate dynamically unstable

Next test:
- DFT relaxation of 4 structurally diverse B-rich candidates
```

## 8. Tool interface

The controller should eventually call a small number of high-level tools rather than shell commands.

Possible interface:

```text
search_candidates(filters)
generate_candidates(spec)
screen_candidates(candidate_ids, model)
submit_dft(candidate_ids, preset)
submit_phonons(candidate_ids, preset)
inspect_failures(calculation_ids)
cluster_candidates(candidate_ids)
get_budget()
get_branch_summary(branch_id)
create_hypothesis(...)
close_branch(...)
```

This makes controller behavior auditable and reduces the risk of arbitrary execution.

## 9. Guardrails

Before execution, a policy validator checks:

- candidate IDs exist;
- requested backend is allowed;
- estimated resource use fits the budget;
- numerical preset exists;
- no forbidden element is introduced;
- maximum batch size is respected;
- duplicate calculations are not launched accidentally.

The controller proposes. The orchestrator validates and executes.

## 10. Evaluation of the AI layer

The controller should be benchmarked independently of the physical models.

Possible metrics:

- best candidate found versus DFT budget;
- number of redundant calculations;
- diversity of sampled space;
- fraction of high-uncertainty regions explored;
- number of successful recovery decisions after numerical failures;
- regret versus a baseline acquisition strategy;
- reproducibility of decisions under repeated runs.

A strong baseline is essential. The AI controller should have to outperform simpler policies such as random sampling, top-N surrogate ranking, and conventional uncertainty sampling.
