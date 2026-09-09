# Research priority: experimentally accessible formulations

## Decision

Make known-component water formulations the primary experimental track. Retain
V5b/V5c/V5d/V5e code and outputs as an archived research track; do not schedule
additional novel-molecule acquisition merely to complete a reference-count target.
No existing experiment or validated water baseline is deleted or relabeled.

The immediate scope includes miscible liquids and dissolved solids. Surfactants
are measurement-only controls until specific data exist. Suspensions are a
separate technology study, not a liquid-mixture fallback. A retained granulate,
porous coating or fixed bed belongs to a surface/hardware experiment rather than
being credited solely to the fluid.

## Last molecular-acquisition report

The supplied V5e-3 run explored 120 of 3,409 eligible candidates. It reported 37
identity-verified probes, 18 local exact-CAS matches, and eight selected targets
(six calibration, two evaluation). These are resolver outcomes, not evidence
that the other structures cannot exist or cannot be synthesized.

Projected post-audit counts were:

| Property | Calibration / evaluation before | Projected after |
| --- | --- | --- |
| Critical temperature | 6 / 3 | 10 / 3 |
| Critical pressure | 6 / 3 | 10 / 3 |
| Vaporization enthalpy | 6 / 2 | 10 / 2 |
| Vapor pressure | 8 / 3 | 12 / 5 |

Only the vapor-pressure planning target of 12/5 was projected closed. No new
reference audit, error benchmark, model fit or rankable promotion was reported.
The adaptive, split-aware acquisition is not an untouched random evaluation
sample. Keep these results as planning artifacts, not a validation upgrade.

## New objective

Find a reproducible formulation and concentration that improves matched cooling
performance relative to water while satisfying preparation, stability, delivery,
pressure, thermal and residue constraints. Include the additive in total consumed
mass. Keep formulation effects separate from altered surface/geometry effects.

The first implementation is [V5m-1](../experiments/formulation-screen-v5m/README.md).
It performs offline cold-liquid and bounded hot-window calculations, declares
missing physics, and creates a gated, randomized test plan plus measurement
analysis. It does not substitute a theoretical property score for a physical
experiment and does not require synthesis of a new molecule.
