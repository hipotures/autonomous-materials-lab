# Project Goal

## Primary research objective

The project is **not** a water-cooled heat-shield project.

Its primary liquid-TPS research objective is:

> Discover, model and rank liquids or liquid mixtures that can protect an atmospheric-entry vehicle while requiring as little working-fluid mass as possible, with any useful braking effect included in the system score.

The working fluid is a **search variable**.

Water is only:

- a reference fluid with well-known properties;
- a calibration case for the evaluator;
- a benchmark denominator for later comparisons.

Water is **not** the intended final working fluid, and optimizing a water-cooled wall is not the project's scientific objective.

## Core question

For a fixed mission and a fixed comparison architecture:

~~~text
Which fluid requires the least mass?
~~~

subject to:

~~~text
thermal constraints satisfied
required vehicle state / deceleration satisfied
fluid can be stored and delivered
flow model remains physically valid
material compatibility constraints satisfied
~~~

The first fluid-only objective is:

~~~text
minimize M_fluid_required
~~~

Later, after hardware models are available:

~~~text
minimize M_total_system
~~~

where the system mass includes fluid, tank, manifold, pressure hardware, porous structure/nozzles and control hardware.

## What is being discovered

The search can progress through increasingly uncertain candidate classes:

~~~text
known pure liquids
        |
        v
known mixtures
        |
        v
new mixture ratios / formulations
        |
        v
hypothetical molecular candidates
        |
        v
AI-generated or simulation-derived candidates
~~~

The evaluator is the common scoring function for all of them.

## Evaluation concept

Each candidate fluid is passed through the same mission evaluator:

~~~text
candidate fluid / mixture
        |
        v
thermodynamic and transport properties
        |
        v
heating + phase-change model
        |
        v
delivery / vapor-generation model
        |
        v
optional directed-jet momentum model
        |
        v
aerodynamic interaction model
        |
        v
trajectory / thermal constraints
        |
        v
required fluid mass
~~~

A candidate can gain value in two ways:

1. absorb or divert more heat per unit mass;
2. create useful net braking per unit mass.

The second effect must include any change in ordinary aerodynamic drag caused by injection.

## Role of water

Water should be evaluated by the same scoring path wherever possible, but its material properties do not need to be discovered.

Its role is:

~~~text
known reference data
        |
        v
same evaluator
        |
        v
M_water_reference
~~~

Later:

~~~text
candidate_score =
    M_candidate_required / M_water_reference
~~~

A result below 1.0 means the candidate requires less working-fluid mass than water under the same model and assumptions.

The water calculation is therefore a **calibration and comparison task**, not the primary research target.

## Fair comparison rule

When comparing fluids, keep geometry and control assumptions fixed unless the experiment explicitly studies co-optimization.

Otherwise the project can accidentally credit a fluid for an improvement actually caused by:

- a better nozzle;
- different injection timing;
- different pressure;
- a different porous surface;
- a different trajectory-control law.

Recommended progression:

~~~text
Experiment A:
fixed architecture + fixed control
compare fluids only

Experiment B:
fluid-specific control optimization
same hardware bounds

Experiment C:
joint fluid + geometry + control optimization
include hardware mass
~~~

## Discovery versus engineering

The project contains two related but separate problems.

### Materials / fluid discovery

Find fluids and mixtures with useful combinations of:

- heat capacity;
- vaporization enthalpy;
- density;
- vapor pressure;
- viscosity;
- thermal conductivity;
- molecular weight of vapor;
- decomposition chemistry;
- residue tendency;
- compatibility;
- storage properties.

### System engineering

Determine how a selected fluid performs in:

- channels;
- porous surfaces;
- vapor chambers;
- micro-nozzles;
- hypersonic flow;
- a complete entry trajectory.

The system evaluator connects these problems, but neither should be confused with the other.

## Role of AI

AI is used to choose what to evaluate next.

It may:

- select chemical families;
- propose mixtures;
- identify missing property calculations;
- choose which candidate deserves MD / ab-initio / CFD;
- balance exploration and exploitation.

It must not declare a fluid superior based on language-model intuition.

The final ranking comes from the evaluator and recorded physical evidence.

## Success criterion

A meaningful project result is not:

> Water cooling works.

A meaningful result is:

> Under a stated mission, architecture and model fidelity, candidate X requires Y% less working-fluid mass than the water reference while satisfying the same constraints.

A stronger later result is:

> Candidate X produces the lowest total system mass among the tested candidates after hardware and control optimization.
