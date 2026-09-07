# Benchmark: Water Reference Case

## 1. Purpose

Water is a **reference and verification fluid**.

It is not the target material of the project, and this document does not define a water-cooled thermal-protection system as the desired design.

Water is used because its thermophysical properties are well characterized. That makes it useful for checking whether the generic fluid evaluator behaves correctly before evaluating less familiar candidates.

## 2. What the water case is for

The water case has three jobs:

1. verify the property-model implementation;
2. verify mass, energy and phase-change accounting;
3. provide a comparison value for candidate fluids.

It should produce a reference quantity such as:

~~~text
M_water_reference
~~~

for a fixed mission, fixed hardware configuration and fixed control-law class.

## 3. What the water case is not for

The project's scientific objective is **not** to optimize a water-cooled wall.

The initial comparison must not give water a specially optimized system while evaluating other fluids under different assumptions.

For the first fluid-ranking experiment, hold fixed:

- geometry;
- porous/nozzle architecture;
- pressure limits;
- control-law form;
- mission;
- numerical model.

Then vary the fluid.

## 4. Comparison metric

For a candidate fluid:

~~~text
mass_ratio =
    M_candidate_required /
    M_water_reference
~~~

Interpretation:

~~~text
mass_ratio < 1.0  -> candidate requires less working-fluid mass
mass_ratio = 1.0  -> equal to water under the same assumptions
mass_ratio > 1.0  -> candidate requires more working-fluid mass
~~~

This ratio is meaningful only for matching evaluator versions, model fidelity and comparison assumptions.

## 5. Verification sequence

### W-PROP — property verification

Compare selected water states against independent reference values.

Check:

- liquid enthalpy;
- vapor enthalpy;
- density;
- heat capacity;
- saturation behavior;
- viscosity / conductivity where used.

### W-ENERGY — thermal accounting verification

Run simple cases with known analytic or reference behavior.

Check:

- heating without phase change;
- heating through phase change;
- conservation residuals;
- fluid depletion.

### W-SYS — reference system evaluation

Run water through the **same generic evaluator** that will later evaluate every other fluid.

Record:

- required loaded mass;
- consumed mass;
- peak temperatures;
- thermal constraint margins;
- pressure / flow state;
- any modeled braking contribution;
- evaluator version;
- hardware and control assumptions.

This creates the water comparison point.

## 6. Generic evaluator requirement

The API should conceptually be:

~~~text
evaluate(candidate_fluid, mission, hardware, control_policy)
~~~

not:

~~~text
evaluate_water(...)
~~~

Water may use a high-quality specialized property backend, but the surrounding thermal, force and trajectory logic must be fluid-agnostic.

See [Fluid Evaluator: Physical and Numerical Contract](fluid-evaluator-contract.md).

## 7. Relation to discovery

After water verification, the next milestone is immediately:

~~~text
water
known fluid A
known fluid B
known fluid C
...
    |
    v
same evaluator
    |
    v
rank by required mass
~~~

The purpose of water is to make those later numbers interpretable.

## 8. Later optimized-water comparison

A separately optimized water design may be useful later, but only when every candidate fluid receives comparable optimization effort.

For example:

~~~text
water + optimized control
candidate A + optimized control
candidate B + optimized control
~~~

and eventually:

~~~text
water + optimized hardware + control
candidate A + optimized hardware + control
~~~

At that stage total system mass should be compared, not fluid mass alone.

## 9. Reference principle

Water is the ruler.

It is not the object being designed.
