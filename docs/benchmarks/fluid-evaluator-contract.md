# Fluid Evaluator: Physical and Numerical Contract

Status: design contract. It applies to **any candidate fluid or mixture**. Water is one verification/reference case.

## 1. Purpose

The evaluator is the scientific scoring function of the liquid-TPS discovery track.

It must answer:

> Given a candidate fluid, a mission definition and a specified delivery/injection architecture, how much fluid is required to satisfy the thermal and trajectory constraints?

The evaluator must not contain water-specific logic except where a water reference dataset is explicitly used for verification.

## 2. Separate model fidelities

Keep distinct experiment classes.

### T0 — material/property screening

No vehicle simulation.

Evaluate candidate properties and reject clearly incompatible states or fluids.

### T1 — prescribed-boundary thermal evaluator

External heating is prescribed.

Solve:

- wall/fluid energy transfer;
- phase change;
- fluid consumption;
- delivery constraints.

No net trajectory claim is allowed.

### J1 — nozzle / directed-discharge component

Add:

- chamber state;
- mass flow;
- nozzle energy balance;
- reaction force.

External aerodynamic conditions remain prescribed.

This provides a force diagnostic, not yet a demonstrated mission-level braking benefit.

### C1 — coupled aerothermal / trajectory evaluator

Couple:

- injection-dependent heat transfer;
- injection-dependent aerodynamic force;
- vehicle dynamics;
- decreasing vehicle mass;
- fluid consumption.

Only this class can evaluate net trajectory benefit.

## 3. Control-volume energy accounting

Use explicit signs and boundaries.

For a wall control volume:

~~~text
dU_wall/dt =
    Q_external_in
  - Q_radiation_out
  - Q_backface_out
  - Q_wall_to_fluid
~~~

Use either incident external heating plus separate radiation, or a net heating boundary. Do not subtract the same loss twice.

For a steady fluid stream between defined stations:

~~~text
Q_wall_to_fluid + P_external_to_fluid
    =
mdot * (
    h_out - h_in
  + (u_out^2 - u_in^2)/2
  + g*(z_out-z_in)
)
~~~

For transient vapor accumulation or changing chamber state, include stored fluid mass and energy.

Do not prescribe vapor generation independently of heat-transfer capacity and chamber pressure.

## 4. Avoid double counting

If the outlet state is defined before an adiabatic nozzle, nozzle kinetic energy comes from a later enthalpy drop.

If the outlet state is the nozzle exit, the stream energy balance must already include exit kinetic energy.

Never count:

~~~text
vaporization enthalpy
+
the same enthalpy again as nozzle kinetic energy
~~~

Blowing-induced reduction of external heating is also a separate system effect. It is not an intrinsic J/kg property of the fluid.

## 5. Mass and phase validity

Reject or explicitly mark unsupported:

- negative mass;
- NaN / infinite states;
- out-of-range property requests;
- impossible pressure/temperature combinations;
- unsupported two-phase nozzle states;
- depletion before mission completion;
- impossible delivery pressure;
- solver timeout or non-convergence.

A numerical failure is not a physical rejection of the candidate.

## 6. Directed-jet force convention

Define the nozzle outward unit vector n.

For a simple single-phase exhaust model:

~~~text
F_vehicle_from_jet =
    -(mdot*u_exit + (p_exit-p_external_local)*A_exit) * n
~~~

Braking force is the component opposing vehicle velocity.

Verification requirements:

- reversing nozzle direction reverses the reaction-force sign;
- zero mass flow produces zero momentum thrust;
- no fictitious pressure thrust remains when the discharge path is physically closed.

## 7. Net braking metric

The relevant system quantity is not nozzle thrust alone.

At matched external state:

~~~text
delta F_braking =
    F_jet_reaction
  + D_with_injection
  - D_without_injection
~~~

For full missions, compare complete trajectories rather than integrating unmatched states as if they were directly equivalent.

A directed jet can reduce ordinary body drag, so positive nozzle reaction does not guarantee positive net deceleration benefit.

## 8. Fluid comparison discipline

For a fluid-only comparison, hold constant:

- mission;
- external model;
- hardware geometry;
- nozzle/porous architecture;
- pressure limits;
- control-law class;
- numerical tolerances.

Then report:

~~~text
best_feasible_fluid_mass_found
~~~

Do not claim a certified global minimum unless the optimization method actually supports such a proof.

Later co-optimization experiments may vary geometry and control, but must include hardware-mass effects before making system-level ranking claims.

## 9. Mixtures

A mixture is not a mass-weighted average of pure-fluid boiling points or latent heats.

Record:

- species;
- mole or mass basis;
- normalized composition;
- interaction model;
- property backend;
- validated state range.

Unsupported mixture states must remain unsupported rather than being silently extrapolated.

## 10. Numerical verification

Before candidate ranking, the evaluator should pass:

| Test | Expected behavior |
| --- | --- |
| Closed adiabatic zero-flow case | mass and total energy conserved to solver tolerance |
| Constant heating, constant heat capacity | analytic temperature rise reproduced |
| Prescribed steady enthalpy uptake | Q = mdot * delta_h when kinetic/potential terms are negligible |
| Time-step refinement | key outputs converge within predefined tolerance |
| Fluid depletion | infeasible, never silently clipped |
| Unsupported phase/property state | explicit unsupported result |
| Nozzle direction reversal | reaction-force sign reverses |
| Failed solver | failure state, never interpreted as a good candidate |

Store residuals and convergence metrics with every accepted result.

## 11. Physical validation

Numerical verification only proves that the equations were solved consistently.

Physical validation separately requires appropriate reference data for:

- thermophysical properties;
- heat transfer;
- phase change;
- porous flow;
- nozzle behavior;
- hypersonic injection interaction.

A higher-complexity model is not automatically more trustworthy.

## 12. Water-specific verification

Water is useful because high-quality property references exist.

Use it to test:

- property implementation;
- saturation transitions;
- energy balance;
- baseline mass calculation.

Passing water verification does not make water the design target.

The same evaluator must then accept non-water candidates without architectural changes.

## 13. Provenance

Every evaluation records at minimum:

~~~text
candidate identity and revision
composition / mixture definition
property backend and version
model class and version
mission definition
geometry/control configuration
solver tolerances
hardware/software versions where relevant
attempt status
mass/energy residuals
constraint margins
raw artifact checksums
~~~

The evaluator is part of the scientific result.
