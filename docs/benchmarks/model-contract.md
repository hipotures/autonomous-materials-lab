# Water Benchmark: Physical and Numerical Contract

Status: proposed implementation contract; no simulation or flight validation has been performed.

## 1. Separate the experiments

- **W1-R (first milestone):** prescribed external heating, fixed geometry, water supply and wall model. Solve thermal response and coolant consumption. No trajectory or braking claim.
- **W2-R:** add an energy-conserving vapor chamber/nozzle model. Report reaction force as a diagnostic at prescribed external conditions, not as a demonstrated trajectory benefit.
- **W1-C / W2-C:** couple flow-dependent heating and aerodynamic forces to variable-mass trajectory dynamics. These require a validated aerodynamic interaction model.
- **W3:** optimize an explicitly named R or C model, with finite design bounds and recorded search budget.

W0 uses the same boundary conditions and passive wall for mechanism comparisons. A separate mission design comparison may resize the passive TPS, but must then include its mass. Never compare differently sized passive structures as though only the coolant changed.

## 2. Required inputs and validity

An executable case must specify time interval, heated area (distinct from aerodynamic reference area), wall thickness/material/heat capacity/conductivity/emissivity, initial wall and backface states, thermal boundary conditions, coolant inlet state, local external pressure, delivery limits, and control parameter bounds. Mission-coupled cases additionally require initial wet and dry masses, geometry, atmosphere, aerodynamic coefficients, reference frames and terminal events.

All internal quantities use SI; pressures are absolute; enthalpy is mass-specific J/kg. Store the property backend, version, reference-state convention and validity ranges. Mixtures specify species, mass or mole fraction basis, normalized composition and validated interaction parameters. Do not treat a mixture as a mass-weighted table of pure-fluid boiling points.

Reject NaN, infinite values, negative masses, missing units, unsorted time grids and out-of-domain property requests. Distinguish invalid model input from a physically infeasible design. Do not clamp an unsupported thermodynamic state into a feasible one.

A lumped wall temperature is allowed only with a justified spatial approximation (for example, a sufficiently small Biot number for the relevant boundary conditions). Otherwise resolve wall conduction and backface temperature. A global average must not certify a local hotspot constraint.

## 3. Energy accounting

For a simple wall control volume, define positive heat flows explicitly:

```text
dU_wall/dt = Q_external_in - Q_radiation_out - Q_backface_out - Q_wall_to_coolant
```

Use either an incident heating boundary with separate radiation losses or a net heating boundary; do not subtract radiation twice. Store signed backface heat transfer when its direction can reverse.

For a steady coolant stream between stated inlet and outlet stations:

```text
Q_wall_to_coolant + P_external_to_fluid
    = mdot * (h_out - h_in + (u_out^2 - u_in^2)/2 + g*(z_out-z_in))
```

For startup, vapor accumulation or a changing chamber state, use the corresponding transient control-volume mass and energy balances, including stored fluid energy. Account for losses and pump/pressurant energy at their actual boundaries. Vapor production cannot be prescribed independently of heat-transfer capacity and chamber pressure.

If h_out is the chamber enthalpy before an adiabatic nozzle, nozzle kinetic energy comes from the subsequent enthalpy drop. If h_out is the nozzle-exit enthalpy, include exit kinetic energy in the stream balance above. Never count the same available energy twice.

Blowing changes Q_external_in through a separately validated closure. It is not an intrinsic coolant enthalpy. A postprocessed effective J/kg benefit may include integrated avoided heating divided by consumed mass only for a named comparison, a nonzero consumed mass and fixed comparison conditions. Disable uncalibrated blowing credit in the first conservative case, and label that case accordingly; it is not a validated transpiration prediction.

## 4. Pressure, flow and nozzle closure

Pressure, mass flow, throat area and chamber temperature are coupled. Select independent control variables and solve for the rest. Enforce delivery pressure losses, pump power, permeability, flow limits and tank depletion. Local external pressure near the injection site is not necessarily free-stream atmospheric pressure.

For single-phase quasi-steady exhaust, with unit vector n pointing out of the vehicle:

```text
F_jet_vector = -(mdot*u_exit + (p_exit-p_local_external)*A_exit) * n
F_braking = -dot(F_jet_vector, velocity_unit_vector)
```

Thus exhaust directed along vehicle motion produces an opposing reaction. Reverse the exhaust direction in a verification test and check the force sign. Zero mass flow also requires a physically closed/no-discharge pressure boundary; do not leave a fictitious pressure thrust active.

Use the ideal-gas choked-nozzle equation only inside its specified dry-vapor, continuum and thermodynamic assumptions, with a checked critical pressure ratio. Wet expansion, flashing, condensation, rarefaction and backflow require another closure or an explicit unsupported status. A property-library phase lookup does not itself solve two-phase nozzle flow.

Define the force integration surface. If CFD already includes the jet momentum/pressure contribution to total vehicle force, do not add that contribution again. If CFD returns only external body pressure and shear, combine it with a consistently defined nozzle reaction.

## 5. Coupled trajectories and comparisons

Use decreasing vehicle mass with discharged coolant and a thrust convention consistent with the momentum balance; do not add another mass-loss momentum term when it is already included in nozzle thrust. Integrate to specified terminal events, not an arbitrary favorable stopping time. Record unachieved terminal conditions and any skipped heating tail.

Distinguish a force difference at the same external state from a difference between complete trajectories. Integrated active/passive force differences on different trajectories are not a pure isolated jet effect. Report final velocity, elapsed time, remaining mass and constraint margins alongside any impulse diagnostic.

Without a validated injection-dependent heat-transfer and drag closure, coupled braking benefit remains unknown. Test sensitivity to adverse drag change as well as favorable change.

## 6. Optimization and mass reporting

The optimizer reports `best_feasible_mass_found`, not a certified global minimum. Retain `M_water_min` only as the benchmark shorthand with that qualification. Store search bounds, initialization/seeds, evaluator hash, stopping reason, simulation count and constraint margins. Bisection is allowed only when feasibility is established to be monotone in the selected scalar variable under fixed other controls.

Report consumed coolant, loaded coolant (including residual/reserve), passive TPS mass, dry cooling hardware and initial wet mass separately. A fluid-only comparison holds the hardware and limits fixed and discloses their mass. If a hardware mass model is unavailable, system-mass ranking is unavailable; use sensitivity bounds, not a zero hardware mass assumption. A zero-water feasible solution makes a fluid mass ratio undefined.

Independent starts and a coarse bounded search check whether the optimizer has exploited solver failures or a poor local optimum. Re-evaluate selected designs with tighter numerics and independent reference data before making comparative claims. Nominal feasibility and feasibility under uncertain heating, permeability or model parameters are separate outputs.

## 7. Verification gates

These are proposed acceptance tolerances for implementation, not observed results or physical uncertainty estimates. Record and justify any revision before comparing designs.

| Test | Required result |
| --- | --- |
| Closed, adiabatic, zero-flow case | Constant total energy and mass to solver tolerance |
| Constant heating with constant heat capacity, no losses | Match analytic temperature rise within 1e-6 relative error away from zero, with stated absolute tolerance near zero |
| Steady prescribed enthalpy uptake | Match Q = mdot * delta_h for negligible kinetic/potential terms |
| Full transient run | Normalized integrated mass residual <= 1e-6 and energy residual <= 1e-4; specify nonzero normalization scales |
| Time-step/tolerance refinement | Peak wall temperature changes <= 0.5 K and consumed mass <= 0.5%; tighten further if feasibility margins are smaller |
| Limits and failure handling | Depletion, unsupported phase, impossible pressure/flow and timeout never return feasible |
| Nozzle component | Match an independent analytic dry-gas case; reversing direction reverses force |
| Property model | Compare independent reference states over the actual pressure/temperature range, including saturation boundaries |
| Optimizer | Recover a manufactured problem with known feasible optimum; re-evaluate the winning design independently |

Use published IAPWS check values to verify a water-property implementation. Agreement between two libraries using the same equation of state checks implementation consistency, not independent physical accuracy; assess the underlying data and stated formulation uncertainty as well.

Numerical verification establishes equation-solving correctness. Physical validation additionally needs independent property and heat-transfer data for the relevant regime. A deterministic rerun alone establishes neither physical validity nor global optimality.

## Sources

- [IAPWS-95 release and check values](https://iapws.org/technical-guidance/release/IAPWS-95): reference water thermodynamic formulation.

- [CoolProp high-level interface](https://coolprop.org/coolprop/HighLevelAPI.html): phase-aware property inputs and enthalpy conventions.
- [CoolProp mixture documentation](https://coolprop.org/fluid_properties/Mixtures.html): mixture model and interaction data requirements.
- [NASA rocket thrust equations](https://www.grc.nasa.gov/www/k-12/BGP/rktthsum.html): momentum and pressure contributions.
- [NASA transpiration study](https://ntrs.nasa.gov/api/citations/20000012950/downloads/20000012950.pdf): regime-dependent flow/chemistry coupling; not a calibration for this water vehicle.
