# Benchmark: Water Reference Case

## 1. Purpose

Water is the first reference fluid for the liquid thermal-protection use case.

The benchmark does **not** assume that water is globally optimal. It establishes a reproducible system-level reference against which every later fluid, mixture and active-cooling design can be compared.

## 2. Why water is useful

For water, the early materials-discovery pipeline is unnecessary.

We do not need to:

- generate a molecule;
- predict basic chemical plausibility;
- estimate basic thermodynamic properties from an ML model;
- rank it among unknown candidates before simulation.

Instead, trusted thermodynamic and transport-property data can be supplied directly to the system evaluator.

This isolates and validates the difficult system physics:

~~~text
trajectory
+
heat transfer
+
phase change
+
fluid delivery
+
jet / transpiration behavior
+
vehicle dynamics
~~~

## 3. Baseline configurations

### W0 — no active liquid cooling

Reference vehicle / passive TPS configuration.

Purpose:

- establish baseline heat load;
- establish baseline drag and trajectory;
- provide comparison for active concepts.

### W1 — water transpiration only

Water passes through a porous or microperforated surface without deliberate opposing-nozzle momentum recovery.

Purpose:

- measure the benefit of thermal absorption and blowing alone.

### W2 — water with directed microjets

Water is heated internally, converted to vapor and discharged through directed nozzles.

Purpose:

- measure the combined effect of cooling, flow-field modification and jet reaction.

### W3 — optimized water system

Optimize:

- water mass;
- flow schedule;
- pressure schedule;
- nozzle distribution;
- nozzle geometry within selected bounds.

Purpose:

- establish the strongest water baseline before comparing new fluids.

## 4. Primary benchmark result

The central benchmark number is:

~~~text
M_water_min
~~~

defined as the minimum water mass required to satisfy all thermal and trajectory constraints under the benchmark mission definition.

Later candidates are compared using:

~~~text
mass_ratio = M_candidate_min / M_water_min
~~~

Interpretation:

~~~text
mass_ratio < 1.0  -> candidate requires less fluid than water
mass_ratio = 1.0  -> equal to water baseline
mass_ratio > 1.0  -> candidate is worse on fluid mass
~~~

Example:

~~~text
mass_ratio = 0.72
~~~

means the candidate requires 28% less fluid mass than the optimized water benchmark under the same mission assumptions.

## 5. Later system-level metric

After the fluid-only benchmark is stable, replace the metric with:

~~~text
system_mass_ratio =
    M_total_candidate_system /
    M_total_water_system
~~~

This includes:

- coolant;
- tank;
- pressure/pump hardware;
- manifold;
- porous structure;
- nozzles;
- additional control hardware.

A chemically superior coolant should not win if its storage and delivery hardware makes the total system heavier.

## 6. Required mission definition

The benchmark is meaningless unless the mission is fixed.

The configuration should eventually include:

~~~yaml
vehicle:
  mass_kg: TBD
  geometry: TBD
  reference_area_m2: TBD
  wall_temperature_limit_K: TBD

entry:
  initial_altitude_m: TBD
  initial_velocity_m_s: TBD
  flight_path_angle_deg: TBD
  atmosphere_model: TBD
  angle_of_attack_profile: TBD

target:
  final_altitude_m: TBD
  maximum_final_velocity_m_s: TBD

cooling_system:
  maximum_pressure_Pa: TBD
  available_volume_m3: TBD
  maximum_local_mass_flux_kg_m2_s: TBD
~~~

The first benchmark should use a deliberately simple, documented trajectory rather than attempting to reproduce a complete flight mission immediately.

## 7. Water property inputs

Required property functions include:

- density as a function of temperature and pressure;
- heat capacity;
- enthalpy;
- vapor pressure;
- phase-change enthalpy;
- viscosity;
- thermal conductivity;
- surface tension;
- vapor thermodynamic properties.

Every property source and version must be recorded in provenance.

## 8. Reduced-model acceptance test

Before CFD or molecular simulation is introduced, the first evaluator should be able to answer:

~~~text
Given:
- fixed entry heat-load history
- water initial state
- wall temperature limit
- simple nozzle/transpiration model

Compute:
- water mass consumed
- vapor generation rate
- wall-temperature history
- approximate jet momentum
- approximate net braking contribution
~~~

This model does not need to be high fidelity. It needs to be transparent, numerically stable and useful as a baseline.

## 9. Full benchmark outputs

Eventually record:

~~~text
M_water_min
peak wall temperature
peak backface temperature
water mass flow versus time
remaining water mass versus time
vapor pressure versus time
jet thrust versus time
body drag versus time
net braking force versus time
vehicle velocity versus time
vehicle altitude versus time
heat flux versus time
total absorbed heat
total braking impulse attributable to active system
maximum system pressure
~~~

## 10. Comparison ladder

Report gains separately:

~~~text
Passive baseline:
TPS mass-equivalent metric = ...

Water transpiration:
fluid mass = ...

Water directed microjets:
fluid mass = ...
net braking gain = ...

Optimized water:
fluid mass = ...

Candidate X:
fluid mass = ...
mass ratio vs optimized water = ...
~~~

This prevents a new fluid from receiving credit for improvements actually caused by nozzle geometry or control policy.

## 11. Benchmark philosophy

Water is not the scientific conclusion.

It is the calibration point.

A new candidate is interesting only if it beats water after both are given comparable optimization effort and evaluated with the same physical model.
