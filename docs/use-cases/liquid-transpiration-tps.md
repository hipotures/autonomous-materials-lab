# Use Case: Working-Fluid Discovery for Transpiration / Microjet Thermal Protection

## 1. Purpose

This use case is a concrete inverse-design problem whose primary search variable is the **working fluid**. The project is not intended to design a water-cooled heat shield; water appears only as a known reference fluid.

The objective is to discover, model and rank liquids or liquid mixtures that minimize the mass required to protect and, where useful, additionally decelerate an atmospheric-entry vehicle.

The candidate fluid may contribute through two coupled mechanisms:

1. **thermal protection** — sensible heating, boiling, vaporization, superheating, decomposition or other endothermic processes absorb heat;
2. **momentum exchange** — directed vapor injection can produce braking impulse and modify the external hypersonic flow field.

The system should optimize the complete coupled effect rather than any single fluid property in isolation.

## 2. Primary objective

Initial objective:

~~~text
minimize M_fluid
~~~

subject to:

~~~text
T_wall(x,t) <= T_limit
v(t_end) <= v_target
p_system(t) <= p_limit
all operational constraints satisfied
~~~

The later system-level objective should become:

~~~text
minimize M_total_TPS

M_total_TPS =
    M_fluid
  + M_tank
  + M_manifold
  + M_nozzle_or_porous_structure
  + M_pressure_system
  + M_control_hardware
~~~

This prevents the optimizer from selecting a fluid that appears mass-efficient but requires an impractically heavy delivery system.

## 3. Why mass is the top-level metric

Latent heat, heat capacity, viscosity, exhaust velocity and vapor pressure are intermediate properties, not the mission objective.

A fluid with excellent latent heat may still lose if it:

- requires excessive storage pressure;
- produces poor directed momentum;
- reduces aerodynamic drag too strongly;
- freezes;
- clogs microchannels;
- decomposes into residue;
- attacks structural materials;
- requires a heavy tank or pump.

The real question is:

> How much total system mass is required for a specified entry mission?

## 4. Two principal effectiveness measures

### Thermal effectiveness

Use the control-volume accounting defined in the [generic fluid evaluator contract](../benchmarks/fluid-evaluator-contract.md).

Intrinsic fluid energy uptake may include:

~~~text
sensible heating
+ phase-change enthalpy
+ vapor superheating
+ consistently defined endothermic chemistry
~~~

Blowing-induced reduction of external heat transfer is a **system effect**, not an intrinsic fluid J/kg property. It should be calculated separately and must not be counted again if the wall heat-flux boundary has already been reduced by the blowing model.

Likewise, nozzle kinetic energy must come from a consistent enthalpy drop; it cannot be added as free energy on top of the same fluid enthalpy budget.

### Braking effectiveness

Define net braking impulse generated per unit fluid mass:

~~~text
j_eff = delta braking impulse / M_fluid
~~~

The net impulse must include both:

- direct reaction force from directed exhaust;
- the change in aerodynamic drag caused by injected flow.

Therefore:

~~~text
net braking effect =
    jet reaction
  + drag_with_injection
  - baseline_drag
~~~

A jet that creates retrothrust but strongly reduces body drag may have less net braking benefit than expected.

## 5. Coupled trajectory effect

Cooling and braking cannot be optimized independently.

~~~text
mass flow
   |
   v
jet / flow-field interaction
   |
   v
vehicle acceleration
   |
   v
velocity history
   |
   v
future aerodynamic heating
   |
   v
future coolant demand
~~~

The final evaluator must eventually simulate the complete entry trajectory rather than score candidates at one fixed heat flux.

## 6. Candidate design variables

### Fluid

- pure liquid or mixture;
- molecular composition;
- density;
- specific heat;
- vaporization enthalpy;
- vapor pressure;
- boiling curve;
- critical properties;
- viscosity;
- thermal conductivity;
- surface tension;
- wetting behavior;
- freezing point;
- decomposition chemistry;
- residue/coking tendency;
- corrosion/material compatibility;
- vapor molecular weight;
- gas heat capacity ratio;
- high-temperature chemistry.

### Delivery architecture

- storage pressure;
- pump or pressure-feed strategy;
- manifold topology;
- channel depth;
- porous-medium permeability;
- nozzle throat diameter;
- nozzle expansion geometry;
- nozzle angle;
- nozzle spatial distribution;
- number of nozzles;
- local mass flux.

### Control law

- mass flow as a function of time;
- pressure as a function of time;
- spatial flow distribution;
- switching between transpiration and directed injection;
- optional mixture ratio as a function of time.

## 7. Hybrid surface concept

The system should not assume identical holes across the entire heat shield.

A plausible design family is:

~~~text
central stagnation region
    -> directed opposing microjets

mid-radius region
    -> mixed microjet + transpiration cooling

outer region
    -> low-rate transpiration cooling
~~~

The optimizer should be allowed to decide whether such specialization is beneficial.

## 8. Simulation layers

### Level 0 — tabulated fluid properties

Known fluids are evaluated from trusted thermodynamic and transport-property data.

### Level 1 — reduced-order thermal model

Calculate:

- wall heat load;
- fluid heating;
- boiling/vaporization;
- required mass flow;
- approximate porous/nozzle heat transfer.

### Level 2 — reduced nozzle and force model

Calculate:

- chamber state;
- nozzle mass flow;
- exhaust velocity;
- pressure thrust;
- first-order braking impulse.

### Level 3 — aerodynamic interaction

Model:

- shock displacement;
- boundary-layer modification;
- heat-flux reduction;
- aerodynamic drag change;
- interaction between multiple jets.

This may initially use reduced-order correlations and later CFD.

### Level 4 — molecular / chemical validation

Use higher-fidelity tools for selected fluids:

- classical or ML molecular dynamics;
- ab-initio MD;
- high-temperature reaction chemistry;
- decomposition and residue analysis.

### Level 5 — coupled trajectory simulation

Integrate:

- vehicle dynamics;
- atmospheric model;
- aerodynamic heating;
- fluid consumption;
- nozzle forces;
- drag modification;
- wall temperature.

The final score is computed here.

## 9. Proposed software additions

The use case may require adding:

- **Cantera** — thermodynamics, transport and reaction kinetics;
- **LAMMPS** — large-scale molecular dynamics and ML-potential MD;
- **CP2K** — ab-initio molecular dynamics and higher-fidelity liquid chemistry;
- **SU2 or another open CFD solver** — hypersonic and reacting-flow calculations;
- a custom reduced-order thermal / trajectory solver;
- a custom porous-wall / microjet boundary model.

Quantum ESPRESSO remains useful for selected electronic-structure questions but is not the central solver for this fluid use case.

## 10. Discovery progression

Water verification comes first only to validate the evaluator. The scientific search begins immediately afterward.

~~~text
Evaluator verification:
water reference

Discovery Stage 1:
known pure liquids

Discovery Stage 2:
known mixtures and mixture-ratio search

Discovery Stage 3:
hypothetical / predicted molecular candidates

Discovery Stage 4:
high-fidelity validation of promising candidates

Discovery Stage 5:
joint fluid + geometry + control optimization
~~~

Water is not a discovery stage and does not need the early chemistry-generation or property-prediction flow.

## 11. Validation hierarchy

~~~text
tabulated-property screening
        <
reduced thermal/nozzle model
        <
validated CFD / MD
        <
coupled trajectory model
        <
experimental subsystem data
        <
flight data
~~~

The AI controller selects what to test next. It does not override higher-fidelity evidence.

## 12. Open scientific questions

The project should explicitly track uncertainties such as:

- whether directed injection increases or decreases net deceleration for a given regime;
- optimal balance between transpiration and opposing jets;
- stability of two-phase flow in embedded microchannels;
- survivability of micro-nozzles under severe heating;
- whether vapor can be generated at useful chamber pressure from recovered heat alone;
- sensitivity to angle of attack and atmospheric uncertainty;
- susceptibility to nozzle blockage;
- effect of coolant chemistry on catalytic surface reactions;
- optimal division between passive TPS and active fluid mass.

These are research questions, not assumptions.
