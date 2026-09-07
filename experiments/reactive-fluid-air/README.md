# Reactive Working Fluid + Air Sanity Check

This experiment is the chemical companion to the thermodynamic enthalpy-window screen.

It asks a deliberately narrow question:

> If a hot candidate gas mixes with oxygen-containing air, how large is its thermodynamic oxidation potential?

It does **not** model where the reaction occurs, how fast ignition happens, or how much reaction heat returns to a vehicle wall.

## Why this experiment exists

The previous screening experiment showed that liquid hydrogen can have a very large idealized enthalpy window when heated from a cryogenic liquid state to a hot gas.

That alone does not establish net cooling performance.

Hydrogen, methane and ammonia can react exothermically with oxygen. This experiment keeps the two effects separate:

~~~text
CoolProp:
storage state -> hot gas
thermal energy absorbed

Cantera:
hot gas + air -> equilibrium products
chemical energy potentially released
~~~

The first quantity is favorable for cooling.

The second is a chemical-risk scale, not automatically heat returned to the wall.

## Model

The default Cantera mechanism is:

~~~text
gri30.yaml
~~~

The script uses chemical **equilibrium**, not finite-rate combustion kinetics.

For each candidate, hot-gas temperature and equivalence ratio:

1. construct a fuel/air mixture;
2. record the initial mixture enthalpy;
3. equilibrate at constant temperature and pressure (TP);
4. calculate the equilibrium chemical enthalpy decrease;
5. reset the initial state;
6. equilibrate adiabatically at constant enthalpy and pressure (HP);
7. report the adiabatic equilibrium temperature.

Cantera documents `equilibrate("TP")` and `equilibrate("HP")` for these equilibrium constraints.

## Fuels

Initial candidates:

- H2;
- CH4;
- NH3.

The oxidizer is dry simplified air:

~~~text
O2:1, N2:3.76
~~~

## Equivalence ratio

The equivalence ratio `phi` controls fuel richness.

Conceptually:

~~~text
phi < 1  -> lean, excess oxygen
phi = 1  -> stoichiometric
phi > 1  -> fuel-rich, oxygen-limited
~~~

The default sweep is:

~~~text
0.5  1.0  2.0  4.0
~~~

This is useful for the hydrogen question because a hydrogen-rich near-wall layer and a fully mixed stoichiometric hydrogen/air region are chemically very different situations.

## Temperatures

Default incoming hot-gas temperatures:

~~~text
500 K
750 K
1000 K
1250 K
1500 K
~~~

These are screening states, not a reentry boundary-layer model.

## Reported quantities

For every case:

~~~text
fuel
T_initial
pressure
phi
fuel mass fraction in mixture
thermal_q_from_storage_MJ_per_kg_fuel
equilibrium_T_HP
chemical_release_TP_MJ_per_kg_fuel
chemical_to_thermal_ratio
selected equilibrium product mole fractions
~~~

### Thermal quantity

When a storage state is configured and CoolProp supports the requested state:

~~~text
q_thermal =
    h_fuel(T_hot, P)
  - h_fuel(T_storage, P_storage)
~~~

This is the same idealized enthalpy-window concept used in the previous experiment.

### Chemical quantity

At the same initial hot-mixture temperature and pressure:

~~~text
q_chemical =
    (h_initial_mixture - h_equilibrium_TP)
    / initial_fuel_mass_fraction
~~~

Units:

~~~text
MJ per kg of initially supplied fuel
~~~

A positive value means the equilibrium products are lower in enthalpy at the same T and P, so oxidation/reaction can release chemical energy.

This is a theoretical equilibrium scale.

It does **not** mean that this energy is released instantly or that it is deposited into the heat shield.

## Installation

Recommended environment:

~~~text
Python 3.12
Cantera 3.2.0
CoolProp 8.0.0
~~~

Create it with:

~~~bash
python3.12 -m venv .venv
source .venv/bin/activate

python -m pip install --upgrade pip
python -m pip install -r requirements.txt
~~~

Cantera 3.2.0 provides CPython 3.12 Linux wheels.

## Run

~~~bash
python reactive_air_check.py
~~~

Write the full result table:

~~~bash
python reactive_air_check.py --output results.csv
~~~

Use a finer richness sweep:

~~~bash
python reactive_air_check.py \
  --phi 0.5 1 1.5 2 3 4 6 \
  --output results.csv
~~~

Change pressure:

~~~bash
python reactive_air_check.py --pressure-pa 50000
~~~

## Interpretation

Do not rank fluids by simply subtracting the chemical number from the thermal number.

That would imply that all reaction energy returns to the protected wall, which is not known.

Instead, use the output as two separate axes:

~~~text
thermal heat-sink potential
vs.
chemical oxidation potential
~~~

For example, hydrogen may have excellent thermal performance and simultaneously very large oxidation potential.

Whether that is acceptable depends later on:

- mixing distance from the wall;
- oxygen availability;
- residence time;
- ignition delay;
- radical pool;
- local pressure;
- wall catalysis;
- hypersonic transport.

Those are later reacting-flow questions.

## Mechanism limitation

GRI-Mech 3.0 contains H2, CH4, NH3 and the relevant H/C/N/O species, so it is sufficient for this **equilibrium sanity check**.

It should not be treated as a validated modern ammonia ignition mechanism or as a reentry-air mechanism.

Equilibrium uses the available species thermochemistry and elemental constraints. Finite-rate ignition calculations require a mechanism selected and validated for the specific fuel, pressure, temperature and composition range.

## Sources

- Cantera 3.2 equilibrium / Python documentation: https://cantera.org/stable/
- GRI-Mech 3.0 mechanism distributed with Cantera: https://www.cantera.org/stable/examples/input/gri30.html
