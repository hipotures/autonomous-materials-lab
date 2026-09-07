# Fluid Enthalpy-Window Screening

This experiment extends the first ambient-liquid screening to **cryogenic and pressurized liquids**.

The key change is simple:

> Each candidate has its own initial storage state `T0, P0`.

That allows water at room temperature to be compared with liquid methane, oxygen, nitrogen, hydrogen and other fluids without pretending that all of them are liquids at 293 K.

## What this experiment computes

For each configured fluid:

~~~text
initial liquid state: (T0, P0)
        |
        +--> saturated vapor at the same pressure
        |
        +--> gas/fluid state at 500 K, same pressure
        |
        +--> gas/fluid state at 1000 K, same pressure
~~~

The corresponding idealized heat-uptake metrics are:

~~~text
q_to_vapor   = h_sat_vapor(P0) - h_initial(T0, P0)
q_to_500K    = h(500 K, P0)     - h_initial(T0, P0)
q_to_1000K   = h(1000 K, P0)    - h_initial(T0, P0)
~~~

Units are MJ/kg.

The script also reports volumetric values in MJ/L using the initial liquid density.

## Important interpretation

This is still a **thermodynamic screening calculation**, not a capsule simulation.

It does not model:

- nozzles;
- porous flow;
- vehicle geometry;
- hypersonic aerodynamics;
- braking;
- chemical decomposition;
- combustion;
- dissociation;
- tank mass;
- insulation mass.

For reactive fluids, a CoolProp enthalpy value at high temperature does **not** prove that the same chemical species remains intact or inert in the external atmosphere. High-temperature columns are therefore exploratory and must later be filtered by chemistry/stability models.

This is especially important for hydrogen and hydrocarbons. A large sensible-enthalpy window can coexist with a very large exothermic oxidation potential after the hot gas mixes with oxygen-containing air. Chemical heat release is deliberately **not** subtracted in this experiment; it belongs in a separate reacting-flow/chemistry stage.

## Why this is useful

The first experiment answered:

> Which CoolProp fluids that are already liquid at 293.15 K and 1 atm absorb the most heat before vaporization?

This experiment asks a broader question:

> If each candidate starts from a physically meaningful liquid storage state, how much enthalpy headroom does one kilogram provide?

This admits cryogenic fluids to the comparison.

## Files

~~~text
compare_fluids.py   screening script
fluids.csv          editable candidate/storage-state list
requirements.txt    pinned Python dependency
~~~

## Environment

Recommended:

~~~text
Python 3.12
CoolProp 8.0.0
~~~

Create the environment:

~~~bash
python3.12 -m venv .venv
source .venv/bin/activate

python -m pip install --upgrade pip
python -m pip install -r requirements.txt
~~~

## Run

~~~bash
python compare_fluids.py
~~~

Write CSV output:

~~~bash
python compare_fluids.py --output results.csv
~~~

Use another candidate file:

~~~bash
python compare_fluids.py --fluids my-fluids.csv
~~~

Change the high-temperature checkpoints:

~~~bash
python compare_fluids.py --targets-k 400 600 800
~~~

## Temperature sweep and water crossover

The script also performs a default sweep from 300 K to 1500 K in 25 K steps.

It writes:

~~~text
enthalpy_sweep.csv
~~~

with:

~~~text
temperature_K
fluid label
q_MJ_kg
mass_ratio_vs_water
~~~

and prints an approximate crossover temperature where each candidate first reaches or exceeds water in MJ/kg.

Example:

~~~bash
python compare_fluids.py \
  --sweep-start-k 300 \
  --sweep-end-k 1500 \
  --sweep-step-k 10 \
  --sweep-output enthalpy_sweep.csv
~~~

A missing crossover means either that the candidate did not beat water in the searched interval or that CoolProp stopped supporting one of the required states before a crossover could be established.

## Candidate file format

`fluids.csv` contains:

~~~text
label,coolprop_name,T0_K,P0_Pa,notes
~~~

Example:

~~~text
Water,Water,293.15,101325,room-temperature reference
Liquid methane,Methane,110.0,101325,illustrative subcooled liquid state
~~~

The included states are **screening reference states**, not certified spacecraft tank conditions.

Change them when a concrete vehicle/storage architecture is defined.

## Water normalization

For every metric that can be calculated for both water and the candidate:

~~~text
mass_ratio_vs_water = q_water / q_candidate
~~~

Therefore:

~~~text
0.80 -> idealized required mass is 20% lower than water
1.00 -> equal to water
1.20 -> idealized required mass is 20% higher than water
~~~

This comparison inherits the chosen initial states.

## Unsupported states

CoolProp does not cover every fluid over every temperature range.

If a final state is outside the equation-of-state validity region, the script prints a blank value instead of treating the fluid as bad.

That distinction matters:

~~~text
unsupported by this property model != physically poor candidate
~~~

## Suggested next calculation

Run the included candidate set first.

Then edit `fluids.csv` and add:

- other cryogenic liquids;
- pressurized liquids;
- refrigerants;
- candidate working fluids you want to inspect.

Only after the simple thermodynamic ranking becomes informative should we add chemical stability, mixtures, or molecular simulation.
