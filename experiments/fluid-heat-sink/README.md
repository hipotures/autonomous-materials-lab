# Fluid Heat-Sink Screening Experiment

This is intentionally the **simplest useful calculation** in the liquid-working-fluid track.

It does not model:

- capsule geometry;
- nozzles;
- porous flow;
- hypersonic aerodynamics;
- trajectory;
- braking;
- chemical decomposition;
- hardware mass.

It asks only:

> At a fixed pressure, how much heat can one kilogram of a liquid absorb while going from its initial liquid state to saturated vapor?

## Metric

For each fluid:

~~~text
q_total =
    h_saturated_vapor(P)
  - h_initial_liquid(T0, P)
~~~

This naturally includes:

~~~text
sensible heating to the boiling point
+
latent heat of vaporization
~~~

The primary unit is:

~~~text
MJ/kg
~~~

For a fixed heat load, the idealized required fluid mass is inversely proportional to this value.

Water is used only as a reference:

~~~text
mass_ratio_vs_water =
    q_water / q_candidate
~~~

Interpretation:

~~~text
< 1.0  -> less fluid mass than water in this simplified model
= 1.0  -> same as water
> 1.0  -> more fluid mass than water
~~~

A secondary volumetric metric is also reported:

~~~text
MJ/L
~~~

because storage volume may matter later.

## Assumptions

Default reference state:

~~~text
initial temperature = 293.15 K
pressure            = 101325 Pa
~~~

A candidate is included only if CoolProp can represent it as a liquid at the initial state and can compute a saturation state at the selected pressure.

This is **not** a prediction of real atmospheric-entry performance. It is a first screening calculation and a way to start using real thermophysical libraries with a transparent physical meaning.

## Install

Use Python 3.12.

~~~bash
python3.12 -m venv .venv
source .venv/bin/activate

python -m pip install --upgrade pip
python -m pip install -r requirements.txt
~~~

## Run

Rank all eligible CoolProp fluids:

~~~bash
python rank_fluids.py
~~~

Show only the top 30:

~~~bash
python rank_fluids.py --top 30
~~~

Change initial temperature and pressure:

~~~bash
python rank_fluids.py \
  --temperature-k 300 \
  --pressure-pa 200000
~~~

Write the full table to CSV:

~~~bash
python rank_fluids.py --csv results.csv
~~~

## What to look for

The first questions are deliberately simple:

1. Where does water rank by MJ/kg?
2. Are any ordinary non-cryogenic liquids clearly better?
3. Does the ranking change strongly with pressure?
4. Which candidates are good by mass but bad by volume?
5. Which results are obviously impractical because the fluid is toxic, corrosive, highly flammable, unstable or otherwise unsuitable?

Those engineering filters come **after** the thermodynamic ranking, not before it.

## Next step

If this experiment is useful, the next extension should still remain simple:

~~~text
sweep pressure
+
sweep initial temperature
+
add basic engineering filters
~~~

Only after that should the project consider molecular simulation, mixtures, AI-guided candidate generation or aerodynamic effects.
