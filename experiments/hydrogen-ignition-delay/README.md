# Hydrogen / Air Ignition-Delay Sweep

This experiment asks one narrow kinetic question:

> After hot hydrogen mixes with oxygen-containing air, how long does a homogeneous gas parcel take to undergo a substantial ignition event?

It follows the previous thermodynamic and equilibrium screens but does **not** model a heat shield, boundary layer, injection geometry, diffusion, or CFD.

## Why this is different from the equilibrium experiment

The equilibrium calculation answered:

> How much chemical energy could be released if the H2 / air mixture reaches equilibrium?

It did not answer:

> How quickly can that reaction happen?

This experiment adds finite-rate chemical kinetics.

The output is therefore a characteristic homogeneous ignition time, not a wall-reaction time and not a complete reentry prediction.

## Mechanism

Default:

~~~text
h2o2.yaml
~~~

This is the H2 / O2 mechanism distributed with Cantera and includes N2 as a diluent.

Cantera explicitly describes its distributed mechanisms as example/convenience data. The results are appropriate for this first kinetics sanity check, not yet a validated reentry-chemistry model.

## Reactor model

Default:

~~~text
adiabatic
constant pressure
homogeneous
zero-dimensional
~~~

The script uses:

~~~python
ct.IdealGasConstPressureReactor(...)
~~~

A constant-volume option is also available:

~~~bash
python ignition_delay.py --reactor constant-volume
~~~

The two models represent different idealizations and should not be interpreted as a boundary-layer calculation.

## Initial composition

Simplified dry air:

~~~text
O2:1, N2:3.76
~~~

Hydrogen richness is controlled by equivalence ratio:

~~~text
phi < 1   lean
phi = 1   stoichiometric
phi > 1   hydrogen-rich / oxygen-limited
~~~

Default values:

~~~text
0.5  1  2  4  8  16
~~~

The high-phi cases are important because the conceptual near-wall layer may contain much more hydrogen than locally available oxygen.

## Default sweep

Temperatures:

~~~text
500 750 1000 1250 1500 K
~~~

Pressures:

~~~text
0.1 0.3 1 3 bar
~~~

Equivalence ratios:

~~~text
0.5 1 2 4 8 16
~~~

Total:

~~~text
5 x 4 x 6 = 120 cases
~~~

## Ignition definition

The script records the reactor trajectory while limiting the temperature change of an integration advance.

A case is classified as a substantial ignition event only if:

~~~text
Tmax - Tinitial >= 200 K
~~~

by default.

For an ignited case it reports two timing markers:

### tau_dTdt

Time corresponding to the largest sampled temperature-rise rate:

~~~text
max(dT/dt)
~~~

This is the primary ignition-delay metric in this experiment.

### tau_OH_peak

Time of the largest sampled OH mole fraction within the simulated trajectory.

This is a secondary marker.

If the required temperature rise is not reached before the configured integration horizon, the result is:

~~~text
no_ignition_within_limit
~~~

This means exactly that. It does not mean that the mixture can never ignite.

## Time horizon

Default:

~~~text
10 s
~~~

Change it with:

~~~bash
python ignition_delay.py --max-time-s 1
python ignition_delay.py --max-time-s 100
~~~

The low-temperature cases may require a much longer horizon than the high-temperature cases.

## Installation

~~~bash
python3.12 -m venv .venv
source .venv/bin/activate

python -m pip install --upgrade pip
python -m pip install -r requirements.txt
~~~

## Run

~~~bash
python ignition_delay.py
~~~

Write all cases to CSV:

~~~bash
python ignition_delay.py --output ignition_delay.csv
~~~

Use a finer temperature sweep:

~~~bash
python ignition_delay.py \
  --temperatures-k 600 650 700 750 800 850 900 950 1000 \
  --pressures-bar 0.1 0.3 1 3 \
  --phi 1 2 4 8 16 \
  --output ignition_delay.csv
~~~

Change the ignition threshold:

~~~bash
python ignition_delay.py --ignition-rise-k 100
~~~

## Output

The main columns are:

~~~text
T0_K
P0_bar
phi
tau_dTdt_s
tau_OH_peak_s
Tmax_K
delta_T_K
OH_peak
status
~~~

For convenience the table also prints the primary delay in milliseconds and microseconds.

## Interpretation

The quantity we ultimately care about later is a comparison between chemical and transport timescales:

~~~text
tau_chem
vs.
tau_residence / tau_mixing
~~~

If chemistry is much slower than transport out of the near-wall region, substantial oxidation may occur farther from the wall.

If chemistry is much faster than transport, a hydrogen-rich cooling layer may react close to the protected surface.

This script computes only the chemical side of that comparison.

## Important limitations

This experiment assumes a perfectly homogeneous premixed gas parcel.

It does not include:

- diffusion of H2 into air;
- finite mixing time;
- hot atomic oxygen or other reentry radicals;
- shock-layer nonequilibrium chemistry;
- wall catalysis;
- surface reactions;
- heat transfer to a wall;
- pressure/temperature histories along a trajectory;
- radiation;
- CFD.

Therefore the result should be called a **homogeneous ignition-delay screening value**, not a reentry ignition prediction.

## References

Cantera 3.2 documentation:

- zero-dimensional reactor networks;
- constant-pressure ideal-gas reactor;
- ignition-delay examples;
- distributed `h2o2.yaml` mechanism.


## V4 dense lookup table

The entry-evaluator V4 chemistry gate uses a denser table and log-space
interpolation across temperature and pressure. Generate the default lookup with:

```bash
python generate_v4_table.py --workers 16
```

The generated `ignition_delay_v4.csv` and manifest are local run artifacts and
are intentionally ignored by Git. The default table spans 600-1100 K,
0.001-30 bar and phi = 2, 4, 8, 16. See
`../entry-evaluator/V4.md` for how the table is consumed.
