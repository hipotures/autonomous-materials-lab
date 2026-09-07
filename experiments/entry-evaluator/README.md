# Low-Fidelity Earth Entry Evaluator (V3)

This experiment is the first end-to-end evaluator for the working-fluid discovery track.

Instead of manually selecting one temperature and pressure, it propagates a complete atmospheric-entry trajectory and evaluates the fluid along the evolving environment.

## Current scope

V1 separates local stagnation heating from integrated forebody power and verifies time-step/grid convergence. V2 adds deterministic physical-sensitivity screening and energy-weighted validity diagnostics. V3 adds a selectable Brandis-Johnston 2014 Earth-entry heating backend and a matched cross-model comparison against the legacy Sutton-Graves / Tauber-Sutton backend.

See [V1](V1.md) for numerical verification, [V2](V2.md) for physical-sensitivity screening and [V3](V3.md) for the heating-backend comparison. Existing configurations without a `surface` section use the uniform V0 area model; the supplied `config.yaml` selects a cosine profile. The default single-run angle stays at -8 degrees, while the V1 verification benchmark uses -11 degrees.

## Primary question

> For the same vehicle and entry trajectory, how much working fluid is required to keep a representative stagnation-region wall below its thermal limit?

The default hard case starts at:

```text
altitude = 200 km
velocity = 15 km/s
flight-path angle = -8 deg
```

All of these values are configuration parameters.

## Architecture

```text
entry state
    |
    v
NRLMSIS 2.1 atmosphere table
    |
    v
planar spherical-Earth trajectory
    |
    +--> drag / dynamic pressure / surface-pressure estimate
    |
    +--> selectable stagnation heating backend
         |-- legacy: Sutton-Graves + Tauber-Sutton
         |-- V3: Brandis-Johnston 2014
    |
    v
angular wall zones + zero-area stagnation probe
    |
    v
minimum required local coolant heat removal
    |
    v
CoolProp enthalpy window
    |
    +--> optional ignition-delay lookup / fixed chemistry limit
    |
    v
coolant mass flow + total consumed mass
```

A single trajectory is sequential because state at time `t + dt` depends on the state at `t`.

Independent cases are parallelized with `ProcessPoolExecutor` in `run_batch.py`.

## Environment

Recommended:

```text
Python 3.12
CoolProp 8.0.0
pymsis 0.12.0
```

Install:

```bash
python3.12 -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
python -m pip install -r requirements.txt
```

`pymsis` is configured with explicit F10.7 / F10.7a / Ap values in `config.yaml`, so the evaluator does not need to fetch current space-weather inputs during a run.

## Single run

```bash
python run_entry.py
```

Outputs:

```text
entry_history.csv
entry_summary.json
```

The history contains the full time series. The summary contains peak values, minimum altitude, minimum Knudsen number, terminal/skip status, and total coolant consumed.

For hyperbolic or shallow skip trajectories the run terminates with `status: atmospheric_exit` after the vehicle descends and then climbs back through the initial altitude. It no longer integrates thousands of seconds into space.

## Parallel batch

On a 16-core CPU:

```bash
python run_batch.py --workers 16
```

The default batch compares hydrogen, water, methane, and ammonia and applies a four-value entry-angle matrix, producing 16 independent trajectories. Each worker executes one complete sequential trajectory. No nested parallelism is used.

A generic Cartesian matrix can be added with dotted configuration paths, for example:

```yaml
matrix:
  entry.flight_path_angle_deg: [-5, -7, -9, -11]
  entry.velocity_km_s: [12, 15]
```

That example produces every angle / velocity combination for every fluid case.

## Atmosphere and rarefied-flow handling

NRLMSIS 2.1 is precomputed on an altitude grid for each process and then interpolated during integration.

The evaluator estimates mean free path and Knudsen number using the vehicle nose radius as the characteristic length.

V0 applies a screening transition:

```text
Kn <= 0.01   continuum correlations fully active
Kn >= 0.1    continuum heating correlations disabled
between      log-linear transition
```

This is not DSMC and is not a validated transitional-flow model. It only prevents continuum stagnation-heating correlations from being applied at 200 km as if the gas were dense continuum flow.

## Trajectory model

The state contains altitude, velocity, flight-path angle, and downrange. The equations are a planar, spherical-Earth entry model with drag, optional lift-to-drag ratio, and gravity variation with radius. RK4 uses a configurable fixed time step.

V0 does not model Earth rotation, winds, guidance, angle-of-attack schedules, 6-DOF attitude dynamics, or skip-entry guidance.

## Convective heating

The stagnation-point Sutton-Graves Earth correlation is used:

```text
q_conv = k * sqrt(rho / Rn) * V^3
k = 1.74153e-4
```

with SI inputs and output in W/m2.

## Radiative heating

V0 implements the Earth Tauber-Sutton form:

```text
q_rad = 4.736e8 * Rn^a * rho^1.22 * f(V)   [W/m2]
a = 1.072e6 * V^-1.88 * rho^-0.325
```

with the published Earth `f(V)` table from 9 to 16 km/s.

The summary reports both the fraction of heating steps and the fraction of integrated radiative energy inside the nominal Tauber-Sutton envelope:

```text
V:   10-16 km/s
rho: 6.66e-5 .. 6.31e-4 kg/m3
Rn:  0.3 .. 3.0 m
```

The V0 result is not high-fidelity radiative shock-layer analysis.

## Wall model

Each angular zone and the stagnation probe have an independent lumped wall state per unit area:

```text
C_areal * dTwall/dt =
    q_external - q_reradiation - q_backface - q_coolant
```

The controller requests the minimum coolant heat flux required to prevent the next thermal step from crossing the configured wall-temperature setpoint.

## Coolant model

For each active zone, CoolProp computes (exact state cache, without rounded first-visitor values):

```text
delta_h = h(T_exit, P_surface) - h(T_storage, P_storage)
mass_flux = q_coolant / delta_h
mass_flow_zone = mass_flux_zone * zone_area
mass_flow_total = sum(mass_flow_zone)
```

By default V0 uses a **fixed trajectory mass** while integrating the required coolant mass. This is deliberate: every candidate fluid is first scored on the same vehicle trajectory, so a fluid is not rewarded or penalized by changing the trajectory while it is being compared.

Set `vehicle.couple_coolant_mass_to_trajectory: true` only for later system-level studies where the carried coolant inventory is known and should change vehicle dynamics.

Likewise, `coolant.available_mass_kg: null` means "score the required mass without an artificial tank-cap failure." Give it a finite value only when testing a specific carried inventory.

V0 deliberately gives no aerodynamic blowing credit to the coolant. External heating is calculated as if the transpiration film did not reduce incoming heat flux.

## Chemistry

Two modes are supported.

### fixed_limit

Use a configured maximum outlet temperature. The default H2 case uses 750 K as a conservative first ceiling motivated by the separate H2/air kinetics sweep.

### ignition_csv

Point `chemistry.ignition_csv` at a CSV generated by:

```text
../hydrogen-ignition-delay/ignition_delay.py
```

The evaluator selects the nearest pressure / equivalence-ratio curve and limits outlet temperature using:

```text
ignition_delay >= residence_time * ignition_safety_factor
```

No Cantera integration is performed inside the trajectory loop.

## What V1 reports

Important outputs include:

```text
coolant_used_kg
peak_total_heat_flux
peak_radiative_heat_flux
peak_dynamic_pressure
peak_deceleration
peak_wall_temperature
peak_coolant_flow
peak_required_injection_pressure
minimum_ignition_delay
status (terminal_velocity / terminal_altitude / atmospheric_exit / failed / max_time)
vehicle_incident_heat_mj
peak_vehicle_heating_power_w
vehicle_coolant_heat_mj
vehicle_radiative_incident_heat_mj
vehicle_radiative_valid_heat_mj
vehicle_radiative_energy_valid_fraction
vehicle_radiative_energy_fraction
wall_energy_relative_residual
rarefied_heating_disabled_time_s
```

## What V0 cannot prove

Terminal status is not a physical-validity flag or flight qualification. The model does not yet contain CFD, DSMC, a reacting boundary layer, porous-media flow, film-cooling effectiveness, multilayer wall conduction, ablation, structural stresses, pump/tank mass, real injector geometry, radiation-flow coupling, or catalytic wall chemistry.

A strong candidate should be escalated to higher-fidelity models rather than accepted as a final design.

## References used for V0 correlations

- Sutton, K. and Graves, R. A., *A General Stagnation-Point Convective-Heating Equation for Arbitrary Gas Mixtures*, NASA TR R-376.
- Tauber, M. E. and Sutton, K., *Stagnation-Point Radiative Heating Relations for Earth and Mars Entries*, Journal of Spacecraft and Rockets 28(1), 1991.
- NASA entry sizing studies using Sutton-Graves and Tauber-Sutton engineering correlations.
- `pymsis` / NRLMSIS 2.1 documentation for density, temperature, and species output.


## V2 physical-sensitivity study

Run the paired 64-sample screening study on 16 CPU workers:

```bash
python run_physical_sensitivity.py --workers 16 --output-dir physical-v2
```

The default run evaluates 65 matched physical scenarios (64 Latin-hypercube samples plus the nominal point) for hydrogen, water, methane and ammonia: 260 trajectories total. See [V2.md](V2.md) for ranges, outputs and interpretation limits.


## V3 heating-model comparison

The evaluator now accepts:

```yaml
heating:
  backend: legacy
```

or:

```yaml
heating:
  backend: brandis_johnston_2014
```

To replay the exact V2 physical scenarios under both backends:

```bash
python run_heating_model_comparison.py --workers 16 --output-dir heating-model-v3
```

The default comparison executes 65 scenarios × 4 fluids × 2 backends = 520 trajectories and reports model-to-model changes in heat load, peak heat flux, coolant mass, H2/water ranking and correlation-validity coverage. See [V3.md](V3.md).
