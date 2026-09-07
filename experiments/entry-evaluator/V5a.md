# V5a: Known binary liquid-mixture property layer

## Purpose

V5a is the first stage of the actual liquid-candidate track.

The objective is narrow:

> Make the entry evaluator accept a real binary mixture through the same
> coolant interface as a pure fluid, preserve pure-component endpoints, and
> score composition-dependent thermodynamics without trajectory-specific
> mixture code.

V5a is a property-pipeline benchmark, not yet an optimizer over hypothetical
molecules.

## Property-provider contract

The coolant model consumes a common provider:

    enthalpy_j_kg(T, P)
    phase(T, P)
    state(T, P)
    saturation_at_pressure(P)
    metadata()

Legacy pure-fluid configurations remain valid:

    coolant:
      coolprop_name: Water

A mixture can use:

    coolant:
      property_provider:
        type: coolprop
        backend: HEOS
        components: [Water, Ethanol]
        composition_basis: mole
        fractions: [0.5, 0.5]

The entry evaluator does not need to know whether the candidate is a pure fluid
or mixture.

## Pure endpoints

Zero-fraction components are removed before the CoolProp state is created.

For Water/Ethanol:

    x_water = 0.0 -> HEOS:Ethanol
    x_water = 1.0 -> HEOS:Water

The endpoints are therefore true pure-fluid reference calculations.

## First benchmark: Water / Ethanol

The default study is:

    mixture-v5a.yaml

It evaluates:

    x_water = 0.0, 0.1, ... 0.9, 1.0

Water/Ethanol is used as a known binary-mixture benchmark. This is not a claim
that ethanol is a promising entry coolant.

All compositions use the same nominal physical entry assumptions and the
Brandis-Johnston heating backend. Mixture chemistry/decomposition is disabled
because V5a does not yet contain a validated composition-dependent chemistry
model.

## Required CoolProp backend

Do not run the Water/Ethanol benchmark with released CoolProp 8.0.0.

CoolProp issue #1900 is specifically about incorrect or failed Water/Ethanol
mixture calculations. The maintainer reports the fix in the development branch.

V5a pins the exact upstream revision:

    d5b0cfb51cd9a9343284cc5af8ebd6a8bd0eecc0

Install it with:

    python -m pip install -r requirements-v5a.txt

The ordinary:

    requirements.txt

remains pinned to CoolProp 8.0.0 so V1-V4 keep their original dependency
contract.

The V5a runner checks CoolProp's embedded git revision before launching any
property or trajectory calculations. The manifest records both version and git
revision.

References:

- https://github.com/CoolProp/CoolProp/issues/1900
- https://github.com/CoolProp/CoolProp/issues/3243

## PT and saturation evaluation

V5a uses two separate CoolProp AbstractState objects:

    property_state
        -> PT property evaluation

    saturation_state
        -> PQ bubble/dew evaluation

For mixtures, the saturation state is used only to determine whether a PT point
is unambiguously outside the two-phase envelope at the same pressure:

    T < min(Tbubble, Tdew)
        -> specify liquid on property_state

    T > max(Tbubble, Tdew)
        -> specify gas on property_state

    otherwise
        -> leave phase unspecified and use the full mixture PT flash

The PT state and saturation state are separate because CoolProp issue #3243
documents stale VLE state/cache corrupting later imposed-phase PT root
selection when both operations share one AbstractState.

This guard is used only with the pinned development revision where the
imposed-phase PT path has the upstream fix. V5a does not invent a phase boundary
or interpolate pure-component boiling points; the phase decision comes from the
mixture bubble/dew calculation at the same composition and pressure.

## Regression target for Water/Ethanol

CoolProp 8.0.0 produced pathological states for some mid-range compositions,
including dense metastable roots and enthalpies of order -1e8 to -1e9 J/kg.

V5a regression tests therefore check physical quantities rather than relying on
CoolProp's textual phase label.

For:

    x_water = 0.4, 0.5, 0.6
    T = 500 K
    P = 25 kPa

the test requires:

    density < 10 kg/m3
    h(500 K, 25 kPa) - h(293.15 K, 1 atm) > 0
    abs(h(500 K, 25 kPa)) < 10 MJ/kg

These are deliberately broad guards. Their purpose is to reject the known
wrong dense-liquid root without encoding a new thermodynamic model in this
repository.

## Property validation

Before trajectory scoring, each composition is evaluated at the storage state
and on a small T/P grid.

Recorded properties include:

    phase
    density
    specific enthalpy
    cp
    viscosity
    thermal conductivity
    bubble temperature
    dew temperature

Missing transport properties are reported as coverage gaps and are not turned
into favorable candidate scores.

Continuity between adjacent composition grid points is diagnostic only.
Non-ideal mixtures may be non-monotonic.

## Entry coupling

The quantity directly coupled into the coolant mass balance is:

    delta_h =
        h_candidate(T_exit, P_surface)
        - h_candidate(T_storage, P_storage)

and:

    mdot = Qdot_coolant / delta_h

Thus composition-dependent enthalpy already affects every active coolant step.

The following are recorded but not yet coupled to system physics:

- density -> tank volume / tank mass;
- viscosity -> porous pressure drop;
- thermal conductivity -> internal heat transfer;
- surface tension;
- decomposition / oxidation chemistry;
- mixture-dependent blowing effectiveness.

V5a is therefore a thermophysical screening result, not a complete system
ranking.

## Run

From experiments/entry-evaluator:

    python -m pip install -r requirements-v5a.txt

    python -m unittest discover -s tests -v

    python run_mixture_sweep.py \
      --stability-algorithm 1 \
      --workers 16 \
      --output-dir mixture-v5a-results

The default benchmark executes 11 entry trajectories after property validation.

Algorithm 1 is the CoolProp Michelsen path and is the V5a default. Algorithm 0
is retained only as an explicit diagnostic comparison; a result that succeeds
only under the legacy path is solver-sensitive and is not promoted as robust.

## Outputs

    mixture-v5a-results/
      manifest.json
      property_validation.csv
      property_state_grid.csv
      continuity_report.json
      entry_results.csv
      summary.json

The summary reports:

- pure-endpoint validation;
- liquid-storage coverage;
- property-grid coverage;
- transport-property coverage;
- coolant mass for every composition;
- mass ratio to the configured reference;
- failed candidates with explicit reasons;
- best evaluated composition.

"Best" means best among the discrete evaluated grid points under V5a
assumptions. It is not a continuous global optimum.

## V5b

V5b should introduce the first property model that is not simply a direct
CoolProp lookup.

Recommended sequence:

    known pure fluids / mixtures
        -> hold out selected compositions
        -> predict properties
        -> compare against trusted data/backend
        -> attach uncertainty
        -> allow optimizer-driven candidate proposals

The same PropertyProvider interface is intended to accept that surrogate
without modifying the entry evaluator.
