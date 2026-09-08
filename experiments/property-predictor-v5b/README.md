# V5b-1: structure-derived property prediction

V5b-1 is the first project stage that evaluates a coolant candidate without
reading its known thermophysical properties during prediction.

The prediction path is:

    SMILES
      -> Sauer group decomposition
      -> Rehner 2023 heterosegmented GC-PC-SAFT
      + Joback/Reid ideal-gas heat capacity
      -> FeOS Helmholtz EOS
      -> PropertyProvider
      -> entry evaluator

CoolProp is not imported or queried by the FeOS provider. In the holdout
benchmark it is opened only after the prediction artifact has been written, as
an independent known-fluid reference.

## Why the ideal-gas model is included

GC-PC-SAFT supplies the residual Helmholtz contribution. Total enthalpy and heat
capacity also require an ideal-gas contribution. V5b-1 therefore composes the
FeOS GC-PC-SAFT residual model with Joback/Reid heat-capacity groups derived
from the same SMILES.

This prevents the entry evaluator from ranking candidates using residual
enthalpy alone.

## Pinned model data

The repository vendors the exact FeOS v0.10.1 group files used by V5b-1:

    ../entry-evaluator/parameters/v5b/sauer2014_smarts.json
    ../entry-evaluator/parameters/v5b/rehner2023_hetero.json
    ../entry-evaluator/parameters/v5b/joback1987.json

Each run records SHA-256 hashes of these files.

Upstream sources:

- FeOS v0.10.1 parameters/pcsaft/sauer2014_smarts.json
- FeOS v0.10.1 parameters/pcsaft/rehner2023_hetero.json
- FeOS v0.10.1 parameters/ideal_gas/joback1987.json

The first file defines structural group matching. The second provides the
GC-PC-SAFT residual parameters. The third supplies Joback ideal-gas heat
capacity coefficients.

## Holdout protocol

The first holdouts are ethanol and acetone.

For each candidate the predictor receives:

    SMILES
    storage T/P
    outlet-temperature limit

It does not receive the CoolProp fluid name or reference property values.

The runner:

1. predicts the storage state, a T/P property grid, saturation temperature and
   nominal-entry coolant mass;
2. writes predictions_before_reference.json;
3. only then evaluates the known CoolProp reference;
4. compares predicted and reference results.

Absolute enthalpies from different backends are not compared because their
reference zeros may differ. The caloric comparison uses:

    delta_h(T,P)
      = h(T,P) - h(T_storage,P_storage)

within each backend and compares those heat-uptake differences.

The entry-level metric compares required coolant mass for the same candidate,
mission and terminal condition.

## Run

From the repository root:

    cd experiments/entry-evaluator
    python -m pip install -r requirements-v5b.txt
    python -m unittest discover -s tests -v

Then:

    cd ../property-predictor-v5b
    python run_holdout.py --output-dir property-v5b1-results

Outputs:

    property-v5b1-results/
      predictions_before_reference.json
      property_comparison.csv
      entry_comparison.csv
      summary.json
      manifest.json

## Acceptance criterion

V5b-1 succeeds when at least one held-out known liquid completes:

    SMILES
      -> structure decomposition
      -> predicted thermodynamic surface
      -> PropertyProvider
      -> entry evaluator

and the runner reports the error against the subsequently revealed reference.

study_complete=true means all configured holdouts produced comparable terminal
entry results. It does not mean the predictor is accurate enough for material
selection.

## Limits

V5b-1 is deliberately not yet:

- a molecular generator;
- a trained ML model;
- a calibrated uncertainty model;
- a mixture predictor;
- a transport-property predictor;
- a decomposition or reactivity model;
- a system-level tank/feed model.

The error distribution from the holdouts is empirical validation evidence.
V5b-2 should broaden the holdout set, define an applicability domain and turn
validation error into calibrated candidate-specific uncertainty before inverse
design begins.
