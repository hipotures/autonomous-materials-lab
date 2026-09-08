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


# V5b-2: applicability domain and calibrated uncertainty

V5b-2 keeps the V5b-1 blind-prediction contract and adds a larger structural
holdout set plus a reusable uncertainty calibration model.

The benchmark is:

    benchmark-v5b2.yaml

It currently contains 29 cases:

- 25 intended in-model holdouts spanning alcohol, ketone, linear/branched
  hydrocarbons, aromatics, cyclic hydrocarbons, ether and alkenes;
- 4 deliberate model-domain probes: methanol, water, ammonia and cyclopropane.

The domain probes are not favorable or unfavorable coolant examples. Their
purpose is to verify that the predictor can refuse structures that cannot be
decomposed by the pinned Sauer/Rehner/Joback group set.

## Structural applicability domain

V5b-2 uses RDKit Morgan fingerprints:

    radius = 2
    fpSize = 2048

with Tanimoto similarity.

For every candidate, structural support is measured only against the
calibration split. The default screening thresholds are:

    in-domain: nearest similarity >= 0.45
               and at least 2 neighbors >= 0.25

    edge:      nearest similarity >= 0.25
               but insufficient calibrated neighborhood

    out-of-domain:
               nearest similarity < 0.25
               or the GC parameterization cannot represent the structure

Only in-domain candidates receive a screening-supported uncertainty flag.
Numerical uncertainty values for edge/out-of-domain cases are diagnostics only.

## Calibration split

Candidate IDs are deterministically hashed with a fixed seed and split into:

    70% calibration
    30% evaluation

The reference backend is still not evaluated until all prediction artifacts
have been written.

The calibration set supplies observed errors for:

    required coolant mass
    storage density
    median delta_h error over the property grid
    median Cp error over the property grid
    saturation-temperature error

## Local uncertainty scale

For a target structure, V5b-2 selects the k most similar calibrated molecules
(default k=5) and forms a similarity-weighted local absolute-error scale.

For every calibration molecule, this scale is calculated leave-one-out. The
ratio:

    observed absolute relative error / local structural error scale

is the nonconformity score.

A finite-sample split-conformal order statistic at the configured target
coverage (default 90%) becomes a multiplicative calibration factor.

For a new in-domain candidate:

    relative uncertainty
      = conformal factor
        * local structural error scale

This gives candidate-specific uncertainty that grows when nearby known
structures were difficult for the predictor.

This is intentionally simpler than a learned graph model. V5b-2 first
establishes whether empirical structural neighborhoods can support useful
uncertainty calibration before adding a trainable model.

## Run

From experiments/entry-evaluator:

    uv pip install -r requirements-v5b.txt
    python -m unittest discover -s tests -v

Then:

    cd ../property-predictor-v5b

    python run_calibration.py \
      --workers 16 \
      --output-dir property-v5b2-results

Outputs include:

    predictions_before_reference.json
    property_comparison.csv
    entry_comparison.csv
    candidate_uncertainty.csv
    uncertainty_model.json
    summary.json
    manifest.json

The key summary fields are:

    study_complete
    entry_comparable_count
    in_domain_evaluation_count

    applicability_domain:
      in_domain_count
      edge_count
      out_of_domain_count

    uncertainty:
      metric_calibrations
      evaluation_in_domain_observed_coverage

## V5b-2 acceptance gate

The default run succeeds only if:

- at least 12 candidates produce comparable entry results;
- an entry-mass uncertainty calibration factor can be constructed;
- at least 3 candidates remain in the independent evaluation split;
- at least 2 evaluation candidates are structurally in-domain.

study_complete does not mean 90% coverage has been scientifically established.
The evaluation coverage is reported explicitly and must be inspected.

The calibration set is still small and chemically nonuniform. V5b-2 therefore
provides a screening uncertainty model, not a metrological guarantee.

## What comes after V5b-2

If the observed evaluation coverage is useful and the error bounds separate
promising from non-promising candidates, the next step is V5c inverse search.

If coverage is poor, V5b should be expanded first with either:

- more reference liquids in sparse structural regions;
- a learned residual model;
- a better descriptor / distance metric;
- chemistry-family-specific uncertainty calibration.

A molecular generator should not be allowed to exploit regions marked
out-of-domain.


# V5b-2.1: robustness before inverse search

V5b-2.1 strengthens the V5b-2 uncertainty result without rerunning the expensive
property and trajectory calculations.

It adds three checks:

    failure taxonomy
    32 deterministic calibration/evaluation splits
    leave-one-family-out structural challenges

## Failure taxonomy

Prediction outcomes are classified separately as:

    success
    structure_out_of_domain
    storage_state_infeasible
    property_model_failure
    entry_evaluator_failure
    prediction_pipeline_failure

This distinction matters for probes such as methanol. A structure can be
representable by GC-PC-SAFT while still violating the configured liquid-storage
contract. That is an application-state failure, not structural OOD.

The benchmark now declares an expected prediction outcome for every deliberate
probe. A mismatch blocks study_complete.

## Repeated split robustness

The original V5b-2 split is only one deterministic 70/30 partition. V5b-2.1
recomputes the domain and uncertainty model over 32 deterministic seeds using
the same revealed holdout errors.

For every metric it reports:

    per-split coverage min / median / max
    pooled eligible count
    pooled covered count
    pooled coverage
    distribution of in-domain evaluation counts
    valid split fraction

This is a split-sensitivity test, not 32 new independent experiments. The same
finite holdout population is reused.

## Family challenge

The benchmark assigns successful holdouts to broad challenge families:

    alkanes
    alkenes
    aromatics
    cyclic_hydrocarbons
    oxygenated

Each family is removed entirely from calibration and evaluated against the
remaining families.

If the applicability-domain model marks a held-out family as edge or
out-of-domain, that is acceptable behavior. If it marks members in-domain, the
reported uncertainty is checked against their observed entry error.

## V5c readiness

V5b-2.1 reports v5c_ready separately from study_complete.

Default readiness gates require:

    >= 75% repeated splits with enough in-domain evaluation cases
    pooled entry-error coverage >= 80%
    median per-split entry coverage >= 80%
    family challenges with >=2 in-domain cases to meet the same entry coverage
    zero failure-taxonomy mismatches

These thresholds are screening gates, not proof of 90% statistical coverage.

## Run without recomputing trajectories

If property-v5b2-results already exists:

    python run_robustness.py \
      --input-dir property-v5b2-results \
      --output-dir property-v5b21-results

This reads the already revealed V5b-2 errors and produces:

    property-v5b21-results/
      summary.json
      repeated_splits.json
      family_challenge.json

To regenerate the blind benchmark with the new failure taxonomy first, run:

    python run_calibration.py \
      --workers 16 \
      --output-dir property-v5b2-results-new

then point run_robustness.py at that directory.
