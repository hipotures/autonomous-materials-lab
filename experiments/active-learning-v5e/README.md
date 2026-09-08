# V5e-1: Active-learning calibration expansion

V5e-1 turns the V5d exploratory population into a concrete calibration plan.

The problem after V5d-1 is no longer a lack of generated molecules. The search
produces many thermodynamically interesting structures, but most of them are
outside the empirical V5b calibration domain.

V5e-1 therefore selects a compact set of structures for which independent
reference property data should be acquired next.

The selection objective combines:

    thermodynamic property merit
    novelty relative to the current V5b calibration
    structural diversity relative to already selected targets
    a small bonus for candidates that also received favorable nominal V5d entry scores

The output is a calibration acquisition set, not a new coolant ranking.

## Inputs

By default V5e reads:

    ../molecular-search-v5d/search-v5d1-results/full_prescreen_results.json
    ../molecular-search-v5d/search-v5d1-results/entry_results.json
    ../property-predictor-v5b/property-v5b2-results/summary.json

Only V5d candidates that passed the full thermodynamic prescreen and are in the
exploratory lane are considered by default.

## Acquisition score

For each candidate V5e calculates:

    property merit percentile
    nearest Tanimoto similarity to V5b calibration
    calibration novelty = 1 - nearest similarity
    optional nominal entry percentile

The static acquisition score is:

    property_merit_weight * property_percentile
      + calibration_novelty_weight * calibration_novelty
      + entry_merit_weight * entry_percentile

Selection is greedy. At each step an additional diversity term is added:

    selected_diversity_weight
      * minimum Tanimoto distance to already selected targets

The default weights are:

    property merit       0.38
    calibration novelty  0.32
    entry merit          0.10
    selected diversity   0.20

The entry term is intentionally smaller because only a small subset of V5d
exploratory candidates received full entry evaluation and those predictions are
not calibrated in the expanded domain.

## Family balance

The default acquisition set contains 50 targets.

To prevent the very large mixed-functional population from dominating:

    minimum_per_family = 4
    maximum_family_fraction = 0.30

The minimum is applied when enough candidates from that family exist.

## Novelty bands

Targets are labeled relative to the current V5b structural domain:

    near_domain  nearest similarity >= 0.45
    bridge       0.25 <= nearest similarity < 0.45
    frontier     nearest similarity < 0.25

These labels describe structural distance only. They do not imply calibrated
uncertainty.

## Acquisition roles

Every selected target is also tagged as:

    exploit
    bridge
    frontier

This helps separate candidates selected mainly for strong property merit from
those selected mainly to expand structural coverage.

## Reference resolution is explicit

V5e does not silently treat a FeOS prediction as ground truth.

Every selected target carries canonical SMILES, InChIKey when available, molecular formula and molecular weight to support identity resolution. It initially has:

    reference_status: unresolved

The generated:

    reference-resolution.yaml

must be filled only after:

    molecular identity is independently verified
    a trusted reference property source is identified
    the reference source is independent of the predictor being calibrated

Preferred reference sources are currently:

    CoolProp
    NIST
    peer-reviewed EOS / high-quality property data

A candidate may be removed from the acquisition set if no sufficiently reliable
reference data can be obtained.

## Structural domain what-if

V5e reports a hypothetical structural-domain expansion:

    before
    after_if_all_selected_are_validated

For the complete exploratory pool it calculates the fraction whose nearest
reference similarity would exceed:

    0.25
    0.45

after adding all selected targets as hypothetical calibration anchors.

This is only a planning calculation.

Selecting a target does not expand the calibrated uncertainty domain. The domain
changes only after trusted reference properties are acquired and the blind V5b
validation/calibration pipeline is rerun.

## Outputs

    active-learning-v5e1-results/
      acquisition-pool.csv
      calibration-targets.csv
      acquisition-set.json
      domain-expansion-whatif.json
      reference-resolution.yaml
      v5b-expansion-template.yaml
      summary.json
      manifest.json

`calibration-targets.csv` is the compact review table.

`acquisition-set.json` preserves all selection metadata.

`reference-resolution.yaml` is the handoff for source resolution.

`v5b-expansion-template.yaml` is deliberately marked:

    executable: false

It is not a valid benchmark until independent reference identities and property
sources have been resolved.

## Run

After updating main:

    cd experiments/entry-evaluator
    source .venv/bin/activate
    python -m unittest discover -s tests -v

Then:

    cd ../active-learning-v5e

    python run_acquisition.py \
      --output-dir active-learning-v5e1-results

The V5d and V5b input paths can also be overridden explicitly.

## Interpretation

V5e-1 answers:

> Which small set of molecules should receive the next expensive reference-data
> effort so that the calibrated predictor expands toward the most promising and
> structurally underrepresented regions found by V5d?

The next stage after V5e-1 is not another search run. It is reference resolution
and V5b calibration expansion using the resolved subset of these targets.
