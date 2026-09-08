# V5c-1: constrained inverse molecular search

V5c-1 is the first project stage that generates molecular structures that were
not supplied as known coolant identities and ranks them through the full
structure-derived entry pipeline.

It is deliberately constrained. It is not an unrestricted generative-chemistry
model.

The search path is:

    deterministic molecular templates
      -> canonical SMILES
      -> remove known V5b reference structures
      -> descriptor constraints
      -> GC-PC-SAFT + Joback representability
      -> liquid storage-state check
      -> saturation support
      -> V5b structural applicability domain
      -> screening entry-mass uncertainty
      -> nominal entry evaluator
      -> conservative ranking against Water

## Required inputs

V5c-1 requires the outputs of the completed V5b validation stages:

    ../property-predictor-v5b/property-v5b2-results/summary.json
    ../property-predictor-v5b/property-v5b21-results/summary.json

The V5b-2 summary supplies observed entry-mass errors for known reference
molecules. V5c rebuilds a production screening uncertainty model using all of
those revealed errors.

The V5b-2.1 summary is a hard authorization gate. By default V5c refuses to run
unless:

    v5c_ready == true

This keeps molecular generation downstream of the uncertainty-robustness check.

## Candidate generation

The generator is deterministic and records provenance for every candidate.

Current operations include:

    linear alkane chain length
    single-methyl alkane branching
    terminal-alkene chain length
    primary-alcohol chain length
    ketone chain length / carbonyl position
    ether carbon partition
    linear alkyl substitution on benzene
    cycloalkane ring size
    methyl-cycloalkane substitution

These templates intentionally remain close to chemistry already represented in
the pinned Sauer/Rehner/Joback parameterization and the V5b calibration set.

A generated structure is not automatically considered valid. The property
provider and applicability-domain gates remain authoritative.

## V5c-1 search domain

The default search restricts candidates to:

    neutral molecules
    C/O elements only
    40 <= molecular weight <= 180 g/mol
    3 <= heavy atoms <= 12
    at most one ring

The storage contract is intentionally stricter than the broad V5b benchmark:

    T_storage = 293.15 K
    P_storage = 101325 Pa
    predicted phase must be liquid

V5c-1 does not allow an unpenalized pressure vessel to make a volatile fluid
appear favorable. Pressurized-storage design belongs in a later system model.

## Known-structure exclusion

The canonical SMILES in the V5b calibration summary are removed before search.

Therefore a V5c ranked candidate is novel relative to the V5b reference
population used to calibrate the current predictor. This is not a claim of
novelty in chemical literature or patent space.

## Applicability domain

A candidate must pass the same Morgan/Tanimoto structural-domain logic used in
V5b:

    in_domain_similarity = 0.45
    edge_similarity = 0.25
    minimum_neighbors = 2

By default:

    edge          -> reject
    out_of_domain -> reject

The generator is not allowed to gain score by moving into chemistry where the
current error model has no support.

## Production screening uncertainty

After V5b-2.1 validates split robustness, V5c uses all available known-fluid
entry errors to rebuild the screening entry-error model.

For a candidate:

    u_entry
      = conformal factor
        * local structural error scale

V5b defines relative error against the revealed reference:

    abs(predicted - reference) / reference <= u_entry

Therefore the finite conservative upper bound is:

    conservative_coolant_kg
      = predicted_coolant_kg
        / (1 - multiplier * u_entry)

for effective uncertainty below 1. If the effective uncertainty is at least 1,
the current error definition gives no finite upper bound and the candidate is
rejected from ranking.

The default multiplier is 1.

Candidates are sorted by conservative_coolant_kg, not by the raw predicted
coolant mass.

## Water reference

Water is recalculated in the same V5c run using CoolProp and the same mission,
cooling area, storage T/P and 500 K outlet ceiling.

A candidate has:

    beats_water_predicted = true

when the raw prediction is below the Water mass.

The stronger screening result is:

    beats_water_conservative = true

which requires the uncertainty-penalized score to remain below Water.

A conservative Water win is only a trigger for higher-fidelity validation. It
is not flight-level evidence and does not include tank/feed-system penalties,
transport validation, decomposition chemistry, porous-flow coupling or blowing.

## Run

The V5b environment already contains the required dependencies.

From the repository after updating main:

    cd experiments/entry-evaluator
    source .venv/bin/activate
    python -m unittest discover -s tests -v

Then:

    cd ../molecular-search-v5c

    python run_search.py \
      --workers 16 \
      --output-dir search-v5c1-results

The default input paths point to the result directories produced in the
adjacent property-predictor-v5b experiment.

## Outputs

    search-v5c1-results/
      generated_candidates.json
      ranking.csv
      rejections.csv
      summary.json
      water_reference.json
      production_uncertainty_model.json
      manifest.json

The stdout and summary report:

    water_reference_kg
    raw / novel generation counts
    ranked_candidate_count
    rejection_counts
    predicted_water_winner_count
    conservative_water_winner_count
    best_candidate

For each ranked candidate the result includes:

    SMILES
    generation provenance
    molecular descriptors
    structural-domain metadata
    screening entry uncertainty
    predicted coolant mass
    conservative coolant mass
    Water ratios
    storage density / implied coolant volume
    complete entry summary
    property-model provenance

## Interpretation

V5c-1 answers a narrow question:

> Within the chemistry that the current structure-derived model can represent
> and for which V5b provides structural error support, can a deterministic
> inverse search produce a previously untested room-temperature liquid whose
> uncertainty-penalized working-fluid mass beats the Water reference?

If the answer is yes, the next stage should validate the top candidates against
independent property data and add chemistry/system penalties before expanding
the molecular search space.
