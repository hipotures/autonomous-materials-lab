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


# V5c-2: property-targeted expansion of the solution space

V5c-2 keeps the V5c-1 uncertainty and Water-ranking contracts but changes the
search strategy in two ways:

1. it expands the molecular templates into additional groups that are present
   in the pinned Sauer/Rehner/Joback parameterization;
2. it evaluates a cheap thermodynamic property grid before the full entry
   trajectory.

The objective is not to generate more molecules indiscriminately. The objective
is to spend full trajectory evaluations only on structures with physically
useful heat-uptake behavior.

## Expanded molecular space

V5c-2 includes the V5c-1 templates plus:

    linear aldehydes
    terminal alkynes
    esters (formates / alkanoates with varied alkoxy partition)
    branched primary alcohols
    primary amines

These directions come directly from the pinned group files:

    CH=O
    C=O
    HCOO
    COO
    C#CH
    OH
    NH2

Ring generation is restricted to 5- and 6-member saturated carbon rings because
the pinned residual model contains dedicated CH/CH2 parameters for r5 and r6,
not a generic arbitrary-ring parameterization.

## Two search lanes

V5c-2 explicitly separates:

    rankable
    exploratory_domain_expansion

A candidate is rankable only when all of the following hold:

    its elements are already represented in successful V5b calibration liquids
    its chemical family is in the configured calibrated-family list
    its structural applicability status is in_domain
    its entry-error uncertainty has a finite conservative upper bound
    it passes the property target gate

With the current V5b calibration, the successful calibrated elements are C/O.
The configured rankable families are:

    alkanes
    alkenes
    aromatics
    cyclic_hydrocarbons
    oxygenated

The newly introduced families:

    aldehydes
    alkynes
    esters
    amines

are exploratory by default.

This is deliberate. A new functional family does not become statistically
rankable merely because a Morgan fingerprint is close to an existing molecule.

Primary amines are an even stronger expansion case: NH2 exists in the pinned
GC-PC-SAFT/Joback model, but nitrogen is not represented by a successful V5b
reference liquid. Their predictions are therefore diagnostic only.

## Thermodynamic pre-screen

Every model-representable room-temperature liquid is evaluated on:

    T = [350, 400, 450, 500] K

and:

    P = [0.025, 0.101325, 0.3, 1, 3, 10] MPa

relative to the common storage state:

    293.15 K
    101325 Pa

For each grid point:

    delta_h(T,P)
      = h(T,P) - h(storage)

V5c-2 reports:

    positive-delta-h fraction
    minimum delta_h
    q25 delta_h
    median delta_h
    maximum delta_h
    median Cp
    storage density
    normal boiling temperature
    latent-enthalpy proxy around the boiling point

The same property grid is calculated for Water with CoolProp.

The primary pre-entry priority metric is:

    candidate q25(delta_h) / Water q25(delta_h)

The lower quartile is used instead of a single favorable state so that a
candidate with weak heat uptake over a material fraction of the T/P envelope is
penalized before trajectory evaluation.

The default hard property gate requires:

    >= 90% of requested grid states with positive cross-pressure delta_h
    q25/median enthalpy merit >= 0.25 of Water

This threshold is intentionally permissive. It removes clearly uncompetitive
or pathological structures without assuming the cheap property metric is an
entry-mass surrogate accurate enough to replace the trajectory.

## Avoiding repeated V5c-1 work

If this artifact exists:

    search-v5c1-results/generated_candidates.json

V5c-2 automatically removes every canonical structure already generated in
V5c-1.

The file is optional. A clean checkout without prior local search results still
runs correctly, but may revisit part of the V5c-1 space.

## Full-entry budget

After the property pre-screen, V5c-2 sends only the best candidates to the full
entry evaluator:

    up to 60 rankable candidates
    up to 12 exploratory candidates

Both limits are configurable.

Rankable candidates retain the V5c-1 conservative score:

    reference upper bound
      = predicted / (1 - u)

for effective u < 1.

Exploratory candidates receive a nominal predicted entry mass only. They never
increment conservative_water_winner_count and never receive a conservative
rank.

## What an exploratory result means

An exploratory candidate below Water is not a discovery claim.

It means:

> this functional family is sufficiently promising under the current
> structure-derived property model that adding independent reference compounds
> for that family should be considered.

The output includes exploratory_family_summary with, for every exploratory
family:

    prescreen pass count
    best property-priority score
    number sent to entry
    best nominal coolant mass
    best nominal ratio to Water

This directly identifies which chemistry should be added to the next V5b
calibration expansion.

## Run

After updating main:

    cd experiments/entry-evaluator
    source .venv/bin/activate
    python -m unittest discover -s tests -v

Then:

    cd ../molecular-search-v5c

    python run_search_v5c2.py \
      --workers 16 \
      --output-dir search-v5c2-results

Default inputs reuse:

    ../property-predictor-v5b/property-v5b2-results/summary.json
    ../property-predictor-v5b/property-v5b21-results/summary.json

and, when present:

    search-v5c1-results/generated_candidates.json

## V5c-2 outputs

    search-v5c2-results/
      generated_candidates.json
      water_property_reference.json
      water_reference.json
      prescreen_results.json
      prescreen.csv
      ranking.csv
      exploratory.csv
      rejections.csv
      production_uncertainty_model.json
      summary.json
      manifest.json

The main stdout fields are:

    water_reference_kg
    generation
    prescreen
    rejection_counts
    ranked_candidate_count
    conservative_water_winner_count
    exploratory_evaluated_count
    exploratory_predicted_below_water_count
    exploratory_family_summary
    best_rankable_candidate
    best_exploratory_candidate

A conservative rankable Water winner remains the strongest screening result.
A strong exploratory result instead tells us where the validation domain should
be expanded next.
