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

is the machine-readable handoff into V5e-2. V5e-2 now automates the normal
resolution path. A human only needs to inspect records that remain unresolved,
have incomplete thermophysical support or fail cross-source consistency checks.

The same reference gate still applies:

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


# V5e-2: Automated reference resolution

V5e-2 removes the manual dataset-building step from the active-learning loop.

It consumes the V5e-1 `reference-resolution.yaml` acquisition handoff and
attempts to resolve every target automatically:

    V5e-1 calibration targets
      -> PubChem exact identity lookup
      -> CAS / InChIKey consistency
      -> exact CoolProp identity match
      -> NIST Chemistry WebBook corroboration
      -> storage-state property extraction
      -> source-consistency gate
      -> V5b-3 calibration handoff

The expected workflow is therefore not to search for 50 compounds by hand.
Human review is restricted to the residual unresolved/review queue.

## Identity contract

The locally generated InChIKey remains the primary identity key.

PubChem PUG REST is queried by exact InChIKey. The returned record must agree
with the acquisition target on:

    exact InChIKey
    molecular formula
    molecular weight within configured tolerance

PubChem identity or computed metadata alone is never accepted as a
thermophysical calibration reference.

CAS Registry Numbers are extracted from PubChem synonyms only when the CAS
checksum is valid.

## Independent property reference

The predictor being calibrated is:

    FeOS GC-PC-SAFT + Joback

V5e-2 explicitly forbids using that prediction path as its own reference.

The current automatic calibration-grade backend is CoolProp. A CoolProp fluid
is accepted only through an exact persistent-identifier match:

    exact InChIKey
    or
    exact CAS Registry Number

Name/fuzzy matching is intentionally not used.

For an accepted fluid V5e-2 records:

    critical temperature
    critical pressure
    acentric factor
    normal boiling temperature
    storage density
    storage Cp
    storage enthalpy
    outlet enthalpy
    enthalpy window
    latent heat of vaporization

For low-boiling fluids the resolver automatically raises storage pressure above
the saturation pressure using the configured margin, subject to the configured
maximum storage pressure. A target is not automatically ingested if the
configured storage temperature is at/above the critical temperature, the
required pressure exceeds the limit, the resolved phase is not liquid, or
required reference properties are missing.

## NIST Chemistry WebBook

NIST Chemistry WebBook SRD 69 is queried by verified CAS number.

V5e-2 records the phase-change page when available and also attempts the NIST
fluid-system isobaric endpoint at the selected storage state. Comparable
CoolProp/NIST values are checked against:

    maximum_cross_source_relative_difference

The default is 5%.

A target exceeding that tolerance is routed to:

    resolved_conflict_review

instead of being silently accepted.

NIST-only records are retained as useful reference evidence but are not yet
automatically converted into V5b holdouts because the current V5b reference
runner consumes CoolProp identities. A later property-provider abstraction may
promote sufficiently complete NIST/literature-only targets without weakening
the blind-reference contract.

## Reproducible source cache

Every HTTP response used by the run is cached under:

    reference-v5e2-results/source-cache/

The cache stores the exact response bytes plus URL/SHA-256 metadata. Re-running
the same URL from the same output tree therefore does not silently substitute
new remote content.

The runner also records package versions, input hashes and the source contract
in the manifest.

## Outputs

    reference-v5e2-results/
      resolved-targets.csv
      unresolved-targets.csv
      reference-audit.csv
      reference-properties.json
      provenance.json
      calibration-expansion.yaml
      v5b3-calibration-input.yaml
      summary.json
      manifest.json
      source-cache/

`resolved-targets.csv` contains only candidates that pass the automatic
calibration gate.

`unresolved-targets.csv` is the human-review queue.

`reference-audit.csv` is a compact audit table covering every attempted
candidate.

`reference-properties.json` preserves normalized identity, property,
cross-source and failure evidence.

`v5b3-calibration-input.yaml` contains only exact-identity, independent
CoolProp-backed targets. It is marked executable only when at least the
configured minimum number of targets were resolved.

## Run

First generate the V5e-1 acquisition set if it does not already exist:

    cd experiments/active-learning-v5e

    python run_acquisition.py \
      --output-dir active-learning-v5e1-results

Then run automated reference resolution:

    python run_reference_resolution.py \
      --output-dir reference-v5e2-results

A small network smoke test can be run with:

    python run_reference_resolution.py \
      --limit 3 \
      --output-dir reference-v5e2-smoke-results

The normal full run should not use `--limit`.

## Interpretation

V5e-2 answers:

> Which V5e-1 acquisition targets can be independently identified and supplied
> with trustworthy reference-property support automatically, and which residual
> cases actually require human review?

A V5e-2 target marked `usable_for_v5b_calibration=true` is still not a
rankable coolant candidate. It is only eligible to enter the next calibration
experiment.

The next stage is V5b-3:

    original V5b calibration
      + V5e-2 resolved reference subset
      -> blind prediction
      -> expanded-domain validation
      -> uncertainty recalibration
      -> applicability-domain update

Only after that loop is rerun may V5d exploratory candidates become rankable.


# V5e-2.1: Multi-source empirical reference resolution

V5e-2 showed that the active-learning targets are mostly identifiable molecules
but are outside CoolProp's compact fluid catalog. V5e-2.1 broadens the
thermophysical reference path without relaxing the independence contract.

It consumes:

    reference-v5e2-results/reference-properties.json

and reuses the exact PubChem identity/CAS work already completed by V5e-2. It
does not repeat the slow network crawl.

The local resolution path is:

    exact V5e-2 identity
      -> verified CAS candidates
      -> chemicals 1.5.2 constant-source inventory
      -> thermo 0.6.1 T-dependent source inventory
      -> explicit empirical allowlist
      -> explicit predicted-method deny gate
      -> property coverage / source-family scoring
      -> V5b-3 property-calibration handoff

## Reference policy

V5e-2.1 never trusts the default method chosen by `chemicals` or `thermo`.
It enumerates all available methods and selects only methods explicitly listed
in `config-v5e21.yaml`.

Accepted constant sources include reviewed or compiled reference data such as:

    HEOS / REFPROP-class data
    IUPAC
    Matthews
    CRC
    NIST WebBook
    CAS Common Chemistry

Accepted T-dependent sources include selected tabulations and correlations
fitted to reference data, for example:

    HEOS_FIT
    Zabransky liquid Cp
    NIST WebBook Shomate
    VDI tabular / PPDS
    DIPPR / Perry
    Antoine / Wagner reference correlations
    CRC / Poling reference constants

Structure-only and generalized estimators are not accepted as calibration
evidence. The hard deny gate includes Joback, Wilson-Jasperson, Dadgostar-Shaw,
Rowlinson, Lastovka, corresponding-state vaporization estimators, COSTALD,
Rackett-family methods, generic EOS, Lee-Kesler, Ambrose-Walton and related
prediction methods.

Every available-but-rejected method remains visible in:

    source-method-audit.csv

so the gate is auditable rather than implicit.

## Property packet

For each CAS candidate V5e-2.1 attempts to build:

    critical_temperature_k
    critical_pressure_pa
    normal_boiling_temperature_k
    vapor_pressure_pa(T)
    liquid_density_kg_m3(T)
    liquid_cp_j_kg_k(T)
    latent_heat_vaporization_j_kg(T)

The temperature grid is shared with V5b:

    293.15, 350, 400, 450, 500 K

No extrapolation is requested from the thermo property objects. A method only
contributes grid points where the library reports it as valid and evaluation
succeeds.

For molecules carrying multiple valid CAS aliases, all CAS candidates are
evaluated and the strongest independently supported packet is selected
deterministically.

## Calibration gate

The default property-calibration gate requires:

    Tc present
    Pc present
    at least 2 dynamic properties
    at least 4 total properties
    at least 2 independent source families
    no forbidden selected method

A higher-quality packet requires:

    at least 3 dynamic properties
    at least 6 total properties
    at least 3 source families

This is a property-calibration gate, not an entry-trajectory reference gate.

A molecule can therefore be:

    ready_for_v5b3_property_calibration = true
    entry_reference_ready = false

That distinction is intentional. V5b-3 can use these empirical anchors to
expand the property/applicability model while retaining the original
CoolProp-backed holdouts for entry-level error calibration.

## Outputs

    reference-v5e21-results/
      calibration-ready.csv
      partial-or-unresolved.csv
      reference-audit.csv
      source-method-audit.csv
      empirical-reference-properties.json
      v5b3-property-calibration-input.yaml
      summary.json
      manifest.json

`empirical-reference-properties.json` is the complete machine-readable
reference packet.

`source-method-audit.csv` records both accepted and rejected methods.

`v5b3-property-calibration-input.yaml` contains only candidates that pass the
property-reference gate. It is executable once the configured minimum number of
anchors is reached.

## Run

V5e-2 must already have completed successfully enough to create
`reference-properties.json`.

Then:

    cd experiments/entry-evaluator
    source .venv/bin/activate
    python -m unittest discover -s tests -v

    cd ../active-learning-v5e

    python run_reference_resolution_v5e21.py \
      --output-dir reference-v5e21-results

This stage is local data lookup and correlation evaluation. It should not have
the long HTTP wait profile of V5e-2.

## Interpretation

V5e-2.1 answers:

> How many identity-resolved active-learning targets have enough independent,
> non-predictive thermophysical evidence to become empirical property anchors
> for the next calibration-domain expansion?

The next stage is V5b-3. It must consume the property packets without pretending
that a property-only anchor is a complete entry-level reference fluid.
