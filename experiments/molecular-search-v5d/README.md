# V5d-1: high-throughput evolutionary molecular search

V5d-1 replaces the hand-written molecular template enumeration used in V5c
with deterministic graph mutations and a hierarchical computational funnel.

The purpose is to search tens to hundreds of thousands of unique molecular
graphs without wasting full entry trajectories on every structure.

## Search hierarchy

The default V5d-1 funnel is:

    V5b calibration + best V5c-2 prescreen structures
      -> deterministic RDKit graph mutations
      -> sanitize / canonicalize / deduplicate
      -> descriptor and explicit family policy
      -> stratified structural sampling
      -> coarse GC-PC-SAFT property screen
      -> score + diversity beam
      -> repeat for multiple generations
      -> full V5c-2 T/P property prescreen
      -> calibrated-domain / exploratory lane
      -> limited full entry evaluations
      -> uncertainty-conservative ranking against Water

The expensive part of the pipeline therefore gets progressively smaller.

## Default scale

The default configuration is intentionally much larger than V5c:

    generations                         6
    children per parent                24
    beam width                       1800
    unique target per generation    22000
    total unique target            120000
    coarse FeOS budget/generation    9000
    full T/P prescreen budget        3500
    rankable entry budget             300
    exploratory entry budget           60

These are upper budgets, not promises that every stage will fill. Invalid,
duplicate, GC-unsupported, non-liquid and property-poor structures reduce the
population.

The run is expected to be materially more expensive than V5c-1/V5c-2. The
correct measure of search seriousness is the funnel size and chemical diversity,
not elapsed time by itself.

## Graph mutations

Every child has explicit parent/operator/generation provenance.

Current mutation operators are:

    add_terminal
    delete_terminal
    change_atom
    change_bond
    insert_atom
    close_ring_5_6
    remove_ring_bond

Allowed mutation elements are:

    C
    N
    O

RDKit sanitization is authoritative for graph valence validity. Disconnected
products are rejected.

The graph generator is deterministic. Each mutation is seeded from:

    global search seed
    parent canonical SMILES
    generation number
    mutation index

The same repository state and configuration therefore produce the same
structural mutation sequence.

## Chemical family policy

V5d does not infer statistical validity solely from elemental composition.

The current rankable families remain:

    alkanes
    alkenes
    aromatics
    cyclic_hydrocarbons
    oxygenated

The exploratory families are:

    aldehydes
    alkynes
    esters
    amines
    mixed_functional

mixed_functional is deliberately broad. A molecule combining an established
functional group with a ring, C=C, aromatic system, or a second distinct
functional motif is exploratory unless a future validation stage explicitly
calibrates that combination.

This prevents a graph mutation such as a cyclic alcohol from inheriting the
ordinary alcohol uncertainty model by accident.

## Initial population

V5d seeds its first generation from:

    successful V5b calibration structures

and, when available:

    ../molecular-search-v5c/search-v5c2-results/prescreen_results.json

The best V5c-2 prescreen structures are preferred as parents.

Structures that have already appeared in:

    search-v5c1-results/generated_candidates.json
    search-v5c2-results/generated_candidates.json

are added to the global seen set when those optional artifacts exist.

## Coarse evolutionary screen

Each generation may produce up to 22k novel valid molecular graphs.

Before FeOS, candidates are sampled across structural buckets based on:

    family
    heavy atom count
    hetero atom count
    ring count

This prevents a large alkane population from consuming the entire EOS budget.

The coarse screen uses six T/P states:

    T = [400, 500] K
    P = [0.101325, 1, 10] MPa

relative to the common 293.15 K / 1 atm storage state.

The coarse property gate is intentionally permissive:

    positive delta_h fraction >= 0.66
    enthalpy merit >= 0.15 of Water

Its job is to remove obvious failures, not replace the full prescreen.

## Evolutionary beam

Candidates that pass the coarse property screen become possible parents for the
next generation.

The beam width is 1800.

By default:

    65% of the beam is selected by thermodynamic property score
    35% is reserved for structural diversity

The diversity portion uses deterministic round-robin sampling across structural
buckets. This avoids collapsing the search to a single homologous series.

## Full prescreen and entry

After all generations, the best/diverse coarse candidates are reduced to up to:

    3500 full property prescreens

using the complete V5c-2 grid:

    T = [350, 400, 450, 500] K
    P = [0.025, 0.101325, 0.3, 1, 3, 10] MPa

The existing V5c-2 rankable/exploratory rules are then applied.

Finally:

    up to 300 rankable candidates
    up to 60 exploratory candidates

receive full current entry evaluations.

Only rankable candidates can produce:

    beats_water_conservative = true

Exploratory candidates remain domain-expansion hypotheses.

## Outputs

V5d writes:

    search-v5d1-results/
      structural_candidates.jsonl.gz
      progress.json
      generation_report.json
      coarse_pass.csv
      full_prescreen.csv
      full_prescreen_results.json
      ranking.csv
      exploratory.csv
      entry_results.json
      water_reference.json
      summary.json
      manifest.json

The compressed JSONL catalog stores the complete accepted structural search
population without turning summary.json into a hundred-megabyte artifact.

## Search-space diagnostics

summary.json reports:

    family distribution
    mutation-operator distribution
    molecular-weight distribution
    approximate Morgan/Tanimoto diversity
    nearest-similarity novelty relative to V5b calibration
    generation-by-generation funnel counts
    coarse rejection counts
    full/entry rejection counts

The diversity statistic uses a deterministic sample of up to 400 searched
structures and reports mean/median pairwise Tanimoto distance.

The same deterministic sample is compared with the V5b calibration population.
The report includes the nearest-reference similarity distribution and fractions
below the 0.45 and 0.25 structural-domain thresholds.

## Run

Update main and run tests first:

    cd ~/DEV/autonomous-materials-lab
    git fetch origin
    git reset --hard origin/main

    cd experiments/entry-evaluator
    source .venv/bin/activate
    python -m unittest discover -s tests -v

Then:

    cd ../molecular-search-v5d

    python run_search_v5d.py \
      --workers 16 \
      --output-dir search-v5d1-results

The local V5b/V5c result directories are reused automatically when present.

During the run, progress.json is overwritten after every completed generation
with the current beam and funnel state. This is a diagnostic checkpoint, not a
resume mechanism.

One compact JSON progress line is also printed after every generation, for
example:

    {"generation": 3, "new_structures": 22000, ...}

The final stdout contains the complete search funnel and the best rankable and
exploratory candidates.

## Interpretation

V5d-1 is the first stage intended to be a genuine high-throughput molecular
search rather than a template-enumeration proof of infrastructure.

It is still a reduced-order discovery stage. Even a conservative Water winner
must be independently validated for properties, chemistry, transport,
materials compatibility and higher-fidelity entry physics before it can be
treated as a credible coolant candidate.
