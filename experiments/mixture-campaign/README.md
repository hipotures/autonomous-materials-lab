# V5m-3: Repeatable mixture campaigns

One command now discovers pairs, freezes each numerical model, evaluates the
common baseline, automatically refines composition/temperature/pressure edges,
runs every required study, compares revisions and publishes compressed results.
There is no acetone-specific path and no manual selection of the next pair.

This is a campaign-engine refactor with an explicit first physics adapter, not
an all-materials simulator or a newly trained ML model. The initial adapter uses
the existing V5m-2 binary-water ChemSep NRTL and Dortmund UNIFAC implementations.
All earlier experiments, their files and their restart contracts remain intact.

## Run

From the repository root with the existing environment:

```bash
git pull --ff-only
cd experiments/entry-evaluator
source .venv/bin/activate
python -m unittest discover -s tests -p 'test_campaign.py' -v
cd ../mixture-campaign
python run_campaign.py
```

The default paths are relative to this script, independent of your shell's
working directory. The runner neither installs packages nor downloads data.
The existing numerical dependencies remain thermo 0.6.1, chemicals 1.5.2,
CoolProp 8.0.0 and RDKit 2026.3.6. Python/library fingerprints are recorded.

Run exactly the same command after interruption, a new study, a configuration
change, or a bug fix. It automatically reuses compatible cached tasks and
executes new or invalidated tasks. Every invocation publishes a NEW immutable
snapshot rather than refusing an old output directory or overwriting it.
There is no separate `--resume` flag to remember.

```bash
python run_campaign.py
```

Retry recorded infrastructure failures explicitly:

```bash
python run_campaign.py --retry-failures
```

For a deliberate full numerical replay, retaining earlier outcomes:

```bash
python run_campaign.py --recompute
```

`--recompute` bypasses persistent task references but still deduplicates work
inside one invocation. Same-key changed outcomes are reported, not silently
claimed reproducible. Unchanged tasks normally produce identical scientific
hashes; timestamps, cache statistics and publication paths are excluded from that
hash. A separate `--cache-dir` provides an independent cold rebuild as well.

A smoke test can process three catalog pairs:

```bash
python run_campaign.py --limit-pairs 3
```

Smoke snapshots are published but NEVER replace the `latest.json` pointer for
full campaigns. `--workers 8` changes scheduling, not numerical task identity;
each numerical child uses one BLAS/OpenMP thread. Different process counts are
not evidence of physical accuracy. All paths may be overridden explicitly.

## What the first campaign covers

The same automatically enumerated catalog policy as V5m-2 is the default.
Previously encountered CAS identifiers are read from the last published registry
and rechecked, including earlier rejections. New discovery is never restricted
to prior finalists. `additional_cas` can request other known pairs; all normal
identity, model-parameter and scope checks still apply.

The default composition grid has 24 additive mass fractions between 0.0001 and
0.30, outlet temperatures 350/400/450 K, and pressures 101325/200000 Pa.
Inlet temperature is 293.15 K. That is 144 coordinates and 288 model states per
pair. Zero-additive water is an explicit normalization competitor, not a new
measurement. This is an EXPLORATORY fixed-composition thermodynamic scenario,
not a justified TPS temperature/pressure requirement or a physical-test recipe.

Every catalog-eligible pair gets the same baseline and every required study.
Adaptive work is assigned by the same rules to all pairs, not a hard-coded name
or manually chosen finalist list. Defaults allow five rounds, 36 new paired
coordinates per round, and at most 180 new coordinates per pair. With 26 pairs,
this bounds baseline plus adaptive work at 16,848 model states, excluding any
separately requested reference comparisons. Limits are per pair, not contingent
on execution order or which parallel worker finishes first.

### Automatic escalation

Adjacent, collinear sampled edges are inspected along additive mass fraction,
outlet temperature and pressure. Additional midpoints are requested for:

- a change in phases or model feasibility;
- crossing the paired gain-reporting threshold;
- a large difference between model predictions;
- a sufficiently steep gradient in a model-supported gain region.

Logarithmic midpoints are used for pressure and dilute composition, arithmetic
midpoints elsewhere. Round-robin allocation across axes prevents the entire
budget being spent on composition. Existing coordinates are never recomputed
merely because they participate in another study. New points stay inside the
predeclared campaign scope; a best point at its boundary is explicitly flagged.
The engine does not silently extend the operating requirements to obtain a win.

Stops are reported as `no_triggered_sampled_edges`, `point_budget`, `max_rounds`
or a planner failure. The first means only that the defined rules do not request
more sampled-edge work. It is NOT a proof of global convergence, absence of
unsampled phase islands, or identification of a thermodynamic binodal. A budget
stop remains unresolved. Invalid phase states and failed numerical jobs remain
in denominators; no water/ideal-mixture fallback fills missing predictions.

## Required studies

The standard suite currently contains:

| Study | Output |
| --- | --- |
| `regimes` | Local paired gains, per-T/P reports, baseline completeness and worst-case ratio, best-point boundary flags |
| `phase_boundaries` | Remaining triggered brackets and unresolved edge count, not certified phase boundaries |
| `model_disagreement` | Comparisons at the SAME composition, including phase fractions and stored caloric evidence |
| `reference_audit` | Independent source-declared enthalpy comparisons when supplied; explicit gaps otherwise |

Phase balance and pure-water endpoint guards from V5m-2 run in the numerical
adapter. No mixture is a physical winner or production-rankable as a result of
this suite. It still does not simulate pore-scale transport, contact angle,
critical heat flux, sequential evaporative fractionation, chemical reactions,
precipitation, suspensions, tank mass or a flight trajectory. Those mechanisms
need appropriate new adapters/studies and data. They are not invented by the
orchestrator. Global Pareto uses only the common BASELINE; uneven adaptive grids
are not used to give one pair an easier global score than another.

## Dependency graph and cache correctness

Each task key hashes its kind, effective inputs, implementation and environment,
plus dependency task keys AND outcome hashes. All cached JSON is normalized and
rejects NaN/Infinity. Objects are compressed and content-verified on reuse.

```text
catalog (policy + current/historical identifiers + library data)
  -> pair-specific parameter projection
     -> frozen model (parameters + fixed correlation-selection context)
        -> exact state (model + mass fraction + temperature + pressure)
           -> observation bundle
              -> adaptive plan -> additional exact states -> new bundle
              -> declared studies and their dependencies
                 -> qualifications and immutable publication
```

Pair-specific projection intentionally does not depend on the entire catalog
artifact key. Adding another compound therefore does not invalidate existing
pair states. Both parent key and outcome hash matter: a failed task that succeeds
on retry must invalidate analyses that consumed the failure, even though its
original request key is unchanged.

Changes to a model implementation or numerical library invalidate its dependent
states. Changes to analysis or plugin studies do not invalidate solver states.
Changing composition-grid resolution keeps already computed matching states.
Changing the fixed correlation context creates a NEW frozen model and dependent
states. Numerical library-tree hashing is conservative: an installation/data
change can invalidate more than the minimum theoretically necessary subset.
This prioritizes correctness over maximum cache reuse. Edits to README files,
Git publication commits and execution worker counts do not change task keys.

The fixed `model.correlation_temperatures_k` and declared temperature envelope
choose pure correlations ONCE per model context, independently of coarse or
adaptive samples. Every child verifies those selected methods before evaluation.
A wider temperature scope can change that context and correctly require a new
numerical series. Coarse and refined points never mix different correlation
choices without the change entering their dependency identity.

The store records immutable outcome objects and append-only invocation events.
Small task refs are mutable cache pointers; retry does not erase the old outcome.
A lock prevents concurrent campaigns from writing the same cache or result root.
Numerical subprocess timeouts kill the child; each missing requested state is
recorded as a failure. A crashed publication has no completed snapshot or latest
pointer; restarting reuses completed numerical tasks. Cache is local and can grow;
no automatic cache deletion or destructive cleanup is implemented.

Old V5m-2 rows are NOT automatically imported into this cache. Their run-level
signature is not an equivalent per-state/frozen-context proof. The FIRST campaign
therefore recalculates the catalog deliberately; later campaigns reuse the new
cache. All previously published reports stay available for comparison.

## Add a study without changing the orchestrator

Create a trusted local Python module exposing `register_studies(registry)`:

```python
from campaign_studies import Study


def run_new_check(dependencies, settings):
    rows = dependencies['observations']['rows']
    # Real implementation goes here; never label predictions as measurements.
    return {'examined_state_count': len(rows), 'physical_validation': False}


def register_studies(registry):
    registry['new_check'] = Study(
        name='new_check',
        version='1',
        requires=('observations',),
        config_keys=(),
        run=run_new_check,
    )
```

Add its module name to `plugins` and `new_check` to `studies` in the campaign
YAML, then use the same command. Every eligible pair, including historical ones,
gets the new required study. Dependencies are topologically ordered and cycles
are rejected. Failed required studies prevent a complete-suite qualification.
Plugins cannot replace existing registered studies silently.

Study callbacks must be deterministic functions of declared dependencies and
settings. They must not read hidden files, clocks, networks or mutable globals.
Their source file is hashed automatically; declare helper source/data paths in
`source_files` and additional third-party distributions in `packages`. Configure
all relevant settings through `config_keys`. Package versions AND installed-file
hashes are recorded for declared plugin packages. Undeclared dependencies cannot
be inferred reliably by this engine; test every new module's invalidation rules.
Plugins are explicitly trusted local code, not a sandbox for downloaded scripts.

## Optional measured enthalpy references

`references` may point to a reviewed JSON file with schema
`mixture-campaign-measurements-v1` and an `observations` array. Each record needs
`id`, `source`, `evidence_kind: measured`, `water_cas: 7732-18-5`, `additive_cas`,
`composition_basis: mass`, `additive_mass_fraction`, `inlet_temperature_k`,
`outlet_temperature_k`, `pressure_pa`, `property: delta_h_j_kg`, `unit: J/kg` and
positive `value`. Optional nonnegative `standard_uncertainty` is retained.

These must describe the SAME fixed overall composition at inlet and outlet.
Within-scope coordinates at the campaign inlet are evaluated automatically in a
separate reference-prediction bundle, not inserted into the discovery score.
Outside-scope or different-inlet observations remain explicitly unmatched.
Reference values never enter a numerical prediction request. No model fitting,
calibrated uncertainty claim, independent-source authentication or training-set
overlap claim is made. No file means `missing_external_references`, not a pass.
This stage does not implement an automatic literature/measurement downloader.

## Publish and review revisions

Each completed run goes to `campaign-v5m3-results/run-<UTC>-<id>/`.
`latest.json` identifies the latest full campaign. Snapshots contain:

- summary, manifest, registry, changes and a complete publication index;
- immutable scientific results, task graph and per-pair dossiers in gzip;
- small UTF-8 CSV shards with all sampled coordinates and statuses;
- common-baseline Pareto with explicit water, and task execution/cache counts.

Gzip is round-trip checked; text shards are at most 32 KiB by default. The GitHub
text connector is not assumed to decode archives. Original/raw observations and
prior snapshots are never rewritten. A published snapshot is not a local cache.
The revision report compares every pair, including additions and pairs removed
by the current scope, and lists changed qualification/result fields. Changes to
required-suite implementations are visible even when numerical values coincide.

From this directory after completion:

```bash
git add -- campaign-v5m3-results/
git diff --cached --stat -- campaign-v5m3-results/
git commit --only -m "Add mixture campaign results" -- campaign-v5m3-results/
git push origin HEAD
```

Only approved report/archive filenames are admitted by `.gitignore`; task blobs,
workers, logs, locks and partial-build trees stay local. Alternate output paths
need an explicit ignore-policy update. Do not use forced add or `git clean -fdx`.
Exit 1 means numerical failures or incomplete required studies; exit 2 denotes
invalid configuration/dependencies or a fatal integrity error. Budget-limited
adaptive studies may finish numerically with exit 0 but remain explicitly
unresolved in their reports. Physical safety approval is never inferred.

## Tests and implementation sources

Architecture tests use synthetic fixtures and test clean versus incremental
rebuilds, cache corruption, dependency invalidation, retry, extension backfill,
immutable revisions, budget limits, compressed output and Git ignore behavior.
Two optional integration tests exercise the actual installed thermodynamics
stack. Missing libraries are reported as skips, never numerical verification.

The first adapter preserves the model assumptions documented in
`../formulation-screen-v5m/V5m-2.md`. Primary API references:

- https://thermo.readthedocs.io/thermo.flash.html
- https://thermo.readthedocs.io/thermo.phases.html
- https://docs.python.org/3/library/subprocess.html

No universal synthesis, hardware control or scientifically exhaustive conclusion
is claimed by a completed computational campaign.
