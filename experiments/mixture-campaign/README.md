# Reusable mixture campaigns

## Default: V5m-4.1 full-inventory continuation

```bash
python run_campaign.py
```

The default now processes every eligible locally discovered water/additive pair
in checkpointed batches. It preserves historical candidates, applies the same
required studies to every pair, diagnoses incomplete coverage, and evaluates a
zero-additive control with each binary model. Matching V5m-4 numerical cache
entries remain reusable. No new thermochemical or pore-flow model is implied.

Read [V5m-4.1.md](V5m-4.1.md) for migration, tests, partial-run budgets,
new-study backfill, scientific limitations and compressed publication.
The default configuration is `mission-campaign-v41.yaml`; results are written to
`mission-v5m41-results`. This is distinct from the previous result directory.

```bash
python run_campaign.py --status
python run_campaign.py --max-batches 1
```

The first command only reads the latest published summary. The second processes
at most one batch of pending pairs; it does not change the scientific inventory.
Repeat the normal command to continue. Do not remove `.campaign-cache`.

## Preserved V5m-4 single-run protocol

```bash
python run_campaign.py --mode mission-v4
```

The following original V5m-4 notes refer to that explicit mode.

The preserved mission-v4 entry point runs **V5m-4: broad known-mixture search on a shared
atmospheric-entry heat-demand profile**, rather than a 26-pair property grid.

```bash
python run_campaign.py --mode mission-v4
```

Read [V5m-4.md](V5m-4.md) for configuration, assumptions, tests and publication.
The default configuration is `mission-campaign.yaml`, with up to 500 model-backed
pairs and 40 initial mass compositions per pair. The actual eligible population
is discovered locally and reported. Neither a 500-pair count nor a physical
cooling improvement is assumed.

The pipeline performs broad catalog discovery, fixed-entry-demand screening,
automatic composition refinement and full discrete-load replay for finalists.
It reuses the common content-addressed Store and writes immutable compressed
snapshots with small connector-readable CSV reports. Repeat the same command to
reuse compatible calculations; use `--retry-failures` for failed tasks.

The score is a frozen-species, fixed-load working-fluid mass estimate. It is not
a coupled reactive, porous-flow or ablative TPS prediction. Missing physical
models, unsupported states and model/reference discrepancies remain explicit.
No experimental winner or production rankable fluid is declared.

## Preserved property-only protocol

V5m-3.1 remains available without changing its configuration, numerical backend
or old results:

```bash
python run_campaign.py --mode properties
```

Its original runner is preserved in `run_property_campaign.py`. See
[V5m-3-properties.md](V5m-3-properties.md) for the original protocol, using the
`--mode properties` switch with its documented commands, and
[V5m-3.1.md](V5m-3.1.md) for water controls and evidence-aware refinement.

Mission results use `mission-v5m4-results`; property results continue to use
`campaign-v5m3-results`. Both share `.campaign-cache` with separate task kinds.
