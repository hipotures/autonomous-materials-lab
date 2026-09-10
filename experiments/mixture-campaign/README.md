# Reusable mixture campaigns

The default entry point now runs **V5m-4: broad known-mixture search on a shared
atmospheric-entry heat-demand profile**, rather than a 26-pair property grid.

```bash
python run_campaign.py
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
