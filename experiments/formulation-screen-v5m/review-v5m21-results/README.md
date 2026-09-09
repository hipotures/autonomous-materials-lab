# V5m-2.1: frozen-grid review

Start with summary.json and publication-index.json. No thermodynamic model was rerun.

report-regime-summary-* contains each pair's best sampled composition per T/P.
report-model-regimes-* keeps models separate; report-disagreements-* shows paired
phase fractions at the largest discrepancy. report-pair-* retains every paired
composition/condition, including incomplete comparisons. report-pareto-water-*
adds pure water explicitly to the existing descriptive Pareto objectives.

All report tables and all allowlisted original CSV/JSON files have verified gzip
copies. Original artifacts and their hashes were not changed. The compressed
archives preserve original bytes; these text shards are derived reports.

A 0.1% default reporting margin is only a sensitivity threshold, not calibrated
uncertainty. Per-regime gains remain retrospective model hypotheses; selecting
a favorable regime after seeing results is not new validation. No wetting,
pore-flow, handling, phase or caloric accuracy has been validated by this step.
