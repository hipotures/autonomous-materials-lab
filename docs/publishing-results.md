# Publishing computational results

## V5m-2

The formulation experiment uses an exact-file allowlist for directories named
`predictive-v5m2*-results` directly under `experiments/formulation-screen-v5m`.
The default `predictive-v5m2-results` directory is covered.

Publish completed scientific outputs: the summary, provenance manifest, artifact
hashes, catalog, grid/refinement plans, prediction CSV/JSON, tradeoff tables,
Pareto tables and validation shortlist. `predictions.json` is retained because
it contains model metadata in addition to the CSV's state rows.

Everything else inside these run directories stays local, including `workers/`,
`checkpoints/`, `source-cache/`, `progress.json`, `.run.lock`, temporary files and
unexpected new outputs. Earlier experiment result directories remain ignored.
The root `.gitignore` also excludes local environments, `.env` files and common
runtime/test caches. This policy is not a secret-content scanner or a size limit.
Inspect artifacts before publishing them, especially in a public repository.

From the formulation experiment directory, after the run has finished:

```bash
git pull --ff-only
git add -- predictive-v5m2-results/
git diff --cached --stat -- predictive-v5m2-results/
git diff --cached --name-only
```

After checking the staged files:

```bash
git commit --only -m "Add V5m-2 predictive grid results" -- predictive-v5m2-results/
git push origin HEAD
```

`--only` scopes this commit to the result directory and leaves unrelated staged
changes out of it. Use the actual branch name when requesting a review after
publishing on a branch other than `main`. No forced add, forced push, index reset,
archive creation or experiment rerun is required. Do not use `git clean -fdx`:
it would remove ignored local results, environments and restart caches.

Before a first commit, `git status --short --untracked-files=all --
predictive-v5m2-results/` lists the visible output files. `git status --short
--ignored -- predictive-v5m2-results/` also shows ignored paths. Files already
tracked remain tracked regardless of ignore rules, so always inspect staged paths.

The run's `manifest.json` records the computation's original code revision. Do
not rewrite it to the later publication commit. `artifact-hashes.json` describes
the original local output set; it may mention local-only files such as
`progress.json`. That does not imply those files were uploaded. A clone of these
published reports is not a restart checkpoint and cannot by itself resume a run.

A future output filename or a different directory naming scheme needs an explicit
ignore-policy update before publication. Renamed runs should get new directories;
keep previously published results immutable rather than recomputing over them.
