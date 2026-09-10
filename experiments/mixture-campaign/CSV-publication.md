# CSV publication and repacking

All CSV report shards and full tables are published as `.csv.gz`. Shard byte
limits apply to the uncompressed content. Compression is deterministic within the
same Python/zlib environment and every gzip output is verified by decompression.
Each index entry carries compressed and uncompressed byte counts and SHA-256.

`preview-*.md` contains at most the first eight input-order records and 12 KiB,
with explicit truncation labels. A preview is not a selected scientific ranking
and does not substitute for the complete compressed dataset. The GitHub text
connector is not assumed to decode gzip. Small JSON summaries, manifests, and
publication indexes remain plain text; existing full JSON archives stay gzip.

## Existing results: no numerical rerun

Run these commands only after this change has been applied to the repository:

```bash
cd /home/user/DEV/autonomous-materials-lab
python -m unittest discover -s experiments/entry-evaluator/tests -p 'test_csv_publication.py' -v
cd experiments/mixture-campaign
python repack_csv_reports.py --results-dir mission-v5m4-results --check
python repack_csv_reports.py --results-dir mission-v5m4-results
```

The latest published snapshot is selected from `latest.json`. Use `--all-runs`
only to explicitly migrate all run snapshots in that result directory. A result
lock prevents concurrent publication. All source/destination preflight checks run
before any migration writes. Unindexed CSV files and conflicting destinations
are refused. Source report bytes are checked against their published hashes;
original CSVs are removed only after their gzip replacements have been verified
and the updated index has been atomically committed. Rerunning after an interrupted
cleanup safely removes remaining verified originals. Repeating a completed
migration is a no-op.

This is an explicit packaging revision, not an in-place change to scientific
results. The original index is retained as
`publication-index.before-csv-repack.json.gz`. The new index records original CSV
hashes and a `csv_repack` record. Scientific JSON, the result digest, `latest.json`,
models and caches remain unchanged. Git history is not rewritten.

Review and stage only the result path after repacking:

```bash
cd /home/user/DEV/autonomous-materials-lab
git add -A -- experiments/mixture-campaign/mission-v5m4-results/
git diff --cached --stat -- experiments/mixture-campaign/mission-v5m4-results/
```

`git add -A` stages verified original-file removals as well as new gzip files.
Do not use `git clean`, force-add ignored files, or rerun the numerical campaign
just to change publication formatting. Future campaign runs use the gzip writer
automatically through the existing shared `campaign_publication.table` interface.
