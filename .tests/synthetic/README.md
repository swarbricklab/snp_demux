# Synthetic-data harnesses for PR review

Small, self-contained generators and harnesses used to test PRs where the OneK1K
fixture (`.tests/integration/`) can't exercise the change: an edge case it doesn't
contain, or a component that should be compared old-vs-new. None of these run in
CI. Run them from inside the repo.

| script | exercises | needs |
|---|---|---|
| `dvc_template_deps.sh REF` | `.dvc/template.yaml`: renders it into a throwaway DVC repo, perturbs each input, and checks `dvc status` | dvc (env `dt3`) |
| `make_vireo_outputs.py OUTDIR` | generates vireo outputs (donor_ids / prob_* / prop_ambient / summary) for 5 scenarios: `missing1`, `all_known`, `zero_cells`, `learned_zero`, `missing2` | python3 stdlib |
| `relabel_vireo_equivalence.sh OLD NEW` | `relabel_vireo.R` at two refs, run under a mocked `snakemake` object on the synthetic outputs; diffs every output | singularity + the `r` image |
| `make_doublet_captures.py OUTDIR` | simulates pooled captures with a known doublet rate `d` and donor proportions `p`, by random pairing; only cross-donor doublets are visible. Writes `cell_assignment.tsv` (current + legacy NA-status format) and vireo-like `donor_ids.tsv` | python3 stdlib |
| `doublet_stack_check.py WORKTREE [--reps N]` | runs `doublet_posteriors.py` + `doublet_summary.py` (#59–#61) via their `main()` on N simulated replicates; reports f_x, d_total bias and 95% coverage of `d_total_sd` against the truth | pandas (env `snakemake_7.32.4`) |
| `compare_calls.py GOLDEN FRESH [capture]` | per-method and consensus call agreement between two result dirs (for PRs expected to DIFFER from the golden) | pandas |
| `bcftools_shims.sh DIR` | PATH shims running the cached bcftools image, for bcftools-backed unit tests (#62) and offline `scenario.sh` runs (#64) | singularity |
