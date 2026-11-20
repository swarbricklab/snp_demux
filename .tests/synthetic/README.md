# Synthetic-data harnesses for PR review

Small, self-contained generators and harnesses used to test PRs where the OneK1K
fixture (`.tests/integration/`) can't exercise the change: an edge case it doesn't
contain, or a component that should be compared old-vs-new. None of these run in
CI. Run them from inside the repo.

| script | exercises | needs |
|---|---|---|
| `dvc_template_deps.sh REF` | `.dvc/template.yaml`: renders it into a throwaway DVC repo, perturbs each input, and checks `dvc status` | dvc (env `dt3`) |
| `make_vireo_outputs.py OUTDIR` | generates vireo outputs (donor_ids / prob_* / prop_ambient / summary) for 4 scenarios: `missing1`, `all_known`, `zero_cells`, `missing2` | python3 stdlib |
| `relabel_vireo_equivalence.sh OLD NEW` | `relabel_vireo.R` at two refs, run under a mocked `snakemake` object on the synthetic outputs; diffs every output | singularity + the `r` image |
