# Job costs: `job-costs`

`job-costs` reads the PBS "Resource Usage" block at the end of job logs and
writes CPU time, memory and service units (SU) per job, per rule and per
capture. It needs numpy and pandas only (`jobcosts.py`).

## What it reads

Any file whose text holds the block PBS appends at the end of a job:

```
======================================================================================
                  Resource Usage on 2026-01-02 03:04:05:
   Job Id:             123456.gadi-pbs
   Project:            p00
   Exit Status:        0
   Service Units:      0.64
   NCPUs Requested:    32                  CPU Time Used: 00:07:23
   Memory Requested:   31.25GB               Memory Used: 6.13GB
   Walltime Requested: 24:00:00            Walltime Used: 00:00:36
   JobFS Requested:    100.0MB                JobFS Used: 0B
======================================================================================
```

That covers Snakemake job logs written through a PBS profile and qxub's
`.log` files.

- **Each block is one job.** A file with two blocks gives two rows.
- **The Snakemake header** (`rule NAME:`, `input:`, `output:`,
  `wildcards:`, `reason:`, `resources:`) gives the rule, the capture (the
  wildcard `capture_id` or `capture`, set by `--capture-keys`), the first
  input and output, and the requested threads, memory and runtime. It is read
  from the log, or from the `.err` file beside it when the log has none.
- **Two jobs under one file name.** A later job that writes the same log
  name normally replaces the file (see Limits). When two such jobs end close
  together, the second job's block can land in the first job's file instead,
  under the first job's header. The header's `tmpdir=/jobfs/<job id>` names
  its own job, so that block takes the header and any other is marked
  `header_lost`, with no capture.
- **A job killed by PBS** (`PBS: job killed: walltime ... exceeded limit`)
  carries the reason in `pbs_kill`.
- **A block missing fields** is kept, with the missing fields empty and
  named in `partial_block`. A block cut before its closing line is kept too,
  marked `truncated`. Such a block gets no `charge_basis` and no
  `su_expected`, and `rule_costs.tsv` and `capture_costs.tsv` count partial
  blocks per group, so a missing figure is never added as 0. A total over a
  group with no known value is empty.
- A file without a header takes its file name, up to the first `.`, as its
  rule.

## How SU are charged

SU = queue rate × the larger of the CPUs requested and (memory requested /
a node's memory × a node's CPUs) × the walltime used. Unused walltime is not
charged. On `normal` (48 CPUs and 192 GB per node) one CPU carries 4 GB at no
charge: a 1-CPU job asking for 16 GB is charged as 4 CPUs.

| Queue | SU per CPU hour | Node |
|---|---|---|
| `normal` | 2.0 | 48 CPUs, 192 GB |
| `normalsr` | 2.0 | 104 CPUs, 500 GB |
| `normalsl` | 1.5 | not recorded: give `--node CPUS:GB` |
| `normalbw` | 1.25 | not recorded: give `--node CPUS:GB` |
| `express` | 6.0 | not recorded: give `--node CPUS:GB` |

The log does not name the queue. `su` is what PBS charged, as printed;
`su_expected` recomputes it for `--queue` (default `normal`) as a check, and
`su_rate` is `su` over walltime times the charged CPUs.

## Running it

```bash
python -I tools/demux_tools job-costs --source cohort=logs/joblogs --out costs/
python -I tools/demux_tools job-costs --source run1=logs/joblogs --source qx=$HOME/qxub \
    --glob '**/*.log' --group-by rule,label --label 'data/snp/(pileup|experiments/[^/]+)/' --out costs/
```

- `--source NAME=FOLDER` is repeatable; `--glob` selects files below each
  folder (default `**/*.log`).
- `--group-by` takes jobs.tsv columns after `source` (default `rule`).
- `--label` is a regular expression searched in the first output (else the
  first input); its groups, joined by `/`, make the `label` column. It
  separates runs of one rule that wrote to different folders.

## Outputs

| File | One row per | Columns |
|---|---|---|
| `jobs.tsv` | PBS block | `source`, `log`, `block`, `n_blocks`, `header_lost`, `rule`, `capture`, `wildcards`, `first_input`, `first_output`, `label`, `reason`, `job_id`, `finished`, `exit_status`, `pbs_kill`, `ncpus`, `mem_req_gb`, `walltime_req_h`, `walltime_h`, `cpu_h`, `mem_used_gb`, `jobfs_req_gb`, `jobfs_used_gb`, `su`, `cpu_share`, `charged_cpus`, `charge_basis`, `su_rate`, `su_expected`, `partial_block`, `smk_threads`, `smk_mem_mb`, `smk_runtime_min`, `smk_resources` |
| `rule_costs.tsv` | source and group | `jobs`, `failed`, `failed_on_walltime`, `captures`, ranges of `ncpus`, `mem_req_gb` and `walltime_req_h`, `median_walltime_min`, `max_walltime_min`, `total_walltime_h`, `total_cpu_h`, `median_cpu_share`, `median_mem_used_gb`, `p90_mem_used_gb`, `max_mem_used_gb`, `median_su`, `total_su`, `failed_su`, `su_expected`, `partial_blocks`, `first`, `last` |
| `capture_costs.tsv` | source and capture | `jobs`, `su`, `cpu_h`, `walltime_h`, `failed_su`, `partial_blocks` |
| `provenance.json` | | Inputs, versions, totals |

- `cpu_share` is CPU time over walltime times CPUs requested.
- A job failed when its exit status is not 0. It failed on walltime when PBS
  says so or it used 98 per cent of its walltime.
- `mem_used_gb` is PBS's figure, which includes the page cache of files the
  job read, so it overstates the memory a process needed. Measure the
  process itself (`/usr/bin/time -v`) before cutting a request.

## Limits

- A log file is normally replaced when a later job writes the same name, so
  a folder holds the last attempt for each name. Failed attempts that a retry
  overwrote are missing, and their SU with them. A Snakemake controller's own
  log keeps every submission and can bound what is missing.
- A capture is known only from a wildcard. Jobs that span captures have none.
