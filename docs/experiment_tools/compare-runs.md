# Compare runs: `calls`, `compare-runs` and `sex-check`

These tools read the calls of any demultiplexing output into one schema,
compare two runs cell by cell, and check donor labels against per-cell sex
calls. They need numpy and pandas only (`calls.py`, `compare.py`).

## The calls schema

Every reader returns one row per barcode, all columns text:

| Column | Holds |
|---|---|
| `barcode` | The cell barcode as written |
| `status` | `singlet`, `doublet`, `unassigned` or `recovered`; any other word is kept as written. hier-v2's `low_cell` becomes `unassigned` |
| `status_raw` | The status word as the file wrote it |
| `donor` | The donor of a singlet or recovered cell; empty otherwise |
| `pair` | A doublet's two donors in string order joined by a comma, such as `d0,d2`; empty when the file names none |
| `confidence` | The file's own score: vireo's `prob_max`, the consensus `label_confidence`, the Scorer's `margin`; empty for the others |
| `cluster` | souporcell's cluster (`3`, or `1/4` for a doublet) where the file has one |
| `source` | The kind the rows were read as |

vireo writes a best pair for every cell, so its rows carry a pair whatever
the call. A cell whose pair probabilities all underflow to 0 carries the
name of pair 0, as vireo writes it; `vireo-pairs` (Core) blanks those.

## Kinds of calls

| Kind | File | Notes |
|---|---|---|
| `vireo` | `donor_ids.tsv` | The pair is corrected with the header of `prob_doublet.tsv.gz`, read from beside the file unless given; vireoSNP 0.5.7 writes the wrong name for most pairs (Core) |
| `souporcell` | `clusters.tsv` | Clusters become donors through a cluster map; without one a singlet's donor is its cluster number |
| `hier` | `hierarchical_consensus.tsv` | hier-v2's consensus, with the `recovered` tier |
| `assignment` | `cell_assignment.tsv` | A status other than singlet, doublet, unassigned, recovered or `low_cell` (main writes `NA`) is split with `--sba` into doublet and unassigned, from `Consensus_DropletType` or else `AtLeastHalfSinglet_DropletType` |
| `demuxafy:TOOL` | `sample_barcode_assignments.tsv` | `TOOL_DropletType` and `TOOL_Individual_Assignment`, e.g. `demuxafy:Vireo` |
| `scorer` | The Scorer's per-cell table | The pair is `best_pair`, the confidence `margin` |
| `generic:barcode=COL,status=COL[,donor=COL,pair=COL,confidence=COL,cluster=COL]` | Any table | Values as written |

Without a kind, the reader picks one from the header. A sample barcode file
needs `demuxafy:TOOL`.

Labels become donors in three steps:

1. A label map (`--label-map`), a table with columns `label` and `donor`, or
   `sample_id` and `donor_id` (a sample sheet works as it is).
2. A regular expression (`--label-pattern`) whose first group is the donor,
   for labels the map does not hold.
3. Otherwise the label itself.

A cluster map is a table with columns `cluster` and `donor`, or
`sba:PATH` to take Demuxafy's own naming of the souporcell singlets in a
`sample_barcode_assignments.tsv`. In `compare-runs`, `majority` names run B's
clusters by the donor most of their singlets carry in run A, counting only
cells that both runs call singlet.

## `calls`

```bash
python -I tools/demux_tools calls demux/vireo/CAP/donor_ids.tsv --out calls.tsv.gz
```

Writes the normalised table and its provenance sidecar.

## `compare-runs`

Joins two runs on the barcode, over every barcode of either run.

```bash
python -I tools/demux_tools compare-runs \
    --a hier:demux/assignment/CAP/hierarchical_consensus.tsv \
    --b assignment:demux/assignment/CAP/cell_assignment.tsv --b-label-map samples.csv \
    --capture CAP --sex-calls CAP/sex_calls.csv --male-donors D6 --out cmp/CAP
```

Many captures go in a manifest with columns `capture`, `a_path` and
`b_path`, and optionally `a_kind`, `b_kind`, `sex_calls` and, per side,
`a_prob_doublet`, `a_label_map`, `a_label_pattern`, `a_cluster_map` and
`a_sba` (`b_` likewise). The command-line options per side are the defaults
for the manifest's empty fields.

| Output | Columns |
|---|---|
| `PREFIX.cells.tsv.gz` | `capture`, `barcode`, `status_a`, `status_b`, `donor_a`, `donor_b`, `pair_a`, `pair_b`, `in_a`, `in_b`, `same_status`, `same_donor`, `same_pair`, `same_call`, `transition`; sorted by barcode within a capture |
| `PREFIX.transitions.tsv` | `capture`, `status_a`, `status_b`, `cells`, `same_donor`, `different_donor` |
| `PREFIX.summary.tsv` | `capture`, `cells_a`, `cells_b`, `cells_both`, `only_a`, `only_b`, `same_call`, `same_call_pct`, `singlet_donor_changes`, `doublets_both`, `doublet_pairs_compared`, `doublet_pairs_same` |
| `PREFIX.sex_check.tsv` | Written with `--sex-calls` and donor sex: the sex check of both runs |
| `PREFIX.provenance.json` | Inputs, versions, totals |

Definitions, for a barcode in both runs:

- `same_status`: the statuses are equal.
- `same_donor`: the donor fields are equal. They are empty for doublets and
  unassigned cells, so two doublets have the same donor.
- `same_pair`: `true` or `false` when both runs name a pair, empty when either
  names none. The consensus, `cell_assignment.tsv` and Demuxafy name no pair;
  vireo names one for every cell.
- `same_call`: same status and same donor, and for two doublets whose pairs
  both runs name, the same pair. A doublet in a run that names no pair is
  not counted as a different call.
- `transition`: `status_a -> status_b`. A barcode missing from a run has the
  status `absent` there.
- In the transitions table, `same_donor` and `different_donor` count only
  cells where both runs name a donor.
- In the summary, `doublets_both` counts cells both runs call doublet,
  `doublet_pairs_compared` those with a pair on both sides and
  `doublet_pairs_same` those with the same pair.

A worked example with four cells:

| Barcode | Run A | Run B | Transition | Counted as |
|---|---|---|---|---|
| 1 | singlet `d0` | singlet `d0` | `singlet -> singlet` | same donor, same call |
| 2 | singlet `d1` | singlet `d2` | `singlet -> singlet` | different donor; a singlet donor change |
| 3 | recovered `d2` | singlet `d2` | `recovered -> singlet` | same donor |
| 4 | doublet `d0,d1` | absent | `doublet -> absent` | only in A |

## `sex-check`

Counts the cells whose donor's sex disagrees with the cell's own sex call,
per run, status and label sex.

```bash
python -I tools/demux_tools sex-check --calls hier:CAP/hierarchical_consensus.tsv \
    --sex-calls CAP/sex_calls.csv --male-donors D6,D7 --run consensus --capture CAP --out sex.tsv
```

- Donor sex comes from `--donors`, a table with columns `donor` and `sex`,
  or from `--male-donors`, which makes every other donor female.
- A donor in neither has sex `unknown` and is not counted.
- `with_sex_call` counts cells whose call is `male` or `female`, in any
  case. `discordant` counts those whose call is the other sex.
- A manifest takes columns `capture`, `run`, `path` and `sex_calls`, and
  optionally `kind`, `prob_doublet`, `label_map`, `label_pattern`,
  `cluster_map` and `sba`.

Output columns: `capture`, `run`, `status`, `label_sex`, `cells`,
`with_sex_call`, `discordant`, `discordant_share`.

The share of female calls among cells labelled a male donor is an upper
bound on wrong-donor calls. A male cell low in RNA reads female when it has
an ambient XIST UMI and no chrY UMI.
