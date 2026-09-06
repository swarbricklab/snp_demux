# Hierarchical, sex-aware, WGS+micro demux — design

**Branch:** `hierarchical-sex-consensus` (module `swarbricklab/snp_demux`).
**Status:** design + implementation in progress. This doc is the PR narrative.

## Motivation — the failure taxonomy

An investigation into the previous **micro-only** production run (see the consuming
workspace `bcarc_chromium/data/snp/experiments/overnight/ANALYSIS_HANDOFF.md`) sorted
every hard capture into three failure modes:

1. **MAF-culling of minority-donor sites.** A donor's private site has pooled AF ≈
   (cell-fraction × 0.5), so `minMAF 0.1` systematically culled the distinctive sites of
   spike-ins / small donors (keep@0.1: 4262 = 1 %, PBMC = 2–10 % vs majority donors 80 %).
   *Fixable:* pile up at `minMAF 0`; add WGS density (15–27× more surviving sites).
2. **Genuinely too few cells** (4262/Pool_5d, 3937): no signal at any MAF, negative de-novo
   silhouette. *Not fixable by genotyping — flag.*
3. **Pool-wide ambient contamination** (4714/P1b, Pool_9a): 16 % of souporcell singlets carry
   a high doublet-posterior, and confident vs uncertain singlets show *no* separation in vireo
   `prob_doublet` (both ≈ 0.45). *Not fixable by genotyping — flag as low-quality.*

**The structural lesson:** the old Demuxafy consensus is a *symmetric* vote
(`Combine_Results.R --method AtLeastHalfSinglet`). When souporcell de-novo **recovers** a
minority donor as a clean cluster but genotype-based vireo **fails to label** it, the two
"disagree" and the cell is dropped. Recovery and labelling must be decoupled:

> **De-novo clustering decides *which cells group together* (recovery); genotypes decide
> *which donor a group is* (labelling). Keep recovered cells even when the label is weak,
> and report recovery-confidence and label-confidence separately.**

Success is **not** "assign every cell" — it is "assign confidently where the data allows and
flag honestly where it does not."

## What we can and cannot get

Incoming WGS (a few weeks): the upcoming sample list + **PBMC-6**. That upgrades **36/138**
pools to full-WGS (PBMC-6 alone flips +24, being the shared spike-in). The remaining difficult
donors (**4870, 4034, 4713, 3589, 4262**) **cannot** be WGS'd (insufficient material), so their
pools stay permanently mixed (WGS majority + micro-only minority). The method must therefore
resolve a **micro-only minority donor against a WGS majority** — and upgrade automatically when
the WGS lands, with no code change.

## Strategy

### A. Genotype panel — principled WGS+micro, adaptive per pool (UPSTREAM)

The module already consumes a per-pool panel at `panel_root/<capture_id>/informative.vcf.gz`
(the `deps.vcf` key is dead). Panel construction stays **upstream** in the consuming
workspace's `build_informative_panels`; this module is source-agnostic. The panel is built as:

- **Tier 1 — discriminative core:** sites where *every* pool donor has a real call. Full-WGS
  pool → WGS intersection (dense); mixed pool → WGS ∩ micro (every donor genuinely typed, **no
  `./.` confound**). This is the "one common complete site set" and drives distinctiveness +
  the primary genotype match.
- **Tier 2 — WGS boost (missingness-aware, A/B):** WGS-private sites present only for WGS
  donors, micro donors explicit `./.`. Consumed *only* by vireo (which models missing GTs),
  never by the distinctiveness metric. Sharpens the WGS majority so a micro-only minority
  falls out by elimination. Shipped as a panel A/B, not the default.

Pileup at `minMAF 0` (retain minority sites) — already the cohort default.

### B0. Configurable tool set

`params.tools` (subset of `{souporcell, vireo, demuxlet, freemuxlet}`, default all four)
selects which per-cell demux tools run. `rule all`, the per-tool aggregates, and the Demuxafy
consensus adapt: only enabled tools are required/passed to `Combine_Results.R`, and `--ref` is
chosen from the enabled set (Vireo > Souporcell). Demuxafy needs ≥2 tools or it fails fast with
an actionable message. This lets a workspace without popscle pileups run `[souporcell, vireo]`
end-to-end (both consensuses + report) without demuxlet/freemuxlet.

### B. Pipeline — souporcell twice + vireo + hierarchical consensus (ALONGSIDE Demuxafy)

- **souporcell de-novo** (existing `souporcell.smk`) — the recovery backbone. Keep.
- **souporcell known-GT** (NEW second mode) — `--known_genotypes` pointed at the real donor
  panel, for a direct genotype assignment to compare against de-novo. Writes to a parallel
  `souporcell/{c}/known/…` path.
- **vireo** genotype-based (existing) — per-cell donor + `prob_doublet`.
- **vireosnp** cluster↔genotype matcher (existing `annotate.py`) — the `GPb_diff` margin.
- **Demuxafy `Combine_Results.R`** — **kept, unchanged**, runs in parallel for head-to-head.
- **NEW hierarchical consensus** (`hierarchical_consensus.smk`) — consumes the same per-tool
  inputs plus vireosnp `clusters_annotated.tsv` and the sex join, and emits
  `assignment/{c}/hierarchical_consensus.tsv` with **separate `recovery_confidence` and
  `label_confidence`** columns. Decision order:
  1. **cluster membership** from de-novo souporcell (high trust);
  2. **cluster→donor label** = genotype margin (vireosnp `GPb_diff` + vireo), **sex-constrained**;
  3. **per-cell confidence** = f(cluster posterior, genotype margin, sex agreement);
  4. **doublet** = souporcell doublet-posterior ∨ vireo `prob_doublet` ∨ sex-conflict.
  Cells de-novo recovers with a weak label are kept as `recovered / low-confidence-label`,
  never silently dropped. Agreement between the two souporcell modes and vireo is a per-cluster
  QC signal, not a gate.

### C. Sex-per-barcode — cross-check *and* integrated prior

Per-barcode sex calls already exist upstream
(`data/cellranger/count/predemux/<capture>/sex_calls.csv`:
`barcode,…,sex_call∈{female,male,no_signal},confidence`). **Donor sex is derived from the
panel itself** (chrX heterozygosity / chrY calls) — no external metadata, fully principled.

- **Constraint:** aggregate high-confidence per-barcode sex per de-novo cluster; a cluster's
  donor must match its sex. This partitions the pool into within-sex subproblems and **breaks
  low-margin ties** — the direct attack on the "flat cluster" problem.
- **QC:** cells whose sex conflicts with the assigned donor's sex, or `no_signal` enrichment,
  flag doublet/ambient. Sex-signal degradation is itself an ambient metric (expected on P1b).

### D. Reports (every capture, every round) — as QC gates, not just diagnostics

Donor **distinctiveness** (Tier-1 unique sites, VCF vs scRNA-detectable, with the MAF-impact
curve), **PCA/UMAP** coloured by de-novo cluster, per-cluster **silhouette**, cluster↔donor
`GPb_diff` margin tables, **sex-concordance**, ambient/doublet fractions. Silhouette +
distinctiveness become gates that decide *trusted* vs *flagged* per pool.

### E. Iterative refinement (the rescue for un-WGS-able minorities)

- **Round 1:** full pipeline → confident cells, flagged clusters, doublet/ambient set.
- **Round 2:** (a) drop flagged doublet/ambient barcodes, re-cluster the clean set (sharper
  separation; marginal help for ambient pools); (b) **build an empirical consensus genotype per
  de-novo cluster from Round-1 confident cells, then match *that* to the micro panel.** A
  pooled-cell empirical genotype is far denser than any single-cell readout — the principled way
  to label a weak micro-only donor (3589/4034/4713/4870). This is where the biggest gains on the
  un-WGS-able donors will come from.
- **Leave-one-out** per pool: drop each majority donor, confirm the minority cluster survives as
  novel/unassigned (distinct ✓) rather than being absorbed (not separable ✗). Also detects swaps.

### F. Forward-compatibility

Panel build is parameterised by per-donor source availability, so when the WGS lands a re-run
upgrades 36 pools mixed→full-WGS automatically; the un-WGS-able difficult donors stay on the
mixed path. Round-1 artifacts (empirical cluster genotypes, margin tables) are durable so Round-2
and the WGS upgrade reuse them.

## Tuning experiments (pilot on the 6-pool shortlist first)

Shortlist spread: 2 full-WGS, 2 mixed-recoverable (Pool_7, Pool15), 2 quality-limited (P1b,
Pool_5d). Settle these knobs before cohort-scale compute:

1. **Panel A/B** — Tier-1 vs Tier-1+2 boost → margin on mixed pools.
2. **k / overclustering** — k=n vs k=n+2-then-merge → does overclustering isolate ambient/doublet?
3. **Sex-prior ablation** — assignment ± sex constraint → Δ margin on low-margin clusters.
4. **Empirical-genotype rescue** — quantify Round-2 boost on the un-WGS-able donors.
5. **Ambient denoise** — cellbender / soup ambient correction on P1b/9a → Δ silhouette.
6. **Min-cell floor** — sweep donor cell-count vs recoverability → auto-flag threshold (4262 anchors it).

## Implementation map (files touched / added in this module)

| change | file(s) | note |
|---|---|---|
| souporcell known-GT second mode | `workflow/rules/souporcell.smk` (+ parallel troublet/consensus → `souporcell/{c}/known/`) | de-novo run stays as-is |
| hierarchical consensus rule | **NEW** `workflow/rules/hierarchical_consensus.smk` + `workflow/scripts/hierarchical_consensus.py` | runs alongside `demuxafy` |
| sex-call input + join | `workflow/rules/common.smk` (`sex_path` helper, `deps.sex_calls`), **NEW** `workflow/scripts/join_sex.py` | CSV keyed on capture |
| report panels (distinctiveness/UMAP/silhouette) | **NEW** `workflow/scripts/report_variant_embedding.py` + `report_distinctiveness.py`; inputs added to `report.smk`/`report_capture.Rmd` | py via `envs/vireosnp.yaml` (pandas/matplotlib) |
| leave-one-out validator | **NEW** `workflow/scripts/leave_one_out.py` + optional rule (diagnostic target) | per-pool QC |
| config keys | `config/template_config.yaml` (+ consuming workspace config) | `deps.sex_calls`, `params.hierarchical.*`, `params.souporcell.known_gt: true`, panel tier flag |
| aggregate/targets | `workflow/rules/aggregate.smk`, `Snakefile` | new `aggregate_hier` sentinel; keep the 5 existing |

Demuxafy's `demuxafy.smk` / `Combine_Results.R` path is **not modified** — the hierarchical
consensus is additive so the PR can show both side by side.
