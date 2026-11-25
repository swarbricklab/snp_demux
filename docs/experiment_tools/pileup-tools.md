# Pileup tools: `pileup-compare`, `pileup-subset` and `pileup-join`

These tools compare, cut and join the pileup folders sc_pileup writes, so
that a site-list or barcode-list variant can be cut from an existing pileup
instead of piling the BAM again, and two captures can be joined into one
in-silico pool. They need numpy and pandas only (`pileup.py`).

## The two engines' folders

| Engine | Sites (one record per matrix row) | Barcodes (one per column) | Matrices |
|---|---|---|---|
| cellsnp | `cellSNP.base.vcf.gz` | `cellSNP.samples.tsv` | `cellSNP.tag.AD.mtx`, `cellSNP.tag.DP.mtx`, `cellSNP.tag.OTH.mtx` |
| souporcell | `variants.vcf` | `barcodes.tsv` | `alt.mtx`, `ref.mtx` |

The engine is read from the files present. A site is the CHROM, POS, REF and
ALT of its record. A site written twice keeps both rows, and the k-th copy
matches the k-th copy in another folder. souporcell's `--known_genotypes`
matches `variants.vcf` records to matrix rows by order, so the tools always
cut the two together.

## What a subset reproduces, and its limits

A site's counts do not depend on which other sites were piled, and a
barcode's counts do not depend on which other barcodes were piled. A subset
of a pileup therefore has the counts that a fresh pileup of the same sites
or barcodes would have. This holds for one BAM, one barcode list, one tool
version and one set of parameters. Two comparisons made before these tools,
with scripts, support it:

- Cut to the sites of a second panel, a pileup of about 6 million panel
  records gives the second panel's own pileup: every shared entry was
  identical for both engines on 17 captures (76.8 million cellsnp entries,
  69.6 million vartrix entries).
- Cut to one barcode list, with cellsnp's minimum count applied again to that
  list's UMIs, a pileup of two barcode lists gives that list's own pileup.

The evidence of this branch's pull request repeats the first on two captures
with these tools, and cuts a pileup to a barcode list byte for byte as a hand
cut did.

Limits:

- A site that was not piled cannot be added. Cut from a pileup of a larger
  panel, never from a smaller one.
- cellsnp drops a site whose UMIs over the piled barcodes fall under
  `--minCOUNT`. A barcode subset keeps every row unless `--min-count`
  applies the same rule to the kept barcodes (DP + OTH summed over them).
  cellsnp's `--minMAF` is not applied again: OTH does not say which other
  allele was seen. sc_pileup pipes with `minMAF 0`.
- souporcell's rows come from a scan of the BAM. sc_pileup scans the whole
  Cell Ranger BAM when it skips souporcell's remapping (`skip_remap: true`,
  its default), so the rows do not depend on the barcodes and are never
  filtered by count. With the remapping, the scanned BAM holds the listed
  barcodes' reads only, and a barcode subset can keep rows a fresh pileup
  of those barcodes would not have.

## Files as the tools write them

- Matrices keep their header lines (the banner and any `%` lines) and the
  separator of their size line and entries: cellsnp writes tabs and an
  `integer` banner, vartrix spaces and a `real` banner. Entries keep their
  order with rows and columns renumbered, so a matrix sorted by row then
  column stays sorted. vartrix's explicit zero entries are kept.
- VCF header lines and records are copied as read. cellsnp's INFO
  (`AD=..;DP=..;OTH=..`, the sums over the barcodes) is recomputed when the
  barcodes change (`--info auto`). A gzip input is written as BGZF.
- Barcodes keep the input's order.
- A file whose content does not change is copied byte for byte.

## `pileup-subset`

```bash
python -I tools/demux_tools pileup-subset IN OUT --sites panel.vcf.gz
python -I tools/demux_tools pileup-subset IN OUT --barcodes cells.tsv.gz --min-count 2
```

| Option | Does |
|---|---|
| `--sites F` | Keeps the rows whose site is in F: a VCF, or a table of CHROM, POS, REF and ALT (header optional) |
| `--barcodes F` | Keeps the barcodes in F, one per line (gzip allowed), in the input's order |
| `--min-count N` | cellsnp: drops rows whose DP + OTH over the kept barcodes is under N |
| `--strip-prefix P` | Removes P from every kept barcode (each must start with it) |
| `--info auto\|recompute\|keep` | cellsnp INFO: recompute when the barcodes change (default), always, or never |
| `--allow-empty` | Writes a subset with no site or no barcode; refused by default |
| `--force` | Writes into a folder that holds files; see below |

- The command prints, and records in `provenance.json`, how many sites and
  barcodes were requested and how many of them the pileup holds. A list
  whose names do not match the pileup (`chr1` against `1`, a barcode suffix)
  finds none, and the empty result is refused.
- An output folder that holds files is refused unless `--force` is given.
  `--force` first removes the pileup files this run will not write again,
  such as `cellSNP.tag.OTH.mtx` of an earlier run, so none is left stale.
  Every file is written under a temporary name and renamed into place, so a
  link in the folder (e.g. to a DVC cache object) is replaced, never written
  through.

The folder gets a `provenance.json` beside the pileup files.

## `pileup-compare`

```bash
python -I tools/demux_tools pileup-compare A B --out compare.tsv [--check]
```

Compares two folders of one engine as sets of (site, barcode, count),
ignoring row and column order. The table has one row for sites, one for
barcodes and one per matrix, with columns `layer`, `entries_a`,
`entries_b`, `in_both`, `equal`, `only_a`, `only_b`. The last line printed
is the verdict: `same sites and counts` when the two hold the same sites,
the same barcodes and equal entries in every matrix, or `differ` followed by
what differs. `--check` exits with status 1 when they differ.

- A site is the CHROM, POS, REF and ALT of its record. Other VCF columns and
  the header lines are not compared.
- A site written twice in one pileup is matched copy by copy: the first copy
  with the first, the second with the second.

| Option | Does |
|---|---|
| `--sites F` | Compares only at these sites, e.g. the sites two panels share |
| `--columns-of-b N` | cellsnp: first cuts A to B's barcodes and to rows with DP + OTH of N or more over them |
| `--layers` | Compares only these matrices |
| `--keep-zeros` | Counts an explicit zero as an entry; by default a zero and a missing entry are the same count |

## `pileup-join`

```bash
python -I tools/demux_tools pileup-join --in capA --prefix A_ --in capB --prefix B_ --sites union pool/
```

Joins folders of one engine into one in-silico pool. Barcodes take each
input's prefix and must then be distinct.

- With `--sites union`, the cells of one input have no counts at a site only
  another input piled, although a pileup of a real pool would count their
  reads there. Which sites a cell has counts at then tells which capture it
  came from, so such a pool is easier to demultiplex than a real one. Pile
  both captures with one site list, or join with `--sites intersection`. The
  command prints, and records in `provenance.json`, each input's rows, the
  rows every input has and the rows only one input has.
- Every input must hold every matrix joined. `--layers` names the matrices
  to join when an input lacks one.
- An output folder that holds files is refused unless `--force`, as for
  `pileup-subset`.

- `--sites union` (default) keeps every site of any input, ordered by
  chromosome (in order of first appearance) and position. Records at one
  position keep the order of the input that lists them first, as a panel
  can list two alleles of one position in either order.
  `--sites intersection` keeps sites in every input, in the first input's
  order.
- Entries are sorted by row then column. A cellsnp record's INFO is summed
  over all the joined barcodes.
- A souporcell `variants.vcf` keeps its first eight columns only, because
  the inputs' genotype columns belong to different donors. Build the pool's
  genotype VCF separately (e.g. `bcftools merge` of the two panels) and cut
  it to these records before passing it to souporcell.

Splitting a joined folder back gives each input:

```bash
python -I tools/demux_tools pileup-subset pool/ backA --barcodes A_barcodes.txt --strip-prefix A_ \
    --sites capA/cellSNP.base.vcf.gz
python -I tools/demux_tools pileup-compare capA backA --check --out check.tsv   # same sites and counts
```
