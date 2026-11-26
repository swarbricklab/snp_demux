# Core: the package, shared readers, `synth`, `vireo-pairs` and `version`

This layer holds the dispatcher, the readers and writers every other tool
uses, a synthetic capture generator for tests, and one command that corrects
vireo's doublet pair names. It needs numpy and pandas only.

## Layout

```
tools/
  README.md
  demux_tools/
    __init__.py      the version; imports nothing
    __main__.py      entry point for `python -I tools/demux_tools`
    cli.py           the COMMANDS table, the dispatcher, the version command
    _io.py           text and tables, deterministic gzip, provenance sidecars
    cellsnp.py       cellsnp-lite pileups: read, filter rows, write
    vireo_files.py   vireo's output files, true doublet pairs, vireo-pairs
    synth.py         synthetic captures, the synth command
  tests/
    conftest.py      fixtures and markers
    test_tools_*.py  one file per module
    envs/image-pins.yaml
```

Each function's docstring states its exact contract. The summary below lists
what later tools build on.

## Shared modules

### `_io`

| Function | Does |
|---|---|
| `read_table(path, sep="\t", usecols=None, header="infer", names=None, comment=None)` | Reads a table with every column as text; "" and "NA" stay as written; gzip by magic bytes |
| `write_table(df, path, float_formats=None, columns=None)` | Writes a tab-separated table with a header, atomically; `float_formats` maps a column to a %-format; `.gz` gives mtime 0 |
| `open_text(path, mode)` | `"rt"` returns a file (gzip by magic); `"wt"` is a context manager that writes atomically |
| `atomic_path(path)`, `open_binary(path)`, `write_text(path, text)`, `write_json(obj, path)`, `read_lines(path)` | Smaller writers and readers with the same rules |
| `fingerprint(paths)` | sha1 over each file's real path, size and modification time |
| `check_output_path(out, inputs)` | Raises `UsageError` when the output is one of the inputs |
| `start_clock()`, `provenance(argv, inputs, clock, extra)`, `write_provenance(out, ...)`, `sidecar_path(out)` | The provenance record and where it goes |
| `git_sha(path=None)` | The commit checked out where the tools live, read from `.git` without running git |
| `library_versions()` | Python, numpy, pandas, scipy and vireoSNP versions; never imports scipy or vireoSNP |
| `UsageError` | A usage problem: the dispatcher prints it and exits with status 2 |

### `cellsnp`

| Name | Does |
|---|---|
| `Pileup` | `chrom`, `pos`, `ref`, `alt`, `barcodes`, `layers` (`{name: (rows, cols, vals)}`, 0-based int32), `rows` (each variant's row in its source), `variants` (vireo ids), `to_dense`, `to_csc` (imports scipy), `subset(rows, cols)` |
| `resolve_cellsnp(cellsnp=None, base=None, samples=None, ad=None, dp=None, oth=None)` | The five file paths of a pileup; explicit files override the folder, so DVC objects are read in place |
| `read_mtx(path, keep_rows=None, chunksize=5000000)` | `(shape, i, j, v)` with pandas' C reader, keeping only the rows asked for while it reads; checks the entry count |
| `read_cellsnp(paths, layers=("AD", "DP"), keep_rows=None)` | A `Pileup`; with `keep_rows` only those rows, in file order, with `rows` mapping back |
| `read_base(path)`, `read_samples(path)`, `mtx_header(path)` | The parts on their own |
| `write_cellsnp(pileup, out_dir, layers=None)`, `write_mtx(...)` | Files as cellsnp-lite writes them: entries sorted by row then column, zeros left out |
| `variant_id(chrom, pos, ref, alt)` | vireo's id, `CHROM_POS_REF_ALT` |

A site listed twice in `cellSNP.base.vcf.gz` keeps both rows.

### `vireo_files`

| Name | Does |
|---|---|
| `read_donor_ids(path)` | `donor_ids.tsv` as text; accepts vireoSNP 0.5.10's `prob_max_doublet` as `prob_doublet` |
| `read_prob_header(path)`, `donor_order(prob_singlet_path)` | Header names; vireo's donor order |
| `pair_naming(prob_doublet_path, donors=None)` | Whether the pair header uses permutations (0.5.7) or combinations (0.5.9), the donor order, and the true pair of each value column |
| `true_pairs(prob_doublet_path, written, donors=None)` | The pair vireo scored for each `best_doublet` it wrote |
| `best_pairs_from_matrix(prob_doublet_path, donors=None)` | The pair with the largest written probability per cell, with underflow and tie flags, as a cross-check |
| `vireo_pairs(donor_ids, prob_doublet, prob_singlet=None)` | The table that `vireo-pairs` writes, and its counts |
| `pair_key(a, b)` | The two names in string order joined by a comma |

### `synth`

| Name | Does |
|---|---|
| `make_capture(n_donors=3, n_sites=20, n_cells=10, seed=0, theta=(0.01, 0.5, 0.99), depth=(1, 4), coverage=0.6, doublets=((0, 1, 0.5),), blank_rate=0.0, donor_names=None, oth_rate=0.0, n_chroms=2)` | A `SynthCapture`: donors, a `Pileup` with AD, DP and OTH, true genotypes, blanks and a truth table |
| `write_capture(sc, out_dir, tags=("GT",), gp_noise=0.0, vcf_rows=None)` | `cellsnp/`, `donors.vcf.gz` and `truth.tsv` |
| `write_vcf(sc, path, tags, gp_noise, rows)` | The donor VCF alone, with GT, GP or PL |
| `build_on_disk(out_dir, tags, gp_noise, vcf_rows, **make_kwargs)` | Both steps; the result also holds the capture |

## Commands

### `version`

Prints the package version, the git commit of the checkout it runs from, and
the library versions.

### `synth`

```bash
python -I tools/demux_tools synth --out DIR [--donors 3] [--sites 20] [--cells 10] [--seed 0]
    [--doublet A:B:SHARE ...] [--no-doublets] [--blank-rate 0] [--coverage 0.6] [--depth 1:4]
    [--theta 0.01,0.5,0.99] [--oth-rate 0] [--tags GT[,GP,PL]] [--gp-noise 0] [--donor-names a,b,c]
```

It writes `DIR/cellsnp/` (the five cellsnp-lite files), `DIR/donors.vcf.gz`,
`DIR/truth.tsv` and `DIR/provenance.json`. The same options give the same
files, with numpy 1.21 to 2.x, because the draws come from numpy's legacy
`RandomState`.

Reads follow vireo's model. At a site, each UMI of a cell is alternative with
probability `theta[g]`, where g is the cell's donor's genotype (0, 1 or 2
alternative alleles). A doublet draws each UMI from its first donor with
probability `share`. For example, at a site where `d0` is 0/0 and `d1` is 0/1,
with the default theta:

| Cell | Alternative probability per UMI |
|---|---|
| Singlet of `d0` | 0.01 |
| Singlet of `d1` | 0.5 |
| 50:50 doublet of `d0` and `d1` | 0.5 × 0.01 + 0.5 × 0.5 = 0.255 |

Singlet c belongs to donor c mod n. Doublets follow the singlets. Every site is
polymorphic among the donors. A blank genotype is written `./.` in the VCF,
while the reads still follow the true genotype.

`truth.tsv` has one row per barcode, in pileup order:

| Column | Singlet | Doublet |
|---|---|---|
| `barcode` | The barcode | The barcode |
| `kind` | `singlet` | `doublet` |
| `donor` | The donor | "" |
| `pair` | "" | The two donors, `a,b` in string order |
| `donor_a` | The donor | The donor listed first in `--doublet` |
| `share_a` | `1` | That donor's share of UMIs |

### `vireo-pairs`

```bash
python -I tools/demux_tools vireo-pairs --vireo-dir DIR --out pairs.tsv.gz
python -I tools/demux_tools vireo-pairs --donor-ids F --prob-doublet F [--prob-singlet F] --out pairs.tsv.gz
```

vireo scores doublets over every pair of donors taken in combination order,
but vireoSNP 0.5.7 names the pairs in permutation order. The name written for
pair k is the k-th permutation, so only pairs that include the first donor are
named right. With donors `d0`, `d1` and `d2`:

| k | Pair scored | Name written by 0.5.7 |
|---|---|---|
| 0 | `d0,d1` | `d0,d1` |
| 1 | `d0,d2` | `d0,d2` |
| 2 | `d1,d2` | `d1,d0` |

A cell with pair probabilities 0.001, 0.002 and 0.95 is best explained by
`d1` and `d2`, and 0.5.7 writes `best_doublet` as `d1,d0`. The header of
`prob_doublet.tsv.gz` lists the six permutation names over three value
columns. The command reads vireo's donor order from that header and replaces
each written name by the pair at its position, which here is `d1,d2`.
vireoSNP 0.5.9 names pairs by combination, and the command then returns the
written pair in string order.

When every pair probability of a cell underflows to 0 (`prob_doublet` written
`0.00e+00`), vireo's choice is pair 0 whatever the data, and the command
leaves `best_pair` empty.

The output holds vireo's eight columns as written (`prob_max_doublet` renamed
`prob_doublet`) and three more:

| Column | Meaning |
|---|---|
| `best_pair` | The pair vireo scored highest, as `a,b` in string order; "" when every pair probability is 0 |
| `pair_prob_underflow` | `true` when every pair probability is 0, else `false` |
| `vireo_naming` | `permutations` (0.5.7) or `combinations` (0.5.9 and later) |

The command also checks each pair against the largest written value in
`prob_doublet.tsv.gz`. It stops with an error if a cell's largest value is a
clear winner and names another pair. Cells whose largest value is tied at the
written precision are counted, not checked. The counts go to the provenance
sidecar under `results`.

## Tests

`tools/tests/conftest.py` gives every layer's tests:

- `tools/` on `sys.path`, so tests `import demux_tools`;
- the markers `scipy` and `vireo`, whose tests skip where scipy or vireoSNP
  does not import;
- `run_tools(*args, cwd=None, python=None, blocked=(), package=None)`, which
  runs `python -I tools/demux_tools` in a subprocess, from a fresh folder by
  default, optionally with modules such as scipy made unimportable;
- `synthetic_capture(name="capture", tags=("GT",), gp_noise=0.0, vcf_rows=None, **make_capture_kwargs)`,
  which writes a capture under the test's temporary folder and returns its
  paths (`dir`, `cellsnp`, `vcf`, `truth`), its `donors` and the `capture`.

The core tests check, among other things:

- a gzip table written twice gives identical bytes, with mtime 0;
- `read_mtx` with `keep_rows` equals a full read followed by a filter, at
  several chunk sizes;
- a repeated site keeps both rows through a write and a read;
- the pair naming on hand-built headers (three donors, unsorted donors,
  combination headers, two donors, names with commas) and on files written by
  vireoSNP's own writer;
- every command's `--help` runs with scipy and vireoSNP unimportable;
- `python -I tools/demux_tools version` and `synth`, run from another folder,
  leave no `__pycache__` in the checkout.
