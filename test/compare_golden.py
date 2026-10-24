#!/usr/bin/env python3
"""Compare a fresh snp_demux test run against the DVC-tracked golden.

The donor assignments and consensus are deterministic, but gzip wrappers embed an
mtime and VCF headers carry volatile ``##`` lines, so we compare *content*, not
raw bytes:

  * .vcf / .vcf.gz : compare the data records only (lines not starting with '#')
  * other .gz      : compare the DECOMPRESSED bytes (ignores the gzip mtime)
  * everything else: compare raw bytes

Deliberately NOT compared (non-deterministic or non-call QC; see
``docs/reproducibility.md``):

  * ``*.png`` / ``*.html`` — plots and reports embed timestamps and are
    visualisations, not calls.
  * vireo's ambient-RNA QC (``prop_ambient.tsv``, ``ambient_summary.csv``,
    ``ambient_stats.csv``) — ``--callAmbientRNAs`` uses non-associative parallel
    float reductions, so these wobble run-to-run. The donor assignments
    (``donor_ids.tsv``) and everything downstream are byte-stable.
  * souporcell's ``ambient_rna.txt`` — likewise an ambient estimate, not a call.
  * index files (``.bai``/``.csi``/``.tbi``/``.fai``).

Exits non-zero (and lists them) if any golden file is missing from the fresh
output or differs in content.

Usage:  compare_golden.py GOLDEN_DIR FRESH_DIR
"""
import gzip
import hashlib
import sys
from pathlib import Path

IGNORE_SUFFIXES = (".bai", ".csi", ".tbi", ".fai", ".png", ".html", ".pdf")
# Non-deterministic QC outputs excluded by basename (see module docstring).
IGNORE_NAMES = frozenset(
    {
        "prop_ambient.tsv",
        "ambient_summary.csv",
        "ambient_stats.csv",
        "ambient_rna.txt",
    }
)
# Whole subtrees excluded by relative-path prefix. The report/ tree is rendered
# summaries + plots: the HTML/PNG embed timestamps, and report/metrics/*.tsv
# carries the non-deterministic ambient_souporcell / ambient_vireo_median columns
# (every *call* column in it is byte-stable). The underlying calls are checked in
# assignment/, vireo/, souporcell/, freemuxlet/, demuxlet/.
IGNORE_DIR_PREFIXES = ("report/",)


def _digest(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def content_key(path: Path) -> str:
    """A content fingerprint that ignores gzip mtime and VCF header volatility.

    Falls back to a raw-byte digest for files carrying a ``.gz``/``.vcf.gz`` name
    that are not actually gzip — e.g. vireo writes a plain-text placeholder
    ("No donor vcf file produced by vireo") as ``GT_donors.vireo.vcf.gz`` when a
    donor panel is supplied, so there is nothing to decompress.
    """
    name = path.name
    try:
        if name.endswith(".vcf.gz"):
            with gzip.open(path, "rt") as fh:
                records = [ln for ln in fh if not ln.startswith("#")]
            return _digest("".join(records).encode())
        if name.endswith(".vcf"):
            with open(path, "rt") as fh:
                records = [ln for ln in fh if not ln.startswith("#")]
            return _digest("".join(records).encode())
        if name.endswith(".gz"):
            with gzip.open(path, "rb") as fh:
                return _digest(fh.read())
    except (gzip.BadGzipFile, OSError):
        pass  # not actually gzip — fall through to a raw-byte compare
    return _digest(path.read_bytes())


def _skip(rel: str, name: str) -> bool:
    return (
        name.endswith(IGNORE_SUFFIXES)
        or name in IGNORE_NAMES
        or any(rel.startswith(pre) for pre in IGNORE_DIR_PREFIXES)
    )


def collect(root: Path):
    out = {}
    for p in sorted(root.rglob("*")):
        if not p.is_file():
            continue
        rel = p.relative_to(root).as_posix()
        if not _skip(rel, p.name):
            out[rel] = p
    return out


def main() -> int:
    if len(sys.argv) != 3:
        print(__doc__)
        return 2
    golden_dir, fresh_dir = Path(sys.argv[1]), Path(sys.argv[2])
    golden, fresh = collect(golden_dir), collect(fresh_dir)

    missing, differ = [], []
    for rel, gpath in golden.items():
        fpath = fresh.get(rel)
        if fpath is None:
            missing.append(rel)
        elif content_key(gpath) != content_key(fpath):
            differ.append(rel)

    extra = sorted(set(fresh) - set(golden))

    print(f"golden files: {len(golden)}   fresh files: {len(fresh)}")
    print(f"OK:      {len(golden) - len(missing) - len(differ)}")
    for rel in missing:
        print(f"MISSING: {rel}")
    for rel in differ:
        print(f"DIFFER:  {rel}")
    for rel in extra:
        print(f"extra (ignored): {rel}")

    if missing or differ:
        print("\nFAIL: fresh output does not match the golden.")
        return 1
    print("\nPASS: fresh output matches the golden.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
