# Experiment tools

`tools/demux_tools` is a Python package of helpers for experiments on this
workflow's inputs and outputs: pileups, donor VCFs and the calls each method
writes. It sits outside `workflow/` and nothing in the workflow imports it, so
the tools change no workflow output. A consumer that lists
`modules/snp_demux/workflow` as a DVC dependency sees no change when the tools
change.

## Running the tools

Run a tool by path from any checkout:

```bash
python -I <checkout>/tools/demux_tools --help             # the commands
python -I <checkout>/tools/demux_tools <command> --help   # one command's options
python -I <checkout>/tools/demux_tools <command> [options]
```

- `-I` keeps the user's site-packages and `PYTHON*` variables out of the run.
- The entry point turns bytecode writing off before it imports anything, so a
  run leaves no `__pycache__` in the checkout. Python caches the entry point
  itself before that line runs, and the entry point deletes that one file.
  `python -I -B` avoids even that write.
- Each command imports only its own module, so a command's `--help` needs only
  the libraries that command uses.
- The tools read outputs of either snp_demux line and import nothing from
  `workflow/`.

To use the modules as a library:

```python
import sys
sys.dont_write_bytecode = True
sys.path.insert(0, "<checkout>/tools")
from demux_tools import cellsnp, vireo_files
```

## Environments

| Tools | Needs | Runs in |
|---|---|---|
| Every command except the Scorer's | Python 3.7 or later, numpy 1.21 to 2.x, pandas 1.3 to 2.2 | The workflow's Snakemake environment (no scipy), an analysis environment with numpy 2, or the vireo image |
| The Scorer | scipy as well, and vireoSNP 0.5.7 or 0.5.9 for its vireo back end | The vireo image (`snp_demux-vireo:0.5.7`, Python 3.7) for exact reproduction |

The code keeps to Python 3.7 so that it runs in the vireo image. Only the
Scorer's modules may import scipy, and never at module level.

## Conventions shared by every tool

- **Outputs go only where the command line names them.** A tool never writes
  beside an input, and it refuses an output path that is one of its inputs.
  Every file is written to a temporary name and renamed when complete.
- **Tables are tab-separated text with a header.** They are read with every
  column as text, so a donor id such as `0012` or `NA` keeps its spelling.
- **Gzip output is deterministic.** A `.gz` output has mtime 0 and no file
  name in its header, so the same content gives the same bytes in one
  environment. Gzip input is recognised by its first bytes, not its name.
- **Every command records its provenance.** Most write a sidecar
  `<output>.provenance.json`, or `provenance.json` inside an output folder.
  `score` keeps the record inside `PREFIX.fit.json`, and `validate-vireo`
  writes one only with `--out`. The record holds the tool version
  and git commit, the tools' path, the command line, the Python, numpy,
  pandas, scipy and vireoSNP versions, each input's path, size and
  modification time, the time and peak memory used, and the command's own
  counts under `results`.
- **A donor pair is written as the two donor names in string order, joined by
  a comma**, for example `d0,d2`.
- **The repository holds no cohort data.** Tests build synthetic captures.
  Manifests that name real captures, donors or paths stay with the
  experiment that uses them.

## Adding a command

1. Add one line to `COMMANDS` in `tools/demux_tools/cli.py`, in alphabetical
   order: `"name": ("module", "one-line help"),`.
2. In that module, define `configure_<name>(parser)`, which adds arguments to
   an `argparse` parser, and `run_<name>(args)`, which returns an exit status.
   A `-` in the command name is written `_` in the function names. Raise
   `demux_tools._io.UsageError` for a usage problem.
3. Put the tests in `tools/tests/test_tools_<module>.py`. The fixtures and the
   `scipy` and `vireo` markers are described in `tools/tests/conftest.py`.
4. Describe the command in its own file under `docs/experiment_tools/` and
   change its line below from planned to available.

## Tests

```bash
pytest -q -p no:cacheprovider tools/tests
```

`.github/workflows/tools.yml` runs the tests in three environments on every
pull request that touches the tools or these pages:

| Leg | Environment | Skips |
|---|---|---|
| `image-pins` | The vireo image's pins (`tools/tests/envs/image-pins.yaml`): Python 3.7, numpy 1.21, pandas 1.3, scipy 1.7, vireoSNP 0.5.7 | None |
| `modern` | Python 3.12, numpy 2, pandas 2.2, scipy, vireoSNP 0.5.9 | None |
| `pandas-only` | Python 3.12, numpy 1.26, pandas 2.1, like the Snakemake environment | Tests marked `scipy` or `vireo` |

## Sections

One page per group of tools:

- [Core](experiment_tools/core.md) (available): the package, the shared readers and writers, `synth`, `vireo-pairs` and `version`.

- [Scorer](experiment_tools/scorer.md) (available): vireo's fixed-mode calls reproduced from the pileup and the donor VCF, with a margin per cell.

- [Scorer diagnostics](experiment_tools/scorer-diagnostics.md) (available): native fit, site masks and groups, genotype overrides and per-site detail.

- [Compare runs](experiment_tools/compare-runs.md) (planned): call readers for both snp_demux lines, run-to-run comparison and the sex check.

- [Scorecard](experiment_tools/scorecard.md) (planned): per-cell tables, scorecards against truth, calibration and doublet excess.

- [Job costs](experiment_tools/job-costs.md) (planned): service units per job, rule and capture from PBS logs.

- [Pileup tools](experiment_tools/pileup-tools.md) (planned): compare, subset and join pileups.
