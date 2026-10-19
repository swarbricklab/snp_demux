import pandas as pd
from pathlib import Path
from snakemake.utils import validate

# Fail fast with a clear message if the config is missing/misnamed keys, rather
# than deep inside a rule. Schema in workflow/schemas/config.schema.yaml.
validate(config, "../schemas/config.schema.yaml")

samples = pd.read_csv(config["deps"]["samples"], dtype=str, index_col="sample_id")
samples.index = samples.index.astype(str)
captures = pd.read_csv(config["deps"]["captures"], dtype=str, index_col="capture_id")
captures.index = captures.index.astype(str)

# --- sc_pileup integration -------------------------------------------------
# Per-cell REF/ALT pileups are produced upstream by the sc_pileup module and
# consumed here (this repo no longer runs the pileups itself). Layout:
#   <pileup_root>/<engine>/<platform>/<capture_id>/<files>
# See docs (sc_pileup: docs/snp_demux_integration.md) for the contract.
_pileup_root = Path(config["deps"]["pileup_root"])
# Per-pool informative panels (donor genotypes). Built by an upstream stage
# (e.g. bcarc_chromium build_informative_panels). Defaults to <pileup_root>/panel
# for standalone sc_pileup runs, but is normally set to the orchestration
# panel dir (e.g. data/snp/vcf/panel).
_panel_root = Path(config["deps"].get("panel_root", str(_pileup_root / "panel")))
_default_platform = config["deps"].get("default_platform", "short")


def _capture_platform(capture_id):
    if "platform" in captures.columns:
        v = captures.loc[capture_id, "platform"]
        if isinstance(v, str) and v and v.lower() != "nan":
            return v
    return _default_platform


def pileup_path(engine, *parts):
    """Return an input-function giving the sc_pileup output path for
    <engine>/<platform>/<capture_id>/<parts...>. NB: pathlib join (no f-strings
    in .smk — snakemake mangles `{}` inside them)."""

    def _inner(wildcards):
        plat = _capture_platform(wildcards.capture_id)
        return str(_pileup_root.joinpath(engine, plat, wildcards.capture_id, *parts))

    return _inner


def panel_path(*parts):
    """Return an input-function giving the per-pool informative panel path:
    <panel_root>/<capture_id>/<parts...> (no engine/platform level).
    This panel carries the pool's donor genotypes, already subset to the pool's
    donors and filtered to informative (donor-discriminative) sites by sc_pileup,
    so it replaces the old cohort `deps.vcf` as the donor genotype source. Using
    the SAME file sc_pileup piled against keeps pileup sites and demux genotypes
    identical by construction."""

    def _inner(wildcards):
        return str(_panel_root.joinpath(wildcards.capture_id, *parts))

    return _inner


def lookup(df, field):
    """
    Look up field in the given dataframe
    """

    def handle_wildcards(wildcards):
        """
        Internal function for processing wildcards
        """
        result = df.loc[wildcards.capture_id, field]
        return result.split(";")

    return handle_wildcards


def get_donors_in_capture(wildcards):
    """
    Return a comma-separated list of the donors in wildcards.capture_id
    """
    capture = samples.loc[samples["capture_id"] == wildcards.capture_id]
    donor_list = ",".join(capture["donor_id"].tolist())
    return donor_list


def get_number_of_samples_in_capture_pool(wildcards):
    """
    Return the number of samples in wildcards.capture_id
    """
    return len(samples.loc[samples["capture_id"] == wildcards.capture_id])


def get_samples_in_capture(wildcards):
    """
    Return the sample_ids declared for wildcards.capture_id, comma-separated.
    """
    capture = samples.loc[samples["capture_id"] == wildcards.capture_id]
    return ",".join(capture.index.tolist())


def get_pool_for_capture(wildcards):
    """
    Return the pool_id wildcards.capture_id belongs to, or "" if unrecorded.

    pool_id is optional and used only for display and grouping in the reports;
    no other rule reads it. A blank cell in captures.csv reaches us as a pandas
    NaN, which Snakemake serialises into the generated R script as a bare `nan`
    -- an undefined symbol in R, so the job dies with "object 'nan' not found"
    before the script runs. Coerce to a string so a missing optional field
    cannot fail a rule.
    """
    if "pool_id" not in captures.columns:
        return ""
    pool = captures.loc[wildcards.capture_id, "pool_id"]
    if pool is None or (isinstance(pool, float) and pd.isna(pool)):
        return ""
    return str(pool)
