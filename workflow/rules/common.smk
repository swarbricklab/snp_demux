import pandas as pd
from pathlib import Path

samples = pd.read_csv(config["deps"]["samples"], dtype=str, index_col='sample_id')
samples.index = samples.index.astype(str)
captures= pd.read_csv(config["deps"]["captures"], dtype=str, index_col='capture_id')
captures.index = captures.index.astype(str)

# --- sc_pileup integration -------------------------------------------------
# Per-cell REF/ALT pileups are produced upstream by the sc_pileup module and
# consumed here (this repo no longer runs the pileups itself). Layout:
#   <pileup_root>/<engine>/<platform>/<capture_id>/<files>
# See docs (sc_pileup: docs/snp_demux_integration.md) for the contract.
_pileup_root = Path(config["deps"]["pileup_root"])
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

def lookup(df,field):
    """
    Look up field in the given dataframe
    """
    def handle_wildcards(wildcards):
        """
        Internal function for processing wildcards
        """
        result=df.loc[wildcards.capture_id, field]
        return result.split(";")
    return handle_wildcards


def get_donors_in_capture(wildcards):
    """
    Return a comma-separated list of the donors in wildcards.capture_id
    """
    capture=samples.loc[samples['capture_id']==wildcards.capture_id]
    donor_list=",".join(capture['donor_id'].tolist())
    return donor_list

def get_number_of_samples_in_capture_pool(wildcards):
    """
    Return the number of samples in wildcards.capture_id
    """
    return len(samples.loc[samples['capture_id']==wildcards.capture_id])