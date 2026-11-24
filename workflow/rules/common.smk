import pandas as pd
from pathlib import Path
from snakemake.exceptions import WorkflowError

# --- selectable per-cell demux tools ---------------------------------------
# The pipeline can run/aggregate/consensus any subset of these four per-cell
# demux tools. Selection is via config['params']['tools']; default = all four
# (fully backward compatible). demuxlet + freemuxlet need popscle pileups which
# are not always available, so a run may enable only {souporcell, vireo}.
KNOWN_TOOLS = ["souporcell", "vireo", "demuxlet", "freemuxlet"]


def enabled_tools():
    """Return the list of enabled per-cell demux tools (config['params']['tools'],
    defaulting to all four). Validated against KNOWN_TOOLS; an unknown entry raises
    a clear WorkflowError. Importable/callable at parse time."""
    tools = config.get('params', {}).get('tools', list(KNOWN_TOOLS))
    if isinstance(tools, str):
        tools = [tools]
    unknown = [t for t in tools if t not in KNOWN_TOOLS]
    if unknown:
        raise WorkflowError(
            "config params.tools contains unknown tool(s) %s; "
            "valid values are %s" % (unknown, KNOWN_TOOLS)
        )
    return list(tools)


def tool_enabled(name):
    """True if per-cell demux tool `name` is enabled for this run."""
    return name in enabled_tools()


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
# Per-pool informative panels (donor genotypes). Built by an upstream stage
# (e.g. bcarc_chromium build_informative_panels). Defaults to <pileup_root>/panel
# for standalone sc_pileup runs, but is normally set to the orchestration
# panel dir (e.g. data/snp/vcf/panel).
_panel_root = Path(config["deps"].get("panel_root", str(_pileup_root / "panel")))
# Per-pool panel VCF filename under <panel_root>/<capture>/. Defaults to the
# historical informative.vcf.gz; set deps.panel_file (e.g. panel.vcf.gz) for
# panels that aren't filtered to informative sites. The index is <file>.csi.
_panel_file = config["deps"].get("panel_file", "informative.vcf.gz")
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


def panel_vcf(suffix=""):
    """Input-function for the per-pool panel VCF (deps.panel_file), or its index
    with suffix=".csi"."""
    return panel_path(_panel_file + suffix)


# Per-capture sex-call CSV (barcode,sex_call,confidence), produced upstream by
# cellranger predemux. Configured as a Python format-string keyed on {capture_id}
# (config deps.sex_calls). May be None if the key is absent.
_sex_calls_template = config["deps"].get("sex_calls")


def sex_path(*parts):
    """Return an input-function giving the per-capture sex-call CSV path, resolved
    from the deps.sex_calls format-string on {capture_id} (mirrors pileup_path in
    spirit; the template already points at the file so *parts is normally unused,
    but any extra parts are appended for forward flexibility). Returns "" when the
    key is unset so the hierarchical join can degrade to empty sex columns.

    NB: consume this via `params:` (not `input:`) — the sex CSV is an optional
    upstream source with no producing rule, so making it a hard input would break
    the DAG whenever it is absent. join_sex.py / hierarchical_consensus.py check
    for existence at run time."""
    def _inner(wildcards):
        if not _sex_calls_template:
            return ""
        path = _sex_calls_template.format(capture_id=wildcards.capture_id)
        if parts:
            path = str(Path(path).joinpath(*parts))
        return path
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