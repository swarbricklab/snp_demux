import pandas as pd
from pathlib import Path
from snakemake.exceptions import WorkflowError
from snakemake.utils import validate

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
    tools = config.get("params", {}).get("tools", list(KNOWN_TOOLS))
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


# --- hierarchical consensus parameters ---------------------------------------
# Optional block params.hierarchical; each key falls back to the documented
# default (config/template_config.yaml, workflow/schemas/config.schema.yaml).
HIERARCHICAL_DEFAULTS = {
    "margin_min": 0.15,
    "sex_weight": 1.0,
    "min_cells_floor": 20,
    "low_conf_keep": True,
}


def hier_param(key):
    """Return params.hierarchical[key], or its default."""
    return config["params"].get("hierarchical", {}).get(key, HIERARCHICAL_DEFAULTS[key])


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


# Doublet rate per 1,000 cells recovered (here, barcodes demultiplexed), by chip,
# from 10x's multiplet-rate table.
# Overridable via params.doublet.rate_per_1k; a capture picks its chip from an
# optional `chip` column in captures.csv (matched case-insensitively).
_DOUBLET_RATE_DEFAULTS = {"standard": 0.008, "ht": 0.004}


def get_doublet_rate_per_1k(wildcards):
    """
    Return the expected doublet rate per 1,000 cells for wildcards.capture_id's chip.
    """
    cfg = config["params"].get("doublet", {})
    rates = {
        str(k).lower(): float(v)
        for k, v in cfg.get("rate_per_1k", _DOUBLET_RATE_DEFAULTS).items()
    }
    chip = str(cfg.get("default_chip", "standard"))
    if "chip" in captures.columns:
        v = captures.loc[wildcards.capture_id, "chip"]
        if isinstance(v, str) and v and v.lower() != "nan":
            chip = v
    if chip.lower() not in rates:
        raise ValueError(
            "capture "
            + wildcards.capture_id
            + ": chip '"
            + chip
            + "' has no params.doublet.rate_per_1k entry (known: "
            + ", ".join(sorted(rates))
            + ")"
        )
    return rates[chip.lower()]


# --- demuxafy consensus over the enabled tools ---------------------------------
# Per-tool outputs handed to Demuxafy's Combine_Results.R. Kept as a single
# source of truth so the input files (below) and the shell --<tool> <dir> flags
# stay in lock-step. NB: Combine_Results.R takes a DIRECTORY per tool (the rule
# passes dirname(<file>)); flag names are literally --demuxlet/--freemuxlet/
# --souporcell/--vireo.
def _demuxafy_tool_outputs():
    """Map tool -> its per-capture output pattern. A function, not a dict, because
    the rules it names are defined after common.smk is parsed; it is only ever
    called from input/params functions, i.e. after the whole workflow is parsed."""
    return {
        "demuxlet": rules.demuxlet.output.best,
        "freemuxlet": rules.freemuxlet.output.samples,
        "souporcell": rules.souporcell_troublet.output.clusters,
        "vireo": rules.relabel_vireo.output.donor_ids,
    }


def _demuxafy_require_two_tools():
    """Demuxafy consensus needs >=2 per-cell tools to combine. With fewer, it
    cannot form a consensus; fail early with a clear message rather than crashing
    inside Combine_Results.R. (Default is 4 tools, so this never triggers in the
    backward-compatible case.)"""
    if len(enabled_tools()) < 2:
        raise WorkflowError(
            "Demuxafy consensus (rule demuxafy) requires at least 2 enabled "
            "per-cell demux tools in config params.tools, but only %s "
            "is/are enabled. Enable a second tool or drop the demuxafy targets."
            % enabled_tools()
        )


def demuxafy_inputs(wildcards):
    """Only the ENABLED tools' outputs are declared as inputs, so disabled tools
    (e.g. demuxlet/freemuxlet without popscle pileups) are never demanded."""
    _demuxafy_require_two_tools()
    c = wildcards.capture_id
    return {
        tool: _demuxafy_tool_outputs()[tool].format(capture_id=c)
        for tool in enabled_tools()
    }


def demuxafy_combine_flags(wildcards):
    """Build the per-tool Combine_Results.R flags for just the enabled tools,
    preserving the original flag order/values so the invocation is byte-for-byte
    identical when all four tools are enabled."""
    _demuxafy_require_two_tools()
    c = wildcards.capture_id
    p = config["params"]["demuxafy"]
    outs = _demuxafy_tool_outputs()

    def tool_dir(tool):
        return str(Path(outs[tool].format(capture_id=c)).parent)

    parts = []
    if tool_enabled("demuxlet"):
        parts += ["--demuxlet", tool_dir("demuxlet")]
    if tool_enabled("freemuxlet"):
        parts += ["--freemuxlet", tool_dir("freemuxlet")]
        parts += [
            "--freemuxlet_correlation_limit",
            str(p["freemuxlet_correlation_limit"]),
        ]
    if tool_enabled("souporcell"):
        parts += ["--souporcell", tool_dir("souporcell")]
        parts += [
            "--souporcell_correlation_limit",
            str(p["souporcell_correlation_limit"]),
        ]
    if tool_enabled("vireo"):
        parts += ["--vireo", tool_dir("vireo")]
    return " ".join(parts)


def demuxafy_ref():
    """Combine_Results.R --ref must name an ENABLED tool. Prefer Vireo, else
    Souporcell, else the first enabled tool (capitalised, e.g. Demuxlet)."""
    if tool_enabled("vireo"):
        return "Vireo"
    if tool_enabled("souporcell"):
        return "Souporcell"
    return enabled_tools()[0].capitalize()


# --- report helpers --------------------------------------------------------------
def freemuxlet_corr_fig(wildcards):
    """Path to the freemuxlet reference-cluster correlation plot for the capture
    page's method-diagnostics section.

    report_capture.Rmd hard-codes this figure in its `figs` list, so the input
    key must always resolve to an EXISTING file (a NULL/absent key would crash
    the Rmd on normalizePath). When freemuxlet is enabled we point at its plot;
    when it is disabled (no popscle pileup) we fall back to the already-required,
    always-present per-capture metrics TSV, so the report still builds and the
    freemuxlet panel simply degrades to a placeholder. The .Rmd is intentionally
    left untouched."""
    c = wildcards.capture_id
    if tool_enabled("freemuxlet"):
        return str(out_dir / "freemuxlet" / c / "ref_clust_pearson_correlation.png")
    return str(report_dir / "metrics" / "{}.tsv".format(c))


# --- rule all ----------------------------------------------------------------------
# Build the `rule all` target list, requiring each per-tool aggregate sentinel
# ONLY if that tool is enabled (config params.tools). The demuxafy/hierarchical/
# diagnostics sentinels and the reports are always required (demuxafy adapts to
# the enabled subset internally). This is the single gate that keeps disabled
# tools out of the DAG. Called from `rule all` in the Snakefile, after every
# rule is defined, so referencing rules.* here is safe.
def rule_all_targets():
    targets = []
    _per_tool_aggregates = {
        "demuxlet": rules.aggregate_demuxlet.output.checkpoint_file,
        "freemuxlet": rules.aggregate_freemuxlet.output.checkpoint_file,
        "vireo": rules.aggregate_vireo.output.checkpoint_file,
        "souporcell": rules.aggregate_souporcell.output.checkpoint_file,
    }
    for tool in ["souporcell", "vireo", "demuxlet", "freemuxlet"]:
        if tool_enabled(tool):
            targets.append(_per_tool_aggregates[tool])
    targets += [
        rules.aggregate_demuxafy.output.checkpoint_file,
        rules.aggregate_hier.output.checkpoint_file,
        rules.aggregate_diagnostics.output.checkpoint_file,
        rules.aggregate_present_absent.output.checkpoint_file,
        rules.aggregate_negative_test.output.checkpoint_file,
        out_dir / "assignment/report/all_singlets_stats_full.tsv",
        out_dir / "assignment/report/all_cells_stats_full.tsv",
        out_dir / "assignment/report/combined_results.tsv",
        out_dir / "assignment/report/doublet_summary.csv",
        report_dir / "index.html",
    ]
    return targets
