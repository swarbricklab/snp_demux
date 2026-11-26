"""Service units per job, rule and capture from PBS job logs; the `job-costs` command.

A PBS job's log ends with a "Resource Usage" block:

    ======================================================================================
                      Resource Usage on 2026-01-02 03:04:05:
       Job Id:             123456.gadi-pbs
       Exit Status:        0
       Service Units:      0.64
       NCPUs Requested:    32                  CPU Time Used: 00:07:23
       Memory Requested:   31.25GB               Memory Used: 6.13GB
       Walltime Requested: 24:00:00            Walltime Used: 00:00:36
       JobFS Requested:    100.0MB                JobFS Used: 0B
    ======================================================================================

Any file with this block is read: a Snakemake job log written through a PBS
profile, or a qxub log. Each block is one job. A Snakemake log also holds the
job's header (`rule NAME:`, `input:`, `output:`, `wildcards:`, `reason:`,
`resources:`), read from the log or else from the `.err` file beside it.

A later job that writes the same log file name normally replaces the file,
so the earlier job's record is lost. When two such jobs end close together,
the second job's block can land in the first job's file instead, under the
first job's header. The header's `tmpdir=/jobfs/<job id>` names its own job,
so that block takes the header and any other block is marked header_lost (no
rule fields beyond the file name).

A block cut before its closing line is kept and marked "truncated" in
partial_block; a block missing fields lists them there. Such a block has no
charge basis and no su_expected, and the per-rule and per-capture tables
count partial blocks, so a missing figure is never added as 0.

The charge: SU = queue rate x max(CPUs requested, memory requested / a node's
memory x a node's CPUs) x walltime used. PBS's own "Service Units" is kept as
`su`; `su_expected` recomputes it from the queue given, as a check.
"""
import glob
import os
import re
import warnings

import numpy as np
import pandas as pd

from demux_tools._io import UsageError, check_output_path, start_clock, write_provenance, write_table

# queue: (SU per CPU hour, CPUs per node, GB per node); None where the node size is not recorded here
QUEUES = {
    "normal": (2.0, 48, 192.0),
    "normalsr": (2.0, 104, 500.0),
    "normalsl": (1.5, None, None),
    "normalbw": (1.25, None, None),
    "express": (6.0, None, None),
}
# A block runs to its closing line of "=", or to the end of the text when the log was cut short.
BLOCK = re.compile(r"Resource Usage on (\S+ \S+):(.*?)(={10}|\Z)", re.S)
KILL = re.compile(r"PBS: job killed: ([^\n]+)")
RULE = re.compile(r"^(?:local)?rule (\S+):\s*$", re.M)
FIELDS = ("Job Id", "Exit Status", "Service Units", "NCPUs Requested", "CPU Time Used", "Memory Requested",
          "Memory Used", "Walltime Requested", "Walltime Used", "JobFS Requested", "JobFS Used")
JOB_COLUMNS = ["source", "log", "block", "n_blocks", "header_lost", "rule", "capture", "wildcards", "first_input",
               "first_output", "label", "reason", "job_id", "finished", "exit_status", "pbs_kill", "ncpus",
               "mem_req_gb", "walltime_req_h", "walltime_h", "cpu_h", "mem_used_gb", "jobfs_req_gb",
               "jobfs_used_gb", "su", "cpu_share", "charged_cpus", "charge_basis", "su_rate", "su_expected",
               "partial_block", "smk_threads", "smk_mem_mb", "smk_runtime_min", "smk_resources"]


def hours(text):
    """"HH:MM:SS" (hours may exceed 24) in hours; NaN for None."""
    if text is None:
        return float("nan")
    h, m, s = (int(x) for x in text.split(":"))
    return h + m / 60.0 + s / 3600.0


def gigabytes(text):
    """"31.25GB", "100.0MB", "0B", ... in GB (1 GB = 1024 MB); NaN for None."""
    if text is None:
        return float("nan")
    m = re.match(r"([\d.]+)([KMGT]?B)$", text)
    if not m:
        raise ValueError("not a PBS size: %r" % text)
    scale = {"B": 1.0 / 1024 ** 3, "KB": 1.0 / 1024 ** 2, "MB": 1.0 / 1024, "GB": 1.0, "TB": 1024.0}
    return float(m.group(1)) * scale[m.group(2)]


def charge(queue, ncpus, mem_gb, walltime_h, node=None):
    """(su, charged_cpus, basis) for a request on `queue`; basis is "cpus" or "memory".

    node: (CPUs, GB) of a node, overriding the queue table. Raises ValueError
    for a queue whose node size is unknown and no node given.
    """
    if queue not in QUEUES:
        raise ValueError("unknown queue %r (known: %s)" % (queue, ", ".join(sorted(QUEUES))))
    rate, cpus, gb = QUEUES[queue]
    if node is not None:
        cpus, gb = node
    if cpus is None:
        raise ValueError("the node size of queue %s is not recorded; give --node CPUS:GB" % queue)
    mem_cpus = mem_gb / (gb / float(cpus))
    charged = max(float(ncpus), mem_cpus)
    return rate * charged * walltime_h, charged, ("memory" if mem_cpus > ncpus else "cpus")


def _header(text, key):
    m = re.search(r"\n\s+" + key + r": ([^\n]+)", text)
    return m.group(1).strip() if m else ""


def _int(pattern, text):
    m = re.search(pattern, text)
    return int(m.group(1)) if m else None


def parse_header(text):
    """The Snakemake job header of a log: rule, input, output, wildcards, reason, resources (text)."""
    rules = RULE.findall(text)
    return {"rule": rules[-1] if rules else "", "input": _header(text, "input"), "output": _header(text, "output"),
            "wildcards": _header(text, "wildcards"), "reason": _header(text, "reason"),
            "resources": _header(text, "resources")}


def parse_blocks(text):
    """[(finished, {field: text})] for every Resource Usage block, in file order.

    A block cut before its closing line is kept, with the key "_truncated" set.
    """
    out = []
    for finished, body, end in BLOCK.findall(text):
        f = {}
        for name in FIELDS:
            m = re.search(re.escape(name) + r":\s+(\S+)", body)
            f[name] = m.group(1) if m else None
        f["_truncated"] = not end
        out.append((finished, f))
    return out


def parse_log(text, err_text="", capture_keys=("capture_id", "capture")):
    """One dict per PBS block of a log text (empty list when it has none). See the module docstring."""
    blocks = parse_blocks(text)
    if not blocks:
        return []
    head = parse_header(text)
    if err_text and not (head["wildcards"] or head["input"] or head["output"]):
        head = dict(parse_header(err_text), rule=head["rule"] or parse_header(err_text)["rule"])
    kills = sorted(set(k.strip() for k in KILL.findall(text) + KILL.findall(err_text or "")))
    res = head["resources"]
    tmp_job = re.search(r"tmpdir=/jobfs/(\d+)", res)
    cap = None
    for k in capture_keys:
        cap = re.search(r"(?:^|[,\s])" + re.escape(k) + r"=([^,\s]+)", head["wildcards"])
        if cap:
            break
    rows = []
    for i, (finished, f) in enumerate(blocks):
        job = f["Job Id"] or ""
        own = tmp_job is None or job.split(".")[0] == tmp_job.group(1)
        missing = [k for k in FIELDS if f[k] is None and not k.startswith("JobFS")]
        ncpus = int(f["NCPUs Requested"]) if f["NCPUs Requested"] else None
        wall = hours(f["Walltime Used"])
        cpu = hours(f["CPU Time Used"])
        rows.append({
            "block": i + 1, "n_blocks": len(blocks), "header_lost": not own,
            "rule": head["rule"],
            "capture": (cap.group(1) if cap else "") if own else "",
            "wildcards": head["wildcards"] if own else "",
            "first_input": head["input"].split(",")[0].strip() if own else "",
            "first_output": head["output"].split(",")[0].strip() if own else "",
            "reason": head["reason"] if own else "",
            "job_id": job, "finished": finished,
            "exit_status": int(f["Exit Status"]) if f["Exit Status"] is not None else None,
            "pbs_kill": " | ".join(kills) if own else "",
            "ncpus": ncpus,
            "mem_req_gb": gigabytes(f["Memory Requested"]),
            "walltime_req_h": hours(f["Walltime Requested"]),
            "walltime_h": wall, "cpu_h": cpu,
            "mem_used_gb": gigabytes(f["Memory Used"]),
            "jobfs_req_gb": gigabytes(f["JobFS Requested"]), "jobfs_used_gb": gigabytes(f["JobFS Used"]),
            "su": float(f["Service Units"]) if f["Service Units"] is not None else float("nan"),
            "partial_block": ",".join((["truncated"] if f["_truncated"] else []) + missing),
            "smk_threads": _int(r"threads=(\d+)", res) if own else None,
            "smk_mem_mb": _int(r"mem_mb=(\d+)", res) if own else None,
            "smk_runtime_min": _int(r"runtime=(\d+)", res) if own else None,
            "smk_resources": res if own else "",
        })
    return rows


def read_logs(sources, pattern="**/*.log", capture_keys=("capture_id", "capture"), label=None):
    """Every PBS block under each source folder, as a jobs DataFrame (JOB_COLUMNS, no charge columns yet).

    sources: [(name, folder)]. pattern: glob below each folder (recursive "**").
    label: a regular expression searched in the first output (else the first
    input); the label is its groups joined by "/", or the whole match.
    """
    lab = re.compile(label) if label else None
    rows = []
    for name, root in sources:
        for path in sorted(glob.glob(os.path.join(root, pattern), recursive=True)):
            if not os.path.isfile(path):
                continue
            with open(path, errors="replace") as f:
                text = f.read()
            err = os.path.splitext(path)[0] + ".err"
            err_text = ""
            if os.path.exists(err) and err != path:
                with open(err, errors="replace") as f:
                    err_text = f.read()
            for r in parse_log(text, err_text, capture_keys):
                if not r["rule"]:
                    r["rule"] = os.path.basename(path).split(".")[0]
                r["source"] = name
                r["log"] = os.path.relpath(path, root)
                r["label"] = ""
                if lab is not None:
                    m = lab.search(r["first_output"] or r["first_input"])
                    if m:
                        groups = [g for g in m.groups() if g]
                        r["label"] = "/".join(groups) if m.groups() else m.group(0)
                rows.append(r)
    return pd.DataFrame(rows, columns=[c for c in JOB_COLUMNS if c not in
                                       ("cpu_share", "charged_cpus", "charge_basis", "su_rate", "su_expected")])


def add_charges(jobs, queue="normal", node=None):
    """Add cpu_share, charged_cpus, charge_basis, su_rate and su_expected; round as the tables print them."""
    j = jobs.copy()
    j["cpu_share"] = (j.cpu_h / (j.walltime_h * j.ncpus)).where(j.walltime_h > 0)
    rate, cpus, gb = QUEUES[queue]
    if node is not None:
        cpus, gb = node
    if cpus is None:
        raise UsageError("the node size of queue %s is not recorded; give --node CPUS:GB" % queue)
    mem_cpus = j.mem_req_gb / (gb / float(cpus))
    j["charged_cpus"] = np.maximum(j.ncpus.astype(float), mem_cpus)
    known = j.ncpus.notna() & j.mem_req_gb.notna()
    j["charge_basis"] = np.where(~known, "", np.where(mem_cpus > j.ncpus, "memory", "cpus"))
    j["su_rate"] = (j.su / (j.walltime_h * j.charged_cpus)).where(j.walltime_h > 0)
    j["su_expected"] = rate * j.charged_cpus * j.walltime_h
    for c, n in (("mem_req_gb", 3), ("walltime_req_h", 4), ("walltime_h", 5), ("cpu_h", 5), ("mem_used_gb", 3),
                 ("jobfs_req_gb", 3), ("jobfs_used_gb", 3), ("cpu_share", 3), ("charged_cpus", 3), ("su_rate", 3),
                 ("su_expected", 3)):
        j[c] = j[c].round(n)
    for c in ("exit_status", "ncpus", "smk_threads", "smk_mem_mb", "smk_runtime_min"):
        j[c] = pd.to_numeric(j[c]).astype("Int64")
    return j[JOB_COLUMNS]


def _range(s):
    s = s.dropna()
    if s.empty:
        return ""
    return "%g" % s.min() if s.min() == s.max() else "%g-%g" % (s.min(), s.max())


def summarise(jobs, by=("source", "rule")):
    """Per group: jobs, failures, request ranges, walltime, CPU, memory and SU.

    A statistic over a group whose values are all missing (e.g. memory used of
    partial blocks) is NaN.
    """
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", message="Mean of empty slice")
        warnings.filterwarnings("ignore", message="All-NaN slice encountered")
        return _summarise(jobs, by)


def _summarise(jobs, by):
    j = jobs.copy()
    j["failed"] = j.exit_status.fillna(-1) != 0
    j["hit_walltime"] = j.pbs_kill.fillna("").str.contains("walltime") | (
        j.failed & (j.walltime_h >= 0.98 * j.walltime_req_h))
    rows = []
    for key, g in j.groupby(list(by), sort=True):
        key = key if isinstance(key, tuple) else (key,)
        caps = g.capture[g.capture != ""]
        rows.append(dict(zip(by, key), **{
            "jobs": len(g), "failed": int(g.failed.sum()), "failed_on_walltime": int(g.hit_walltime.sum()),
            "captures": int(caps.nunique()),
            "ncpus": _range(g.ncpus.astype(float)), "mem_req_gb": _range(g.mem_req_gb.round(2)),
            "walltime_req_h": _range(g.walltime_req_h.round(2)),
            "median_walltime_min": round(g.walltime_h.median() * 60, 1),
            "max_walltime_min": round(g.walltime_h.max() * 60, 1),
            "total_walltime_h": round(g.walltime_h.sum(), 2), "total_cpu_h": round(g.cpu_h.sum(), 2),
            "median_cpu_share": round(g.cpu_share.median(), 3),
            "median_mem_used_gb": round(g.mem_used_gb.median(), 2),
            "p90_mem_used_gb": round(g.mem_used_gb.quantile(0.9), 2),
            "max_mem_used_gb": round(g.mem_used_gb.max(), 2),
            "median_su": round(g.su.median(), 2), "total_su": round(g.su.sum(min_count=1), 2),
            "failed_su": round(g.su[g.failed].sum(), 2), "su_expected": round(g.su_expected.sum(min_count=1), 2),
            "partial_blocks": int((g.partial_block != "").sum()),
            "first": g.finished.min()[:10], "last": g.finished.max()[:10]}))
    return pd.DataFrame(rows)


def capture_costs(jobs, by=("source", "capture")):
    """Per source and capture (jobs with a capture only): jobs, SU, CPU hours, walltime, failed SU, partial blocks."""
    j = jobs[jobs.capture != ""].copy()
    j["failed"] = j.exit_status.fillna(-1) != 0
    rows = []
    for key, g in j.groupby(list(by), sort=True):
        key = key if isinstance(key, tuple) else (key,)
        rows.append(dict(zip(by, key), jobs=len(g), su=round(g.su.sum(min_count=1), 2),
                         cpu_h=round(g.cpu_h.sum(min_count=1), 4), walltime_h=round(g.walltime_h.sum(min_count=1), 4),
                         failed_su=round(g.su[g.failed].sum(), 2), partial_blocks=int((g.partial_block != "").sum())))
    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------
# Command: job-costs


def configure_job_costs(parser):
    parser.description = ("Read every PBS Resource Usage block in job logs (Snakemake or qxub) and write "
                          "OUT/jobs.tsv (one row per job), OUT/rule_costs.tsv (per source and group), "
                          "OUT/capture_costs.tsv (per source and capture) and OUT/provenance.json.")
    parser.add_argument("--source", action="append", required=True, metavar="NAME=FOLDER",
                        help="a folder of job logs and its name; repeatable")
    parser.add_argument("--glob", default="**/*.log", help="log files below each folder (default **/*.log)")
    parser.add_argument("--group-by", default="rule",
                        help="comma-separated jobs.tsv columns for rule_costs.tsv, after source (default rule)")
    parser.add_argument("--capture-keys", default="capture_id,capture",
                        help="wildcards that name the capture, first found wins (default capture_id,capture)")
    parser.add_argument("--label", default=None,
                        help="regular expression on the first output (else input); its groups make the label column")
    parser.add_argument("--queue", default="normal", help="queue for su_expected (%s)" % ", ".join(sorted(QUEUES)))
    parser.add_argument("--node", default=None, metavar="CPUS:GB", help="node size, overriding the queue's")
    parser.add_argument("--out", required=True, help="output folder (created if needed)")


def run_job_costs(args):
    clock = start_clock()
    sources = []
    for s in args.source:
        if "=" not in s:
            raise UsageError("--source needs NAME=FOLDER, not %r" % s)
        name, root = s.split("=", 1)
        if not os.path.isdir(root):
            raise UsageError("no such folder: %s" % root)
        sources.append((name, root))
    if args.queue not in QUEUES:
        raise UsageError("unknown queue %r (known: %s)" % (args.queue, ", ".join(sorted(QUEUES))))
    node = None
    if args.node:
        c, g = args.node.split(":")
        node = (int(c), float(g))
    by = ["source"] + [c for c in args.group_by.split(",") if c]
    bad = [c for c in by if c not in JOB_COLUMNS]
    if bad:
        raise UsageError("--group-by: no column %s in jobs.tsv" % ", ".join(bad))
    jobs = read_logs(sources, args.glob, tuple(k for k in args.capture_keys.split(",") if k), args.label)
    if jobs.empty:
        raise UsageError("no PBS Resource Usage block found under %s" % ", ".join(r for _, r in sources))
    jobs = add_charges(jobs, args.queue, node)
    os.makedirs(args.out, exist_ok=True)
    outs = [os.path.join(args.out, n) for n in ("jobs.tsv", "rule_costs.tsv", "capture_costs.tsv")]
    for o in outs:
        check_output_path(o, [r for _, r in sources])
    write_table(jobs, outs[0])
    rules = summarise(jobs, by)
    write_table(rules, outs[1])
    write_table(capture_costs(jobs), outs[2])
    totals = {"jobs": len(jobs), "logs": int(jobs[["source", "log"]].drop_duplicates().shape[0]),
              "header_lost": int(jobs.header_lost.sum()), "partial_blocks": int((jobs.partial_block != "").sum()),
              "truncated_blocks": int(jobs.partial_block.str.startswith("truncated").sum()),
              "su": round(float(jobs.su.sum()), 2), "cpu_h": round(float(jobs.cpu_h.sum()), 2)}
    write_provenance(args.out, args.argv, [r for _, r in sources], clock=clock,
                     extra={"totals": totals, "outputs": outs, "queue": args.queue})
    print("job-costs: %(jobs)d PBS blocks from %(logs)d logs, %(su).2f SU, %(cpu_h).2f CPU hours; "
          "%(header_lost)d blocks without their header, %(partial_blocks)d partial (%(truncated_blocks)d cut "
          "before their closing line)" % totals)
    return 0
