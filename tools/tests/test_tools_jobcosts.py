"""demux_tools.jobcosts: the PBS Resource Usage parser, the charge rule and the job-costs command."""
import json

import pandas as pd
import pytest

from demux_tools import jobcosts as J

BAR = "=" * 86


def _block(job, su, ncpus, cpu, mem_req, mem_used, wall_req, wall, exit_status=0, when="2026-01-02 03:04:05",
           drop=()):
    lines = [BAR, "                  Resource Usage on %s:" % when, "   Job Id:             %s.gadi-pbs" % job,
             "   Project:            p00", "   Exit Status:        %d" % exit_status,
             "   Service Units:      %s" % su,
             "   NCPUs Requested:    %d                  CPU Time Used: %s        " % (ncpus, cpu),
             "   Memory Requested:   %s               Memory Used: %s          " % (mem_req, mem_used),
             "   Walltime Requested: %s            Walltime Used: %s        " % (wall_req, wall),
             "   JobFS Requested:    100.0MB                JobFS Used: 0B              ", BAR]
    return "\n".join(l for l in lines if not any(d in l for d in drop)) + "\n"


def _header(rule, cap, job, threads=4, mem_mb=16000):
    return ("\n[Mon Jan  1 00:00:00 2026]\nrule %s:\n    input: data/bam/%s.bam, data/panel/%s.vcf.gz\n"
            "    output: results/%s/%s/out.tsv\n    jobid: 0\n    reason: Missing output files\n"
            "    wildcards: platform=short, capture_id=%s\n"
            "    resources: mem_mb=%d, tmpdir=/jobfs/%s.gadi-pbs, runtime=360, threads=%d\n\n"
            % (rule, cap, cap, rule, cap, cap, mem_mb, job, threads))


def test_one_full_block():
    text = _header("vartrix", "capA", "111") + _block("111", "8.53", 16, "06:24:00", "31.25GB", "18.5GB",
                                                       "12:00:00", "00:24:00")
    (r,) = J.parse_log(text)
    assert (r["rule"], r["capture"], r["job_id"], r["exit_status"], r["ncpus"]) == ("vartrix", "capA", "111.gadi-pbs", 0, 16)
    assert r["cpu_h"] == pytest.approx(6.4) and r["walltime_h"] == pytest.approx(0.4)
    assert r["mem_req_gb"] == pytest.approx(31.25) and r["jobfs_req_gb"] == pytest.approx(100 / 1024.0)
    assert r["su"] == 8.53 and r["first_input"] == "data/bam/capA.bam" and r["first_output"] == "results/vartrix/capA/out.tsv"
    assert (r["smk_threads"], r["smk_mem_mb"], r["smk_runtime_min"]) == (4, 16000, 360)
    assert not r["header_lost"] and r["partial_block"] == ""


def test_partial_block_and_kill():
    text = (_header("scan", "capB", "222") + "PBS: job killed: walltime 21650 exceeded limit 21600\n"
            + _block("222", "48.1", 1, "06:00:50", "15.62GB", "15.62GB", "06:00:00", "06:00:50", exit_status=271,
                     drop=("Memory Used",)))
    (r,) = J.parse_log(text)
    assert r["pbs_kill"] == "walltime 21650 exceeded limit 21600" and r["exit_status"] == 271
    assert r["partial_block"] == "Memory Requested,Memory Used"   # the line holds both fields
    j = J.add_charges(pd.DataFrame([dict(r, source="s", log="x.log", label="")]))
    s = J.summarise(j)
    assert s.failed_on_walltime.tolist() == [1] and s.failed_su.tolist() == [48.1]


def test_two_blocks_one_header():
    # PBS appended a second job's block; the header names job 333 in its tmpdir
    text = (_header("cellsnp", "capC", "333") + _block("444", "1.0", 1, "00:10:00", "4.0GB", "1GB", "01:00:00",
                                                       "00:15:00") + _block("333", "2.0", 1, "00:20:00", "4.0GB",
                                                                            "1GB", "01:00:00", "00:30:00"))
    a, b = J.parse_log(text)
    assert (a["job_id"], a["header_lost"], a["capture"]) == ("444.gadi-pbs", True, "")
    assert (b["job_id"], b["header_lost"], b["capture"]) == ("333.gadi-pbs", False, "capC")
    assert a["rule"] == b["rule"] == "cellsnp" and (a["block"], b["block"], b["n_blocks"]) == (1, 2, 2)


def test_err_file_gives_the_header(tmp_path):
    d = tmp_path / "logs" / "souporcell"
    d.mkdir(parents=True)
    (d / "souporcell.7.log").write_text(_block("555", "0.5", 2, "00:20:00", "8.0GB", "2GB", "01:00:00", "00:15:00"))
    (d / "souporcell.7.err").write_text(_header("souporcell", "capD", "555") + "PBS: job killed: mem 9GB\n")
    jobs = J.read_logs([("run", str(tmp_path / "logs"))])
    r = jobs.iloc[0]
    assert (r.rule, r.capture, r.log, r.pbs_kill) == ("souporcell", "capD", "souporcell/souporcell.7.log", "mem 9GB")


def test_qxub_log_without_header_takes_the_file_stem(tmp_path):
    (tmp_path / "myjob_20260101_000000.log").write_text("some output\n" + _block("666", "0.02", 1, "00:00:21",
                                                         "4.0GB", "0.1GB", "00:30:00", "00:00:40"))
    (tmp_path / "notes.log").write_text("no block here\n")
    jobs = J.read_logs([("qx", str(tmp_path))], pattern="*.log")
    assert len(jobs) == 1 and jobs.rule.tolist() == ["myjob_20260101_000000"] and jobs.capture.tolist() == [""]


def test_charge_rule():
    su, charged, basis = J.charge("normal", 1, 16.0, 2.0)
    assert (su, charged, basis) == (16.0, 4.0, "memory")    # 16 GB on normal is 4 CPUs' worth
    su, charged, basis = J.charge("normal", 8, 8.0, 0.5)
    assert (su, charged, basis) == (8.0, 8.0, "cpus")
    assert J.charge("normalsr", 1, 4.8, 1.0)[0] == 2.0      # 4.8 GB fits one normalsr CPU (500 GB / 104)
    su, charged, _ = J.charge("normalsr", 1, 9.6, 1.0)
    assert charged == pytest.approx(9.6 / (500 / 104.0)) and su == pytest.approx(2.0 * charged)
    with pytest.raises(ValueError):
        J.charge("normalbw", 1, 4.0, 1.0)
    assert J.charge("normalbw", 2, 4.0, 1.0, node=(28, 128.0))[0] == pytest.approx(2.5)
    assert J.gigabytes("1024MB") == 1.0 and J.hours("26:30:00") == 26.5


def test_job_costs_command(tmp_path, run_tools):
    root = tmp_path / "joblogs"
    for rule, jobs in (("cellsnp", [("c1", "11", "1.5"), ("c2", "12", "2.5")]), ("vartrix", [("c1", "21", "4.0")])):
        (root / rule).mkdir(parents=True)
        for n, (cap, job, su) in enumerate(jobs):
            (root / rule / ("%s.%d.log" % (rule, n))).write_text(
                _header(rule, cap, job) + _block(job, su, 4, "01:00:00", "16.0GB", "3GB", "02:00:00", "00:30:00"))
    out = tmp_path / "costs"
    res = run_tools("job-costs", "--source", "run=%s" % root, "--label", r"results/([^/]+)/", "--out", out)
    assert res.returncode == 0, res.stderr
    jobs = pd.read_csv(str(out / "jobs.tsv"), sep="\t")
    assert len(jobs) == 3 and set(jobs.label) == {"cellsnp", "vartrix"}
    assert jobs.su_expected.tolist() == [4.0, 4.0, 4.0]   # 2 SU x 4 CPUs x 0.5 h
    rules = pd.read_csv(str(out / "rule_costs.tsv"), sep="\t").set_index("rule")
    assert rules.at["cellsnp", "total_su"] == 4.0 and rules.at["cellsnp", "captures"] == 2
    caps = pd.read_csv(str(out / "capture_costs.tsv"), sep="\t").set_index("capture")
    assert caps.at["c1", "su"] == 5.5 and caps.at["c1", "jobs"] == 2
    prov = json.load(open(str(out / "provenance.json")))
    assert prov["results"]["totals"]["jobs"] == 3
    res = run_tools("job-costs", "--source", "run=%s" % root, "--group-by", "rule,capture", "--out", tmp_path / "c2")
    assert res.returncode == 0, res.stderr
    assert len(pd.read_csv(str(tmp_path / "c2" / "rule_costs.tsv"), sep="\t")) == 3


def test_truncated_and_partial_blocks_are_counted_not_charged_as_zero(tmp_path, run_tools):
    d = tmp_path / "logs" / "r"
    d.mkdir(parents=True)
    (d / "full.1.log").write_text(_header("full", "c1", "11") + _block("11", "2.0", 1, "00:30:00", "4.0GB", "1GB",
                                                                       "01:00:00", "01:00:00"))
    cut = _block("22", "0.5", 4, "01:00:00", "16.0GB", "2GB", "01:00:00", "01:00:05", exit_status=271)
    cut = cut[:cut.index("   Service Units")]           # the log ends before the block does
    (d / "cut.2.log").write_text(_header("cut", "c2", "22") + "PBS: job killed: walltime 3605 exceeded limit 3600\n"
                                 + cut)
    part = _block("33", "0.5", 4, "01:00:00", "16.0GB", "2GB", "01:00:00", "01:00:05", drop=("Memory Used",))
    (d / "part.3.log").write_text(_header("part", "c3", "33") + part)
    (d / "none.4.log").write_text("no block here\n")
    jobs = J.add_charges(J.read_logs([("s", str(tmp_path / "logs"))]))
    j = jobs.set_index("rule")
    assert len(jobs) == 3
    assert j.at["cut", "partial_block"].startswith("truncated")
    assert j.at["cut", "pbs_kill"] == "walltime 3605 exceeded limit 3600" and j.at["cut", "exit_status"] == 271
    assert j.at["part", "charge_basis"] == "" and pd.isna(j.at["part", "su_expected"])
    assert j.at["full", "charge_basis"] == "cpus" and j.at["full", "su_expected"] == 2.0
    rules = J.summarise(jobs).set_index("rule")
    assert pd.isna(rules.at["part", "su_expected"]) and rules.at["part", "partial_blocks"] == 1
    assert pd.isna(rules.at["cut", "total_su"]) and rules.at["full", "partial_blocks"] == 0
    res = run_tools("job-costs", "--source", "s=%s" % (tmp_path / "logs"), "--out", tmp_path / "out")
    assert res.returncode == 0, res.stderr
    assert "2 partial (1 cut before their closing line)" in res.stdout
    prov = json.load(open(str(tmp_path / "out" / "provenance.json")))
    assert prov["results"]["totals"]["truncated_blocks"] == 1
