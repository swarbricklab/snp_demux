import io

import pytest

from conftest import load_script

ap = load_script("add_posteriors")

SAMPLES = """sample_id,donor_id,capture_id
007,D7,cap1
S2,D2,cap1
S3,D3,cap1
X9,D9,cap2
"""

# The consensus table (donor space). Souporcell clusters: 0 -> D7, 1 -> D2;
# cluster 2 is never named; cluster 3 is named inconsistently (duplicate rows).
TABLE = """Barcode\tVireo_Individual_Assignment\tSouporcell_Cluster\tSouporcell_Individual_Assignment\tSouporcell_DropletType
AAA-1\tD7\t0\tD7\tsinglet
CCC-1\tdoublet\tdoublet\tdoublet\tdoublet
GGG-1\tunassigned\t2\tdoublet\tdoublet
TTT-1\tD2\t1\tD2\tsinglet
NNN-1\tdonor5\tunassigned\tunassigned\tunassigned
DUP-1\tD3\t3\tD3\tsinglet
DUP-1\tD3\t3\tD2\tsinglet
"""

VIREO = """cell\tdonor_id\tprob_max\tprob_doublet\tn_vars\tbest_singlet\tbest_doublet\tdoublet_logLikRatio
AAA-1\tD7\t9.99e-01\t1.2e-05\t40\tD7\tD7,D2\t-3.1
CCC-1\tdoublet\t0.02\t0.97\t35\tD2\tD2,D3\t4.0
GGG-1\tunassigned\t0.4\t0.3\t2\tD3\tD7,D3\t0.1
TTT-1\tD2\t1\t0\t50\tD2\tD2,D3\t-5
NNN-1\tunassigned\t0.5\t0.1\t3\tdonor5\tD2,D3\t-1
"""

# troublet 2.5 layout. NNN-1 is missing; GGG-1's best cluster is 2 (unnamed).
SOUPORCELL = """barcode\tstatus\tassignment\tsinglet_posterior\tdoublet_posterior\tlog_prob_singleton\tlog_prob_doublet\tcluster0\tcluster1\tcluster2\tcluster3
AAA-1\tsinglet\t0\t1\t2.5e-30\t-100\t-170\t-100\t-300\t-310\t-320
CCC-1\tdoublet\t0/1\t0.62\t1\t-400\t-250\t-400\t-410\t-500\t-600
GGG-1\tunassigned\t2\t0.55\t0.4\t-90\t-91\t-200\t-210\t-90\t-300
TTT-1\tsinglet\t1\t0.999\t1e-12\t-50\t-80\t-120\t-50\t-130\t-140
DUP-1\tsinglet\t3\t0.98\t0.001\t-60\t-90\t-200\t-210\t-220\t-60
"""


def write(tmp_path, name, text):
    p = tmp_path / name
    p.write_text(text)
    return p


def run(tmp_path, souporcell=SOUPORCELL, table=TABLE):
    out = tmp_path / "out.tsv"
    log = io.StringIO()
    ap.run(
        write(tmp_path, "table.tsv", table),
        write(tmp_path, "donor_ids.tsv", VIREO),
        write(tmp_path, "clusters.tsv", souporcell),
        write(tmp_path, "samples.csv", SAMPLES),
        "cap1",
        out,
        log,
    )
    header, rows = ap.read_table(out)
    return header, rows, log.getvalue()


def by_barcode(rows):
    return {r["Barcode"]: r for r in rows}


def test_input_rows_and_columns_pass_through(tmp_path):
    header, rows, _ = run(tmp_path)
    in_header, in_rows = ap.read_table(write(tmp_path, "t.tsv", TABLE))
    assert header == in_header + ap.COLUMNS
    assert len(rows) == len(in_rows)
    for r, i in zip(rows, in_rows):
        assert {k: r[k] for k in in_header} == i


def test_vireo_values_verbatim(tmp_path):
    rows = by_barcode(run(tmp_path)[1])
    assert rows["AAA-1"]["Vireo_prob_max"] == "9.99e-01"
    assert rows["AAA-1"]["Vireo_prob_doublet"] == "1.2e-05"
    assert rows["CCC-1"]["Vireo_prob_doublet"] == "0.97"
    assert rows["AAA-1"]["Vireo_best_singlet"] == "D7"


def test_vireo_best_singlet_not_a_capture_donor_is_na(tmp_path):
    rows = by_barcode(run(tmp_path)[1])
    assert rows["NNN-1"]["Vireo_best_singlet"] == "NA"
    assert rows["NNN-1"]["Vireo_prob_max"] == "0.5"


def test_souporcell_values_verbatim(tmp_path):
    rows = by_barcode(run(tmp_path)[1])
    assert rows["AAA-1"]["Souporcell_singlet_posterior"] == "1"
    assert rows["AAA-1"]["Souporcell_doublet_posterior"] == "2.5e-30"
    assert rows["TTT-1"]["Souporcell_doublet_posterior"] == "1e-12"


def test_souporcell_best_singlet_named_from_table(tmp_path):
    rows = by_barcode(run(tmp_path)[1])
    assert rows["AAA-1"]["Souporcell_best_singlet"] == "D7"
    assert rows["TTT-1"]["Souporcell_best_singlet"] == "D2"
    # A doublet's best single cluster is still reported (cluster 0 -> D7).
    assert rows["CCC-1"]["Souporcell_best_singlet"] == "D7"


def test_souporcell_unnamed_cluster_is_na_not_a_cluster_number(tmp_path):
    rows = by_barcode(run(tmp_path)[1])
    assert rows["GGG-1"]["Souporcell_best_singlet"] == "NA"


def test_souporcell_inconsistently_named_cluster_is_na(tmp_path):
    rows = run(tmp_path)[1]
    assert [r["Souporcell_best_singlet"] for r in rows if r["Barcode"] == "DUP-1"] == [
        "NA",
        "NA",
    ]


def test_missing_barcode_gets_na_and_is_logged(tmp_path):
    _, rows, log = run(tmp_path)
    r = by_barcode(rows)["NNN-1"]
    assert [r[c] for c in ap.SOUPORCELL_COLUMNS] == ["NA", "NA", "NA"]
    assert "1 barcodes have no souporcell row" in log


def test_cluster_names_read_from_table():
    rows = [
        {"Souporcell_Cluster": "0", "Souporcell_Individual_Assignment": "D7"},
        {"Souporcell_Cluster": "0", "Souporcell_Individual_Assignment": "D7"},
        {"Souporcell_Cluster": "1", "Souporcell_Individual_Assignment": "doublet"},
        {"Souporcell_Cluster": "2", "Souporcell_Individual_Assignment": "D2"},
        {"Souporcell_Cluster": "2", "Souporcell_Individual_Assignment": "D3"},
    ]
    assert ap.souporcell_cluster_names(rows, {"D7", "D2", "D3"}) == {"0": "D7"}


def test_old_troublet_without_posteriors_fails_loudly(tmp_path):
    old = "\n".join(
        "\t".join(f for i, f in enumerate(line.split("\t")) if i not in (3, 4))
        for line in SOUPORCELL.splitlines()
    )
    with pytest.raises(ValueError, match="singlet_posterior"):
        run(tmp_path, souporcell=old + "\n")


def test_existing_posterior_columns_rejected(tmp_path):
    lines = TABLE.splitlines()
    lines[0] += "\tVireo_prob_max"
    lines[1:] = [line + "\tx" for line in lines[1:]]
    with pytest.raises(ValueError, match="already has"):
        run(tmp_path, table="\n".join(lines) + "\n")
