"""Unit tests for workflow/scripts/consensus.py (the Demuxafy-equivalent port).

Output-level equivalence with Demuxafy is checked by the golden comparison and
by .tests/synthetic/consensus_vs_demuxafy.py; these pin the individual rules.
"""
from conftest import load_script

c = load_script("consensus")


def popscle(rows):
    return [{"BARCODE": b, "DROPLET.TYPE": t, "BEST.GUESS": g} for b, t, g in rows]


def soup(rows):
    return [{"barcode": b, "status": s, "assignment": a} for b, s, a in rows]


def vireo(rows):
    return [{"cell": b, "donor_id": d} for b, d in rows]


def test_popscle_parsing():
    dem = c.read_popscle(
        popscle([("A", "SNG", "d1,d1,0.00"), ("B", "DBL", "d1,d2,0.50"), ("C", "AMB", "d1,d2,0.00")]),
        "Individual_Assignment",
    )
    assert dem["A"] == {"DropletType": "singlet", "Individual_Assignment": "d1"}
    assert dem["B"] == {"DropletType": "doublet", "Individual_Assignment": "doublet"}
    assert dem["C"] == {"DropletType": "unassigned", "Individual_Assignment": "unassigned"}
    fmx = c.read_popscle(popscle([("A", "SNG", "3,3"), ("B", "DBL", "3,1")]), "Cluster")
    assert fmx["A"] == {"DropletType": "singlet", "Cluster": "3"}
    assert fmx["B"] == {"DropletType": "doublet", "Cluster": "doublet"}


def test_vireo_and_souporcell_parsing():
    v = c.read_vireo(vireo([("A", "d1"), ("B", "doublet"), ("C", "unassigned")]))
    assert [v[b]["DropletType"] for b in "ABC"] == ["singlet", "doublet", "unassigned"]
    s = c.read_souporcell(soup([("A", "singlet", "0"), ("B", "doublet", "0/1"), ("C", "unassigned", "1")]))
    assert [s[b]["Cluster"] for b in "ABC"] == ["0", "doublet", "unassigned"]


def test_cluster_map_majority_and_threshold():
    bcs = list("ABCDE")
    clusters = dict(zip(bcs, ["0", "0", "0", "1", "1"]))
    ref = dict(zip(bcs, ["d1", "d1", "d2", "d2", "d3"]))
    m = c.cluster_map(bcs, clusters, ref, 0.6)
    assert m == {"0": "d1"}  # cluster 1 splits 50/50 -> below 0.6, unmapped


def test_cluster_map_tie_goes_to_pair_seen_last():
    # At pct 0.5 cluster 1 has two candidates with equal counts; Demuxafy's
    # rev(order(Count)) puts the later-seen pair first and downstream kept it.
    bcs = list("ABCD")
    clusters = dict(zip(bcs, ["1", "1", "1", "1"]))
    ref = dict(zip(bcs, ["d2", "d3", "d3", "d2"]))
    assert c.cluster_map(bcs, clusters, ref, 0.5) == {"1": "d3"}


def test_cluster_map_counts_missing_reference_as_a_value():
    bcs = list("ABC")
    clusters = {"A": "0", "B": "0", "C": "0"}
    ref = {"C": "d1"}  # A and B absent from vireo
    assert c.cluster_map(bcs, clusters, ref, 0.6) == {"0": None}


def test_unmapped_cluster_becomes_doublet():
    calls = c.read_souporcell(soup([("A", "singlet", "0"), ("B", "unassigned", "1"), ("C", "singlet", "2")]))
    c.assign_clusters(calls, {"0": "d1", "2": None})
    assert calls["A"] == {"DropletType": "singlet", "Cluster": "0", "Individual_Assignment": "d1"}
    # Demuxafy behaviour preserved here; see #65 for why this is wrong.
    assert calls["B"]["DropletType"] == calls["B"]["Individual_Assignment"] == "doublet"
    assert calls["C"]["DropletType"] == "doublet"


def test_plurality():
    assert c.plurality(["d1", "d1", "d2", "unassigned"]) == "d1"
    assert c.plurality(["d1", "d2", "doublet", "unassigned"]) is None
    assert c.plurality(["doublet", "unassigned"]) is None


def row(types, donors):
    r = {}
    for m, t, d in zip(c.METHODS, types, donors):
        r[m + "_DropletType"], r[m + "_Individual_Assignment"] = t, d
    return r


def test_at_least_half_singlet():
    S, D, U = "singlet", "doublet", "unassigned"
    assert c.at_least_half_singlet(row([S, S, D, D], ["a", "a", D, D])) == (S, "a")
    assert c.at_least_half_singlet(row([S, U, U, D], ["a", U, U, D])) == (D, D)
    assert c.at_least_half_singlet(row([S, S, U, U], ["a", "b", U, U])) == (S, U)


def test_build_end_to_end():
    dem = popscle([("A", "SNG", "d1,d1,0.00"), ("B", "SNG", "d2,d2,0.00"), ("C", "AMB", "d1,d2,0.00")])
    fmx = popscle([("A", "SNG", "0,0"), ("B", "SNG", "1,1")])  # C missing
    sop = soup([("A", "singlet", "5"), ("B", "singlet", "6"), ("C", "unassigned", "5")])
    vir = vireo([("A", "d1"), ("B", "d2"), ("C", "d1")])
    snp = soup([("A", "singlet", "d1"), ("B", "singlet", "d2")])
    out = c.build(dem, fmx, sop, vir, snp, pct_agreement=0.7)
    assert [r["Barcode"] for r in out] == ["A", "B", "C"]
    assert set(out[0]) == set(c.COLUMNS)
    a, b, cc = out
    assert a["AtLeastHalfSinglet_Individual_Assignment"] == "d1"
    assert b["Freemuxlet_Individual_Assignment"] == "d2"
    assert cc["Freemuxlet_DropletType"] == "unassigned"  # missing barcode
    # Demuxafy treats the pseudo-cluster "unassigned" like any cluster: here all
    # of its barcodes are vireo d1, so it "maps" to d1 (droplet type unchanged).
    assert cc["Souporcell_DropletType"] == "unassigned"
    assert cc["Souporcell_Individual_Assignment"] == "d1"
    assert cc["AtLeastHalfSinglet_DropletType"] == "doublet"
    assert cc["Souporcell_vireosnp_DropletType"] == "NA"
