"""Unit tests for workflow/scripts/consensus.py."""
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
    mapping, report = c.cluster_map(bcs, clusters, ref, 0.6)
    assert mapping == {"0": "d1"}  # cluster 1 splits 50/50 -> below 0.6, unnamed
    assert [r[-1] for r in report] == ["d1", "unnamed"]


def test_cluster_map_tie_is_ambiguous_not_first_wins():
    bcs = list("ABCD")
    clusters = dict(zip(bcs, ["1"] * 4))
    ref = dict(zip(bcs, ["d2", "d3", "d3", "d2"]))
    mapping, report = c.cluster_map(bcs, clusters, ref, 0.5)
    assert mapping == {}
    assert report[0][-1] == "ambiguous"


def test_cluster_map_never_names_after_reference_sentinels():
    bcs = list("ABCD")
    clusters = dict(zip(bcs, ["0"] * 4))
    ref = {"A": "unassigned", "B": "unassigned", "C": "unassigned", "D": "d1"}
    assert c.cluster_map(bcs, clusters, ref, 0.7)[0] == {}
    # the denominator still counts them: 1 of 4 is not agreement
    ref = {"A": "doublet", "B": "d1", "C": "d1", "D": "d1"}
    assert c.cluster_map(bcs, clusters, ref, 0.7)[0] == {"0": "d1"}


def test_cluster_map_ignores_pseudo_clusters():
    bcs = list("ABC")
    clusters = {"A": "unassigned", "B": "unassigned", "C": "doublet"}
    ref = {"A": "d1", "B": "d1", "C": "d1"}
    assert c.cluster_map(bcs, clusters, ref, 0.5) == ({}, [])


def test_assign_clusters_keeps_the_methods_own_call():
    calls = c.read_souporcell(soup([("A", "singlet", "0"), ("B", "unassigned", "1"),
                                    ("C", "singlet", "2"), ("D", "doublet", "0/2")]))
    c.assign_clusters(calls, {"0": "d1"})
    assert calls["A"] == {"DropletType": "singlet", "Cluster": "0", "Individual_Assignment": "d1"}
    # #65: Demuxafy made both of these doublets
    assert calls["B"]["DropletType"] == calls["B"]["Individual_Assignment"] == "unassigned"
    # a singlet in an unnamed cluster abstains (Demuxafy: doublet)
    assert calls["C"] == {"DropletType": "unassigned", "Cluster": "2", "Individual_Assignment": "unassigned"}
    assert calls["D"]["DropletType"] == calls["D"]["Individual_Assignment"] == "doublet"


def test_plurality():
    assert c.plurality(["d1", "d1", "d2", "unassigned"]) == "d1"
    assert c.plurality(["d1", "d2", "doublet", "unassigned"]) is None
    assert c.plurality(["doublet", "unassigned"]) is None


def row(types, donors):
    r = {}
    for m, t, d in zip(c.METHODS, types, donors):
        r[m + "_DropletType"], r[m + "_Individual_Assignment"] = t, d
    return r


def test_informative_majority():
    S, D, U = "singlet", "doublet", "unassigned"
    im = c.informative_majority
    assert im(row([S, S, D, U], ["a", "a", D, U])) == (S, "a")
    assert im(row([S, D, D, U], ["a", D, D, U])) == (D, D)
    assert im(row([S, S, D, D], ["a", "a", D, D])) == (U, U)  # tie
    assert im(row([D, U, U, U], [D, U, U, U])) == (U, U)  # one informative call
    assert im(row([D, U, U, U], [D, U, U, U]), min_informative=1) == (D, D)
    assert im(row([U, U, U, U], [U, U, U, U])) == (U, U)
    assert im(row([S, S, U, U], ["a", "b", U, U])) == (U, U)  # identity tie
    # the #65 pattern: one doublet vote among unassigned + one singlet
    assert im(row([U, S, D, U], [U, "a", D, U])) == (U, U)


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
    assert (a["Consensus_DropletType"], a["Consensus_Individual_Assignment"]) == ("singlet", "d1")
    assert (a["Consensus_n_singlet"], a["Consensus_n_doublet"], a["Consensus_n_informative"]) == (4, 0, 4)
    assert b["Freemuxlet_Individual_Assignment"] == "d2"
    assert cc["Freemuxlet_DropletType"] == "unassigned"  # missing barcode
    assert cc["Souporcell_DropletType"] == cc["Souporcell_Individual_Assignment"] == "unassigned"
    # only vireo calls C: AtLeastHalfSinglet says doublet, the consensus abstains
    assert cc["AtLeastHalfSinglet_DropletType"] == "doublet"
    assert cc["Consensus_DropletType"] == "unassigned"
    assert cc["Souporcell_vireosnp_DropletType"] == "NA"
    legacy = c.build(dem, fmx, sop, vir, snp, pct_agreement=0.7, method="at_least_half_singlet")
    assert [r["Consensus_DropletType"] for r in legacy] == [r["AtLeastHalfSinglet_DropletType"] for r in legacy]


def test_build_rejects_unknown_method():
    import pytest
    with pytest.raises(ValueError):
        c.build([], [], [], [], [], 0.7, method="MajoritySinglet")


# ------------------------------------------------------ genotype naming ----
def gpb(matrix):
    """{cluster: {donor: score}} -> rows of a GPb_diff table."""
    donors = sorted({d for row in matrix.values() for d in row})
    return [dict({"known genotype": d}, **{"snp_cluster_" + cl: str(row.get(d, "NA")) for cl, row in matrix.items()})
            for d in donors]


def test_read_scores_drops_na():
    s = c.read_scores(gpb({"0": {"d1": 0.1, "d2": 0.4}, "1": {"d1": "NA", "d2": 0.2}}))
    assert s == {"0": {"d1": 0.1, "d2": 0.4}, "1": {"d2": 0.2}}


def test_genotype_map_names_clear_matches_only():
    bcs = list("ABCD")
    clusters = dict(zip(bcs, ["0", "1", "2", "2"]))
    scores = {
        "0": {"d1": 0.10, "d2": 0.40, "d3": 0.40},  # clear -> d1
        "1": {"d1": 0.38, "d2": 0.40, "d3": 0.41},  # 5% margin -> unnamed
        "2": {"d1": 0.40, "d2": 0.40, "d3": 0.15},  # clear -> d3
    }
    m, report = c.genotype_map(bcs, clusters, scores, 0.1)
    assert m == {"0": "d1", "2": "d3"}
    assert {r[0]: r[4] for r in report} == {"0": "d1", "1": "unnamed", "2": "d3"}


def test_genotype_map_requires_a_mutual_best_match():
    # cluster 1 (e.g. a doublet cluster, or the smaller half of a split donor) is
    # nearest d1, but d1 is nearer cluster 0: only cluster 0 takes d1's name
    clusters = {"A": "0", "B": "1"}
    scores = {"0": {"d1": 0.1, "d2": 0.4}, "1": {"d1": 0.15, "d2": 0.4}}
    m, report = c.genotype_map(list("AB"), clusters, scores, 0.1)
    assert m == {"0": "d1"}
    assert dict((r[0], r[4]) for r in report)["1"].startswith("unnamed (d1 is closer to cluster 0")


def test_genotype_map_ungenotyped_donor_named_by_vireo_agreement():
    # cluster 2 is donor x, who has no genotype: it sits nearest d1 (not mutual)
    clusters = {"A": "0", "B": "1", "C": "2"}
    scores = {"0": {"d1": 0.1, "d2": 0.4}, "1": {"d1": 0.4, "d2": 0.1}, "2": {"d1": 0.2, "d2": 0.4}}
    m, report = c.genotype_map(list("ABC"), clusters, scores, 0.1, ungenotyped={"2": "x"})
    assert m == {"0": "d1", "1": "d2", "2": "x"}
    assert "no genotype" in dict((r[0], r[4]) for r in report)["2"]
    # without vireo's say-so, the leftover cluster is not named by elimination
    assert c.genotype_map(list("ABC"), clusters, scores, 0.1)[0] == {"0": "d1", "1": "d2"}


def test_build_genotype_naming_uses_vireo_for_ungenotyped_donors():
    bcs = list("ABCD")
    dem = popscle([(b, "AMB", "d1,d2,0.00") for b in bcs])
    fmx = popscle([])
    sop = soup([("A", "singlet", "0"), ("B", "singlet", "0"), ("C", "singlet", "1"), ("D", "singlet", "1")])
    vir = vireo([("A", "d1"), ("B", "d1"), ("C", "x"), ("D", "x")])
    scores = {"Freemuxlet": [], "Souporcell": gpb({"0": {"d1": 0.1, "d2": 0.4}, "1": {"d1": 0.2, "d2": 0.4}})}
    out = c.build(dem, fmx, sop, vir, [], 0.7, cluster_naming="genotype", scores=scores, ungenotyped=["x"])
    assert [r["Souporcell_Individual_Assignment"] for r in out] == ["d1", "d1", "x", "x"]
    assert [r["Consensus_Individual_Assignment"] for r in out] == ["d1", "d1", "x", "x"]


def test_genotype_map_needs_two_candidates():
    assert c.genotype_map(["A"], {"A": "0"}, {"0": {"d1": 0.1}}, 0.1)[0] == {}


def test_build_with_genotype_naming():
    dem = popscle([("A", "SNG", "d1,d1,0.00"), ("B", "SNG", "d2,d2,0.00")])
    fmx = popscle([("A", "SNG", "0,0"), ("B", "SNG", "1,1")])
    sop = soup([("A", "singlet", "5"), ("B", "singlet", "6")])
    vir = vireo([("A", "unassigned"), ("B", "unassigned")])  # agreement could name nothing
    snp = soup([])
    scores = {
        "Freemuxlet": gpb({"0": {"d1": 0.1, "d2": 0.4}, "1": {"d1": 0.4, "d2": 0.1}}),
        "Souporcell": gpb({"5": {"d1": 0.1, "d2": 0.4}, "6": {"d1": 0.3, "d2": 0.31}}),
    }
    out = c.build(dem, fmx, sop, vir, snp, 0.7, cluster_naming="genotype", scores=scores)
    a, b = out
    assert (a["Freemuxlet_Individual_Assignment"], a["Souporcell_Individual_Assignment"]) == ("d1", "d1")
    assert a["Consensus_Individual_Assignment"] == "d1"
    assert b["Souporcell_DropletType"] == "unassigned"  # cluster 6 unnamed: abstains
    assert (b["Consensus_DropletType"], b["Consensus_Individual_Assignment"]) == ("singlet", "d2")
    agree = c.build(dem, fmx, sop, vir, snp, 0.7, cluster_naming="vireo_agreement")
    assert agree[0]["Souporcell_DropletType"] == "unassigned"


def test_build_rejects_unknown_naming():
    import pytest
    with pytest.raises(ValueError):
        c.build([], [], [], [], [], 0.7, cluster_naming="correlation")
