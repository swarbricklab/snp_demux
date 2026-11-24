import math

import pandas as pd
import pytest

from conftest import load_script

ds = load_script("doublet_summary")


def _posteriors(assignments, statuses, prob_doublet):
    return pd.DataFrame(
        {
            "capture_id": "cap1",
            "barcode": [f"BC{i}" for i in range(len(statuses))],
            "consensus_status": statuses,
            "consensus_assignment": assignments,
            "vireo_prob_doublet": [None if p is None else str(p) for p in prob_doublet],
            "vireo_prob_max": "0.9",
            "vireo_best_singlet": "S1",
        }
    )


def test_hand_computed_values():
    # 6 singlets (S1 x3, S2 x2, 007 x1), 2 doublets, 2 unassigned; passenger P0.
    a = ["S1", "S1", "S1", "S2", "S2", "007", None, None, None, None]
    s = ["singlet"] * 6 + ["doublet", "doublet", "unassigned", "unassigned"]
    p = [0.0, 0.0, 0.1, 0.0, 0.0, 0.0, 0.9, 1.0, 0.5, None]
    row = ds.summarise(_posteriors(a, s, p), ["S1", "S2", "007", "P0"], "cap1")

    assert row["n_barcodes"] == 10
    assert row["n_singlet"] == 6
    assert row["n_doublet_hard"] == 2
    assert row["n_unassigned"] == 2
    assert row["sum_prob_doublet"] == pytest.approx(2.5)  # NA skipped
    assert row["doublet_rate_cross_donor"] == pytest.approx(0.2)
    assert row["doublet_rate_cross_donor_posterior"] == pytest.approx(0.25)
    assert row["n_samples"] == 4  # passenger counted
    # p = 1/2, 1/3, 1/6, 0  ->  1/4 + 1/9 + 1/36 = 14/36
    assert row["same_donor_pair_prob"] == pytest.approx(14 / 36)
    assert row["sample_props"] == "S1:0.5;S2:0.333333;007:0.166667;P0:0.0"


def test_single_sample_capture_has_sum_p_sq_one():
    a = ["S1"] * 5 + [None]
    s = ["singlet"] * 5 + ["unassigned"]
    row = ds.summarise(_posteriors(a, s, [0] * 6), ["S1"], "cap1")
    assert row["same_donor_pair_prob"] == 1.0
    assert row["n_samples"] == 1
    assert row["doublet_rate_cross_donor"] == 0.0


def test_no_singlets_gives_na_not_crash():
    row = ds.summarise(
        _posteriors([None, None], ["doublet", "unassigned"], [1, 0]), ["S1", "S2"], "c"
    )
    assert math.isnan(row["same_donor_pair_prob"])
    assert math.isnan(row["sample_props"])


def test_assignments_outside_declared_samples_ignored_for_props():
    a = ["S1", "S2", "ZZ"]
    row = ds.summarise(_posteriors(a, ["singlet"] * 3, [0] * 3), ["S1", "S2"], "c")
    assert row["same_donor_pair_prob"] == pytest.approx(0.5)


def test_csv_columns_and_combine(tmp_path):
    a = ["S1", "S2"]
    rows = [
        ds.to_frame(ds.summarise(_posteriors(a, ["singlet"] * 2, [0, 0]), a, cid))
        for cid in ("007", "cap2")
    ]
    paths = []
    for i, df in enumerate(rows):
        paths.append(tmp_path / f"{i}.csv")
        ds.write_csv(df, paths[-1])
    assert paths[0].read_text().splitlines()[0] == ",".join(ds.COLUMNS)

    class SM:
        class input:
            summaries = [str(p) for p in paths]

        class output:
            csv = str(tmp_path / "all.csv")

    ds.main_combine(SM)
    combined = ds.read_str_table(SM.output.csv)
    assert combined["capture_id"].tolist() == ["007", "cap2"]  # bare digits kept
    assert list(combined.columns) == ds.COLUMNS


def test_read_calls_joins_consensus_and_vireo_posterior():
    cell = pd.DataFrame(
        {
            "barcode": ["A-1", "B-1", "C-1"],
            "assignment": ["007", None, None],
            "status": ["singlet", "doublet", "unassigned"],
        }
    )
    sba = pd.DataFrame(
        {
            "Barcode": ["B-1", "A-1", "A-1"],
            "Vireo_prob_doublet": ["0.97", "1.2e-05", "1.2e-05"],
            "Vireo_prob_max": ["0.02", "0.99", "0.99"],
        }
    )
    calls = ds.read_calls(cell, sba)
    assert list(calls["barcode"]) == ["A-1", "B-1", "C-1"]
    assert list(calls["consensus_status"]) == ["singlet", "doublet", "unassigned"]
    assert calls["consensus_assignment"].tolist()[0] == "007"
    assert calls["vireo_prob_doublet"].tolist()[:2] == ["1.2e-05", "0.97"]
    assert pd.isna(calls["vireo_prob_doublet"].tolist()[2])
