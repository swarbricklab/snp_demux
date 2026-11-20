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
    assert row["f_x_hard"] == pytest.approx(0.2)
    assert row["f_x_posterior"] == pytest.approx(0.25)
    assert row["n_samples"] == 4  # passenger counted
    # p = 1/2, 1/3, 1/6, 0  ->  1/4 + 1/9 + 1/36 = 14/36
    assert row["sum_p_sq"] == pytest.approx(14 / 36)
    assert row["sample_props"] == "S1:0.5;S2:0.333333;007:0.166667;P0:0.0"


def test_single_sample_capture_has_sum_p_sq_one():
    a = ["S1"] * 5 + [None]
    s = ["singlet"] * 5 + ["unassigned"]
    row = ds.summarise(_posteriors(a, s, [0] * 6), ["S1"], "cap1")
    assert row["sum_p_sq"] == 1.0
    assert row["n_samples"] == 1
    assert row["f_x_hard"] == 0.0


def test_no_singlets_gives_na_not_crash():
    row = ds.summarise(
        _posteriors([None, None], ["doublet", "unassigned"], [1, 0]), ["S1", "S2"], "c"
    )
    assert math.isnan(row["sum_p_sq"])
    assert math.isnan(row["sample_props"])


def test_assignments_outside_declared_samples_ignored_for_props():
    a = ["S1", "S2", "ZZ"]
    row = ds.summarise(_posteriors(a, ["singlet"] * 3, [0] * 3), ["S1", "S2"], "c")
    assert row["sum_p_sq"] == pytest.approx(0.5)


def test_csv_columns_and_combine(tmp_path):
    a = ["S1", "S2"]
    rows = []
    for cid in ("007", "cap2"):
        row = ds.summarise(_posteriors(a, ["singlet"] * 2, [0, 0]), a, cid)
        row.update(ds.expected_rates(row, 0.008))
        rows.append(ds.to_frame(row))
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


def _equal_four_donor_capture(f_x=0.06, n=10000):
    """4 samples in equal proportion; n barcodes, a fraction f_x hard doublets."""
    n_doublet = round(n * f_x)
    per_sample = (n - n_doublet) // 4
    a = [s for s in ("A", "B", "C", "D") for _ in range(per_sample)]
    statuses = ["singlet"] * len(a) + ["doublet"] * n_doublet
    a += [None] * n_doublet
    return _posteriors(a, statuses, [0] * len(statuses)), ["A", "B", "C", "D"]


def test_expected_rates_equal_four_donor():
    post, samples = _equal_four_donor_capture()
    row = ds.summarise(post, samples, "cap1")
    assert row["sum_p_sq"] == pytest.approx(0.25)
    r = ds.expected_rates(row, 0.008)
    assert r["d_total"] == pytest.approx(0.08)
    assert r["d_within"] == pytest.approx(0.02)
    assert r["d_loading"] == pytest.approx(0.08)  # 10k cells * 0.008 / 1k
    assert r["d_ratio"] == pytest.approx(1.0)
    # binomial SE of f_x, scaled by 1 / (1 - sum_p_sq)
    assert r["d_total_sd"] == pytest.approx(math.sqrt(0.06 * 0.94 / 10000) / 0.75)
    assert r["d_fallback"] == ""


def test_ht_rate_halves_loading_and_doubles_ratio():
    post, samples = _equal_four_donor_capture()
    row = ds.summarise(post, samples, "cap1")
    r = ds.expected_rates(row, 0.004)
    assert r["d_loading"] == pytest.approx(0.04)
    assert r["d_ratio"] == pytest.approx(2.0)


@pytest.mark.parametrize(
    "assignments,statuses,samples,reason",
    [
        (["S1"] * 3, ["singlet"] * 3, ["S1"], "single_sample"),
        ([None, None], ["doublet", "unassigned"], ["S1", "S2"], "no_singlets"),
        (["S1", "S1", None], ["singlet", "singlet", "doublet"], ["S1", "S2"], "one_sample_only"),
    ],
)
def test_fallbacks_give_na_rates(assignments, statuses, samples, reason):
    post = _posteriors(assignments, statuses, [0] * len(statuses))
    row = ds.summarise(post, samples, "cap1")
    r = ds.expected_rates(row, 0.008)
    assert r["d_fallback"] == reason
    for k in ("d_total", "d_total_sd", "d_within", "d_ratio"):
        assert math.isnan(r[k])
    assert r["d_loading"] == pytest.approx(len(statuses) / 1000 * 0.008)


def test_written_row_has_all_columns(tmp_path):
    post, samples = _equal_four_donor_capture()
    row = ds.summarise(post, samples, "cap1")
    row.update(ds.expected_rates(row, 0.008))
    ds.write_csv(ds.to_frame(row), tmp_path / "s.csv")
    back = ds.read_str_table(tmp_path / "s.csv")
    assert list(back.columns) == ds.COLUMNS
    assert back.loc[0, "d_total"] == "0.08"
    assert pd.isna(back.loc[0, "d_fallback"])  # empty -> NA on read


def test_recovered_barcodes_counted_apart_and_left_out_of_p_i():
    a = ["S1", "S2", "S2", None]
    s = ["singlet", "singlet", "recovered", "unassigned"]
    row = ds.summarise(_posteriors(a, s, [0.0, 0.0, 0.0, 0.0]), ["S1", "S2"], "cap1")
    assert (row["n_singlet"], row["n_recovered"], row["n_unassigned"]) == (2, 1, 1)
    assert row["sample_props"] == "S1:0.5;S2:0.5"  # the recovered S2 is not a singlet
