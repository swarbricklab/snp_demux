import io

import pandas as pd
import pytest

from conftest import load_script

dp = load_script("doublet_posteriors")

SAMPLES = """sample_id,donor_id,capture_id
007,D7,cap1
S2,D2,cap1
S3,D3,cap1
X9,D7,cap2
"""

CELL_ASSIGNMENT = """barcode\tassignment\tstatus
AAA-1\t007\tsinglet
CCC-1\tNA\tdoublet
GGG-1\tNA\tunassigned
TTT-1\tS2\tsinglet
NNN-1\tNA\tdoublet
"""

VIREO = """cell\tdonor_id\tprob_max\tprob_doublet\tn_vars\tbest_singlet\tbest_doublet\tdoublet_logLikRatio
AAA-1\tD7\t0.99\t1.2e-05\t40\tD7\tD7,D2\t-3.1
CCC-1\tdoublet\t0.51\t0.97\t35\tD2\tD2,D3\t4.0
GGG-1\tunassigned\t0.4\t0.3\t2\tD3\tD7,D3\t0.1
TTT-1\tD2\t1\t0\t50\tD2\tD2,D3\t-5
"""


def _tables():
    return (
        dp.read_str_table(io.StringIO(SAMPLES), ","),
        dp.read_str_table(io.StringIO(CELL_ASSIGNMENT), "\t"),
        dp.read_str_table(io.StringIO(VIREO), "\t"),
    )


def _build():
    samples, ca, vireo = _tables()
    return dp.build_posteriors(
        ca, vireo, dp.donor_to_sample_map(samples, "cap1"), "cap1"
    )


def test_columns_and_one_row_per_consensus_barcode():
    out = _build()
    assert list(out.columns) == dp.COLUMNS
    assert out["barcode"].tolist() == ["AAA-1", "CCC-1", "GGG-1", "TTT-1", "NNN-1"]
    assert (out["capture_id"] == "cap1").all()


def test_status_counts_match_consensus():
    _, ca, _ = _tables()
    out = _build()
    assert out["consensus_status"].value_counts().to_dict() == (
        ca["status"].value_counts().to_dict()
    )


def test_bare_digit_sample_ids_survive_and_best_singlet_in_sample_space():
    out = _build().set_index("barcode")
    assert out.loc["AAA-1", "consensus_assignment"] == "007"
    assert out.loc["AAA-1", "vireo_best_singlet"] == "007"  # D7 -> 007, not cap2's X9
    assert out.loc["CCC-1", "vireo_best_singlet"] == "S2"


def test_vireo_values_pass_through_verbatim():
    out = _build().set_index("barcode")
    assert out.loc["AAA-1", "vireo_prob_doublet"] == "1.2e-05"
    assert out.loc["TTT-1", "vireo_prob_max"] == "1"


def test_barcode_missing_from_vireo_gets_na():
    out = _build().set_index("barcode")
    assert pd.isna(out.loc["NNN-1", "vireo_prob_doublet"])
    assert pd.isna(out.loc["NNN-1", "consensus_assignment"])


def test_csv_round_trip_writes_na(tmp_path):
    path = tmp_path / "x.csv"
    dp.write_csv(_build(), path)
    text = path.read_text().splitlines()
    assert text[0] == ",".join(dp.COLUMNS)
    assert "cap1,NNN-1,doublet,NA,NA,NA,NA" in text
    assert "cap1,AAA-1,singlet,007,1.2e-05,0.99,007" in text


def test_non_numeric_posterior_fails_loudly():
    samples, ca, vireo = _tables()
    vireo.loc[0, "prob_doublet"] = "oops"
    with pytest.raises(ValueError):
        dp.build_posteriors(ca, vireo, {}, "cap1")
