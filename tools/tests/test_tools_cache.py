"""demux_tools.cache: keys follow the inputs and options; entries round-trip without pickling."""
import os

import numpy as np

from demux_tools import cache


def test_key_changes_with_inputs_and_options(tmp_path):
    a = tmp_path / "a.txt"
    a.write_text("x")
    k1 = cache.key("capture", [str(a)], {"tag": "GT"})
    assert k1 == cache.key("capture", [str(a)], {"tag": "GT"})
    assert k1 != cache.key("capture", [str(a)], {"tag": "GP"})
    assert k1 != cache.key("other", [str(a)], {"tag": "GT"})
    st = os.stat(str(a))
    os.utime(str(a), ns=(st.st_atime_ns, st.st_mtime_ns + 10 ** 9))
    assert k1 != cache.key("capture", [str(a)], {"tag": "GT"})


def test_save_and_load_round_trip(tmp_path):
    arrays = {"ints": np.arange(5, dtype=np.int32), "names": np.array(["a", "bb"]), "mask": np.array([True, False])}
    path = cache.save(str(tmp_path / "c"), "k" * 40, arrays, {"shape": [2, 3]})
    assert os.path.basename(path) == "k" * 40 + ".npz"
    assert sorted(os.listdir(str(tmp_path / "c"))) == ["k" * 40 + ".npz"]      # no temporary left
    got, meta = cache.load(str(tmp_path / "c"), "k" * 40)
    assert meta == {"shape": [2, 3]}
    for name, value in arrays.items():
        assert np.array_equal(got[name], value) and got[name].dtype == value.dtype
    assert cache.load(str(tmp_path / "c"), "0" * 40) is None
    (tmp_path / "c" / ("1" * 40 + ".npz")).write_bytes(b"not an archive")
    assert cache.load(str(tmp_path / "c"), "1" * 40) is None
