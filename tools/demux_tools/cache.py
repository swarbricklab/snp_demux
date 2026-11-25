"""A folder of numpy archives that holds loaded inputs between runs.

An entry is one `.npz` file named by a key. The key is a sha1 over the
cache format version, the kind of entry, each input file's real path, size
and modification time, and the options that shaped the arrays, so replacing
or touching an input, or changing an option, gives a new key. Nothing else is
checked: an input rewritten with the same size and modification time would
be served stale, and some filesystems (Lustre among them) keep modification
times to the second. FORMAT_VERSION changes whenever what an entry holds, or
how a load builds it, changes, so older entries are not read.

Entries are written to a temporary name and renamed into place, so a run
that stops part way leaves no entry, and two runs writing the same key leave
one complete entry. Arrays are stored without pickling (strings as numpy
unicode arrays), and the metadata as JSON text inside the archive. An entry
that cannot be read counts as absent, and an entry that cannot be written
gives a warning: the cache never stops a run.
"""
import hashlib
import json
import os
import sys
import zipfile

import numpy as np

from demux_tools._io import atomic_path, file_record

# 2: the Scorer matches repeated donor VCF ids on every SNP record, as vireo does
FORMAT_VERSION = 2
_META = "__meta__"


def key(kind, inputs, options):
    """Cache key (sha1 hex) of `kind` for the input files `inputs` (paths, order kept) and `options` (JSON-able)."""
    records = []
    for p in inputs:
        r = file_record(p)
        records.append([r["realpath"], r["size"], r["mtime_ns"]])
    text = json.dumps({"format": FORMAT_VERSION, "kind": kind, "inputs": records, "options": options},
                      sort_keys=True, default=str)
    return hashlib.sha1(text.encode("utf-8")).hexdigest()


def entry_path(cache_dir, cache_key):
    """Path of the entry for `cache_key` in `cache_dir`."""
    return os.path.join(os.fspath(cache_dir), cache_key + ".npz")


def load(cache_dir, cache_key):
    """(arrays, meta) of an entry, or None when there is none or it cannot be read.

    arrays: {name: numpy array}; meta: the dict saved with it.
    """
    path = entry_path(cache_dir, cache_key)
    if not os.path.isfile(path):
        return None
    try:
        with np.load(path, allow_pickle=False) as z:
            arrays = {name: z[name] for name in z.files if name != _META}
            meta = json.loads(str(z[_META])) if _META in z.files else {}
    except (OSError, ValueError, KeyError, EOFError, zipfile.BadZipFile):
        return None
    return arrays, meta


def save(cache_dir, cache_key, arrays, meta=None):
    """Write an entry atomically and return its path; None, with a warning on stderr, when it cannot be written.

    Arrays must not need pickling.
    """
    path = entry_path(cache_dir, cache_key)
    payload = dict(arrays)
    payload[_META] = np.array(json.dumps(meta or {}, sort_keys=True, default=str))
    try:
        os.makedirs(os.fspath(cache_dir), exist_ok=True)
        with atomic_path(path) as tmp:
            with open(tmp, "wb") as f:
                np.savez(f, **payload)
    except OSError as e:
        sys.stderr.write("demux_tools: warning: cache entry not written (%s); the run carries on\n" % e)
        return None
    return path
