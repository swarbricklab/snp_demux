"""Shared file helpers: text and tables, deterministic gzip, provenance sidecars.

Contracts every command relies on:

- Tables are read with every column as text (`dtype=str`, no NA parsing), so a
  donor id such as "0012" or "NA" survives a round trip unchanged.
- A path ending in ".gz" is written as gzip with mtime 0 and no file name in the
  header, so the same text gives the same bytes in one environment.
- Every write goes to a temporary file beside the target, renamed into place
  when complete; a failed write leaves neither a partial file nor the temporary.
- Gzip input is recognised by its magic bytes, not its name (BGZF is gzip).
- Outputs go only where the caller names them; `check_output_path` refuses an
  output path that is one of the inputs.

numpy and pandas are imported where they are used, so this module imports on
the standard library alone.
"""
import contextlib
import gzip
import hashlib
import io
import json
import os
import platform
import sys
import time
from datetime import datetime, timezone

GZIP_MAGIC = b"\x1f\x8b"


class UsageError(Exception):
    """A user-facing error: the dispatcher prints the message and exits with status 2."""


# ---------------------------------------------------------------------------
# Opening files


def is_gzip(path):
    """True when the file at `path` starts with the gzip magic bytes (gzip or BGZF)."""
    with open(path, "rb") as f:
        return f.read(2) == GZIP_MAGIC


@contextlib.contextmanager
def atomic_path(path):
    """Yield a temporary path beside `path`; rename it to `path` if the block succeeds.

    The temporary is removed if the block raises. Use it for writers that take a
    path (e.g. numpy or a subprocess); `open_text` and `open_binary` use it.
    """
    path = os.fspath(path)
    directory = os.path.dirname(os.path.abspath(path))
    tmp = os.path.join(directory, ".%s.tmp%d" % (os.path.basename(path), os.getpid()))
    try:
        yield tmp
        os.replace(tmp, path)
    finally:
        if os.path.lexists(tmp):
            os.remove(tmp)


@contextlib.contextmanager
def open_binary(path):
    """Open `path` for writing bytes, atomically; gzip (mtime 0, no name) if it ends in ".gz"."""
    path = os.fspath(path)
    with atomic_path(path) as tmp:
        with open(tmp, "wb") as raw:
            if path.endswith(".gz"):
                with gzip.GzipFile(filename="", mode="wb", fileobj=raw, mtime=0) as gz:
                    yield gz
            else:
                yield raw


def open_text(path, mode="rt", encoding="utf-8"):
    """Open a text file.

    mode "rt": returns a file object; gzip is detected by magic bytes.
    mode "wt": returns a context manager (use it in a `with` block) that writes
    atomically, as gzip with mtime 0 when `path` ends in ".gz", with "\\n" line
    ends whatever the platform.
    """
    path = os.fspath(path)
    if mode == "rt":
        if is_gzip(path):
            return gzip.open(path, "rt", encoding=encoding)
        return open(path, "rt", encoding=encoding)
    if mode == "wt":
        return _open_text_write(path, encoding)
    raise ValueError("open_text mode must be 'rt' or 'wt', not %r" % (mode,))


@contextlib.contextmanager
def _open_text_write(path, encoding):
    with open_binary(path) as raw:
        text = io.TextIOWrapper(raw, encoding=encoding, newline="\n")
        try:
            yield text
            text.flush()
        finally:
            text.detach()


def read_lines(path):
    """All lines of a text file (gzip by magic), without line ends."""
    with open_text(path) as f:
        return [line.rstrip("\r\n") for line in f]


# ---------------------------------------------------------------------------
# Tables


def read_table(path, sep="\t", usecols=None, header="infer", names=None, comment=None):
    """Read a delimited table with every column as text.

    Empty fields stay "" and words such as "NA" stay as written (no NA parsing).
    Gzip is detected by magic bytes. Arguments are pandas.read_csv's.
    Returns a pandas DataFrame with a default RangeIndex.
    """
    import pandas as pd

    path = os.fspath(path)
    return pd.read_csv(
        path,
        sep=sep,
        usecols=usecols,
        header=header,
        names=names,
        comment=comment,
        dtype=str,
        keep_default_na=False,
        na_filter=False,
        compression="gzip" if is_gzip(path) else None,
        engine="c",
    )


def format_column(values, fmt):
    """Format each value with a %-format (e.g. "%.2e"); missing values (None, NaN, NA) become ""."""
    import pandas as pd

    return ["" if pd.isna(v) else fmt % v for v in values]


def write_table(df, path, float_formats=None, columns=None):
    """Write `df` as a tab-separated table with a header and no index.

    float_formats: {column: "%-format"} applied to those columns first (NaN -> "").
    Other values are written as pandas writes them; missing values as "".
    A path ending in ".gz" gives deterministic gzip (mtime 0). The write is
    atomic. Returns `path`.
    """
    path = os.fspath(path)
    if columns is not None:
        df = df.loc[:, list(columns)]
    if float_formats:
        df = df.copy()
        for col, fmt in float_formats.items():
            df[col] = format_column(df[col].tolist(), fmt)
    with open_text(path, "wt") as f:
        df.to_csv(f, sep="\t", index=False, na_rep="")
    return path


def write_text(path, text):
    """Write a string to `path` atomically (gzip with mtime 0 for ".gz"). Returns `path`."""
    with open_text(path, "wt") as f:
        f.write(text)
    return os.fspath(path)


# ---------------------------------------------------------------------------
# Inputs, outputs and provenance


def file_record(path):
    """{path, realpath, size, mtime_ns} of an existing file (symbolic links followed)."""
    path = os.fspath(path)
    st = os.stat(path)
    return {"path": path, "realpath": os.path.realpath(path), "size": st.st_size, "mtime_ns": st.st_mtime_ns}


def fingerprint(paths):
    """sha1 hex digest over (realpath, size, mtime_ns) of each path, in the order given.

    Cheap (no file content is read); changes when a file is replaced or touched.
    A directory contributes its own record only, not its contents.
    """
    h = hashlib.sha1()
    for p in paths:
        r = file_record(p)
        h.update(("%s\t%d\t%d\n" % (r["realpath"], r["size"], r["mtime_ns"])).encode("utf-8"))
    return h.hexdigest()


def check_output_path(out, inputs):
    """Raise UsageError if `out` resolves to one of `inputs` (never overwrite an input)."""
    target = os.path.realpath(os.fspath(out))
    for p in inputs:
        if p is not None and os.path.realpath(os.fspath(p)) == target:
            raise UsageError("output %s is also an input" % (out,))


def _read_first_line(path):
    with open(path, "rt") as f:
        return f.readline().strip()


def _is_sha(text):
    return len(text) in (40, 64) and all(c in "0123456789abcdef" for c in text)


def git_sha(path=None):
    """Commit checked out in the git repository that holds `path` (default: this package).

    Reads .git files directly (git need not be installed, nothing is written).
    Handles worktrees, loose and packed refs, and a detached HEAD. Returns None
    outside a repository or when the ref cannot be resolved. It does not say
    whether the working tree has uncommitted changes.
    """
    d = os.path.dirname(os.path.abspath(path or __file__))
    gitdir = None
    while True:
        dotgit = os.path.join(d, ".git")
        if os.path.isdir(dotgit):
            gitdir = dotgit
            break
        if os.path.isfile(dotgit):
            line = _read_first_line(dotgit)
            if not line.startswith("gitdir:"):
                return None
            gitdir = os.path.normpath(os.path.join(d, line[len("gitdir:"):].strip()))
            break
        parent = os.path.dirname(d)
        if parent == d:
            return None
        d = parent
    try:
        head = _read_first_line(os.path.join(gitdir, "HEAD"))
        if not head.startswith("ref:"):
            return head if _is_sha(head) else None
        ref = head[len("ref:"):].strip()
        commondir = gitdir
        if os.path.isfile(os.path.join(gitdir, "commondir")):
            commondir = os.path.normpath(os.path.join(gitdir, _read_first_line(os.path.join(gitdir, "commondir"))))
        for base in (gitdir, commondir):
            loose = os.path.join(base, ref)
            if os.path.isfile(loose):
                sha = _read_first_line(loose)
                return sha if _is_sha(sha) else None
        packed = os.path.join(commondir, "packed-refs")
        if os.path.isfile(packed):
            with open(packed, "rt") as f:
                for line in f:
                    parts = line.split()
                    if len(parts) == 2 and parts[1] == ref and _is_sha(parts[0]):
                        return parts[0]
    except OSError:
        return None
    return None


def _distribution_version(dist):
    try:
        from importlib import metadata  # Python 3.8+
    except ImportError:
        return None
    try:
        return metadata.version(dist)
    except Exception:
        return None


def library_versions():
    """Versions of Python and the libraries the tools use.

    numpy and pandas are imported. scipy and vireoSNP are never imported here:
    their version is reported when already imported, otherwise from installed
    package metadata (Python 3.8+); "installed" when the package is found but
    its version cannot be read without importing it (Python 3.7); None when absent.
    """
    import importlib.util

    import numpy
    import pandas

    out = {
        "python": platform.python_version(),
        "numpy": numpy.__version__,
        "pandas": pandas.__version__,
    }
    for key, module, dist in (("scipy", "scipy", "scipy"), ("vireoSNP", "vireoSNP", "vireoSNP")):
        mod = sys.modules.get(module)
        version = getattr(mod, "__version__", None) if mod is not None else None
        version = version or _distribution_version(dist)
        if version is None:
            try:
                version = "installed" if importlib.util.find_spec(module) is not None else None
            except (ImportError, ValueError):
                version = None
        out[key] = version
    return out


def start_clock():
    """Start times for `provenance`: {"started", "wall", "cpu"}."""
    return {
        "started": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "wall": time.time(),
        "cpu": time.process_time(),
    }


def peak_rss_mb():
    """Peak resident memory of this process so far, in MB (None where unavailable)."""
    try:
        import resource
    except ImportError:
        return None
    rss = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    return round(rss / (1024.0 * 1024.0) if sys.platform == "darwin" else rss / 1024.0, 1)


def provenance(argv, inputs, clock=None, extra=None):
    """The provenance record written beside every output.

    argv: the command line (list of str). inputs: paths read (files must exist).
    clock: from `start_clock()`, for timings. extra: command-specific values,
    stored under "results". Keys: tool, version, git_sha, tools_path, argv, cwd, versions,
    inputs (file_record each), inputs_fingerprint, started, seconds_wall,
    seconds_cpu, peak_rss_mb, results.
    """
    from demux_tools import __version__

    inputs = [os.fspath(p) for p in inputs if p is not None]
    rec = {
        "tool": "demux_tools",
        "version": __version__,
        "git_sha": git_sha(),
        "tools_path": os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
        "argv": list(argv),
        "cwd": os.getcwd(),
        "versions": library_versions(),
        "inputs": [file_record(p) for p in inputs],
        "inputs_fingerprint": fingerprint(inputs),
    }
    if clock is not None:
        rec["started"] = clock["started"]
        rec["seconds_wall"] = round(time.time() - clock["wall"], 3)
        rec["seconds_cpu"] = round(time.process_time() - clock["cpu"], 3)
    rec["peak_rss_mb"] = peak_rss_mb()
    if extra is not None:
        rec["results"] = extra
    return rec


def sidecar_path(out):
    """Where the provenance of `out` goes: <dir>/provenance.json for a directory, else <out>.provenance.json."""
    out = os.fspath(out)
    if os.path.isdir(out):
        return os.path.join(out, "provenance.json")
    return out + ".provenance.json"


def _json_default(obj):
    import numpy as np

    if isinstance(obj, np.integer):
        return int(obj)
    if isinstance(obj, np.floating):
        return float(obj)
    if isinstance(obj, np.bool_):
        return bool(obj)
    if isinstance(obj, np.ndarray):
        return obj.tolist()
    raise TypeError("not JSON serialisable: %r" % (type(obj),))


def write_json(obj, path):
    """Write `obj` as indented JSON with sorted keys, atomically. numpy scalars and arrays allowed."""
    return write_text(path, json.dumps(obj, indent=1, sort_keys=True, default=_json_default) + "\n")


def write_provenance(out, argv, inputs, clock=None, extra=None):
    """Write the provenance record of output `out` to `sidecar_path(out)`; return that path."""
    path = sidecar_path(out)
    write_json(provenance(argv, inputs, clock=clock, extra=extra), path)
    return path
