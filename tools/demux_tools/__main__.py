"""Entry point for `python -I <checkout>/tools/demux_tools <command> [options]`.

Python writes no bytecode from here on, so running the tools from a checkout
that belongs to someone else leaves no __pycache__ in it.
"""
import os
import sys

sys.dont_write_bytecode = True

_PACKAGE_DIR = os.path.dirname(os.path.abspath(__file__))

# Running a directory compiles this file through the import system, which caches
# its bytecode in __pycache__ before the line above runs. Remove that one file,
# and the folder if it is then empty. `python -I -B` avoids even this write.
_cached = globals().get("__cached__")
if _cached and os.path.dirname(os.path.dirname(os.path.abspath(_cached))) == _PACKAGE_DIR:
    try:
        os.remove(_cached)
        os.rmdir(os.path.dirname(_cached))
    except OSError:
        pass

# Running a directory puts the package directory itself first on sys.path, and
# -I leaves out the folder that holds the package. Swap one for the other, so
# `import demux_tools` works and no module here shadows a top-level name.
sys.path[:] = [p for p in sys.path if os.path.abspath(p or os.curdir) != _PACKAGE_DIR]
_TOOLS_DIR = os.path.dirname(_PACKAGE_DIR)
if _TOOLS_DIR not in sys.path:
    sys.path.insert(0, _TOOLS_DIR)

from demux_tools.cli import main  # noqa: E402

if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
