"""Shared setup for the demux_tools tests (`pytest tools/tests`).

- Puts tools/ on sys.path, so tests `import demux_tools...`.
- Registers two markers and skips their tests where the library is missing:
      @pytest.mark.scipy   needs scipy
      @pytest.mark.vireo   needs an importable vireoSNP (tests compare with it)
- Fixtures:
      tools_dir             Path of tools/
      run_tools(*args, cwd=None, python=None, blocked=(), package=None)
                            runs `python -I <tools>/demux_tools *args` in a
                            subprocess (cwd defaults to a fresh temporary
                            folder; `blocked` names top-level modules made
                            unimportable, e.g. ("scipy",); `package` runs
                            another copy of demux_tools); returns the
                            CompletedProcess with text stdout and stderr
      synthetic_capture(name="capture", tags=("GT",), gp_noise=0.0,
                        vcf_rows=None, **make_capture_kwargs)
                            writes a synthetic capture under tmp_path/<name>
                            with demux_tools.synth.build_on_disk and returns
                            its dict: dir, cellsnp, vcf, truth, donors, capture
Test files are named test_tools_<module>.py, so their basenames stay unique
next to other test folders in this repository.
"""
import importlib
import os
import subprocess
import sys
from pathlib import Path

import pytest

TOOLS_DIR = Path(__file__).resolve().parents[1]
if str(TOOLS_DIR) not in sys.path:
    sys.path.insert(0, str(TOOLS_DIR))
sys.dont_write_bytecode = True

# Run before the tool starts: makes the named top-level modules unimportable.
_BLOCKER = """
import sys
class _Block(object):
    def __init__(self, names):
        self.names = set(names)
    def find_spec(self, name, path=None, target=None):
        if name.split(".")[0] in self.names:
            raise ImportError("blocked for this test: " + name)
        return None
    def find_module(self, name, path=None):
        if name.split(".")[0] in self.names:
            raise ImportError("blocked for this test: " + name)
        return None
sys.meta_path.insert(0, _Block({names!r}))
sys.argv = [{main!r}] + {args!r}
sys.dont_write_bytecode = True
import runpy
runpy.run_path({main!r}, run_name="__main__")
"""


def _importable(name):
    try:
        importlib.import_module(name)
        return True
    except Exception:  # a broken install counts as missing
        return False


def pytest_configure(config):
    config.addinivalue_line("markers", "scipy: needs scipy; skipped where it is not installed")
    config.addinivalue_line("markers", "vireo: compares with vireoSNP; skipped where it does not import")


def pytest_collection_modifyitems(config, items):
    needs = {"scipy": "scipy", "vireo": "vireoSNP"}
    available = {}
    for item in items:
        for marker, module in needs.items():
            if item.get_closest_marker(marker) is None:
                continue
            if module not in available:
                available[module] = _importable(module)
            if not available[module]:
                item.add_marker(pytest.mark.skip(reason="%s is not importable here" % module))


@pytest.fixture
def tools_dir():
    return TOOLS_DIR


@pytest.fixture
def run_tools(tmp_path):
    def run(*args, **kw):
        cwd = kw.pop("cwd", None)
        python = kw.pop("python", None) or sys.executable
        blocked = tuple(kw.pop("blocked", ()))
        package = kw.pop("package", None) or str(TOOLS_DIR / "demux_tools")
        if kw:
            raise TypeError("unknown option(s) %s" % sorted(kw))
        if cwd is None:
            cwd = tmp_path / "cwd"
            cwd.mkdir(exist_ok=True)
        args = [str(a) for a in args]
        if blocked:
            code = _BLOCKER.format(names=list(blocked), main=package, args=args)
            cmd = [python, "-I", "-c", code]
        else:
            cmd = [python, "-I", package] + args
        return subprocess.run(cmd, cwd=str(cwd), stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                              universal_newlines=True, env=dict(os.environ))
    return run


@pytest.fixture
def synthetic_capture(tmp_path):
    from demux_tools.synth import build_on_disk

    def build(name="capture", tags=("GT",), gp_noise=0.0, vcf_rows=None, **make_kwargs):
        return build_on_disk(str(tmp_path / name), tags=tags, gp_noise=gp_noise, vcf_rows=vcf_rows, **make_kwargs)
    return build
