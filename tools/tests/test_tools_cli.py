"""The dispatcher: the command table, --help without scipy, running by path without bytecode."""
import ast
import os
import shutil

import pytest

from demux_tools import cli

# Modules whose commands may need scipy to import. Every other command's module
# must import, and its --help run, with scipy unavailable.
SCIPY_MODULES = set()


def _table_lines():
    """The COMMANDS entries as written in cli.py, one per line."""
    source = open(cli.__file__).read()
    tree = ast.parse(source)
    for node in tree.body:
        if isinstance(node, ast.Assign) and getattr(node.targets[0], "id", None) == "COMMANDS":
            lines = source.splitlines()[node.lineno:node.end_lineno - 1] if hasattr(node, "end_lineno") else None
            return node, lines
    raise AssertionError("no COMMANDS table")


def test_commands_table_is_one_line_per_command_in_order():
    node, lines = _table_lines()
    names = [k.value if hasattr(k, "value") else k.s for k in node.value.keys]
    assert names == sorted(names)
    assert list(cli.COMMANDS) == names
    for name, entry in cli.COMMANDS.items():
        assert isinstance(entry, tuple) and len(entry) == 2 and all(isinstance(x, str) for x in entry), name
    if lines is not None:   # Python 3.8+ records where the table ends
        entries = [l for l in lines if l.strip() and not l.strip().startswith("#")]
        assert len(entries) == len(names)
        for line, name in zip(entries, names):
            assert line.strip().startswith('"%s": (' % name) and line.rstrip().endswith("),")


@pytest.mark.parametrize("name", sorted(cli.COMMANDS))
def test_each_command_has_its_functions(name):
    if cli.COMMANDS[name][0] in SCIPY_MODULES:
        pytest.importorskip("scipy")
    configure, run = cli.command_functions(name)
    assert callable(configure) and callable(run)


@pytest.mark.parametrize("name", sorted(cli.COMMANDS))
def test_help_runs_without_scipy(name, run_tools):
    if cli.COMMANDS[name][0] in SCIPY_MODULES:
        pytest.skip("a Scorer command may need scipy")
    res = run_tools(name, "--help", blocked=("scipy", "vireoSNP"))
    assert res.returncode == 0, res.stderr
    assert res.stdout.startswith("usage: demux_tools %s" % name)


def test_top_level_help_lists_every_command(run_tools):
    res = run_tools("--help")
    assert res.returncode == 0
    for name in cli.COMMANDS:
        assert "\n  %s " % name in res.stdout
    res = run_tools()
    assert res.returncode == 2 and "commands:" in res.stderr


def test_unknown_command(run_tools):
    res = run_tools("vireo-pair")
    assert res.returncode == 2
    assert "did you mean vireo-pairs" in res.stderr


def test_run_by_path_from_another_folder_writes_no_bytecode(tmp_path, run_tools, tools_dir):
    copy = tmp_path / "checkout" / "tools"
    shutil.copytree(str(tools_dir / "demux_tools"), str(copy / "demux_tools"),
                    ignore=shutil.ignore_patterns("__pycache__", "*.pyc"))
    elsewhere = tmp_path / "elsewhere"
    elsewhere.mkdir()
    res = run_tools("version", cwd=elsewhere, package=str(copy / "demux_tools"))
    assert res.returncode == 0, res.stderr
    assert res.stdout.startswith("demux_tools ")
    res = run_tools("synth", "--out", tmp_path / "s", cwd=elsewhere, package=str(copy / "demux_tools"))
    assert res.returncode == 0, res.stderr
    written = [os.path.join(d, f) for d, dirs, files in os.walk(str(copy)) for f in dirs + files]
    assert not [p for p in written if "__pycache__" in p or p.endswith(".pyc")]
    assert os.listdir(str(elsewhere)) == []


def test_main_in_process(capsys):
    assert cli.main(["--help"]) == 0
    assert "vireo-pairs" in capsys.readouterr().out
    assert cli.main(["version"]) == 0
    assert capsys.readouterr().out.startswith("demux_tools ")


def test_command_modules_import_without_scipy(run_tools, tmp_path):
    # Every command's module imports with scipy and vireoSNP blocked (but those in SCIPY_MODULES).
    modules = sorted(set(m for m, _ in cli.COMMANDS.values()) - SCIPY_MODULES
                     | {"_io", "cellsnp", "vireo_files", "synth"})
    script = tmp_path / "imports.py"
    script.write_text("import importlib\nfor m in %r:\n    importlib.import_module('demux_tools.' + m)\n"
                      "print('ok')\n" % modules)
    # run_tools runs a package's __main__; here a plain script stands in for it.
    pkg = tmp_path / "runner"
    pkg.mkdir()
    (pkg / "__main__.py").write_text("import sys\nsys.path.insert(0, %r)\nexec(open(%r).read())\n"
                                     % (str(cli.__file__).rsplit(os.sep, 2)[0], str(script)))
    res = run_tools(blocked=("scipy", "vireoSNP"), package=str(pkg))
    assert res.returncode == 0, res.stderr
    assert res.stdout.strip() == "ok"
