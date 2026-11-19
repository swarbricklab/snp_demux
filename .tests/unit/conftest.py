"""Shared helpers for the unit tests of workflow/scripts/*.py.

Snakemake runs those scripts with an injected global `snakemake`; without it
they only define functions, so the tests import them by path.
"""
import importlib.util
from pathlib import Path

SCRIPTS = Path(__file__).resolve().parents[2] / "workflow" / "scripts"


def load_script(name):
    spec = importlib.util.spec_from_file_location(name, SCRIPTS / f"{name}.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod
