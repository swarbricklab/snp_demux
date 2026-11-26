"""Command table and dispatcher for `python -I <checkout>/tools/demux_tools <command>`.

Only the module of the command that runs is imported, so a command's --help
needs only that module's libraries.

A command's module defines two functions, named after the command with "-"
written "_":

    configure_<name>(parser)   add arguments to an argparse.ArgumentParser
    run_<name>(args)           run with the parsed arguments; return an exit status

`args.argv` holds the command line ("demux_tools", name, options...) for the
provenance sidecar. Raise demux_tools._io.UsageError for a usage problem: the
dispatcher prints its message and exits with status 2.
"""
import argparse
import difflib
import importlib
import sys

# One line per command, in alphabetical order: "name": ("module", "one-line help"),
# A new command is added as one line in alphabetical order.
COMMANDS = {
    "score": ("scorer", "reproduce vireo's fixed-mode calls for one capture, with each cell's evidence"),
    "score-many": ("scorer", "score the captures of a manifest in parallel, with a summary and validation table"),
    "site-terms": ("scorer", "one cell's log-likelihood under each donor, site by site"),
    "sites": ("sites", "site lists for the Scorer: annotate with a GTF and INFO, thin, take from a panel"),
    "synth": ("synth", "write a synthetic capture: cellsnp folder, donor VCF and truth table"),
    "validate-vireo": ("scorer", "compare Scorer calls with vireo's donor_ids.tsv cell by cell"),
    "version": ("cli", "print the package version, its git commit and library versions"),
    "vireo-pairs": ("vireo_files", "add the true doublet pair to vireo's donor_ids.tsv"),
}


def _usage():
    width = max(len(name) for name in COMMANDS)
    lines = ["usage: python -I <checkout>/tools/demux_tools <command> [options]", "",
             "Helpers for experiments on demultiplexing outputs (docs/experiment_tools.md).", "",
             "commands:"]
    for name, (_, text) in COMMANDS.items():
        lines.append("  %-*s  %s" % (width, name, text))
    lines += ["", "Run '<command> --help' for a command's options."]
    return "\n".join(lines)


def command_functions(name):
    """(configure, run) of a registered command; imports its module."""
    module_name, _ = COMMANDS[name]
    module = importlib.import_module("demux_tools." + module_name)
    stem = name.replace("-", "_")
    return getattr(module, "configure_" + stem), getattr(module, "run_" + stem)


def main(argv=None):
    """Dispatch `argv` (without the program name); returns the exit status."""
    from demux_tools._io import UsageError

    argv = list(sys.argv[1:] if argv is None else argv)
    if not argv:
        sys.stderr.write(_usage() + "\n")
        return 2
    if argv[0] in ("-h", "--help"):
        print(_usage())
        return 0
    if argv[0] == "--version":
        argv[0] = "version"
    name = argv[0]
    if name not in COMMANDS:
        close = difflib.get_close_matches(name, list(COMMANDS), n=3)
        hint = " (did you mean %s?)" % ", ".join(close) if close else ""
        sys.stderr.write("demux_tools: unknown command %r%s\n\n%s\n" % (name, hint, _usage()))
        return 2
    configure, run = command_functions(name)
    parser = argparse.ArgumentParser(prog="demux_tools " + name, description=COMMANDS[name][1])
    configure(parser)
    args = parser.parse_args(argv[1:])
    args.argv = ["demux_tools"] + argv
    try:
        status = run(args)
    except UsageError as e:
        sys.stderr.write("demux_tools %s: error: %s\n" % (name, e))
        return 2
    return int(status or 0)


# ---------------------------------------------------------------------------
# Command: version


def configure_version(parser):
    parser.description = "Print the package version, the git commit of its checkout and library versions."


def run_version(args):
    from demux_tools import __version__
    from demux_tools._io import git_sha, library_versions

    print("demux_tools %s" % __version__)
    print("git %s" % (git_sha() or "unknown"))
    for key, value in library_versions().items():
        print("%s %s" % (key, value or "not found"))
    return 0
