# demux_tools

Helpers for experiments on this workflow's pileups, donor VCFs and calls. They
sit outside `workflow/`, and the workflow does not use them.

Run them by path from any checkout:

```bash
python -I <checkout>/tools/demux_tools --help
python -I <checkout>/tools/demux_tools <command> --help
```

The commands need numpy and pandas, and keep to Python 3.7 so that they also
run in the vireo image. A run writes only the outputs its command line names,
plus a provenance sidecar for each, and leaves no `__pycache__` in
the checkout.

Run the tests with:

```bash
pytest -q -p no:cacheprovider tools/tests
```

The documentation is [docs/experiment_tools.md](../docs/experiment_tools.md).
