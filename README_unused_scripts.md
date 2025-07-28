# Unused Scripts Detection

This script (`detect_unused_scripts.py`) identifies scripts in the `workflow/scripts/` directory that are not referenced anywhere in the workflow.

## Usage

Run the script from the repository root:

```bash
python detect_unused_scripts.py
```

## What it does

The script:
1. Scans all files in `workflow/scripts/` to find script files
2. Searches through all workflow files (`.smk`, `Snakefile`, `.Rmd`, `.yaml`, `.yml`) for script references
3. Reports scripts that exist but are not referenced anywhere
4. Provides a summary showing total, referenced, and unused scripts

## Output

The script provides:
- List of all scripts found in `workflow/scripts/`
- List of files where script references were found
- **Main result**: List of unused scripts that can potentially be removed
- Summary statistics
- List of referenced scripts for verification

## Script Reference Patterns Detected

The script looks for various ways scripts can be referenced:
- `script: "../scripts/filename.py"` (Snakemake script directive)
- `scripts/filename.py` (direct references)
- `python workflow/scripts/filename.py` (shell commands)
- `Rscript workflow/scripts/filename.R` (R script commands)
- Include/source statements in R files

## Important Notes

- **Always verify manually** before deleting any scripts identified as unused
- Some scripts might be used in ways not detected by this tool
- The script only scans workflow files, not external documentation or comments
- Scripts marked as unused are safe candidates for cleanup, but manual verification is recommended

## Example Output

```
Found 15 unused scripts:
  - demuxafy_count_cells.py
  - demuxafy_counts.py
  - merge_plp.py
  - merge_umi.py
  - run_vireo.sh
  [...]

Summary:
  Total scripts: 29
  Referenced scripts: 14
  Unused scripts: 15
```