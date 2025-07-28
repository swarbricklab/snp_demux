#!/usr/bin/env python3
"""
Script to detect unused scripts in the workflow/scripts directory.

This script scans all Snakemake workflow files to find script references
and reports any scripts in workflow/scripts that are not being used.
"""

import os
import re
import sys
from pathlib import Path
from typing import Set, List


def find_all_scripts(scripts_dir: Path) -> Set[str]:
    """Find all script files in the scripts directory."""
    if not scripts_dir.exists():
        return set()
    
    scripts = set()
    for file_path in scripts_dir.iterdir():
        if file_path.is_file():
            scripts.add(file_path.name)
    
    return scripts


def scan_file_for_script_references(file_path: Path) -> Set[str]:
    """Scan a single file for script references."""
    script_refs = set()
    
    try:
        with open(file_path, 'r', encoding='utf-8') as f:
            content = f.read()
    except (UnicodeDecodeError, PermissionError):
        # Skip files that can't be read
        return script_refs
    
    # Pattern 1: script: "../scripts/filename" or script: "scripts/filename"
    script_directive_pattern = r'script:\s*["\'].*?scripts/([^"\']+)["\']'
    matches = re.findall(script_directive_pattern, content)
    script_refs.update(matches)
    
    # Pattern 2: Direct references to scripts/ directory in various contexts
    direct_script_pattern = r'(?:scripts/|workflow/scripts/)([a-zA-Z0-9_.-]+\.(py|R|sh|Rmd))'
    matches = re.findall(direct_script_pattern, content)
    script_refs.update([match[0] for match in matches])
    
    # Pattern 3: Shell commands that might reference scripts
    # python workflow/scripts/script.py or Rscript workflow/scripts/script.R
    shell_command_pattern = r'(?:python|Rscript|bash|sh)\s+.*?(?:workflow/)?scripts/([a-zA-Z0-9_.-]+\.(py|R|sh|Rmd))'
    matches = re.findall(shell_command_pattern, content)
    script_refs.update([match[0] for match in matches])
    
    # Pattern 4: Include or source statements that might reference scripts
    include_pattern = r'(?:include|source|import)\s*\(\s*["\'].*?scripts/([a-zA-Z0-9_.-]+\.(py|R|sh|Rmd))["\']'
    matches = re.findall(include_pattern, content)
    script_refs.update([match[0] for match in matches])
    
    # Pattern 5: References without scripts/ prefix but in workflow context
    # Only check this for non-script files to avoid self-references
    if not file_path.name.endswith(('.py', '.R', '.sh', '.Rmd')):
        workflow_script_pattern = r'["\']([a-zA-Z0-9_.-]+\.(py|R|sh|Rmd))["\']'
        potential_matches = re.findall(workflow_script_pattern, content)
        for match in potential_matches:
            script_name = match[0]
            # Only consider it a reference if the script actually exists
            if (Path("workflow/scripts") / script_name).exists():
                script_refs.add(script_name)
    
    return script_refs


def scan_workflow_files(workflow_dir: Path) -> Set[str]:
    """Scan all workflow files for script references."""
    referenced_scripts = set()
    
    # File patterns to scan
    patterns = [
        "**/*.smk",
        "**/Snakefile",
        "**/*.Rmd",
        "**/*.yaml",
        "**/*.yml"
    ]
    
    for pattern in patterns:
        for file_path in workflow_dir.glob(pattern):
            if file_path.is_file():
                refs = scan_file_for_script_references(file_path)
                if refs:
                    print(f"Found script references in {file_path}: {refs}")
                referenced_scripts.update(refs)
    
    return referenced_scripts


def main():
    """Main function to detect unused scripts."""
    # Get the repository root (assuming script is run from repo root)
    repo_root = Path.cwd()
    workflow_dir = repo_root / "workflow"
    scripts_dir = workflow_dir / "scripts"
    
    if not workflow_dir.exists():
        print("Error: workflow directory not found. Please run this script from the repository root.")
        sys.exit(1)
    
    if not scripts_dir.exists():
        print("Error: workflow/scripts directory not found.")
        sys.exit(1)
    
    print("Scanning for scripts in workflow/scripts...")
    all_scripts = find_all_scripts(scripts_dir)
    print(f"Found {len(all_scripts)} scripts in workflow/scripts/")
    for script in sorted(all_scripts):
        print(f"  - {script}")
    
    print("\nScanning workflow files for script references...")
    referenced_scripts = scan_workflow_files(workflow_dir)
    print(f"Found {len(referenced_scripts)} script references in workflow files")
    
    # Find unused scripts
    unused_scripts = all_scripts - referenced_scripts
    
    print("\n" + "="*60)
    print("UNUSED SCRIPTS DETECTION RESULTS")
    print("="*60)
    
    if unused_scripts:
        print(f"\nFound {len(unused_scripts)} unused scripts:")
        for script in sorted(unused_scripts):
            print(f"  - {script}")
        
        print(f"\nThese scripts exist in workflow/scripts/ but are not referenced")
        print(f"in any workflow files (.smk, Snakefile, .Rmd, .yaml, .yml files).")
        print(f"\nIf these scripts are truly unused, they can be safely removed.")
        print(f"However, please verify manually before deletion.")
    else:
        print("\nNo unused scripts found! All scripts in workflow/scripts/ are being used.")
    
    print(f"\nSummary:")
    print(f"  Total scripts: {len(all_scripts)}")
    print(f"  Referenced scripts: {len(referenced_scripts)}")
    print(f"  Unused scripts: {len(unused_scripts)}")
    
    # Show referenced scripts for verification
    if referenced_scripts:
        print(f"\nReferenced scripts:")
        for script in sorted(referenced_scripts):
            print(f"  - {script}")
    
    return len(unused_scripts)


if __name__ == "__main__":
    main()