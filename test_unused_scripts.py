#!/usr/bin/env python3
"""
Simple test script for detect_unused_scripts.py

This validates that the script can run and produces expected output.
"""

import subprocess
import sys
from pathlib import Path


def test_script_runs():
    """Test that the script runs without errors."""
    result = subprocess.run([sys.executable, "detect_unused_scripts.py"], 
                          capture_output=True, text=True)
    
    if result.returncode != 0:
        print(f"Script failed with return code {result.returncode}")
        print("STDERR:", result.stderr)
        return False
    
    output = result.stdout
    
    # Check for expected output sections
    expected_sections = [
        "Scanning for scripts in workflow/scripts",
        "Scanning workflow files for script references",
        "UNUSED SCRIPTS DETECTION RESULTS",
        "Summary:"
    ]
    
    for section in expected_sections:
        if section not in output:
            print(f"Missing expected section: {section}")
            return False
    
    # Check that it found some scripts
    if "Found 0 scripts in workflow/scripts/" in output:
        print("Expected to find some scripts")
        return False
    
    print("✓ Script runs successfully and produces expected output")
    return True


def test_verbose_mode():
    """Test that verbose mode produces additional output."""
    result = subprocess.run([sys.executable, "detect_unused_scripts.py", "-v"], 
                          capture_output=True, text=True)
    
    if result.returncode != 0:
        print(f"Verbose mode failed with return code {result.returncode}")
        return False
    
    # Verbose mode should show "used in:" for referenced scripts
    if "used in:" not in result.stdout:
        print("Verbose mode should show 'used in:' information")
        return False
    
    print("✓ Verbose mode works correctly")
    return True


def main():
    """Run tests for the unused scripts detection tool."""
    print("Testing detect_unused_scripts.py...")
    
    # Change to the repo root
    original_dir = Path.cwd()
    repo_root = Path(__file__).parent
    
    try:
        os.chdir(repo_root)
        
        tests = [
            test_script_runs,
            test_verbose_mode
        ]
        
        passed = 0
        for test in tests:
            if test():
                passed += 1
            else:
                print(f"✗ Test {test.__name__} failed")
        
        print(f"\nResults: {passed}/{len(tests)} tests passed")
        
        if passed == len(tests):
            print("🎉 All tests passed!")
            return 0
        else:
            print("❌ Some tests failed")
            return 1
            
    finally:
        os.chdir(original_dir)


if __name__ == "__main__":
    import os
    sys.exit(main())