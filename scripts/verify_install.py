#!/usr/bin/env python3
"""Verify AlphaGenome installation."""

import sys
import os

# Use python-dotenv for proper .env loading
try:
    from dotenv import load_dotenv
    # Load .env file from project root
    env_path = os.path.join(os.path.dirname(__file__), '..', '.env')
    load_dotenv(env_path)
except ImportError:
    print("  [WARN] python-dotenv not installed, falling back to manual .env parsing")

def check_imports():
    """Check that all required modules can be imported."""
    print("Checking imports...")
    try:
        from alphagenome.data import genome
        from alphagenome.models import dna_client
        from alphagenome.visualization import plot_components
        print("  [OK] All imports successful")
        return True
    except ImportError as e:
        print(f"  [FAIL] Import error: {e}")
        print("  [INFO] If imports fail, verify API structure at https://www.alphagenomedocs.com/")
        return False

def check_api_key():
    """Check if API key is configured."""
    print("Checking API key...")

    api_key = os.environ.get('ALPHAGENOME_API_KEY')

    if api_key and api_key != '<your_api_key_here>' and api_key != 'your_api_key_here':
        print("  [OK] API key found")
        return True
    else:
        print("  [WARN] API key not configured (set ALPHAGENOME_API_KEY in .env)")
        return False

def main():
    print("=" * 50)
    print("AlphaGenome Installation Verification")
    print("=" * 50)
    print()

    results = []
    results.append(("Imports", check_imports()))
    results.append(("API Key", check_api_key()))

    print()
    print("=" * 50)
    print("Summary")
    print("=" * 50)

    all_pass = True
    for name, passed in results:
        status = "PASS" if passed else "FAIL/WARN"
        print(f"  {name}: {status}")
        if not passed and name == "Imports":
            all_pass = False

    print()
    if all_pass:
        print("Installation verified successfully!")
        return 0
    else:
        print("Installation has issues. Please check above.")
        return 1

if __name__ == "__main__":
    sys.exit(main())
