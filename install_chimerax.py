#!/usr/bin/env python3
"""
ISIS ChimeraX Plugin Installer

Run this script with ChimeraX's Python to install:
    /Applications/ChimeraX-1.10.app/Contents/bin/python3.11 install_chimerax.py

Or run from within ChimeraX's Python shell.
"""

import os
import sys
import subprocess
from pathlib import Path


def find_script_dir():
    """Find the directory containing this script."""
    if hasattr(sys, 'frozen'):
        return Path(sys.executable).parent
    return Path(__file__).parent.resolve()


def install_core_library(script_dir):
    """Install the ISIS core library."""
    print("[INFO] Installing ISIS core library...")
    result = subprocess.run(
        # [ml,plot] extras: without them the MHC models cannot load and
        # `isis plot` is unavailable, giving a working-looking but crippled install.
        [sys.executable, "-m", "pip", "install", "--upgrade",
         f"{script_dir}[ml,plot]"],
        capture_output=True,
        text=True
    )
    if result.returncode != 0:
        print(f"[ERROR] Failed to install core library:")
        print(result.stderr)
        return False
    print("[INFO] Core library installed successfully")
    return True


def install_chimerax_bundle(script_dir):
    """Install the ChimeraX bundle."""
    bundle_dir = script_dir / "src" / "isis_chimerax"

    if not bundle_dir.exists():
        print(f"[ERROR] Bundle directory not found: {bundle_dir}")
        return False

    # Try to import ChimeraX
    try:
        from chimerax.core.commands import run
        from chimerax.core import session
        print("[INFO] Running inside ChimeraX, installing bundle directly...")
        run(session, f'devel install "{bundle_dir}"')
        print("[INFO] Bundle installed successfully")
        return True
    except ImportError:
        pass

    # Not running inside ChimeraX, print instructions
    print("[INFO] Not running inside ChimeraX.")
    print("")
    print("To complete installation, run this command in ChimeraX:")
    print(f'    devel install "{bundle_dir}"')
    print("")
    print("Or run this script using ChimeraX's Python:")
    print(f"    /path/to/ChimeraX/bin/python3 {__file__}")
    return True


def main():
    print("")
    print("==================================")
    print("  ISIS ChimeraX Plugin Installer")
    print("==================================")
    print("")
    print(f"Python: {sys.executable}")
    print(f"Version: {sys.version}")
    print("")

    script_dir = find_script_dir()
    print(f"[INFO] Script directory: {script_dir}")

    # Install core library
    if not install_core_library(script_dir):
        sys.exit(1)

    # Install ChimeraX bundle
    install_chimerax_bundle(script_dir)

    print("")
    print("==================================")
    print("[INFO] Installation complete!")
    print("==================================")
    print("")
    print("Next steps:")
    print("  1. Restart ChimeraX")
    print("  2. Open a structure: open 1ubq")
    print("  3. Run prediction:   isis predict #1")
    print("  4. Color by scores:  isis color #1")
    print("")


if __name__ == "__main__":
    main()
