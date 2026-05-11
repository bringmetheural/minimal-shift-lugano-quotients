#!/usr/bin/env python3
"""Minimal entry point: reproduce the quotient bound and partition enumeration."""
import subprocess
import sys
from pathlib import Path

root = Path(__file__).resolve().parents[1]
subprocess.check_call([sys.executable, str(root / "scripts" / "reproduce_all.py"), "--out", str(root / "results_minimal")])
