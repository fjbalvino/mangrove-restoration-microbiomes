#!/usr/bin/env python3
"""Verify locally generated figure-input files before rendering.

SHA256SUMS.json deliberately excludes itself. Integrity is relative to this
manifest; this is not a digital signature or authentication of scientific truth.
"""
import argparse
import hashlib
import json
from pathlib import Path

def verify(root):
    manifest = json.loads((root / "inputs" / "SHA256SUMS.json").read_text())
    errors = []
    for name, expected in manifest.items():
        path = root / "inputs" / name
        if not path.is_file():
            errors.append("MISSING " + name)
            continue
        if hashlib.sha256(path.read_bytes()).hexdigest() != expected:
            errors.append("MODIFIED " + name)
    if errors:
        raise SystemExit("\n".join(errors))
    print(f"PACKAGE PASS | {len(manifest)} archivos SHA-256 verificados")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=Path(__file__).resolve().parent.parent)
    verify(parser.parse_args().root)
