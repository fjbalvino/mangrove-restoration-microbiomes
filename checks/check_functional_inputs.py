#!/usr/bin/env python3
"""Validate the four externally supplied read-based matrices without fitting models."""
import argparse
import csv
import hashlib
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--base", type=Path, default=ROOT)
    args = parser.parse_args()
    failures = []
    with (ROOT / "docs/FUNCTIONAL_INPUTS_REQUIRED.tsv").open() as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    for row in rows:
        path = args.base / row["path"]
        if not path.is_file():
            failures.append("Supply external matrix: " + row["path"])
            continue
        if path.stat().st_size != int(row["bytes"]):
            failures.append("Matrix size differs: " + row["path"])
        elif hashlib.sha256(path.read_bytes()).hexdigest() != row["sha256"]:
            failures.append("Matrix SHA-256 differs: " + row["path"])
    print(json.dumps({"status": "FAIL" if failures else "PASS", "checks": len(rows),
                      "failures": failures, "guide": "docs/FUNCTIONAL_REFERENCE.md"}, indent=2))
    return 1 if failures else 0


if __name__ == "__main__":
    raise SystemExit(main())
