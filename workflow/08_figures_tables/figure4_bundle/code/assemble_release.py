#!/usr/bin/env python3
"""Create release inventory, hash manifest and ZIP from the documented tree.

Run after intentional code/documentation changes and validation. This script
never changes inputs/SHA256SUMS.json; altered scientific inputs still fail the
plotting script. The release manifest excludes itself and its generated archive.
"""
import argparse
import csv
import hashlib
import json
from pathlib import Path
import zipfile

def files(root):
    return sorted(p for p in root.rglob("*") if p.is_file() and not any(
        part in {".venv", ".git", "__pycache__"} for part in p.relative_to(root).parts))

def purpose(path):
    if str(path).startswith("outputs/"):
        return "Salida/copia de procedencia de una ejecución, según README."
    return {"inputs": "Entrada congelada o su manifiesto de integridad.",
        "code": "Código entregado, inspeccionable y ejecutable.",
        "reference": "Fuente histórica preservada; no ejecutada en esta etapa.",
        "docs": "Documentación y trazabilidad.",
        "config": "Configuración gráfica explícita."}.get(path.parts[0], "Guía, entorno, lanzador o manifiesto de distribución.")

def assemble(root, archive):
    inventory = root / "docs" / "FILE_INVENTORY.tsv"
    manifest = root / "SHA256SUMS.json"
    with inventory.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["relative_path", "bytes", "purpose"])
        for path in files(root):
            if path in (inventory, manifest):
                continue
            writer.writerow([str(path.relative_to(root)), path.stat().st_size, purpose(path.relative_to(root))])
    manifest.write_text(json.dumps({str(p.relative_to(root)): hashlib.sha256(p.read_bytes()).hexdigest() for p in files(root) if p != manifest}, indent=2) + "\n")
    if archive:
        archive = archive.resolve()
        if archive.is_relative_to(root):
            raise ValueError("Archive must be outside the package directory")
        with zipfile.ZipFile(archive, "w", compression=zipfile.ZIP_DEFLATED) as z:
            for path in files(root):
                z.write(path, path.relative_to(root.parent))
        digest = hashlib.sha256(archive.read_bytes()).hexdigest()
        archive.with_suffix(archive.suffix + ".sha256").write_text(f"{digest}  {archive.name}\n")
        print(f"Release: {archive}\nSHA256: {digest}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=Path(__file__).resolve().parent.parent)
    parser.add_argument("--archive", type=Path)
    args = parser.parse_args()
    assemble(args.root.resolve(), args.archive)
