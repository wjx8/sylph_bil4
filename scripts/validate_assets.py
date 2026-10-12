#!/usr/bin/env python3
"""Validate the bundled reference FASTA files and Sylph database."""

import argparse
import csv
import hashlib
from pathlib import Path


def md5(path: Path) -> str:
    digest = hashlib.md5()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--root",
        type=Path,
        default=Path(__file__).resolve().parents[1],
        help="method package root (default: package containing this script)",
    )
    parser.add_argument(
        "--skip-database",
        action="store_true",
        help="validate reference FASTA files without requiring BIL4.syldb",
    )
    args = parser.parse_args()
    root = args.root.resolve()
    manifest_path = root / "database/sylph_BIL4/reference_manifest.tsv"
    database_path = root / "database/sylph_BIL4/BIL4.syldb"

    if not manifest_path.is_file():
        raise SystemExit(f"ERROR: missing reference manifest: {manifest_path}")
    if not args.skip_database and (not database_path.is_file() or database_path.stat().st_size == 0):
        raise SystemExit(f"ERROR: missing or empty Sylph database: {database_path}")

    with manifest_path.open(newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))

    expected = {"BI", "BL", "BS", "BX"}
    observed = {row.get("clade", "") for row in rows}
    if observed != expected:
        raise SystemExit(f"ERROR: expected clades {sorted(expected)}, observed {sorted(observed)}")

    for row in rows:
        path = root / row["fasta_path"]
        if not path.is_file() or path.stat().st_size == 0:
            raise SystemExit(f"ERROR: missing or empty reference FASTA: {path}")
        observed_md5 = md5(path)
        if observed_md5 != row["md5"]:
            raise SystemExit(
                f"ERROR: MD5 mismatch for {row['clade']}: "
                f"expected {row['md5']}, observed {observed_md5}"
            )
        print(f"OK {row['clade']} {path.name} md5={observed_md5}")

    if args.skip_database:
        print("OK reference FASTA assets")
    else:
        print(f"OK database {database_path} bytes={database_path.stat().st_size}")


if __name__ == "__main__":
    main()
