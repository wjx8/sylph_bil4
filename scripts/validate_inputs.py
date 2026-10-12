#!/usr/bin/env python3
"""Validate a sample_id/r1/r2 TSV without reading the FASTQ contents."""

import argparse
import csv
import re
from pathlib import Path


REQUIRED_COLUMNS = {"sample_id", "r1", "r2"}
SAMPLE_ID_PATTERN = re.compile(r"^[A-Za-z0-9._-]+$")
MATE_PATTERN = re.compile(
    r"([._-])(?:R)?[12](?:[._-]?\d+)?(?=\.(?:fastq|fq)(?:\.gz)?$)",
    flags=re.IGNORECASE,
)


def resolve_input_path(value: str, manifest_path: Path) -> Path:
    path = Path(value).expanduser()
    if not path.is_absolute():
        path = manifest_path.parent / path
    return path.resolve()


def mate_key(path: Path) -> str:
    return MATE_PATTERN.sub(r"\1M", path.name)


def read_manifest(path: Path):
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        columns = set(reader.fieldnames or [])
        missing = REQUIRED_COLUMNS - columns
        if missing:
            raise SystemExit(
                f"ERROR: manifest is missing columns: {', '.join(sorted(missing))}"
            )
        rows = list(reader)
    if not rows:
        raise SystemExit("ERROR: manifest has no sample rows")
    return rows


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--manifest", required=True, type=Path)
    args = parser.parse_args()
    manifest = args.manifest.expanduser().resolve()
    if not manifest.is_file():
        raise SystemExit(f"ERROR: missing manifest: {manifest}")

    rows = read_manifest(manifest)
    sample_ids = []
    for line_number, row in enumerate(rows, start=2):
        sample_id = (row.get("sample_id") or "").strip()
        if not sample_id or not SAMPLE_ID_PATTERN.fullmatch(sample_id):
            raise SystemExit(
                f"ERROR: invalid sample_id on line {line_number}: {sample_id!r}"
            )
        sample_ids.append(sample_id)

        paths = []
        for column in ("r1", "r2"):
            value = (row.get(column) or "").strip()
            if not value:
                raise SystemExit(f"ERROR: empty {column} for {sample_id}")
            path = resolve_input_path(value, manifest)
            if not path.is_file() or path.stat().st_size == 0:
                raise SystemExit(f"ERROR: missing or empty {column} for {sample_id}: {path}")
            paths.append(path)

        if mate_key(paths[0]) != mate_key(paths[1]):
            raise SystemExit(
                f"ERROR: R1/R2 do not look like a matched pair for {sample_id}: "
                f"{paths[0].name}, {paths[1].name}"
            )

    duplicates = sorted({sample for sample in sample_ids if sample_ids.count(sample) > 1})
    if duplicates:
        raise SystemExit(f"ERROR: duplicate sample_id: {', '.join(duplicates)}")

    print(f"OK: {len(rows)} samples; all FASTQ paths exist and are non-empty")


if __name__ == "__main__":
    main()
