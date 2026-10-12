#!/usr/bin/env python3
"""Run the four-reference Sylph BIL4 pipeline."""

from __future__ import annotations

import argparse
import csv
import datetime as dt
import os
import shutil
import subprocess
import sys
import tempfile
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

SCRIPT_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(SCRIPT_DIR))

from process_results import process_results  # noqa: E402
from validate_assets import main as validate_assets_main  # noqa: E402
from validate_inputs import read_manifest, resolve_input_path  # noqa: E402


def parse_args() -> argparse.Namespace:
    package_root = SCRIPT_DIR.parent
    parser = argparse.ArgumentParser(
        description="Quantify BI/BL/BS/BX from paired-end clean FASTQ with Sylph."
    )
    parser.add_argument("--package-root", type=Path, default=package_root, help=argparse.SUPPRESS)
    parser.add_argument("--manifest", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--threads", type=int, default=8)
    parser.add_argument("--jobs", type=int, default=4)
    parser.add_argument("--sylph", default=os.environ.get("SYLPH_BIN", "sylph"))
    return parser.parse_args()


def run_checked(command, log_path: Path) -> None:
    log_path.parent.mkdir(parents=True, exist_ok=True)
    with log_path.open("w") as log:
        log.write("$ " + " ".join(str(part) for part in command) + "\n\n")
        log.flush()
        subprocess.run(command, stdout=log, stderr=subprocess.STDOUT, check=True)


def resolve_command(command: str) -> str:
    path = Path(command)
    if path.parent != Path("."):
        if not path.is_file() or not os.access(path, os.X_OK):
            raise SystemExit(f"ERROR: Sylph executable is not available: {command}")
        return str(path)
    resolved = shutil.which(command)
    if resolved is None:
        raise SystemExit(f"ERROR: Sylph executable is not in PATH: {command}")
    return resolved


def load_input_rows(manifest_path: Path):
    rows = read_manifest(manifest_path)
    loaded = []
    for row in rows:
        loaded.append(
            {
                "sample_id": row["sample_id"].strip(),
                "r1": str(resolve_input_path(row["r1"].strip(), manifest_path)),
                "r2": str(resolve_input_path(row["r2"].strip(), manifest_path)),
            }
        )
    return loaded


def write_reads_manifest(path: Path, rows) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["sample_id", "r1", "r2"])
        writer.writeheader()
        writer.writerows(rows)


def existing_reads_manifest(path: Path):
    if not path.is_file():
        return None
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle))


def sketch_one(sample, output_dir: Path, sylph: str) -> str:
    sample_id = sample["sample_id"]
    final = output_dir / "sketches" / f"{sample_id}.paired.sylsp"
    done = output_dir / "logs" / f"{sample_id}.sketch.done"
    log = output_dir / "logs" / f"{sample_id}.sketch.log"
    if final.is_file() and final.stat().st_size > 0 and done.is_file():
        return f"reused {sample_id}"

    with tempfile.TemporaryDirectory(prefix=f".{sample_id}.", dir=output_dir / "sketches") as temp:
        command = [
            sylph,
            "sketch",
            "-1",
            sample["r1"],
            "-2",
            sample["r2"],
            "-S",
            sample_id,
            "-d",
            temp,
            "-t",
            "1",
        ]
        run_checked(command, log)
        generated = sorted(Path(temp).glob("*.sylsp"))
        if len(generated) != 1 or generated[0].stat().st_size == 0:
            raise RuntimeError(f"Sylph did not create exactly one sketch for {sample_id}")
        final.unlink(missing_ok=True)
        shutil.move(str(generated[0]), final)
    done.touch()
    return f"completed {sample_id}"


def main() -> None:
    args = parse_args()
    if args.threads < 1 or args.jobs < 1:
        raise SystemExit("ERROR: --threads and --jobs must be positive integers")

    package_root = args.package_root.resolve()
    manifest_path = args.manifest.expanduser().resolve()
    output_dir = args.output.expanduser().resolve()
    database = package_root / "database/sylph_BIL4/BIL4.syldb"
    reference_manifest = package_root / "database/sylph_BIL4/reference_manifest.tsv"
    sylph = resolve_command(args.sylph)

    version = subprocess.check_output([sylph, "--version"], text=True).strip()
    if version != "sylph 0.5.1":
        raise SystemExit(f"ERROR: expected sylph 0.5.1, observed: {version}")

    # Validate the bundled references through the same code used by the shell check.
    old_argv = sys.argv
    sys.argv = ["validate_assets.py", "--root", str(package_root)]
    try:
        validate_assets_main()
    finally:
        sys.argv = old_argv

    rows = load_input_rows(manifest_path)
    for subdir in ("sketches", "raw", "processed", "qc", "logs"):
        (output_dir / subdir).mkdir(parents=True, exist_ok=True)
    reads_manifest = output_dir / "qc/sample_clean_reads.csv"
    previous_rows = existing_reads_manifest(reads_manifest)
    if previous_rows is not None and previous_rows != rows:
        raise SystemExit(
            f"ERROR: {reads_manifest} belongs to a different input manifest; "
            "choose a new --output directory"
        )
    write_reads_manifest(reads_manifest, rows)

    run_info = output_dir / "run_info.txt"
    run_info.write_text(
        "".join(
            [
                "method=sylph_bil4\n",
                f"date={dt.datetime.now().astimezone().isoformat(timespec='seconds')}\n",
                f"sylph_version={version}\n",
                f"threads={args.threads}\n",
                f"jobs={args.jobs}\n",
                "profile_minimum_ani=97\n",
                "BI_BL_BS_postprocess_ani=97\n",
                "BX_postprocess_ani=99\n",
                f"reference_manifest={reference_manifest}\n",
                f"database={database}\n",
            ]
        )
    )

    failures = []
    with ThreadPoolExecutor(max_workers=args.jobs) as executor:
        futures = {
            executor.submit(sketch_one, sample, output_dir, sylph): sample["sample_id"]
            for sample in rows
        }
        for future in as_completed(futures):
            sample_id = futures[future]
            try:
                print(future.result())
            except Exception as exc:  # noqa: BLE001
                failures.append((sample_id, str(exc)))

    if failures:
        for sample_id, error in failures:
            print(f"FAILED {sample_id}: {error}", file=sys.stderr)
        raise SystemExit(1)

    sketch_paths = [output_dir / "sketches" / f"{row['sample_id']}.paired.sylsp" for row in rows]
    if not all(path.is_file() and path.stat().st_size > 0 for path in sketch_paths):
        raise SystemExit("ERROR: one or more expected Sylph sketches are missing")

    profile_inputs = output_dir / "logs/profile_inputs.txt"
    profile_inputs.write_text("".join(f"{path}\n" for path in sketch_paths))
    raw_output = output_dir / "raw/BIL4_sylph_raw.tsv"
    profile_log = output_dir / "logs/profile.log"
    profile_command = [
        sylph,
        "profile",
        str(database),
        "-l",
        str(profile_inputs),
        "-m",
        "97",
        "-t",
        str(args.threads),
        "-o",
        str(raw_output),
    ]
    run_checked(profile_command, profile_log)
    # Sylph exits successfully but writes a zero-byte file when no reference
    # genome passes the profile threshold. That is a valid all-zero result.
    if not raw_output.is_file():
        raise SystemExit(f"ERROR: missing Sylph profile output: {raw_output}")

    process_results(
        raw_path=raw_output,
        reads_path=reads_manifest,
        manifest_path=reference_manifest,
        output_dir=output_dir,
        database_path=database,
    )
    (output_dir / "logs/BIL4.all_done").touch()
    print(f"DONE: {output_dir / 'processed/BIL4_abundance_wide.csv'}")


if __name__ == "__main__":
    main()
