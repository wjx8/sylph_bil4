#!/usr/bin/env python3
"""Convert native Sylph BIL4 output into portable CSV tables and QC."""

from __future__ import annotations

import csv
import hashlib
import pathlib
import statistics
from collections import defaultdict


EXPECTED_CLADES = ("BI", "BL", "BS", "BX")


def normalized(text):
    return "".join(char.lower() for char in str(text) if char.isalnum())


def get_field(row, names):
    normalized_row = {normalized(key): value for key, value in row.items()}
    for name in names:
        if normalized(name) in normalized_row:
            return normalized_row[normalized(name)]
    return ""


def as_float(value):
    try:
        return float(value)
    except (TypeError, ValueError):
        return None


def sample_id_from_sylsp(value):
    name = pathlib.Path(value).name
    for suffix in (".paired.sylsp", ".sylsp"):
        if name.endswith(suffix):
            return name[: -len(suffix)]
    return name


def distribution(values):
    values = [value for value in values if isinstance(value, (int, float))]
    if not values:
        return "NA"
    return f"min={min(values):.6g};median={statistics.median(values):.6g};max={max(values):.6g}"


def file_md5(path: pathlib.Path) -> str:
    digest = hashlib.md5()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_delimited_rows(path: pathlib.Path):
    """Read both the new CSV manifest and the legacy TSV manifest."""
    with path.open(newline="") as handle:
        preview = handle.read(4096)
        handle.seek(0)
        delimiter = "\t" if "\t" in preview and "," not in preview else ","
        if delimiter == "\t":
            rows = list(csv.reader(handle, delimiter="\t"))
            if rows and {"sample_id", "r1", "r2"}.issubset(set(rows[0])):
                return [
                    {key: value for key, value in zip(rows[0], values)}
                    for values in rows[1:]
                ]
            return [
                {"sample_id": values[0], "r1": values[1], "r2": values[2]}
                for values in rows
                if len(values) >= 3
            ]
        return list(csv.DictReader(handle, delimiter=delimiter))


def process_results(
    raw_path: pathlib.Path,
    reads_path: pathlib.Path,
    manifest_path: pathlib.Path,
    output_dir: pathlib.Path,
    database_path: pathlib.Path,
) -> None:
    with manifest_path.open(newline="") as handle:
        references = list(csv.DictReader(handle, delimiter="\t"))
    if tuple(row["clade"] for row in references) != EXPECTED_CLADES:
        raise SystemExit(
            "ERROR: reference manifest must contain clades in order BI, BL, BS, BX"
        )
    for row in references:
        row["ani_threshold"] = float(row["ani_threshold"])
        row["genome_name"] = pathlib.Path(row["fasta_path"]).name

    samples = [
        {
            "sample_id": row["sample_id"],
            "r1": row["r1"],
            "r2": row["r2"],
        }
        for row in read_delimited_rows(reads_path)
    ]

    reference_by_genome = {row["genome_name"]: row for row in references}
    profile_rows = defaultdict(dict)
    unmatched_genomes = set()
    abundance_out_of_range = 0

    with raw_path.open(newline="") as handle:
        for raw in csv.DictReader(handle, delimiter="\t"):
            sample_file = get_field(raw, ["Sample_file", "Sample"])
            genome_file = get_field(raw, ["Genome_file", "Genome"])
            sample_id = sample_id_from_sylsp(sample_file)
            genome_name = pathlib.Path(genome_file).name
            reference = reference_by_genome.get(genome_name)
            if reference is None:
                unmatched_genomes.add(genome_name or genome_file)
                continue
            abundance = as_float(get_field(raw, ["Taxonomic_abundance"]))
            if abundance is not None and not 0 <= abundance <= 100:
                abundance_out_of_range += 1
            profile_rows[sample_id][reference["clade"]] = {
                "sample_file": sample_id,
                "genome_file": genome_name,
                "taxonomic_abundance": abundance,
                "adjusted_ani": as_float(get_field(raw, ["Adjusted_ANI"])),
                "eff_cov": as_float(get_field(raw, ["Eff_cov"])),
                "ani_percentile": get_field(raw, ["ANI_5-95_percentile"]),
                "median_cov": as_float(get_field(raw, ["Median_cov"])),
                "mean_cov_geq1": as_float(get_field(raw, ["Mean_cov_geq1"])),
                "containment_ind": get_field(raw, ["Containment_ind"]),
            }

    long_fields = [
        "sample_id",
        "clade",
        "scientific_name",
        "reference_id",
        "Taxonomic_abundance",
        "Adjusted_ANI",
        "ani_threshold",
        "detected",
        "Eff_cov",
        "ANI_5-95_percentile",
        "Median_cov",
        "Mean_cov_geq1",
        "Containment_ind",
        "Sample_file",
        "Genome_file",
    ]
    long_rows = []
    for sample in samples:
        sample_id = sample["sample_id"]
        for reference in references:
            clade = reference["clade"]
            raw = profile_rows.get(sample_id, {}).get(clade, {})
            ani = raw.get("adjusted_ani")
            detected = int(ani is not None and ani >= reference["ani_threshold"])
            abundance = raw.get("taxonomic_abundance") if detected else 0.0
            long_rows.append(
                {
                    "sample_id": sample_id,
                    "clade": clade,
                    "scientific_name": reference["scientific_name"],
                    "reference_id": reference["reference_id"],
                    "Taxonomic_abundance": 0.0 if abundance is None else abundance,
                    "Adjusted_ANI": "NA" if ani is None else ani,
                    "ani_threshold": reference["ani_threshold"],
                    "detected": detected,
                    "Eff_cov": "NA" if raw.get("eff_cov") is None else raw["eff_cov"],
                    "ANI_5-95_percentile": raw.get("ani_percentile") or "NA",
                    "Median_cov": "NA" if raw.get("median_cov") is None else raw["median_cov"],
                    "Mean_cov_geq1": "NA" if raw.get("mean_cov_geq1") is None else raw["mean_cov_geq1"],
                    "Containment_ind": raw.get("containment_ind") or "NA",
                    "Sample_file": raw.get("sample_file", "NA"),
                    "Genome_file": raw.get("genome_file", "NA"),
                }
            )

    processed_dir = output_dir / "processed"
    qc_dir = output_dir / "qc"
    processed_dir.mkdir(parents=True, exist_ok=True)
    qc_dir.mkdir(parents=True, exist_ok=True)

    long_path = processed_dir / "BIL4_long.csv"
    with long_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=long_fields)
        writer.writeheader()
        writer.writerows(long_rows)

    wide_fields = ["sample_id"]
    for suffix in ("abundance", "ANI", "detected"):
        wide_fields.extend(f"{reference['clade']}_{suffix}" for reference in references)
    wide_fields.append("BIL_total_abundance")
    wide_fields.extend(f"{reference['clade']}_fraction_in_BIL" for reference in references)

    wide_path = processed_dir / "BIL4_abundance_wide.csv"
    with wide_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=wide_fields)
        writer.writeheader()
        for sample in samples:
            sample_rows = {
                row["clade"]: row
                for row in long_rows
                if row["sample_id"] == sample["sample_id"]
            }
            result = {"sample_id": sample["sample_id"]}
            for reference in references:
                clade = reference["clade"]
                result[f"{clade}_abundance"] = sample_rows[clade]["Taxonomic_abundance"]
                result[f"{clade}_ANI"] = sample_rows[clade]["Adjusted_ANI"]
                result[f"{clade}_detected"] = sample_rows[clade]["detected"]
            total = sum(float(result[f"{reference['clade']}_abundance"]) for reference in references)
            result["BIL_total_abundance"] = total
            for reference in references:
                clade = reference["clade"]
                result[f"{clade}_fraction_in_BIL"] = (
                    "NA" if total == 0 else float(result[f"{clade}_abundance"]) / total
                )
            writer.writerow(result)

    sketch_success = sum(
        (output_dir / "sketches" / f"{sample['sample_id']}.paired.sylsp").is_file()
        and (output_dir / "sketches" / f"{sample['sample_id']}.paired.sylsp").stat().st_size > 0
        for sample in samples
    )
    duplicate_count = len(samples) - len({sample["sample_id"] for sample in samples})

    def pair_key(value):
        name = pathlib.Path(value).name
        for token in (".clean_1.fastq.gz", ".clean_2.fastq.gz", ".clean_1.fastq", ".clean_2.fastq"):
            if name.endswith(token):
                return name[: -len(token)]
        return name

    pair_mismatch_count = sum(pair_key(sample["r1"]) != pair_key(sample["r2"]) for sample in samples)
    four_rows_ok = all(
        sum(row["sample_id"] == sample["sample_id"] for row in long_rows) == 4
        for sample in samples
    )

    qc_rows = [
        ("total_samples", len(samples)),
        ("successful_sketch_samples", sketch_success),
        ("successful_profile_samples", sketch_success if raw_path.is_file() else 0),
        ("missing_samples", len(samples) - sketch_success),
        ("duplicate_sample_id_count", duplicate_count),
        ("r1_r2_mismatch_count", pair_mismatch_count),
        ("all_samples_have_four_clades", int(four_rows_ok)),
        ("unmatched_genome_file_count", len(unmatched_genomes)),
        ("taxonomic_abundance_out_of_range_count", abundance_out_of_range),
        (
            "BI_BL_BS_threshold_is_97",
            int(all(row["ani_threshold"] == 97 for row in references if row["clade"] in {"BI", "BL", "BS"})),
        ),
        (
            "BX_threshold_is_99",
            int(next(row["ani_threshold"] for row in references if row["clade"] == "BX") == 99),
        ),
    ]
    for reference in references:
        clade_rows = [row for row in long_rows if row["clade"] == reference["clade"]]
        qc_rows.extend(
            [
                (f"{reference['clade']}_detected_samples", sum(row["detected"] for row in clade_rows)),
                (
                    f"{reference['clade']}_Adjusted_ANI_distribution",
                    distribution([as_float(row["Adjusted_ANI"]) for row in clade_rows]),
                ),
                (
                    f"{reference['clade']}_Eff_cov_distribution",
                    distribution([as_float(row["Eff_cov"]) for row in clade_rows]),
                ),
                (
                    f"{reference['clade']}_Taxonomic_abundance_distribution",
                    distribution([as_float(row["Taxonomic_abundance"]) for row in clade_rows]),
                ),
            ]
        )

    qc_path = qc_dir / "BIL4_summary.csv"
    with qc_path.open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["metric", "value"])
        writer.writerows(qc_rows)

    result_readme = output_dir / "README_results.md"
    result_readme.write_text(
        "\n".join(
            [
                "# Sylph BIL4 results",
                "",
                "- Sylph version: 0.5.1",
                "- References: BI, BL, BS, BX",
                "- ANI thresholds: BI/BL/BS >=97; BX >=99",
                "- Abundance field: Sylph `Taxonomic_abundance`",
                f"- Database MD5: {file_md5(database_path)}",
                "",
                "Primary handoff: `processed/BIL4_abundance_wide.csv`.",
                "",
            ]
        )
    )


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument("--raw", required=True, type=pathlib.Path)
    parser.add_argument("--reads", required=True, type=pathlib.Path)
    parser.add_argument("--manifest", required=True, type=pathlib.Path)
    parser.add_argument("--output", required=True, type=pathlib.Path)
    parser.add_argument("--database", required=True, type=pathlib.Path)
    cli = parser.parse_args()
    process_results(cli.raw, cli.reads, cli.manifest, cli.output, cli.database)
