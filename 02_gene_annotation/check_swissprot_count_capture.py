import csv
import os
from collections import defaultdict

metadata_file = "salmon_sample_metadata.tsv"
mapping_dir = "validated_swissprot_mappings"

sample_report_file = "swissprot_count_capture_by_sample.tsv"
summary_report_file = "swissprot_count_capture_summary.tsv"


def read_mapping(species):
    """Return transcripts with a selected Swiss-Prot best hit."""
    mapping_file = os.path.join(
        mapping_dir,
        f"{species}_transcript_to_swissprot_besthit.tsv"
    )

    if not os.path.exists(mapping_file):
        raise FileNotFoundError(f"Missing mapping file: {mapping_file}")

    mapped_transcripts = set()

    with open(mapping_file, newline="") as infile:
        reader = csv.DictReader(infile, delimiter="\t")

        for row in reader:
            transcript_id = row["transcript_id"].strip()

            if transcript_id:
                mapped_transcripts.add(transcript_id)

    return mapped_transcripts


def quantify_capture(quant_file, mapped_transcripts):
    """Calculate total and Swiss-Prot-annotated Salmon NumReads."""
    total_counts = 0.0
    annotated_counts = 0.0
    expressed_transcripts = 0
    expressed_annotated_transcripts = 0

    with open(quant_file, newline="") as infile:
        reader = csv.DictReader(infile, delimiter="\t")

        for row in reader:
            transcript_id = row["Name"]
            count = float(row["NumReads"])

            total_counts += count

            if count > 0:
                expressed_transcripts += 1

            if transcript_id in mapped_transcripts:
                annotated_counts += count

                if count > 0:
                    expressed_annotated_transcripts += 1

    percent_counts_annotated = (
        100 * annotated_counts / total_counts if total_counts > 0 else 0
    )

    percent_expressed_transcripts_annotated = (
        100 * expressed_annotated_transcripts / expressed_transcripts
        if expressed_transcripts > 0 else 0
    )

    return {
        "total_numreads": total_counts,
        "annotated_numreads": annotated_counts,
        "percent_numreads_annotated": percent_counts_annotated,
        "expressed_transcripts": expressed_transcripts,
        "expressed_annotated_transcripts": expressed_annotated_transcripts,
        "percent_expressed_transcripts_annotated":
            percent_expressed_transcripts_annotated,
    }


# Read sample metadata
samples = []

with open(metadata_file, newline="") as infile:
    reader = csv.DictReader(infile, delimiter="\t")

    for row in reader:
        samples.append(row)

# Cache mappings per species
species_mappings = {}

for species in sorted({row["species"] for row in samples}):
    species_mappings[species] = read_mapping(species)

# Quantify annotation capture for every sample
sample_results = []
group_results = defaultdict(lambda: {"total": 0.0, "annotated": 0.0, "n": 0})

for row in samples:
    species = row["species"]
    sample = row["sample"]
    tissue = row["tissue"]
    sex = row["sex"]
    quant_file = row["quant_file"]

    result = quantify_capture(quant_file, species_mappings[species])

    sample_result = {
        "species": species,
        "sample": sample,
        "tissue": tissue,
        "sex": sex,
        "total_numreads": round(result["total_numreads"], 3),
        "annotated_numreads": round(result["annotated_numreads"], 3),
        "percent_numreads_annotated": round(
            result["percent_numreads_annotated"], 2
        ),
        "expressed_transcripts": result["expressed_transcripts"],
        "expressed_annotated_transcripts":
            result["expressed_annotated_transcripts"],
        "percent_expressed_transcripts_annotated": round(
            result["percent_expressed_transcripts_annotated"], 2
        ),
    }

    sample_results.append(sample_result)

    for group_name in [
        ("species", species),
        ("tissue", tissue),
        ("species_tissue", f"{species}:{tissue}"),
    ]:
        group_results[group_name]["total"] += result["total_numreads"]
        group_results[group_name]["annotated"] += result["annotated_numreads"]
        group_results[group_name]["n"] += 1

# Write per-sample report
with open(sample_report_file, "w", newline="") as outfile:
    fields = [
        "species",
        "sample",
        "tissue",
        "sex",
        "total_numreads",
        "annotated_numreads",
        "percent_numreads_annotated",
        "expressed_transcripts",
        "expressed_annotated_transcripts",
        "percent_expressed_transcripts_annotated",
    ]

    writer = csv.DictWriter(
        outfile,
        delimiter="\t",
        fieldnames=fields,
        lineterminator="\n",
    )

    writer.writeheader()
    writer.writerows(sample_results)

# Write grouped summary
with open(summary_report_file, "w", newline="") as outfile:
    writer = csv.writer(outfile, delimiter="\t", lineterminator="\n")

    writer.writerow([
        "group_type",
        "group",
        "number_of_samples",
        "total_numreads",
        "annotated_numreads",
        "percent_numreads_annotated",
    ])

    for group_key in sorted(group_results):
        group_type, group = group_key
        total = group_results[group_key]["total"]
        annotated = group_results[group_key]["annotated"]
        number_of_samples = group_results[group_key]["n"]

        percent = 100 * annotated / total if total > 0 else 0

        writer.writerow([
            group_type,
            group,
            number_of_samples,
            round(total, 3),
            round(annotated, 3),
            round(percent, 2),
        ])

print("Created:")
print(f"  {sample_report_file}")
print(f"  {summary_report_file}")

print("\nSwiss-Prot annotated NumReads captured by tissue:")
for group_key in sorted(group_results):
    group_type, group = group_key

    if group_type != "tissue":
        continue

    total = group_results[group_key]["total"]
    annotated = group_results[group_key]["annotated"]
    percent = 100 * annotated / total if total > 0 else 0

    print(
        group,
        group_results[group_key]["n"],
        round(percent, 2),
        sep="\t"
    )
