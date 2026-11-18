import csv
import gzip
import os
from collections import defaultdict

species_list = ["as", "at", "cb", "cl", "ga", "gm", "if", "xc", "xr"]
metadata_file = "salmon_sample_metadata.tsv"
report_file = "annotation_mapping_report.tsv"
mapping_dir = "validated_swissprot_mappings"

os.makedirs(mapping_dir, exist_ok=True)


def open_text(path):
    """Open either a plain text file or a gzip-compressed text file."""
    if path.endswith(".gz"):
        return gzip.open(path, "rt", newline="")
    return open(path, "r", newline="")


def remove_orf_suffix(query_id):
    """
    Convert protein/ORF identifier to Salmon transcript identifier:
    NODE_..._i6.p1 -> NODE_..._i6
    """
    if ".p" in query_id:
        return query_id.rsplit(".p", 1)[0]
    return query_id


def get_first_quant_file(species):
    with open(metadata_file, newline="") as infile:
        reader = csv.DictReader(infile, delimiter="\t")
        for row in reader:
            if row["species"] == species:
                return row["quant_file"]
    raise ValueError(f"No Salmon quant file found for species: {species}")


def get_annotation_file(species):
    gz_file = os.path.join(species, f"{species}_anno.tsv.gz")
    plain_file = os.path.join(species, f"{species}_anno.tsv")

    if os.path.exists(gz_file):
        return gz_file
    if os.path.exists(plain_file):
        return plain_file

    raise FileNotFoundError(f"No annotation file found for species: {species}")


def read_salmon_transcripts(quant_file):
    transcript_ids = set()

    with open(quant_file, newline="") as infile:
        reader = csv.DictReader(infile, delimiter="\t")

        for row in reader:
            transcript_ids.add(row["Name"])

    return transcript_ids


def read_annotation_hits(annotation_file):
    """
    Read Swiss-Prot annotation hits and retain information per transcript.
    Multiple .p ORFs from the same Salmon transcript are merged at transcript level.
    """

    hits_by_transcript = defaultdict(list)

    with open_text(annotation_file) as infile:
        reader = csv.DictReader(infile, delimiter="\t")

        required = {"queryID", "targetID", "e-value", "bitscore"}
        missing = required.difference(reader.fieldnames or [])

        if missing:
            raise ValueError(
                f"{annotation_file} is missing required columns: {sorted(missing)}"
            )

        for row in reader:
            query_id = row["queryID"].strip()
            target_id = row["targetID"].strip()

            if not query_id or not target_id:
                continue

            transcript_id = remove_orf_suffix(query_id)

            try:
                bitscore = float(row["bitscore"])
            except (ValueError, TypeError):
                continue

            try:
                evalue = float(row["e-value"])
            except (ValueError, TypeError):
                evalue = float("inf")

            hits_by_transcript[transcript_id].append(
                {
                    "target_id": target_id,
                    "bitscore": bitscore,
                    "evalue": evalue,
                    "query_id": query_id,
                }
            )

    return hits_by_transcript


def choose_best_hit(hits):
    """
    Pick one Swiss-Prot hit per Salmon transcript:
    highest bitscore first, then lowest e-value.
    """
    ranked = sorted(
        hits,
        key=lambda x: (-x["bitscore"], x["evalue"], x["target_id"])
    )
    return ranked[0]


report_rows = []

for species in species_list:
    quant_file = get_first_quant_file(species)
    annotation_file = get_annotation_file(species)

    salmon_transcripts = read_salmon_transcripts(quant_file)
    annotation_hits = read_annotation_hits(annotation_file)

    annotated_transcripts = set(annotation_hits.keys())
    overlapping_transcripts = salmon_transcripts.intersection(annotated_transcripts)

    unique_mapping = {}
    ambiguous_target_count = 0
    multiple_orf_count = 0

    for transcript_id in sorted(overlapping_transcripts):
        hits = annotation_hits[transcript_id]

        targets = {hit["target_id"] for hit in hits}
        orfs = {hit["query_id"] for hit in hits}

        if len(targets) > 1:
            ambiguous_target_count += 1

        if len(orfs) > 1:
            multiple_orf_count += 1

        best_hit = choose_best_hit(hits)
        unique_mapping[transcript_id] = best_hit

    mapping_output = os.path.join(
        mapping_dir,
        f"{species}_transcript_to_swissprot_besthit.tsv"
    )

    with open(mapping_output, "w", newline="") as outfile:
        writer = csv.writer(outfile, delimiter="\t", lineterminator="\n")
        writer.writerow(
            [
                "transcript_id",
                "swissprot_id",
                "bitscore",
                "evalue",
                "selected_query_orf",
            ]
        )

        for transcript_id in sorted(unique_mapping):
            hit = unique_mapping[transcript_id]
            writer.writerow(
                [
                    transcript_id,
                    hit["target_id"],
                    hit["bitscore"],
                    hit["evalue"],
                    hit["query_id"],
                ]
            )

    total_salmon = len(salmon_transcripts)
    mapped = len(overlapping_transcripts)
    percent_mapped = 100 * mapped / total_salmon if total_salmon else 0

    report_rows.append(
        {
            "species": species,
            "salmon_transcripts": total_salmon,
            "annotated_transcripts_in_table": len(annotated_transcripts),
            "mapped_salmon_transcripts": mapped,
            "percent_mapped": round(percent_mapped, 2),
            "transcripts_with_multiple_targets": ambiguous_target_count,
            "transcripts_with_multiple_orfs": multiple_orf_count,
            "annotation_file": annotation_file,
            "mapping_file": mapping_output,
        }
    )


with open(report_file, "w", newline="") as outfile:
    fieldnames = [
        "species",
        "salmon_transcripts",
        "annotated_transcripts_in_table",
        "mapped_salmon_transcripts",
        "percent_mapped",
        "transcripts_with_multiple_targets",
        "transcripts_with_multiple_orfs",
        "annotation_file",
        "mapping_file",
    ]

    writer = csv.DictWriter(
        outfile,
        delimiter="\t",
        fieldnames=fieldnames,
        lineterminator="\n",
    )

    writer.writeheader()
    writer.writerows(report_rows)


print("Created:", report_file)
print()
print(
    "species\tSalmon_transcripts\tmapped\tpercent_mapped\t"
    "multiple_targets\tmultiple_ORFs"
)

for row in report_rows:
    print(
        row["species"],
        row["salmon_transcripts"],
        row["mapped_salmon_transcripts"],
        row["percent_mapped"],
        row["transcripts_with_multiple_targets"],
        row["transcripts_with_multiple_orfs"],
        sep="\t",
    )

print()
print("Candidate best-hit mappings written to:", mapping_dir)
print("Do not use these for final matrices until the mapping report is reviewed.")
