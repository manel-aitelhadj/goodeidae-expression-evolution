import csv
import os
from collections import defaultdict

metadata_file = "salmon_sample_metadata.tsv"
mapping_dir = "validated_swissprot_mappings"
output_dir = "tissue_swissprot_count_matrices"

os.makedirs(output_dir, exist_ok=True)


def read_metadata():
    samples_by_tissue = defaultdict(list)

    with open(metadata_file, newline="") as infile:
        reader = csv.DictReader(infile, delimiter="\t")

        for row in reader:
            if row["tissue"] == "UNKNOWN":
                raise ValueError(f"Unclassified sample found: {row['sample']}")

            samples_by_tissue[row["tissue"]].append(row)

    return samples_by_tissue


def read_mapping(species):
    """
    Read the previously validated best-hit mapping:
    Salmon transcript ID -> Swiss-Prot accession ID.
    """
    mapping_file = os.path.join(
        mapping_dir,
        f"{species}_transcript_to_swissprot_besthit.tsv"
    )

    if not os.path.exists(mapping_file):
        raise FileNotFoundError(f"Missing mapping file: {mapping_file}")

    mapping = {}

    with open(mapping_file, newline="") as infile:
        reader = csv.DictReader(infile, delimiter="\t")

        for row in reader:
            transcript_id = row["transcript_id"].strip()
            swissprot_id = row["swissprot_id"].strip()

            if transcript_id and swissprot_id:
                mapping[transcript_id] = swissprot_id

    return mapping


def aggregate_sample_to_swissprot(quant_file, mapping):
    """
    Sum Salmon NumReads for all transcripts assigned to the same
    Swiss-Prot accession.
    """
    counts = defaultdict(float)
    total_numreads = 0.0
    annotated_numreads = 0.0

    with open(quant_file, newline="") as infile:
        reader = csv.DictReader(infile, delimiter="\t")

        required = {"Name", "NumReads"}
        if not required.issubset(reader.fieldnames or []):
            raise ValueError(f"Unexpected Salmon file format: {quant_file}")

        for row in reader:
            transcript_id = row["Name"]
            numreads = float(row["NumReads"])

            total_numreads += numreads

            if transcript_id in mapping:
                swissprot_id = mapping[transcript_id]
                counts[swissprot_id] += numreads
                annotated_numreads += numreads

    return counts, total_numreads, annotated_numreads


samples_by_tissue = read_metadata()

# Cache mappings once per species
all_species = sorted(
    {
        sample["species"]
        for tissue_samples in samples_by_tissue.values()
        for sample in tissue_samples
    }
)

mapping_by_species = {
    species: read_mapping(species)
    for species in all_species
}

overall_summary = []

for tissue in sorted(samples_by_tissue):
    samples = sorted(
        samples_by_tissue[tissue],
        key=lambda row: (row["species"], row["sample"])
    )

    species_in_tissue = sorted({sample["species"] for sample in samples})

    counts_by_sample = {}
    all_swissprot_ids = set()

    total_tissue_numreads = 0.0
    total_tissue_annotated_numreads = 0.0

    # Aggregate each Salmon sample to Swiss-Prot accessions
    for sample in samples:
        species = sample["species"]
        sample_name = sample["sample"]
        quant_file = sample["quant_file"]

        counts, total_reads, annotated_reads = aggregate_sample_to_swissprot(
            quant_file,
            mapping_by_species[species]
        )

        counts_by_sample[sample_name] = counts
        all_swissprot_ids.update(counts.keys())

        total_tissue_numreads += total_reads
        total_tissue_annotated_numreads += annotated_reads

    all_swissprot_ids = sorted(all_swissprot_ids)

    # Identify which Swiss-Prot accessions are annotated in each species
    species_targets = {
        species: set(mapping_by_species[species].values())
        for species in species_in_tissue
    }

    feature_presence = {}

    for swissprot_id in all_swissprot_ids:
        present_species = [
            species
            for species in species_in_tissue
            if swissprot_id in species_targets[species]
        ]
        feature_presence[swissprot_id] = present_species

    shared_2plus = [
        feature
        for feature in all_swissprot_ids
        if len(feature_presence[feature]) >= 2
    ]

    shared_all_species = [
        feature
        for feature in all_swissprot_ids
        if len(feature_presence[feature]) == len(species_in_tissue)
    ]

    # Write complete annotated matrix: union of Swiss-Prot IDs
    complete_matrix_file = os.path.join(
        output_dir,
        f"{tissue}_swissprot_counts_all_annotated.csv"
    )

    with open(complete_matrix_file, "w", newline="") as outfile:
        writer = csv.writer(outfile, lineterminator="\n")
        writer.writerow(
            ["swissprot_id"] + [sample["sample"] for sample in samples]
        )

        for swissprot_id in all_swissprot_ids:
            writer.writerow(
                [swissprot_id] +
                [
                    round(counts_by_sample[sample["sample"]].get(swissprot_id, 0.0), 6)
                    for sample in samples
                ]
            )

    # Write subset containing IDs represented in at least two species
    shared_matrix_file = os.path.join(
        output_dir,
        f"{tissue}_swissprot_counts_shared_2plus_species.csv"
    )

    with open(shared_matrix_file, "w", newline="") as outfile:
        writer = csv.writer(outfile, lineterminator="\n")
        writer.writerow(
            ["swissprot_id"] + [sample["sample"] for sample in samples]
        )

        for swissprot_id in shared_2plus:
            writer.writerow(
                [swissprot_id] +
                [
                    round(counts_by_sample[sample["sample"]].get(swissprot_id, 0.0), 6)
                    for sample in samples
                ]
            )

    # Write feature presence report
    presence_file = os.path.join(
        output_dir,
        f"{tissue}_swissprot_feature_presence.tsv"
    )

    with open(presence_file, "w", newline="") as outfile:
        writer = csv.writer(outfile, delimiter="\t", lineterminator="\n")

        writer.writerow([
            "swissprot_id",
            "number_of_species_with_annotation",
            "species_with_annotation"
        ])

        for swissprot_id in all_swissprot_ids:
            present_species = feature_presence[swissprot_id]

            writer.writerow([
                swissprot_id,
                len(present_species),
                ",".join(present_species)
            ])

    # Write corresponding sample metadata
    tissue_metadata_file = os.path.join(
        output_dir,
        f"{tissue}_sample_metadata.tsv"
    )

    with open(tissue_metadata_file, "w", newline="") as outfile:
        writer = csv.DictWriter(
            outfile,
            delimiter="\t",
            lineterminator="\n",
            fieldnames=[
                "species",
                "sample",
                "sample_code",
                "tissue",
                "sex",
                "quant_file"
            ]
        )

        writer.writeheader()
        writer.writerows(samples)

    percentage_captured = (
        100 * total_tissue_annotated_numreads / total_tissue_numreads
        if total_tissue_numreads > 0
        else 0
    )

    summary_row = {
        "tissue": tissue,
        "number_of_samples": len(samples),
        "species_present": ",".join(species_in_tissue),
        "number_of_species": len(species_in_tissue),
        "all_annotated_features": len(all_swissprot_ids),
        "features_in_2plus_species": len(shared_2plus),
        "features_in_all_species": len(shared_all_species),
        "percent_numreads_annotated": round(percentage_captured, 2)
    }

    overall_summary.append(summary_row)

    print(
        tissue,
        f"samples={len(samples)}",
        f"species={len(species_in_tissue)}",
        f"features={len(all_swissprot_ids)}",
        f"shared_2plus={len(shared_2plus)}",
        f"shared_all={len(shared_all_species)}",
        f"annotated_reads={round(percentage_captured, 2)}%",
        sep="\t"
    )


summary_file = os.path.join(
    output_dir,
    "tissue_swissprot_matrix_summary.tsv"
)

with open(summary_file, "w", newline="") as outfile:
    writer = csv.DictWriter(
        outfile,
        delimiter="\t",
        lineterminator="\n",
        fieldnames=[
            "tissue",
            "number_of_samples",
            "species_present",
            "number_of_species",
            "all_annotated_features",
            "features_in_2plus_species",
            "features_in_all_species",
            "percent_numreads_annotated"
        ]
    )

    writer.writeheader()
    writer.writerows(overall_summary)

print()
print("Created matrix summary:")
print(summary_file)
