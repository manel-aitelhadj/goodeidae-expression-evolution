import csv
import re
from collections import Counter

input_file = "salmon_sample_inventory.tsv"
output_file = "salmon_sample_metadata.tsv"


def parse_sample(sample):
    """
    Parse tissue and sex from sample names.

    Supported examples:
      As_B_R2          female brain
      Xr_B_R1_F_RNA    female brain
      AtBR5            female brain

      As_Br_M1         male brain

      As_E_R2          embryo
      Gm_E_R8_2        embryo
      Xc_E__R10        embryo
      AtER5            embryo

      As_O_R2          ovary
      Gm_O_R1_RNA      ovary
      AtOR5            ovary

      As_Te_M1         testes

      As_T_R2          trophotaenia
      Xc_T__R10        trophotaenia
      GmTR22           trophotaenia
    """

    # Male brain: test before female brain to preserve B versus Br distinction
    if re.search(r"_Br_M\d+", sample):
        return "Br", "brain", "male"

    # Testes
    if re.search(r"_Te_M\d+", sample):
        return "Te", "testes", "male"

    # Female brain: standard or compact format
    if re.search(r"_B_R\d+", sample) or re.search(r"^[A-Za-z]+BR\d+", sample):
        return "B", "brain", "female"

    # Embryo: allows one or two underscores before R, plus compact format
    if re.search(r"_E_+R\d+", sample) or re.search(r"^[A-Za-z]+ER\d+", sample):
        return "E", "embryo", "undefined"

    # Ovary: allows suffixes such as _RNA, plus compact format
    if re.search(r"_O_R\d+", sample) or re.search(r"^[A-Za-z]+OR\d+", sample):
        return "O", "ovary", "female"

    # Trophotaenia: allows one or two underscores before R, plus compact format
    if re.search(r"_T_+R\d+", sample) or re.search(r"^[A-Za-z]+TR\d+", sample):
        return "T", "trophotaenia", "undefined"

    return "UNKNOWN", "UNKNOWN", "UNKNOWN"


rows = []
counts = Counter()
species_counts = Counter()
unknown = []

with open(input_file, newline="") as infile:
    reader = csv.DictReader(infile, delimiter="\t")

    for row in reader:
        code, tissue, sex = parse_sample(row["sample"])

        output_row = {
            "species": row["species"],
            "sample": row["sample"],
            "sample_code": code,
            "tissue": tissue,
            "sex": sex,
            "quant_file": row["quant_file"],
        }

        rows.append(output_row)
        species_counts[row["species"]] += 1
        counts[(row["species"], tissue, sex)] += 1

        if tissue == "UNKNOWN":
            unknown.append(output_row)

with open(output_file, "w", newline="") as outfile:
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
            "quant_file",
        ],
    )
    writer.writeheader()
    writer.writerows(rows)

print(f"Created: {output_file}")
print(f"Total samples: {len(rows)}")

print("\nSpecies detected:")
for species in sorted(species_counts):
    print(species, species_counts[species], sep="\t")

print("\nSamples still unrecognised:")
if len(unknown) == 0:
    print("None")
else:
    for row in unknown:
        print(row["species"], row["sample"], row["quant_file"], sep="\t")

print("\nSamples per species, tissue and sex:")
for key in sorted(counts):
    species, tissue, sex = key
    print(species, tissue, sex, counts[key], sep="\t")
