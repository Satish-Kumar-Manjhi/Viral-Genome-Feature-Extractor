#!/usr/bin/env python3

"""
Genome Feature Extractor
------------------------
Parses NCBI GenBank (.gb/.gbk) viral genome files and extracts:

- Accession number
- Organism name
- 5' UTR sequence
- Complete CDS
- Individual mature peptides (mat_peptide)
- 3' UTR
- Translated proteome sequence

Author: Satish Kumar Manjhi
"""

import csv
import sys


def read_genbank_file(filepath):
    with open(filepath) as f:
        return f.readlines()


def extract_accession(lines):
    for line in lines:
        if line.startswith("ACCESSION"):
            return line.split()[1]
    return ""


def extract_source(lines):
    for line in lines:
        if line.startswith("SOURCE"):
            return line.split("SOURCE")[1].strip()
    return ""


def extract_sequence(lines):
    sequence = ""
    origin_index = 0

    for i, line in enumerate(lines):
        if line.strip() == "ORIGIN":
            origin_index = i + 1
            break

    for line in lines[origin_index:]:
        if line.startswith("//"):
            break
        parts = line.strip().split()
        sequence += "".join(parts[1:])

    return sequence


def extract_feature_sequence(sequence, location):
    start, end = map(int, location.split(".."))
    return sequence[start - 1:end]


def main():
    if len(sys.argv) != 3:
        print("Usage: python genome_parser.py input.gb output.csv")
        sys.exit(1)

    input_file = sys.argv[1]
    output_file = sys.argv[2]

    lines = read_genbank_file(input_file)

    accession = extract_accession(lines)
    source = extract_source(lines)
    full_sequence = extract_sequence(lines)

    five_utr = ""
    three_utr = ""
    cds_seq = ""
    cds_translation = ""

    peptide_names = []
    peptide_sequences = []
    peptide_proteins = []

    cds_translation_start = 0
    cds_protein_pointer = 0

    for i, line in enumerate(lines):

        if "     5'UTR" in line:
            location = line.split()[-1]
            five_utr = extract_feature_sequence(full_sequence, location)

        if "     3'UTR" in line:
            location = line.split()[-1]
            three_utr = extract_feature_sequence(full_sequence, location)

        if line.strip().startswith("CDS"):
            location = line.strip().split()[1]
            cds_seq = extract_feature_sequence(full_sequence, location)

        if "/translation=" in line:
            cds_translation_start = i
            break

    # Extract CDS translation
    for line in lines[cds_translation_start:]:
        if line.strip().startswith("/"):
            break
        cds_translation += line.strip()

    cds_translation = cds_translation.replace('/translation="', "").replace('"', "")

    # Extract mat_peptides
    for i, line in enumerate(lines):
        if "mat_peptide" in line:
            location = line.strip().split()[1]
            pep_seq = extract_feature_sequence(full_sequence, location)
            peptide_sequences.append(pep_seq)

            pep_len = int(len(pep_seq) / 3)
            peptide_proteins.append(
                cds_translation[cds_protein_pointer:cds_protein_pointer + pep_len]
            )
            cds_protein_pointer += pep_len

        if "/product=" in line:
            name = line.split("=")[1].strip().replace('"', "")
            peptide_names.append(name)

    # Write output
    with open(output_file, "w", newline="") as f:
        writer = csv.writer(f)

        header = ["Organism", "Accession", "5'UTR", "CDS"]
        header.extend(peptide_names)
        header.append("3'UTR")

        writer.writerow(header)

        row = [source, accession, five_utr, cds_seq]
        row.extend(peptide_sequences)
        row.append(three_utr)

        writer.writerow(row)

        protein_row = ["", "Proteome", "", cds_translation]
        protein_row.extend(peptide_proteins)

        writer.writerow(protein_row)


if __name__ == "__main__":
    main()

