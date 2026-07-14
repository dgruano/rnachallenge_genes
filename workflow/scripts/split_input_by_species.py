#!/usr/bin/env python3
import re
import sys
from collections import defaultdict
from pathlib import Path

from Bio import SeqIO


def normalize_species(organism_str):
    """Normalize species name: lowercase, collapse non-alphanumeric to _, strip edges."""
    s = organism_str.strip().lower()
    s = re.sub(r"[^a-z0-9]+", "_", s)
    s = s.strip("_")
    return s


def main():
    resolved_path = snakemake.input.resolved
    input_fastas = snakemake.input.inputs
    output_dir = Path(snakemake.output[0])

    output_dir.mkdir(parents=True, exist_ok=True)

    # Build transcript_id → normalized_species mapping from resolved TSV
    transcript_to_species = {}
    with open(resolved_path) as f:
        f.readline()  # skip header
        for line in f:
            parts = line.rstrip("\n").split("\t")
            transcript_id = parts[0]
            organism = parts[4]  # column 5 (0-indexed: col 4)
            if transcript_id not in transcript_to_species:
                transcript_to_species[transcript_id] = normalize_species(organism)

    # Open file handles per species
    open_files = {}

    # Stream input FASTAs
    unresolved_count = 0
    species_counts = defaultdict(int)

    for fasta_file in input_fastas:
        for record in SeqIO.parse(fasta_file, "fasta"):
            record_id = record.id
            if record_id in transcript_to_species:
                species = transcript_to_species[record_id]
                species_counts[species] += 1

                if species not in open_files:
                    fasta_path = output_dir / f"{species}.fasta"
                    open_files[species] = open(fasta_path, "w")

                SeqIO.write(record, open_files[species], "fasta")
            else:
                unresolved_count += 1
                if "_unresolved" not in open_files:
                    open_files["_unresolved"] = open(
                        output_dir / "_unresolved.fasta", "w"
                    )
                SeqIO.write(record, open_files["_unresolved"], "fasta")

    # Close all files
    for f in open_files.values():
        f.close()

    # Log counts
    total = sum(species_counts.values()) + unresolved_count
    print(f"Total records: {total}")
    print(f"Resolved to species: {sum(species_counts.values())}")
    for species in sorted(species_counts.keys()):
        print(f"  {species}: {species_counts[species]}")
    print(f"Unresolved: {unresolved_count}")


if __name__ == "__main__":
    if "snakemake" not in dir():
        # Self-check: normalization function
        tests = [
            ("Arabidopsis thaliana", "arabidopsis_thaliana"),
            ("arabidopsis_thaliana", "arabidopsis_thaliana"),
            ("Homo sapiens", "homo_sapiens"),
            ("homo sapiens", "homo_sapiens"),
            ("Danio rerio", "danio_rerio"),
            ("danio rerio", "danio_rerio"),
            ("Drosophila yakuba (flies)", "drosophila_yakuba_flies"),
        ]
        for input_str, expected in tests:
            result = normalize_species(input_str)
            assert (
                result == expected
            ), f"normalize_species('{input_str}') = '{result}', expected '{expected}'"
        print("All normalization tests passed.")
    else:
        main()
