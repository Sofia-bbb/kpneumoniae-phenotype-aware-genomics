#!/usr/bin/env python3

from pathlib import Path
import csv
import re

BASE = Path(__file__).resolve().parents[2]

MASTER = BASE / "data" / "GPSKp_genomic_master.tsv"
RESULTS = BASE / "results" / "vfdb"
FASTA = BASE / "data" / "vfdb" / "VFDB_setA_nt.fas"

OUT = BASE / "results" / "vfdb_coverage_FINAL1168_90id_80cov.tsv"

IDENTITY_CUTOFF = 90.0
COVERAGE_CUTOFF = 80.0
EVALUE_CUTOFF = 1e-10


def load_vfdb_reference():
    lengths = {}
    names = {}

    current = None
    current_length = 0

    with open(
        FASTA,
        "r",
        encoding="utf-8",
        errors="replace",
    ) as handle:

        for line in handle:
            line = line.rstrip()

            if line.startswith(">"):

                if current is not None:
                    lengths[current] = current_length

                header = line[1:].strip()
                current = header.split()[0]
                current_length = 0

                match = re.search(r"\[([^\]]+)\]", header)

                if match:
                    names[current] = match.group(1)
                else:
                    names[current] = current

            else:
                current_length += len(line.strip())

        if current is not None:
            lengths[current] = current_length

    return lengths, names


ref_lengths, vf_names = load_vfdb_reference()


with open(
    MASTER,
    "r",
    encoding="utf-8",
    errors="replace",
) as handle:
    reader = csv.DictReader(handle, delimiter="\t")
    final_ids = [
        row["SRR"].strip()
        for row in reader
        if row.get("SRR", "").strip()
    ]


rows = []

missing_files = 0
missing_reference_lengths = 0


for srr in final_ids:

    infile = RESULTS / srr / f"{srr}_vfdb_raw.tsv"

    if not infile.exists():
        missing_files += 1
        continue

    original_factors = set()
    coverage_factors = set()

    with open(
        infile,
        "r",
        encoding="utf-8",
        errors="replace",
    ) as blast:

        for line in blast:

            if not line.strip():
                continue

            cols = line.rstrip().split("\t")

            if len(cols) < 12:
                continue

            subject = cols[1]
            identity = float(cols[2])
            alignment_length = int(cols[3])
            evalue = float(cols[10])

            if (
                identity < IDENTITY_CUTOFF
                or evalue > EVALUE_CUTOFF
            ):
                continue

            factor = vf_names.get(subject, subject)
            original_factors.add(factor)

            ref_len = ref_lengths.get(subject)

            if not ref_len:
                missing_reference_lengths += 1
                continue

            subject_coverage = (
                alignment_length / ref_len
            ) * 100.0

            if subject_coverage >= COVERAGE_CUTOFF:
                coverage_factors.add(factor)

    removed_factors = (
        original_factors
        - coverage_factors
    )

    rows.append(
        {
            "SRR": srr,

            "vf_factor_count_original":
                len(original_factors),

            "vf_factor_list_original":
                ";".join(
                    sorted(original_factors)
                ),

            "vf_factor_count_80cov":
                len(coverage_factors),

            "vf_factor_list_80cov":
                ";".join(
                    sorted(coverage_factors)
                ),

            "removed_factor_count":
                len(removed_factors),

            "removed_factor_list":
                ";".join(
                    sorted(removed_factors)
                ),
        }
    )


with open(
    OUT,
    "w",
    newline="",
    encoding="utf-8",
) as out:

    writer = csv.DictWriter(
        out,
        fieldnames=[
            "SRR",
            "vf_factor_count_original",
            "vf_factor_list_original",
            "vf_factor_count_80cov",
            "vf_factor_list_80cov",
            "removed_factor_count",
            "removed_factor_list",
        ],
        delimiter="\t",
    )

    writer.writeheader()
    writer.writerows(rows)


print("=" * 60)
print("VFDB COVERAGE FILTER - FINAL 1168")
print("=" * 60)
print("Final IDs loaded :", len(final_ids))
print("Samples analyzed :", len(rows))
print("Missing raw files:", missing_files)
print("Missing ref len  :", missing_reference_lengths)
print("Identity cutoff  :", IDENTITY_CUTOFF)
print("Coverage cutoff  :", COVERAGE_CUTOFF)
print("E-value cutoff   :", EVALUE_CUTOFF)
print("Output           :", OUT)
print("=" * 60)


