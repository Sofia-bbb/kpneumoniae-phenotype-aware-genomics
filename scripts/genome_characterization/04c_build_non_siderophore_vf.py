#!/usr/bin/env python3

from pathlib import Path
import csv

BASE = Path(__file__).resolve().parents[2]

INPUT = BASE / "results" / "vfdb_coverage_FINAL1168_90id_80cov.tsv"

OUTPUT = BASE / "results" / "vfdb_non_siderophore_FINAL1168_90id_80cov.tsv"


def is_siderophore_factor(factor):
    """
    Exclude VFDB factors corresponding to the siderophore systems
    used to define the Kleborate siderophore profiles:
    yersiniabactin (ybt), aerobactin (iuc), and salmochelin (iro).

    Matching is based on the VFDB factor names retained in the
    filtered factor list.
    """
    text = factor.lower().strip()

    if "yersiniabactin" in text:
        return True

    if text.startswith("ybt "):
        return True

    if "aerobactin" in text:
        return True

    if "salmochelin" in text:
        return True

    # VFDB Set A uses "Sal (VF0563)" for the salmochelin system.
    if text.startswith("sal (vf0563)"):
        return True

    return False


rows = []

total_factors_before = 0
total_factors_after = 0
total_siderophore_removed = 0

with open(
    INPUT,
    "r",
    newline="",
    encoding="utf-8",
    errors="replace",
) as handle:

    reader = csv.DictReader(handle, delimiter="\t")

    for row in reader:

        srr = row["SRR"]

        factor_list = (
            row.get("vf_factor_list_80cov") or ""
        ).strip()

        if factor_list:
            factors = {
                x.strip()
                for x in factor_list.split(";")
                if x.strip()
            }
        else:
            factors = set()

        siderophore_factors = {
            x for x in factors
            if is_siderophore_factor(x)
        }

        non_siderophore_factors = (
            factors - siderophore_factors
        )

        total_factors_before += len(factors)
        total_factors_after += len(
            non_siderophore_factors
        )
        total_siderophore_removed += len(
            siderophore_factors
        )

        rows.append(
            {
                "SRR": srr,

                "vf_factor_count_80cov":
                    len(factors),

                "vf_factor_list_80cov":
                    ";".join(sorted(factors)),

                "siderophore_vf_count_removed":
                    len(siderophore_factors),

                "siderophore_vf_list_removed":
                    ";".join(
                        sorted(siderophore_factors)
                    ),

                "non_siderophore_vf_count":
                    len(non_siderophore_factors),

                "non_siderophore_vf_list":
                    ";".join(
                        sorted(non_siderophore_factors)
                    ),
            }
        )


with open(
    OUTPUT,
    "w",
    newline="",
    encoding="utf-8",
) as out:

    writer = csv.DictWriter(
        out,
        fieldnames=[
            "SRR",
            "vf_factor_count_80cov",
            "vf_factor_list_80cov",
            "siderophore_vf_count_removed",
            "siderophore_vf_list_removed",
            "non_siderophore_vf_count",
            "non_siderophore_vf_list",
        ],
        delimiter="\t",
    )

    writer.writeheader()
    writer.writerows(rows)


print("=" * 65)
print("NON-SIDEROPHORE VFDB BURDEN - FINAL 1168")
print("=" * 65)
print("Samples processed             :", len(rows))
print("VF factors before exclusion   :", total_factors_before)
print("Siderophore factors removed   :", total_siderophore_removed)
print("Non-siderophore factors remain:", total_factors_after)
print("Output                        :", OUTPUT)
print("=" * 65)



