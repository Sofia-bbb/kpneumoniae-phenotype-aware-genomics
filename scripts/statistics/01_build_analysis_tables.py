#!/usr/bin/env python3

from __future__ import annotations

import csv
from pathlib import Path


# Repository root
BASE = Path(__file__).resolve().parents[2]

GENOMIC = BASE / "data" / "GPSKp_genomic_master.tsv"
PHENOTYPE = BASE / "data" / "GPSKp_genotype_phenotype_long.tsv"

OUTPUT = BASE / "results" / "statistics"
OUTPUT.mkdir(parents=True, exist_ok=True)

MASTER_OUT = OUTPUT / "01_Master_Analysis.tsv"
AMR_OUT = OUTPUT / "02_AMR_Table.tsv"
VIRULENCE_OUT = OUTPUT / "03_Virulence_Table.tsv"
SIDEROPHORE_OUT = OUTPUT / "04_Siderophore_Table.tsv"
CAPSULE_OUT = OUTPUT / "05_Capsule_Table.tsv"
ST_OUT = OUTPUT / "06_ST_Table.tsv"
GENO_PHENO_OUT = OUTPUT / "07_Genotype_Phenotype_Table.tsv"


def clean(value: object) -> str:
    if value is None:
        return ""

    text = str(value).strip()

    if text.lower() in {
        "",
        "-",
        "na",
        "n/a",
        "nan",
        "none",
        "null",
    }:
        return ""

    return text


def first_present(
    row: dict[str, str],
    candidates: list[str],
) -> str:
    for field in candidates:
        value = clean(row.get(field))
        if value:
            return value

    return ""


def locus_present(value: str) -> str:
    value = clean(value).lower()

    if value in {
        "",
        "0",
        "0.0",
        "-",
        "none",
        "na",
        "n/a",
        "nan",
        "null",
        "false",
        "negative",
    }:
        return "0"

    return "1"


def read_tsv(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    with path.open(
        "r",
        encoding="utf-8",
        errors="replace",
        newline="",
    ) as handle:
        reader = csv.DictReader(handle, delimiter="\t")

        fields = list(reader.fieldnames or [])
        rows = [
            {
                key: clean(value)
                for key, value in row.items()
            }
            for row in reader
        ]

    return fields, rows


def write_tsv(
    path: Path,
    fields: list[str],
    rows: list[dict[str, str]],
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)

    with path.open(
        "w",
        encoding="utf-8",
        newline="",
    ) as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=fields,
            delimiter="\t",
            extrasaction="ignore",
        )

        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    if not GENOMIC.is_file():
        raise FileNotFoundError(GENOMIC)

    if not PHENOTYPE.is_file():
        raise FileNotFoundError(PHENOTYPE)

    genomic_fields, genomic_rows = read_tsv(GENOMIC)
    phenotype_fields, phenotype_rows = read_tsv(PHENOTYPE)

    # One-row-per-isolate master analysis table
    write_tsv(
        MASTER_OUT,
        genomic_fields,
        genomic_rows,
    )

    amr_fields = [
        "SRR",
        "species",
        "ST",
        "Country",
        "Collection_Year",
        "Isolation_Source_Metadata",
        "resistance_score",
        "num_resistance_classes",
        "num_resistance_genes",
        "amr_gene_count",
        "amr_gene_list",
        "resistance_class_count",
        "resistance_class_list",
        "point_mutation_count",
        "point_mutation_list",
    ]

    amr_rows = [
        {
            field: clean(row.get(field))
            for field in amr_fields
        }
        for row in genomic_rows
    ]

    virulence_fields = [
        "SRR",
        "species",
        "ST",
        "Country",
        "Collection_Year",
        "Isolation_Source_Metadata",
        "virulence_score",
        "YbST",
        "yersiniabactin",
        "CbST",
        "colibactin",
        "AbST",
        "aerobactin",
        "SmST",
        "salmochelin",
        "RmST",
        "RmpADC",
        "rmpA2",
        "vf_hit_count_raw",
        "vf_factor_count_filtered",
        "vf_factor_list",
    ]

    virulence_rows = [
        {
            field: clean(row.get(field))
            for field in virulence_fields
        }
        for row in genomic_rows
    ]

    siderophore_fields = [
        "SRR",
        "ST",
        "Country",
        "Collection_Year",
        "Isolation_Source_Metadata",
        "virulence_score",
        "resistance_score",
        "yersiniabactin_present",
        "YbST",
        "yersiniabactin",
        "aerobactin_present",
        "AbST",
        "aerobactin",
        "salmochelin_present",
        "SmST",
        "salmochelin",
        "siderophore_count",
        "siderophore_profile",
    ]

    siderophore_rows: list[dict[str, str]] = []

    for row in genomic_rows:
        ybt = first_present(
            row,
            ["yersiniabactin", "YbST"],
        )
        iuc = first_present(
            row,
            ["aerobactin", "AbST"],
        )
        iro = first_present(
            row,
            ["salmochelin", "SmST"],
        )

        ybt_present = int(locus_present(ybt))
        iuc_present = int(locus_present(iuc))
        iro_present = int(locus_present(iro))

        profile_parts: list[str] = []

        if ybt_present:
            profile_parts.append("ybt")
        if iuc_present:
            profile_parts.append("iuc")
        if iro_present:
            profile_parts.append("iro")

        profile = (
            "+".join(profile_parts)
            if profile_parts
            else "none"
        )

        siderophore_rows.append(
            {
                "SRR": clean(row.get("SRR")),
                "ST": clean(row.get("ST")),
                "Country": clean(row.get("Country")),
                "Collection_Year": clean(
                    row.get("Collection_Year")
                ),
                "Isolation_Source_Metadata": clean(
                    row.get("Isolation_Source_Metadata")
                ),
                "virulence_score": clean(
                    row.get("virulence_score")
                ),
                "resistance_score": clean(
                    row.get("resistance_score")
                ),
                "yersiniabactin_present": str(ybt_present),
                "YbST": clean(row.get("YbST")),
                "yersiniabactin": clean(
                    row.get("yersiniabactin")
                ),
                "aerobactin_present": str(iuc_present),
                "AbST": clean(row.get("AbST")),
                "aerobactin": clean(row.get("aerobactin")),
                "salmochelin_present": str(iro_present),
                "SmST": clean(row.get("SmST")),
                "salmochelin": clean(
                    row.get("salmochelin")
                ),
                "siderophore_count": str(
                    ybt_present + iuc_present + iro_present
                ),
                "siderophore_profile": profile,
            }
        )

    capsule_fields = [
        "SRR",
        "ST",
        "Country",
        "Collection_Year",
        "Isolation_Source_Metadata",
        "K_locus",
        "K_type",
        "O_locus",
        "O_type",
        "virulence_score",
        "resistance_score",
    ]

    capsule_rows = [
        {
            field: clean(row.get(field))
            for field in capsule_fields
        }
        for row in genomic_rows
    ]

    st_fields = [
        "SRR",
        "species",
        "ST",
        "Country",
        "Collection_Year",
        "Isolation_Source_Metadata",
        "K_locus",
        "O_locus",
        "virulence_score",
        "resistance_score",
        "amr_gene_count",
        "vf_factor_count_filtered",
        "siderophore_count",
        "siderophore_profile",
    ]

    siderophore_by_srr = {
        row["SRR"]: row
        for row in siderophore_rows
    }

    st_rows: list[dict[str, str]] = []

    for row in genomic_rows:
        srr = clean(row.get("SRR"))
        side = siderophore_by_srr.get(srr, {})

        st_rows.append(
            {
                "SRR": srr,
                "species": clean(row.get("species")),
                "ST": clean(row.get("ST")),
                "Country": clean(row.get("Country")),
                "Collection_Year": clean(
                    row.get("Collection_Year")
                ),
                "Isolation_Source_Metadata": clean(
                    row.get("Isolation_Source_Metadata")
                ),
                "K_locus": clean(row.get("K_locus")),
                "O_locus": clean(row.get("O_locus")),
                "virulence_score": clean(
                    row.get("virulence_score")
                ),
                "resistance_score": clean(
                    row.get("resistance_score")
                ),
                "amr_gene_count": clean(
                    row.get("amr_gene_count")
                ),
                "vf_factor_count_filtered": clean(
                    row.get("vf_factor_count_filtered")
                ),
                "siderophore_count": clean(
                    side.get("siderophore_count")
                ),
                "siderophore_profile": clean(
                    side.get("siderophore_profile")
                ),
            }
        )

    geno_pheno_fields = [
        field
        for field in phenotype_fields
        if field in {
            "SRR",
            "Country",
            "Collection_Year",
            "Isolation_Source_Metadata",
            "ST",
            "K_locus",
            "O_locus",
            "virulence_score",
            "resistance_score",
            "amr_gene_count",
            "amr_gene_list",
            "resistance_class_list",
            "point_mutation_list",
            "yersiniabactin",
            "aerobactin",
            "salmochelin",
            "Antibiotic_Original",
            "Resistance_Phenotype_Original",
            "Measurement_Sign",
            "MIC_Value",
            "Antibiotic",
            "Resistance_Phenotype",
            "MIC_Units",
            "MIC_Reported",
            "Testing_Platform",
            "Method_or_Reagent",
            "Testing_Standard",
        }
    ]

    geno_pheno_rows = [
        {
            field: clean(row.get(field))
            for field in geno_pheno_fields
        }
        for row in phenotype_rows
    ]

    write_tsv(AMR_OUT, amr_fields, amr_rows)
    write_tsv(
        VIRULENCE_OUT,
        virulence_fields,
        virulence_rows,
    )
    write_tsv(
        SIDEROPHORE_OUT,
        siderophore_fields,
        siderophore_rows,
    )
    write_tsv(
        CAPSULE_OUT,
        capsule_fields,
        capsule_rows,
    )
    write_tsv(ST_OUT, st_fields, st_rows)
    write_tsv(
        GENO_PHENO_OUT,
        geno_pheno_fields,
        geno_pheno_rows,
    )

    print("=" * 72)
    print("GPS-Kp Module 7 Analysis Tables")
    print("=" * 72)
    print(f"Master isolates: {len(genomic_rows)}")
    print(f"AMR rows: {len(amr_rows)}")
    print(f"Virulence rows: {len(virulence_rows)}")
    print(f"Siderophore rows: {len(siderophore_rows)}")
    print(f"Capsule rows: {len(capsule_rows)}")
    print(f"ST rows: {len(st_rows)}")
    print(f"Genotype-phenotype rows: {len(geno_pheno_rows)}")
    print(f"Output directory: {OUTPUT}")
    print("=" * 72)


if __name__ == "__main__":
    main()

