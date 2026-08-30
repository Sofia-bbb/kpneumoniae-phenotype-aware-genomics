#!/usr/bin/env python3

import csv
import math
import re
from collections import defaultdict
from pathlib import Path

# Repository root
BASE = Path(__file__).resolve().parents[2]

# Input data
INPUT = BASE / "data" / "Supplementary_Table_S7_Integrated_Genotype_Phenotype_FINAL.tsv"

# Output directory
OUT = BASE / "results" / "genotype_phenotype"
OUT.mkdir(parents=True, exist_ok=True)


# ============================================================
# GENOTYPE PREDICTION RULES
# ============================================================

CARBAPENEMASE_PATTERNS = [
    r"\bblaKPC",
    r"\bblaNDM",
    r"\bblaVIM",
    r"\bblaIMP",
    r"\bblaOXA-48\b",
    r"\bblaOXA-181\b",
    r"\bblaOXA-232\b",
    r"\bblaGES-5\b",
    r"\bblaNMC",
    r"\bblaIMI",
    r"\bblaFRI",
]

FLUOROQUINOLONE_GENE_PATTERNS = [
    r"\bqnrA",
    r"\bqnrB",
    r"\bqnrS",
    r"\bqepA",
    r"aac\(6['’]\)-Ib-cr",
]

FLUOROQUINOLONE_MUTATIONS = {
    "parC_S80I",
    "parC_E84K",
    "gyrA_S83I",
    "gyrA_S83F",
    "gyrA_S83Y",
    "gyrA_S83L",
    "gyrA_D87N",
    "gyrA_D87A",
    "gyrA_D87G",
    "gyrA_D87Y",
}

COLISTIN_GENE_PATTERNS = [
    r"\bmcr-",
    r"\bmcr\d",
]

COLISTIN_MUTATIONS = {
    "pmrB_T157P",
    "pmrB_T140P",
    "pmrB_S85R",
    "pmrA_S42N",
    "pmrA_G53S",
    "pmrA_G53C",
    "phoQ_V24G",
    "phoP_D191Y",
    "mgrB_Q30STOP",
    "mgrB_W47R",
    "mgrB_D31N",
    "crrB_L94M",
    "crrB_G183V",
    "crrB_S195N",
    "crrB_P151L",
}

PORIN_PATTERNS = [
    r"ompK35_",
    r"ompK36_",
]

ESBL_AMP_C_PATTERNS = [
    r"\bblaCTX-M",
    r"\bblaSHV",
    r"\bblaTEM",
    r"\bblaCMY",
    r"\bblaDHA",
    r"\bblaFOX",
    r"\bblaACT",
    r"\bblaMIR",
]


def clean(x):
    if x is None:
        return ""
    x = str(x).strip()
    if x.lower() in {"", "na", "n/a", "nan", "none", "null", "not defined", "-"}:
        return ""
    return x


def split_items(x):
    x = clean(x)
    if not x:
        return []
    return [i.strip() for i in x.split(";") if i.strip()]


def regex_hits(items, patterns):
    hits = []
    text = ";".join(items)

    for pattern in patterns:
        if re.search(pattern, text, flags=re.IGNORECASE):
            hits.append(pattern)

    return hits


def exact_hits(items, allowed):
    return [x for x in items if x in allowed]


def normalize_antibiotic(x):
    x = clean(x).lower()

    mapping = {
        "ciprofloxacin": "Ciprofloxacin",
        "colistin": "Colistin",
        "imipenem": "Imipenem",
        "meropenem": "Meropenem",
    }

    return mapping.get(x, "")


def normalize_phenotype(x):
    x = clean(x).upper()

    if x in {"R", "RESISTANT"}:
        return "R"

    if x in {"S", "SUSCEPTIBLE"}:
        return "S"

    if x in {
        "I",
        "INTERMEDIATE",
        "SDD",
        "SUSCEPTIBLE-DOSE DEPENDENT",
    }:
        return "I"

    return ""


def siderophore_profile(row):
    profile = []

    if clean(row.get("yersiniabactin")):
        profile.append("ybt")

    if clean(row.get("aerobactin")):
        profile.append("iuc")

    if clean(row.get("salmochelin")):
        profile.append("iro")

    return "+".join(profile) if profile else "none"


def genotype_prediction(row, antibiotic):
    genes = split_items(row.get("amr_gene_list"))
    muts = split_items(row.get("point_mutation_list"))

    hits = []

    if antibiotic in {"Imipenem", "Meropenem"}:
        hits = regex_hits(genes, CARBAPENEMASE_PATTERNS)

    elif antibiotic == "Ciprofloxacin":
        hits = regex_hits(
            genes,
            FLUOROQUINOLONE_GENE_PATTERNS,
        )
        hits += exact_hits(
            muts,
            FLUOROQUINOLONE_MUTATIONS,
        )

    elif antibiotic == "Colistin":
        hits = regex_hits(
            genes,
            COLISTIN_GENE_PATTERNS,
        )
        hits += exact_hits(
            muts,
            COLISTIN_MUTATIONS,
        )

    prediction = "R" if hits else "S"

    return prediction, ";".join(sorted(set(hits)))


# ============================================================
# LOAD FINAL S7
# ============================================================

rows = []

with INPUT.open() as handle:
    reader = csv.DictReader(handle, delimiter="\t")

    for row in reader:
        antibiotic = normalize_antibiotic(
            row.get("Antibiotic")
            or row.get("Antibiotic_Original")
        )

        if not antibiotic:
            continue

        phenotype = normalize_phenotype(
            row.get("Resistance_Phenotype")
            or row.get("Resistance_Phenotype_Original")
        )

        row["_antibiotic"] = antibiotic
        row["_phenotype"] = phenotype
        rows.append(row)


print("Input observations:", len(rows))
print("Unique SRR:", len({r["SRR"] for r in rows}))


# ============================================================
# PHENOTYPE SUMMARY
# ============================================================

pheno_summary = defaultdict(
    lambda: {
        "observations": 0,
        "SRR": set(),
        "R": 0,
        "I": 0,
        "S": 0,
        "undefined": 0,
    }
)

for row in rows:
    a = row["_antibiotic"]
    p = row["_phenotype"]

    pheno_summary[a]["observations"] += 1
    pheno_summary[a]["SRR"].add(row["SRR"])

    if p in {"R", "I", "S"}:
        pheno_summary[a][p] += 1
    else:
        pheno_summary[a]["undefined"] += 1


with (OUT / "01_Final_Phenotype_Summary.tsv").open("w", newline="") as f:
    writer = csv.writer(f, delimiter="\t")

    writer.writerow([
        "Antibiotic",
        "unique_isolates",
        "observations",
        "R",
        "I",
        "S",
        "not_defined",
        "defined",
        "R_percent_defined",
    ])

    for a in sorted(pheno_summary):
        d = pheno_summary[a]
        defined = d["R"] + d["I"] + d["S"]

        pct = (
            100 * d["R"] / defined
            if defined
            else float("nan")
        )

        writer.writerow([
            a,
            len(d["SRR"]),
            d["observations"],
            d["R"],
            d["I"],
            d["S"],
            d["undefined"],
            defined,
            f"{pct:.2f}",
        ])


# ============================================================
# DUPLICATE AUDIT + DEDUPLICATION
# ============================================================

groups = defaultdict(list)

for row in rows:
    key = (
        row["SRR"],
        row["_antibiotic"],
    )
    groups[key].append(row)


deduplicated = []
conflicts = []
duplicates = []

for key, group in groups.items():

    if len(group) > 1:
        duplicates.append(
            (
                key[0],
                key[1],
                len(group),
            )
        )

    phenotypes = {
        r["_phenotype"]
        for r in group
        if r["_phenotype"]
    }

    if len(phenotypes) > 1:
        conflicts.append(
            (
                key[0],
                key[1],
                ";".join(sorted(phenotypes)),
                len(group),
            )
        )
        continue

    deduplicated.append(group[0])


with (OUT / "00_Duplicate_Isolate_Antibiotic_Audit.tsv").open(
    "w",
    newline=""
) as f:
    writer = csv.writer(f, delimiter="\t")

    writer.writerow([
        "SRR",
        "Antibiotic",
        "record_count",
    ])

    writer.writerows(duplicates)


with (OUT / "00_Conflicting_Phenotype_Records.tsv").open(
    "w",
    newline=""
) as f:
    writer = csv.writer(f, delimiter="\t")

    writer.writerow([
        "SRR",
        "Antibiotic",
        "Phenotypes",
        "record_count",
    ])

    writer.writerows(conflicts)


# ============================================================
# FINAL CONCORDANCE
# ============================================================

detail = []
confusion = defaultdict(
    lambda: {
        "TP": 0,
        "TN": 0,
        "FP": 0,
        "FN": 0,
        "I": 0,
        "Undefined": 0,
    }
)


for row in deduplicated:

    antibiotic = row["_antibiotic"]
    phenotype = row["_phenotype"]

    if phenotype == "I":
        confusion[antibiotic]["I"] += 1
        continue

    if phenotype not in {"R", "S"}:
        confusion[antibiotic]["Undefined"] += 1
        continue

    prediction, matched = genotype_prediction(
        row,
        antibiotic,
    )

    if phenotype == "R" and prediction == "R":
        category = "TP"

    elif phenotype == "S" and prediction == "S":
        category = "TN"

    elif phenotype == "S" and prediction == "R":
        category = "FP"

    else:
        category = "FN"

    confusion[antibiotic][category] += 1

    detail.append({
        "SRR": row["SRR"],
        "Antibiotic": antibiotic,
        "Phenotype": phenotype,
        "Genotype_Prediction": prediction,
        "Category": category,
        "Matched_Genotype_Rule": matched,
        "ST": clean(row.get("ST")),
        "virulence_score": clean(row.get("virulence_score")),
        "resistance_score": clean(row.get("resistance_score")),
        "siderophore_profile": siderophore_profile(row),
        "amr_gene_list": clean(row.get("amr_gene_list")),
        "point_mutation_list": clean(row.get("point_mutation_list")),
    })


detail_columns = [
    "SRR",
    "Antibiotic",
    "Phenotype",
    "Genotype_Prediction",
    "Category",
    "Matched_Genotype_Rule",
    "ST",
    "virulence_score",
    "resistance_score",
    "siderophore_profile",
    "amr_gene_list",
    "point_mutation_list",
]

with (OUT / "02_Final_Concordance_Detail.tsv").open(
    "w",
    newline=""
) as f:

    writer = csv.DictWriter(
        f,
        delimiter="\t",
        fieldnames=detail_columns,
    )

    writer.writeheader()
    writer.writerows(detail)


def div(a, b):
    return a / b if b else float("nan")


with (OUT / "03_Final_Concordance_Summary.tsv").open(
    "w",
    newline=""
) as f:

    writer = csv.writer(f, delimiter="\t")

    writer.writerow([
        "Antibiotic",
        "TP",
        "TN",
        "FP",
        "FN",
        "Evaluated",
        "Sensitivity",
        "Specificity",
        "PPV",
        "NPV",
        "Accuracy",
        "Excluded_Intermediate",
        "Excluded_Undefined",
    ])

    for a in sorted(confusion):

        d = confusion[a]

        tp = d["TP"]
        tn = d["TN"]
        fp = d["FP"]
        fn = d["FN"]

        evaluated = tp + tn + fp + fn

        writer.writerow([
            a,
            tp,
            tn,
            fp,
            fn,
            evaluated,
            f"{div(tp,tp+fn):.4f}",
            f"{div(tn,tn+fp):.4f}",
            f"{div(tp,tp+fp):.4f}",
            f"{div(tn,tn+fn):.4f}",
            f"{div(tp+tn,evaluated):.4f}",
            d["I"],
            d["Undefined"],
        ])


# ============================================================
# DISCORDANT ISOLATES
# ============================================================

discordant = [
    r for r in detail
    if r["Category"] in {"FP", "FN"}
]

with (OUT / "04_Final_Discordant_Isolates.tsv").open(
    "w",
    newline=""
) as f:

    writer = csv.DictWriter(
        f,
        delimiter="\t",
        fieldnames=detail_columns,
    )

    writer.writeheader()
    writer.writerows(discordant)


# ============================================================
# DISCORDANCE BY ST / SIDEROPHORE
# ============================================================

by_st = defaultdict(
    lambda: defaultdict(int)
)

by_sid = defaultdict(
    lambda: defaultdict(int)
)

for r in detail:

    by_st[
        (r["Antibiotic"], r["ST"])
    ][r["Category"]] += 1

    by_sid[
        (
            r["Antibiotic"],
            r["siderophore_profile"],
        )
    ][r["Category"]] += 1


with (OUT / "05_Final_Discordance_By_ST.tsv").open(
    "w",
    newline=""
) as f:

    writer = csv.writer(f, delimiter="\t")

    writer.writerow([
        "Antibiotic",
        "ST",
        "TP",
        "TN",
        "FP",
        "FN",
        "Evaluated",
        "Discordant",
        "Discordance_percent",
    ])

    for (a, st), d in sorted(by_st.items()):

        total = sum(d.values())
        discord = d["FP"] + d["FN"]

        writer.writerow([
            a,
            st,
            d["TP"],
            d["TN"],
            d["FP"],
            d["FN"],
            total,
            discord,
            f"{100*discord/total:.2f}",
        ])


with (OUT / "06_Final_Discordance_By_Siderophore.tsv").open(
    "w",
    newline=""
) as f:

    writer = csv.writer(f, delimiter="\t")

    writer.writerow([
        "Antibiotic",
        "siderophore_profile",
        "TP",
        "TN",
        "FP",
        "FN",
        "Evaluated",
        "Discordant",
        "Discordance_percent",
    ])

    for (a, sid), d in sorted(by_sid.items()):

        total = sum(d.values())
        discord = d["FP"] + d["FN"]

        writer.writerow([
            a,
            sid,
            d["TP"],
            d["TN"],
            d["FP"],
            d["FN"],
            total,
            discord,
            f"{100*discord/total:.2f}",
        ])


# ============================================================
# CARBAPENEM FALSE-NEGATIVE MECHANISM AUDIT
# ============================================================

with (OUT / "07_Final_Carbapenem_FN_Mechanisms.tsv").open(
    "w",
    newline=""
) as f:

    writer = csv.writer(f, delimiter="\t")

    writer.writerow([
        "SRR",
        "Antibiotic",
        "ST",
        "ESBL_AmpC_detected",
        "Porin_alteration_detected",
        "ESBL_AmpC_hits",
        "Porin_hits",
        "amr_gene_list",
        "point_mutation_list",
    ])

    for r in discordant:

        if (
            r["Category"] != "FN"
            or r["Antibiotic"]
            not in {"Imipenem", "Meropenem"}
        ):
            continue

        genes = split_items(r["amr_gene_list"])
        muts = split_items(r["point_mutation_list"])

        beta_hits = regex_hits(
            genes,
            ESBL_AMP_C_PATTERNS,
        )

        porin_hits = regex_hits(
            muts,
            PORIN_PATTERNS,
        )

        writer.writerow([
            r["SRR"],
            r["Antibiotic"],
            r["ST"],
            "Yes" if beta_hits else "No",
            "Yes" if porin_hits else "No",
            ";".join(beta_hits),
            ";".join(porin_hits),
            r["amr_gene_list"],
            r["point_mutation_list"],
        ])


print("Deduplicated isolate-antibiotic records:", len(deduplicated))
print("Duplicate groups:", len(duplicates))
print("Conflicting phenotype groups:", len(conflicts))
print("Evaluated concordance records:", len(detail))
print("Discordant records:", len(discordant))
print("Final analysis completed.")


