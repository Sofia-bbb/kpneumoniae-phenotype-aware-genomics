# Phenotype-Aware Genomic Surveillance of *Klebsiella pneumoniae*

## Overview

This repository contains the downstream reproducible analysis workflow and supporting analysis-ready datasets for the manuscript:

**Phenotype-Aware Genomic Surveillance Reveals Siderophore-Associated Virulence–Resistance Backgrounds in *Klebsiella pneumoniae***

The study integrates whole-genome sequencing data with antimicrobial susceptibility phenotypes to examine antimicrobial resistance determinants, genotype–phenotype concordance and discordance, siderophore and virulence profiles, sequence types, capsule loci, and phylogenomic structure in *Klebsiella pneumoniae*.

The workflow was designed to support phenotype-aware genomic surveillance by identifying genomic and lineage-specific patterns that may be missed when genomic resistance determinants or phenotypic susceptibility data are considered independently.

---

## Study Workflow

The complete study comprised the following major stages:

```text
1. Data retrieval and cohort construction
        ↓
2. Read processing and quality control
        ↓
3. Genome assembly and assembly quality control
        ↓
4. Species validation and final cohort definition
        ↓
5. Genomic characterization
   ├── Sequence typing
   ├── Capsule and O-antigen loci
   ├── Siderophore profiling
   ├── AMR determinants
   ├── Resistance-associated mutations
   └── Virulence-factor profiling
        ↓
6. Genotype–phenotype integration
        ↓
7. Statistical and concordance analyses
        ↓
8. Core-genome and phylogenomic analyses
        ↓
9. Integrated lineage-specific interpretation
```

The present repository focuses on the **downstream genomic characterization, genotype–phenotype integration, siderophore analyses, statistical analyses, and phylogenomic analyses** used to generate the reported results.

Raw-read retrieval, read trimming, genome assembly, and large intermediate genome-processing files are not redistributed through this repository.

---

## Repository Structure

```text
kpneumoniae-phenotype-aware-genomics/
│
├── data/
│   ├── 01_Stage4_Integrated_Isolate_Master_FINAL_1143.tsv
│   ├── 02_AMR_Table.tsv
│   ├── 02_Eligible_ST_Siderophore_Groups_nGE5.tsv
│   ├── 03_Isolate_Level_Siderophore_AMR_Genotype_Phenotype.tsv
│   ├── 04_Siderophore_Table.tsv
│   ├── 06_ST_Table.tsv
│   ├── 16_concordance_detail_deduplicated.tsv
│   ├── GPSKp_genomic_master.tsv
│   ├── GPSKp_genotype_phenotype_long.tsv
│   ├── KP_Final_1143_core_genome_tree.nwk
│   ├── ST258_Klocus_FINAL1143.tsv
│   ├── Supplementary_Table_S7_Integrated_Genotype_Phenotype_FINAL.tsv
│   ├── vfdb_non_siderophore_FINAL1168_90id_80cov.tsv
│   └── assemblies.example.tsv
│
├── scripts/
│   ├── genome_characterization/
│   ├── genotype_phenotype/
│   ├── siderophore/
│   ├── statistics/
│   └── phylogeny/
│
├── .gitignore
└── README.md
```

---

## Analysis Scripts

### Genome Characterization

`scripts/genome_characterization/`

Contains scripts supporting VFDB-based virulence-factor processing and construction of the non-siderophore virulence-feature dataset.

Key scripts:

- `GPSKp_04_VFDB.sh`
- `04b_vfdb_coverage_filter.py`
- `04c_build_non_siderophore_vf.py`

## Genotype–Phenotype Concordance and Discordance

`scripts/genotype_phenotype/`

Contains the final antimicrobial genotype–phenotype concordance analysisand discordance analyses, including antimicrobial-specific genomic prediction rules and classification of genotype-positive/phenotype-resistant (G+/P-R), genotype-negative/phenotype-susceptible (G−/P-S), genotype-positive/phenotype-susceptible (G+/P-S), and genotype-negative/phenotype-resistant (G−/P-R) observations.

Key script:

- `01_final_concordance.py`

## Siderophore Analyses

`scripts/siderophore/`

Contains scripts supporting siderophore-associated genomic and phenotypic comparisons.

Key scripts:

- `01_siderophore_associations.R`
- `02_adjusted_siderophore_models.R`

## Statistical Analysis

`scripts/statistics/`

Contains scripts used to construct analysis-ready tables for integrated statistical analyses.

Key script:

- `01_build_analysis_tables.py`

## Phylogenomic Analyses

`scripts/phylogeny/`

Contains scripts supporting within-lineage phylogenetic structure, virulence clustering, ST-adjusted carbapenem discordance clustering, and focused ST258 KL106/KL107 analyses.

Key scripts:

- `01_Final1143_Phylogenetic_Structure.R`
- `10b_Stage3_Virulence_Phylogeny_No_Siderophores.R`
- `12_ST_Adjusted_Carbapenem_FINAL.R`
- `13_ST258_KL106_KL107_Integrated_Analysis.R`

---

## Reproducing the Downstream Analyses

The repository is organized so that the deposited analysis-ready datasets in `data/` can be used with the corresponding scripts in `scripts/`.

The general downstream analysis sequence is:

```text
Analysis-ready genomic and phenotype datasets
        ↓
Build integrated analysis tables
        ↓
Genotype–phenotype concordance analysis
        ↓
Siderophore-associated analyses
        ↓
Within-lineage analyses
        ↓
Phylogenetic clustering and structure analyses
```

Scripts use relative paths where applicable to facilitate reuse outside the original high-performance computing environment.

Analyses requiring locally available genome assemblies can use:

`data/assemblies.example.tsv`

as a template for defining local assembly paths.

---

## Software

The study used the following principal software packages:

- FastQC
- Trimmomatic v0.41
- SPAdes
- QUAST v5.0.2
- Kleborate v3.2.4
- AMRFinderPlus v4.0.23
- VFDB Set A
- BLASTN
- Prokka
- Roary v3.13.0
- FastTree v2.1.11
- Python
- R

The AMRFinderPlus database release used in the final analysis was:

**2025-07-16.1**

Exact statistical procedures, filtering thresholds, and antimicrobial-specific genomic prediction rules are described in the manuscript Methods.

---

## Data Availability

The genomic and antimicrobial susceptibility data analysed in this study were obtained from publicly available records in **NCBI Pathogen Detection** and the **Sequence Read Archive (SRA)**.

Accession identifiers for the isolates included in the study are provided in the Supplementary Material accompanying the manuscript.

This repository provides derived, analysis-ready datasets supporting the reported genotype–phenotype, siderophore, virulence, statistical, and phylogenomic analyses.

Raw sequencing reads and genome assemblies are not redistributed through this repository because they remain available through their original public NCBI accessions.

---

## Reproducibility Scope

This repository is intended to support reproduction of the **downstream analyses reported in the manuscript** from the deposited analysis-ready datasets.

The upstream computational stages, including raw SRA retrieval, read processing, genome assembly, assembly quality assessment, and generation of large intermediate files, were performed on high-performance computing infrastructure and are described in the manuscript Methods.

---

## Citation

If you use this workflow or the associated datasets, please cite the corresponding manuscript.

Full citation information will be added following publication.

---

## Authors

**Sofia B. Mohamed**  
Department of Chemistry and Biochemistry  
University of Northern British Columbia (UNBC)  
Prince George, British Columbia, Canada  

National University Biomedical Research Institute (NUBRI)  
National University-Sudan  
Khartoum, Sudan  

**Jacob P. Walsh**  
Department of Chemistry and Biochemistry  
University of Northern British Columbia (UNBC)  
Prince George, British Columbia, Canada

