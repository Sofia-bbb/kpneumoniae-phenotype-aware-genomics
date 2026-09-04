## Phenotype-Aware Genomic Surveillance of Klebsiella pneumoniae

### Overview

This repository contains the reproducible upstream and downstream analysis workflow and supporting analysis-ready datasets for the manuscript:

**Phenotype-Aware Genomic Surveillance Characterizes Siderophore-Associated Virulence–Resistance Backgrounds in *Klebsiella pneumoniae***

The study integrates whole-genome sequencing data with antimicrobial susceptibility phenotypes to examine antimicrobial resistance determinants, genotype–phenotype concordance and discordance, siderophore and virulence profiles, sequence types, capsule loci, and phylogenomic structure in *Klebsiella pneumoniae*.

The workflow was designed to support phenotype-aware genomic surveillance by identifying genomic and lineage-specific patterns that may be missed when genomic resistance determinants or phenotypic susceptibility data are considered independently.

Following assembly-quality assessment and species validation, the final genomic cohort comprised 1,156 *K. pneumoniae* isolates; 12 isolates assigned to other species within the *K. pneumoniae* species complex were excluded. The final phylogenetic cohort comprised 1,143 isolates with valid annotations included in the core-genome analysis. Genotype–phenotype integration included 3,280 initial isolate–antimicrobial observations, 3,271 deduplicated records, and 2,757 observations with evaluable resistant or susceptible classifications.

---

### Study Workflow

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

The repository includes the **upstream data-retrieval, read-processing, assembly, quality-control, species-validation, genomic-characterization, and phylogenomic workflow**, together with the downstream genotype–phenotype, siderophore, statistical, and phylogenetic analyses used to generate the reported results.

Raw sequencing reads, genome assemblies, and large intermediate files are not redistributed; however, the scripts and configuration templates required to reproduce these stages are provided.

---

### Repository Structure

```text
kpneumoniae-phenotype-aware-genomics/
│
├── config/
│   ├── sra_accessions_1181.txt
│   └── assemblies.example.tsv
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
│   └── vfdb_non_siderophore_FINAL1156_90id_80cov.tsv
│
├── scripts/
│   ├── data_retrieval/
│   ├── read_processing/
│   ├── assembly_qc/
│   ├── cohort_construction/
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

### Analysis Scripts

### Data Retrieval and Cohort Construction

`scripts/data_retrieval/` and `scripts/cohort_construction/`

Contains scripts for NCBI record retrieval, metadata harmonization, AST integration, human-clinical filtering, representative-run selection, FASTQ retrieval, and construction of the isolate-level study cohort.

### Read Processing, Assembly, and Quality Control

`scripts/read_processing/` and `scripts/assembly_qc/`

Contains Slurm workflows for FastQC, Trimmomatic, SPAdes, and QUAST, together with generation of the cohort-level assembly-quality summary.

### Genome Characterization

`scripts/genome_characterization/`

Contains scripts for Kleborate species validation and typing, AMRFinderPlus resistance-determinant detection, VFDB-based virulence-factor processing, and construction of the non-siderophore virulence-feature dataset.

Key scripts:

- `01_Kleborate.sh`
- `GPSKp_03_AMRFinder.sh`
- `GPSKp_04_VFDB.sh`
- `04b_vfdb_coverage_filter.py`
- `04c_build_non_siderophore_vf.py`

### Genotype–Phenotype Concordance and Discordance

`scripts/genotype_phenotype/`

Contains the final antimicrobial genotype–phenotype concordance analysisand discordance analyses, including antimicrobial-specific genomic prediction rules and classification of genotype-positive/phenotype-resistant (G+/P-R), genotype-negative/phenotype-susceptible (G−/P-S), genotype-positive/phenotype-susceptible (G+/P-S), and genotype-negative/phenotype-resistant (G−/P-R) observations.

Key script:

- `01_final_concordance.py`

### Siderophore Analyses

`scripts/siderophore/`

Contains scripts supporting siderophore-associated genomic and phenotypic comparisons.

Key scripts:

- `01_siderophore_associations.R`
- `02_adjusted_siderophore_models.R`

### Statistical Analysis

`scripts/statistics/`

Contains scripts used to construct analysis-ready tables for integrated statistical analyses.

Key script:

- `01_build_analysis_tables.py`

### Phylogenomic Analyses

`scripts/phylogeny/`

Contains scripts for Prokka annotation, Roary core-genome alignment, FastTree phylogenetic reconstruction, within-lineage phylogenetic structure, virulence clustering, ST-adjusted carbapenem discordance clustering, and focused ST258 KL106/KL107 analyses.

Key scripts:

- `GPSKp_11_Prokka.sh`
- `run_final_roary.sh`
- `run_FastTree_1143.sh`
- `01_Final1143_Phylogenetic_Structure.R`
- `10b_Stage3_Virulence_Phylogeny_No_Siderophores.R`
- `12_ST_Adjusted_Carbapenem_FINAL.R`
- `13_ST258_KL106_KL107_Integrated_Analysis.R`

---

## Reproducing the Complete Workflow

The repository is organized so that the deposited analysis-ready datasets in `data/` can be used with the corresponding scripts in `scripts/`.

Run scripts from the repository root. Before submitting Slurm jobs, create the local working directories with `mkdir -p logs work results`. Copy `config/assemblies.example.tsv` to `config/assemblies.tsv` and replace the example entries with local assembly paths. Slurm array ranges should be supplied at submission according to the number of records in the relevant manifest; they are not hard-coded in the public scripts.

The default working directory is `work/` and can be overridden with `KP_WORKDIR`. Manifest, container, Roary, and FastTree paths can also be overridden using the environment variables defined in the corresponding scripts. Kleborate and AMRFinderPlus container images are not redistributed through this repository.

The general analysis sequence is:

```text
NCBI data retrieval and cohort construction
        ↓
SRA/FASTQ retrieval and read-quality control
        ↓
Read trimming, genome assembly, and assembly QC
        ↓
Species validation and genomic characterization
        ↓
Prokka annotation, Roary alignment, and FastTree reconstruction
        ↓
Analysis-ready genomic and phenotype datasets
        ↓
Build integrated analysis tables
        ↓
Genotype–phenotype concordance and discordance analyses
        ↓
Siderophore-associated and within-lineage analyses
        ↓
Phylogenetic clustering and integrated interpretation
```

Scripts use relative paths where applicable to facilitate reuse outside the original high-performance computing environment.

Analyses requiring locally available genome assemblies can use `config/assemblies.example.tsv` as a template for defining local assembly paths.

---

### Software

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

### Data Availability

The genomic and antimicrobial susceptibility data analysed in this study were obtained from publicly available records in **NCBI Pathogen Detection** and the **Sequence Read Archive (SRA)**.

Accession identifiers for the isolates included in the study are provided in the Supplementary Material accompanying the manuscript.

This repository provides derived, analysis-ready datasets supporting the reported genotype–phenotype, siderophore, virulence, statistical, and phylogenomic analyses.

Raw sequencing reads and genome assemblies are not redistributed through this repository because they remain available through their original public NCBI accessions.

---

### Reproducibility Scope

This repository provides the scripts used for both the upstream and downstream computational workflow. Upstream stages include public-data retrieval, read processing, genome assembly, assembly-quality assessment, species validation, genomic characterization, annotation, core-genome alignment, and phylogenetic reconstruction.

Because raw reads, assemblies, and large intermediate files are not redistributed, reproducing the upstream stages requires downloading the public SRA records and access to the listed software and suitable computing resources. The deposited analysis-ready datasets permit reproduction of the downstream statistical and phylogenetic analyses without rerunning the upstream workflow.

---

### Citation

If you use this workflow or the associated datasets, please cite the corresponding manuscript.

Full citation information will be added following publication.

---

### Authors

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

