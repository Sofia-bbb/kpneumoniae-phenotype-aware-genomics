## Overview

This repository contains the reproducible analysis workflow and supporting datasets for the manuscript:

**Phenotype-Aware Genomic Surveillance Reveals Siderophore-Associated Virulence–Resistance Backgrounds in *Klebsiella pneumoniae***

The study integrates whole-genome sequencing with antimicrobial susceptibility phenotypes to examine antimicrobial resistance determinants, genotype–phenotype concordance and discordance, siderophore and virulence profiles, sequence types, capsule loci, and phylogenomic structure in Klebsiella pneumoniae.

The workflow was designed to support phenotype-aware genomic surveillance by identifying genomic and lineage-specific patterns that may be missed when genomic resistance determinants or phenotypic susceptibility data are considered independently.

<sub>The workflow was developed by Sofia Ali and Jacob P. Walsh in the Department of Chemistry and Biochemistry, University of Northern British Columbia (UNBC).</sub>
---

## Workflow 

The study followed an integrated genomic–phenotypic workflow comprising the following major stages:

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
   ├── Sequence typing (ST)
   ├── Capsule and O-antigen loci
   ├── Siderophore profiling
   ├── AMR determinants
   ├── Resistance-associated mutations
   └── Virulence-factor profiling
        ↓
6. Genotype–phenotype integration
   ├── Genomic resistance profiles
   ├── Observed AST phenotypes
   └── Isolate–antimicrobial integrated dataset
        ↓
7. Statistical and genotype–phenotype analyses
   ├── Antimicrobial-specific genomic prediction rules
   ├── Concordance analysis
   ├── Discordance analysis
   ├── Diagnostic performance metrics
   └── Siderophore–AMR–phenotype association testing
        ↓
8. Pan-genome and phylogenomic analysis
   ├── Core-genome analysis
   ├── Phylogenetic reconstruction
   ├── Siderophore phylogenetic structuring
   ├── AMR and virulence clustering
   └── Genotype–phenotype discordance clustering
        ↓
9. Integrated surveillance interpretation
   ├── Siderophore–virulence–resistance patterns
   ├── Lineage-specific variation
   ├── Within-lineage analyses
   └── Surveillance-relevant genomic–phenotypic profiles
```

---

## Repository Structure

```text
kpneumoniae-phenotype-aware-genomics/
│
├── data/
│   ├── integrated isolate-level datasets
│   ├── genotype–phenotype datasets
│   ├── AMR and siderophore analysis tables
│   ├── virulence-factor datasets
│   ├── sequence-type and capsule-locus data
│   └── final core-genome phylogenetic tree
│
├── scripts/
│   ├── genome_characterization/
│   ├── statistics/
│   ├── genotype_phenotype/
│   ├── siderophore/
│   └── phylogeny/
│
├── results/
│   └── generated analysis outputs (not tracked)
│
├── .gitignore
└── README.md
```

---

## Reproducibility and Data Availability

The genomic and antimicrobial susceptibility data analysed in this study were obtained from publicly available records in **NCBI Pathogen Detection and the Sequence Read Archive (SRA). Accession identifiers for the isolates included in the study are provided in the Supplementary Material accompanying the manuscript.

This repository contains the analysis scripts and derived, analysis-ready datasets supporting the reported genotype–phenotype, siderophore, statistical, and phylogenomic analyses. Raw sequencing reads and genome assemblies are not redistributed through this repository.

Repository scripts use relative paths where applicable to facilitate reuse outside the original high-performance computing environment. A template assembly manifest (`data/assemblies.example.tsv`) is provided for analyses requiring locally available genome assemblies.

<sub>External resources, including raw sequencing data, the VFDB reference database, and software containers, must be obtained or configured separately where required.</sub>
---

## Authors

Sofia B. Mohamed  
Department of Chemistry and Biochemistry 
University of Northern British Columbia (UNBC)  
Prince George, British Columbia, Canada  

National University Biomedical Research Institute (NUBRI)  
National University-Sudan  
Khartoum, Sudan

Jacob P. Walsh  
Department of Chemistry and Biochemistry  
University of Northern British Columbia (UNBC)  
Prince George, British Columbia, Canada
