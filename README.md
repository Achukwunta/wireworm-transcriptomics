# *De Novo* Transcriptome Assembly and Comparative Analysis of *Hypnoidus abbreviatus*

[![Snakemake](https://img.shields.io/badge/snakemake-≥7.32.4-brightgreen.svg)](https://snakemake.github.io)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

This repository contains a comprehensive pipeline for the *de novo* assembly, quality control, filtering, and annotation of the wireworm (*Hypnoidus abbreviatus*) transcriptome. The project encompasses raw read processing, contamination screening, functional annotation, and comparative analysis against related insect species.

This work is part of my MSc thesis research at Brandon University.

## 📊 Project Overview

Wireworms are soil-dwelling larvae of click beetles (Coleoptera: Elateridae), with a worm-like appearance. They are small (half inch in size), tend to be hard-bodied, segmented and yellowish in colour, with three pairs of legs. Wireworms are among the most challenging agricultural pests to manage worldwide, especially in North America. *Hypnoidus abbreviatus* and *Limonius californicus* are prevalent species in Western Canada, yet little is known about their molecular biology, particularly genes involved in development and insecticide resistanceUnfortunately, due to lack of genomic resources for them. This project aims to construct a high-quality reference transcriptomes to enable future studies on its biology and management.

**Key Goals:**
- Generate a high-quality, contaminant-free *de novo* transcriptome assembly for *H. abbreviatus* and *L. californicus*.
- Annotate putative functions for assembled transcripts.
- Conduct a comparative analysis against other Coleoptera species to identify conserved and unique genes.
- Prepare the transcriptome for submission to public databases (NCBI TSA).

## 📁 Project Structure

```
wireworm-transcriptomics/
├── config/ # Snakemake configuration files
│ └── config.yaml
├── data/ # Raw data (not stored on GitHub)
│ └── raw/ # Place raw FASTQ files here
├── resources/ # Reference databases (e.g., BLAST nt, EggNOG)
├── workflow/ # Snakemake workflow definition
│ └── Snakefile
├── scripts/ # Helper scripts for analysis
├── notebooks/ # Jupyter notebooks for visualization
├── results/ # Pipeline output (assemblies, annotations, plots)
├── logs/ # Log files from pipeline execution
├── .gitignore
├── environment.yml
└── README.md
```


## ⚙️ Installation & Setup

1.  **Clone the repository:**
    ```bash
    git clone git@github.com:Achukwunta/wireworm-transcriptomics.git
    cd wireworm-transcriptomics
    ```

2.  **Download Required Databases:**  
    Key databases (BLAST `nt`, EggNOG, BUSCO `insecta_odb10`) must be downloaded and placed in the `resources/` directory. Paths are specified in `config.yaml`.

3.  **Create and activate the Conda environment:**
    ```bash
    conda env create -f environment.yml
    conda activate wireworm-transcriptomics
    ```

## 🚀 Running the Pipeline

The pipeline is managed with Snakemake. Execute the entire analysis with:

```bash
snakemake --cores all --use-conda
```
### Pipeline Steps
1. - **Quality Control & Trimming**: `FastQC` -> `MultiQC` -> `fastp` (adapter/quality trimming) -> `seqkit rmdup` (duplicate removal).
2. - **De Novo Assembly**: `Trinity` assembly of filtered reads.
3. - **Assembly QC**: `BUSCO` (completeness), `TrinityStats.pl`, `QUAST` (assembly metrics).
4. - **Transcript Filtering**: Length filtering (`seqkit`), expression-based filtering (`Salmon` quantification).
5. - **Contamination Screening**: `BLASTn` vs. `nt` database and `FCS-Adapter` to remove non-metazoan contaminants and vectors.
6. - **Functional Annotation**: `TransDecoder` (ORF prediction) -> `EggNOG-mapper` (functional annotation).
7. - **Comparative Analysis**: `tBLASTx` against transcripts from *Photinus pyralis*, *Agrilus planipennis*, *Onthophagus taurus*, and *Tribolium castaneum*.

## 🔧 Configuration
Project-specific variables are set in config/config.yaml.

Example config.yaml:

- samples: ["Habb_WWA_6", "Habb_WWA_7", "Habb_WWA_8"]
- params:
  trimmomatic: "ILLUMINACLIP:resources/adapters/NexteraPE-PE.fa:2:30:10:2:True LEADING:3 TRAILING:3 SLIDINGWINDOW:5:18 MINLEN:30"
- references:
  nt_db: "resources/blastdb/nt"
  busco_lineage: "resources/busco_downloads/lineages/insecta_odb10"


## 📊 Data Sources
- Raw Reads: Generated in-house (Brandon University, Cassone Lab).

- Comparative Transcriptomes: Downloaded from NCBI for comparative analysis.

- Photinus pyralis (Firefly): GCF_008802855.1

- Agrilus planipennis (Emerald Ash Borer): GCF_000699045.2

- Onthophagus taurus (Dung Beetle): GCF_036711975.1

- Tribolium castaneum (Red Flour Beetle): GCF_031307605.1

## 👨‍💻 Author
Augustine Chukwunta

MSc Candidate in Biology, Brandon University, Manitoba, Canada

Thesis: "Multi-Omics Approaches to Unravel Plant-Pathogen Interactions: Transcriptomics and Microbiome Analysis."

Advisor: Dr. Bryan Cassone

https://img.shields.io/badge/GitHub-Achukwunta-blue?logo=github

## 📜 Citation
This work is part of an ongoing study. If you use this pipeline or results, please cite:

Chukwunta A., & Cassone B. (2025). "De novo assembly and comparative analysis of Hypnoidus abbreviatus and Limonius californicus transcriptomes." (Target Journal: Entomology).

## 🤝 Contributing & Contact
For questions about the analysis, suggestions for improvement, or potential collaboration, please contact Augustine Chukwunta or open an issue on this repository.

Thank you and Happy research!
