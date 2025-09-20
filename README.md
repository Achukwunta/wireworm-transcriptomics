# *De Novo* Transcriptome Assembly and Comparative Analysis of *Hypnoidus abbreviatus*

[![Snakemake](https://img.shields.io/badge/snakemake-≥7.32.4-brightgreen.svg)](https://snakemake.github.io)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

This repository contains a comprehensive pipeline for the *de novo* assembly, quality control, filtering, and annotation of the wireworm (*Hypnoidus abbreviatus*) transcriptome. The project encompasses raw read processing, contamination screening, functional annotation, and comparative analysis against related insect species. This work is part of my MSc thesis research at Brandon University.

## 📊 Project Overview

Wireworms are soil-dwelling larvae of click beetles (Coleoptera: Elateridae), and are among the most challenging agricultural pests to manage worldwide. *Hypnoidus abbreviatus* and *Limonius californicus* are prevalent in Western Canada. Due to lack of genomic resources, little is known about their molecular biology, particularly genes involved in development and fungicide resistance. This project aims to construct a high-quality reference transcriptomes to enable future studies on its biology and management.

**Key Goals:**
- Generate a high-quality, contaminant-free *de novo* transcriptome assembly for *H. abbreviatus*.
- Annotate putative functions for assembled transcripts.
- Conduct a comparative analysis against other Coleoptera species to identify conserved and unique genes.
- Prepare the transcriptome for submission to public databases (NCBI TSA).

## 📁 Project Structure

```
wireworm-transcriptomics/
├── config/            # Snakemake configuration files
│ └── config.yaml
├── data/               # Raw data (not stored on GitHub)
│ └── raw/              # Place raw FASTQ files here
├── resources/          # Reference databases (e.g., BLAST nt, EggNOG)
├── workflow/           # Snakemake workflow definition
│ └── Snakefile
├── scripts/            # Helper scripts for analysis
├── notebooks/          # Jupyter notebooks for visualization
├── results/            # Pipeline output (assemblies, annotations, plots)
├── logs/               # Log files from pipeline execution
├── .gitignore
├── environment.yml
└── README.md
```

## ⚙️ Installation & Setup
1. **Snakemake** are best installed via the [Conda](https://docs.conda.io/en/latest/). It is recommended to install **Conda** via **Miniforge**. To install run:
 ```bash
 conda create -c conda-forge -c bioconda -c nodefaults --name snakemake snakemake
 conda activate snakemake
 ```

2. **Install Git**: If you don't have git, install it first. For Debian-based distributions like Ubuntu, Mint use:
 ```bash
sudo apt-get update
sudo apt-get install git
```

3.  **Clone the repository:**
 ```bash
 git clone https://github.com/Achukwunta/wireworm-transcriptomics.git
 cd wireworm-transcriptomics
 ```

4. **Download Reference Data:**  
Key databases (`BLAST nt`, `EggNOG`, `BUSCO`, `insecta_odb10`) must be downloaded and placed in the `resources/` directory. The `config.yml` file expects them at specific paths. See the **Data Sources** section below.

5.  **Create and activate the Conda environment:**
Install any required dependencies from `environment.yml` if you're not using **conda** to install everything initially:
 ```bash
 conda env create -f environment.yml
 conda activate wireworm-transcriptomics
 ```

## 🚀 Running the Pipeline
The pipeline is managed with Snakemake. You need to first install:
```bash
#run the full analysis with:
snakemake --cores all --use-conda

#run a specific rule, like FastQC:
snakemake fastqc_raw --cores 2 --use-conda
```
Use `--use-conda` to ensure all dependencies are managed in isolated environments
 
### Pipeline Steps
1. **Quality Control & Trimming**: `FastQC` -> `MultiQC` -> `fastp` (adapter/quality trimming) -> `seqkit rmdup` (duplicate removal)
2. **De Novo Assembly**: `Trinity` assembly of filtered reads.
3. **Assembly QC**: `BUSCO` (completeness), `TrinityStats.pl`, `QUAST` (assembly metrics)
4. **Transcript Filtering**: Length filtering (`seqkit`), expression-based filtering (`Salmon` quantification)
5. **Contamination Screening**: `BLASTn` vs. `nt` database and `FCS-Adapter` to remove non-metazoan contaminants and vectors.
6. **Functional Annotation**: `TransDecoder` (ORF prediction) -> `EggNOG-mapper` (functional annotation)
7. **Comparative Analysis**: `tBLASTx` against transcripts from *Photinus pyralis*, *Agrilus planipennis*, *Onthophagus taurus*, and *Tribolium castaneum*

## 🔧 Configuration
Project-specific variables are set in config/config.yml.
Example config.yml:

- **samples**: ["Habb_WWA_6", "Habb_WWA_7", "Habb_WWA_8"]
- **params**:
  trimmomatic: "ILLUMINACLIP:resources/adapters/NexteraPE-PE.fa:2:30:10:2:True LEADING:3 TRAILING:3 SLIDINGWINDOW:5:18 MINLEN:30"
- **references**:
  nt_db: "resources/blastdb/nt"
  busco_lineage: "resources/busco_downloads/lineages/insecta_odb10"

## 📊 Data Sources
- **Raw Reads**: Generated in-house (Brandon University, Cassone Lab).
- **Comparative Transcriptomes**: Downloaded from NCBI for comparative analysis.
- *Photinus pyralis* (Firefly): GCF_008802855.1
- *Agrilus planipennis* (Emerald Ash Borer): GCF_000699045.2
- *Onthophagus taurus* (Dung Beetle): GCF_036711975.1
- *Tribolium castaneum* (Red Flour Beetle): GCF_031307605.1

## 👨‍💻 Author
**Augustine Chukwunta**  
MSc Biology, Brandon University, Manitoba, Canada
Thesis: *"Multi-Omics Approaches to Unravel Plant-Pathogen Interactions: Transcriptomics and Microbiome Analysis."*  
Advisor: Dr. Bryan Cassone  

[![GitHub: Achukwunta](https://img.shields.io/badge/GitHub-Achukwunta-blue?logo=github)](https://github.com/Achukwunta)

## 📜 Citation
This work is part of an ongoing study. If you use this pipeline or results, please cite:
Chukwunta A., & Cassone B. (2025). "De novo assembly and comparative analysis of Hypnoidus abbreviatus and Limonius californicus transcriptomes." (Target Journal: Entomology).

## 🤝 Contributing & Contact
For questions about the analysis, suggestions for improvement, or potential collaboration, please contact Augustine Chukwunta or open an issue on this repository.

Happy research!


