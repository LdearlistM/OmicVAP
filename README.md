<div align="center">

<img src="https://github.com/user-attachments/assets/e94f0aee-41f2-4e23-9298-dfef220df643" width="180" alt="xMetaVar Logo">

# xMetaVar: A Systematic Platform for Strain-Level Microbial Variants

[![Docker Image Version](https://img.shields.io/badge/docker-v1.0.0-blue)](https://github.com/ldearlistm/xmetavar/pkgs/container/xmetavar)
[![Website](https://img.shields.io/badge/Web_Server-Available-green)](https://www.biosino.org/iMAC/xmetavar)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

</div>

---

## 📖 Introduction

**xMetaVar** is a systematic one-stop platform for the integrated calling and biological interpretation of multiple types of **strain-level microbial variants** (SNPs, short INDELs, and structural variants) directly from metagenomic sequencing data.

While we provide a user-friendly [**Web Server**](https://www.biosino.org/iMAC/xmetavar) for immediate analysis without installation, this repository hosts the documentation for the **Docker-based command-line workflow**, enabling high-throughput analysis on local servers or HPC clusters.

The pipeline integrates five optimized callers (GT-Pro, MIDAS, QuickVariant, SGVFinder, PhaseFinder) to seamlessly cover:
*   **SNPs**: Single Nucleotide Polymorphisms
*   **INDELs**: Short Insertions and Deletions
*   **SVs**: Structural Variants (dSVs, vSVs, Inversions, CNVs)

## ⚡ Quick Links
- [Web Server (No installation required)](https://www.biosino.org/iMAC/xmetavar)
- [Web Tutorial & Visualization Guide](https://www.biosino.org/iMAC/xmetavar/tutorial)
- [Docker Hub / Packages](https://github.com/ldearlistm/xmetavar/pkgs/container/xmetavar)

---

## 🛠️ Installation

### Hardware Requirements
The pipeline is designed to run on standard Linux servers. 
- **CPU**: 16+ cores recommended for parallel processing.
- **RAM**: 32GB+ (64GB+ recommended for large cohorts).
- **Disk**: Sufficient storage for raw FASTQ files and reference databases (~50GB for the core database).

### Docker Image
We provide a pre-built Docker image containing all dependencies (Snakemake, Python envs, Bioinformatics tools). You do not need to install the tools individually.

```bash
# Pull the latest image
docker pull ghcr.io/ldearlistm/xmetavar:1.0.0
```

---

## 🚀 Usage (Local / Docker)

To run xMetaVar locally, you need to prepare your input files and mount them into the Docker container.

### 1. Data Preparation
Organize your project directory as follows:

```text
my_project/
├── rawdata/               # Directory containing your .fastq.gz files
├── database/              # Directory containing the reference panel
├── output/                # Directory for results (will be created)
├── samples.tsv            # Sample list
└── config.yaml            # Pipeline configuration
```

#### The `samples.tsv`
A tab-separated file listing your sample IDs (one per line).
```text
SampleA
SampleB
SampleC
```

#### The `config.yaml`
Standard Snakemake configuration file. Ensure parameters match your analysis needs (e.g., read length, trim settings).

### 2. Custom Reference Database
By default, the pipeline requires a reference panel (e.g., the 43 core human-gut species panel provided by xMetaVar). If you wish to use a custom database or the standard xMetaVar database locally, ensure the `database/` folder contains the necessary indexed files for Bowtie2, GT-Pro, and MIDAS.

> **Note:** You can download our pre-compiled core database from the [Web Server Resources](https://www.biosino.org/iMAC/xmetavar/download) (if applicable) or construct your own following the tool-specific indexing requirements.

### 3. Running the Pipeline
Use `docker run` to execute the workflow. You must mount your local directories to the specific paths inside the container:

| Local Path | Container Path | Description |
| :--- | :--- | :--- |
| `/path/to/rawdata` | `/pipeline/rawdata` | Input FASTQ files (Read Only) |
| `/path/to/output` | `/pipeline/results` | Output directory (Read/Write) |
| `/path/to/database` | `/pipeline/database` | Reference Database (Read/Write) |
| `/path/to/config.yaml` | `/pipeline/config.yaml` | Configuration file |
| `/path/to/samples.tsv` | `/pipeline/samples.tsv` | Sample list |

#### Example Command (Run All Modules)
```bash
docker run --rm \
  -v "$(pwd)/rawdata:/pipeline/rawdata:ro" \
  -v "$(pwd)/output:/pipeline/results:rw" \
  -v "$(pwd)/config.yaml:/pipeline/config.yaml:ro" \
  -v "$(pwd)/samples.tsv:/pipeline/samples.tsv:ro" \
  -v "$(pwd)/database:/pipeline/database:rw" \
  ghcr.io/ldearlistm/xmetavar:1.0.0 all --cores 16
```

### 4. Modular Analysis
You can run specific modules by replacing `all` with specific target keywords:

| Keyword | Description | Tools Used |
| :--- | :--- | :--- |
| `snp_gtpro` | Rapid SNP genotyping | GT-Pro |
| `snp_midas` | High-precision SNP calling | MIDAS v3 |
| `indel` | Short INDEL detection | QuickVariants |
| `sv_sgvfinder` | Deletion & Variable SVs | SGVFinder |
| `sv_inversion` | Inversion detection | PhaseFinder |
| `sv_midas` | CNV (Copy Number Variation) | MIDAS v3 |
| **`all`** | **Run full pipeline** | **All of above** |

**Example: Running only INDEL and SV analysis**
```bash
docker run --rm [mount arguments...] ghcr.io/ldearlistm/xmetavar:1.0.0 indel sv_sgvfinder --cores 8
```

---

## 📂 Output Files

Upon successful execution, the results will be organized in the `output/02-variant-calling` directory (mapped from container). Key output files include:

*   **SNPs**: 
    *   `SNP/GT-Pro/across-sample/snp.tsv`
    *   `SNP/MIDAS/across_sample/snps/snps_summary.tsv`
*   **INDELs**:
    *   `INDEL/across-sample/indel_anno.tsv` (Annotation)
    *   `INDEL/across-sample/indels.tsv` (Presence/Absence Matrix)
*   **Structural Variants**:
    *   `SV/dSV-vSV/across-sample/dsgv_anno.tsv` (SGVFinder Results)
    *   `SV/Inversion/across-sample/inversion_anno.tsv` (Inversions)
    *   `SV/CNV/across_sample/merge.genes_copynum.tsv` (Copy Number)

For detailed visualization of these results (Heatmaps, PCoA, Genome Browser), we highly recommend uploading these output matrices to the **[xMetaVar Web Server Visualization Module](https://www.biosino.org/iMAC/xmetavar)**.

---

## 📚 Citation

If you use xMetaVar in your research, please cite:

> *[Authors]. xMetaVar: A Systematic Web Platform for Integrated Analysis of Strain-Level Microbial Variants. [Journal Name]. [Year];[Volume]:[Pages]. DOI: ...*

(Citation will be updated upon publication)

## 📧 Contact & Support

*   **Issues:** Please report bugs via the [GitHub Issues](https://github.com/ldearlistm/xmetavar/issues) page.
*   **Web Server:** https://www.biosino.org/iMAC/xmetavar
