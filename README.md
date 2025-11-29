<div align="center">

<!-- 1. Logo 和 标题区域 -->
<img src="https://github.com/user-attachments/assets/e94f0aee-41f2-4e23-9298-dfef220df643" width="200" alt="xMetaVar Logo">

# xMetaVar

<!-- 2. 副标题（我们之前定的那个） -->
### A Systematic Web Platform for Integrated Analysis of Strain-Level Microbial Variants

<!-- 3. 徽章区域 (Badges) - 显得专业 -->
[![Docker Image Version](https://img.shields.io/badge/docker-v1.0.0-blue)](https://github.com/ldearlistm/xmetavar/pkgs/container/xmetavar)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Website](https://img.shields.io/badge/Website-Online-green)](https://www.biosino.org/iMAC/xmetavar)
[![PRs Welcome](https://img.shields.io/badge/PRs-welcome-brightgreen.svg)](http://makeapullrequest.com)

[**Web Server**](https://www.biosino.org/iMAC/xmetavar) | [**Docker Hub**](https://github.com/ldearlistm/xmetavar/pkgs/container/xmetavar) | [**Documentation**](https://github.com/ldearlistm/xmetavar/wiki)

</div>

---

<!-- 4. 简介 (使用之前润色的第一段) -->
## 📖 Introduction

**xMetaVar** is the first systematic web-based one-stop platform that provides integrated calling and biological interpretation of multiple types of **strain-level microbial variants** (SNPs, short INDELs, and diverse structural variants) directly from metagenomic sequencing data.

The platform is built upon a **unified reference panel of 43 core human-gut species**, systematically curated from >13,000 metagenomes to ensure broad prevalence and high precision. No upstream variant files or complex installation are required, making it accessible to both computational and experimental microbiologists.

<!-- 5. 核心特性 (列点展示) -->
## ✨ Key Features

- **🌐 Fully Web-based**: No installation required. Analyze 30–50 samples (~300GB) in 4–8 hours using free cloud resources.
- **🧬 Multi-type Variant Calling**: Simultaneous detection of SNPs, Indels, and SVs (dSVs, vSVs, Inversions, CNVs) using 5 optimized callers.
- **📊 Interactive Visualization**: Powered by modules for variant-based population ordinations, stratified compositions, and phenotype-genotype association landscapes.
- **🤖 Machine Learning Integration**: Automated biomarker selection and classifier construction with interpretable outputs (SHAP/ROC).
- **🔒 Privacy First**: No tracking cookies, no login requirement, and full data privacy.

---

<!-- 6. 快速开始 (提供 Docker 用法) -->
## 🚀 Quick Start (Docker)

If you prefer running the pipeline locally instead of using the [Web Server](https://www.biosino.org/iMAC/xmetavar), you can pull our pre-built Docker image:

```bash
# Pull the image from GitHub Container Registry
docker pull ghcr.io/ldearlistm/xmetavar:1.0.0

# Run the container (Example)
docker run -it ghcr.io/ldearlistm/xmetavar:1.0.0 /bin/bash
