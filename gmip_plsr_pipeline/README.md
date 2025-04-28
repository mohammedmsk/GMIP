# GMIP-PLSR Pipeline

The **GMIP-PLSR Pipeline** is a modular, reproducible Nextflow workflow designed to perform post-GWAS gene prioritization by integrating GWAS summary statistics with multi-omics data. This pipeline leverages Partial Least Squares Regression (PLSR) to manage multicollinearity among features, providing robust gene reprioritization results.

This pipeline corresponds to the final optimized method presented in Chapter 3 of the dissertation and is the main focus of the GMIP manuscript.

---

## 📊 Workflow Overview

The GMIP-PLSR pipeline consists of four main modules:

1. **SNP2Gene Mapping**
   - Converts GWAS SNP-level p-values to gene-level scores using MAGMA.

2. **Machine Learning Method**
   - Fits a Partial Least Squares Regression model to prioritize genes.

3. **Cross-Validation Strategy**
   - Supports LOCO (Leave-One-Chromosome-Out) cross-validation to prevent information leakage.

4. **Evaluation**
   - Evaluates reprioritized gene lists using Benchmarker methodology:
     - Tau (heritability enrichment) scores

---

## 🚀 Quick Start

### 1. Install Nextflow
```bash
curl -s https://get.nextflow.io | bash
```

### 2. Navigate to the GMIP-PLSR Pipeline Directory
```bash
cd gmip_plsr_pipeline
```

### 3. Run the Pipeline
```bash
nextflow run workflows/main.nf -profile standard
```
*(Adjust the `-profile` according to your environment. Available profiles are defined in `/conf/`.)*

### 4. Required Inputs
- GWAS summary statistics (in required MAGMA format)
- Multi-omics feature matrix (pre-generated or using provided scripts)
- Configuration parameters (see `/conf/` examples)

### 5. Outputs
- Reprioritized gene lists (ranked by PLSR scores)
- Normalized Tau heritability scores

---

## 📚 Folder Structure

| Folder | Description |
|:---|:---|
| `/bin/` | Utility scripts for processing inputs and outputs. |
| `/modules/` | Modular Nextflow processes for MAGMA, feature preparation, model fitting, and evaluation. |
| `/workflows/` | Main pipeline definitions (`main.nf`). |
| `/conf/` | Configuration files for different execution profiles (local, HPC, AWS, etc.). |
| `/misc/` | Example data files and templates. |

---

## 📜 Citation

If you use GMIP-PLSR, please cite:

> Mohammed Shabbir Kanchwala, et al.  
> **GMIP-PLSR: A Reproducible Nextflow Pipeline for Post-GWAS Gene Prioritization via Multi-Omics Integration and Partial Least Squares Regression** (Manuscript in preparation, 2025).

GitHub Repository: [https://github.com/mohammedmsk/GMIP](https://github.com/mohammedmsk/GMIP)

---

## 🛠️ Environment and Reproducibility

- Pipeline developed with **Nextflow**.
- Compatible with Linux, MacOS, HPC clusters, and cloud platforms.
- Recommended Nextflow version: `>=21.04.0`
- Docker/Singularity profiles can be added if containerized execution is desired.

---

## 💎 Highlights

- **Handles multicollinearity** in multi-omics features using PLSR.
- **Fully reproducible** with modular and portable Nextflow design.
- **Benchmarking included** using heritability enrichment.
- **Scalable** to hundreds of GWAS datasets.

---

## 📈 Related Work

The GMIP-PLSR pipeline builds upon lessons learned from flexible benchmarking performed in the original GMIP Framework (see `/gmip_framework/`).

---

## 📬 Contact

For questions, suggestions, or issues, please contact:

**Mohammed Shabbir Kanchwala**  
Email: [mohammedmsk@gmail.com]

---

# ✅ Happy Prioritizing!