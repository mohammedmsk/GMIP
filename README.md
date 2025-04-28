# GMIP: GWAS and Multi-Omics Integration Pipeline

**GMIP** is an open-source, modular framework for post-GWAS gene prioritization. It integrates GWAS results with multi-omics data—such as gene expression, protein-protein interaction networks, and biological pathways—using scalable, reproducible Nextflow pipelines.

This repository provides two related components:

- **GMIP-PLSR Pipeline:**  
  The final, recommended pipeline that uses Partial Least Squares Regression (PLSR) to handle multicollinearity and improve gene prioritization performance.

- **GMIP Framework:**  
  A flexible, earlier framework developed to benchmark different multi-omics feature sets, machine learning models, and cross-validation strategies. Insights from this framework informed the development of GMIP-PLSR.

---

## 📚 Repository Structure

| Folder | Description |
|:---|:---|
| `/gmip_plsr_pipeline/` | GMIP-PLSR Pipeline (Final recommended workflow for reproducible gene prioritization). |
| `/gmip_framework/` | Flexible GMIP Framework (Benchmarking of alternative strategies). |
| `/bin/`, `/modules/`, `/workflows/`, `/conf/` | Nextflow scripts, modules, workflows, and configurations. |

---

## 🚀 Quick Start (Using GMIP-PLSR Pipeline)

### 1. Install Nextflow
```bash
curl -s https://get.nextflow.io | bash
```

### 2. Clone this repository
```bash
git clone https://github.com/mohammedmsk/GMIP.git
cd GMIP/gmip_plsr_pipeline/
```

### 3. Run the main pipeline
```bash
nextflow run workflows/main.nf -profile standard
```
*(Adjust the profile according to your system: local, HPC, AWS, etc.)*

### 4. Inputs Required
- GWAS summary statistics
- Multi-omics feature files (e.g., PoPS-derived, scRNA-seq features)
- Configuration file (example in `/conf/`)

### 5. Outputs
- Re-prioritized gene lists
- Normalized heritability (Tau) scores

---

## 📂 Component Details

### GMIP-PLSR Pipeline
- Focused on handling feature multicollinearity using PLSR.
- Uses PoPS features mainly for gene prioritization but other can also be used.
- Validated using LOCO cross-validation and Benchmarker metrics.

### GMIP Framework
- Benchmarked multiple combinations of features, (PoPS-features, NAGA-features), ML models (NAGA, PoPS) and cross-validation strategies (noCV, k-Fold, LOCO).
- Used to derive key insights leading to GMIP-PLSR optimization.

---

## 📜 Citation

If you use GMIP or GMIP-PLSR in your work, please cite:

> Mohammed Shabbir Kanchwala, et al.  
> **GMIP-PLSR: A Reproducible Nextflow Pipeline for Post-GWAS Gene Prioritization via Multi-Omics Integration and Partial Least Squares Regression** (Manuscript in preparation, 2025).

GitHub repository: [https://github.com/mohammedmsk/GMIP](https://github.com/mohammedmsk/GMIP)

---

## 🛠️ Reproducibility and Environments

- Written in **Nextflow**.
- Compatible with local execution, HPC clusters, AWS Batch, and Azure.
- Example configurations provided in `/conf/`.
- Scripts and workflows use modular, reproducible design patterns.

---

## 📬 Contact

For questions, issues, or contributions, please contact:

**Mohammed Shabbir Kanchwala**  
Email: mohammedmsk@gmail.com

---

# ✅ Summary

- Main GMIP-PLSR pipeline: `/gmip_plsr_pipeline/`
- Benchmarking GMIP Framework: `/gmip_framework/`
- Fully reproducible, open-source, scalable pipelines for post-GWAS gene prioritization.

---
