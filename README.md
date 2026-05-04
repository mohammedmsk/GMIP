# GMIP-PLSR: Post-GWAS Gene Prioritization via Multi-Omics Integration and Partial Least Squares Regression

GMIP-PLSR is a reproducible [Nextflow](https://www.nextflow.io/) pipeline for post-GWAS gene prioritization. It integrates GWAS summary statistics with multi-omics features (gene expression, protein-protein interactions, biological pathways) and applies Partial Least Squares Regression (PLSR) to handle feature multicollinearity and improve prioritization performance.

This repository provides two components:

- **`gmip_plsr_pipeline/`** — The recommended pipeline. Runs MAGMA → PoPS → GMIP-PLSR → Benchmarker.
- **`gmip_framework/`** — The earlier benchmarking framework used to evaluate feature sets, ML models, and cross-validation strategies that informed GMIP-PLSR design.

**Preprint:** [doi:10.64898/2026.04.06.716845](https://doi.org/10.64898/2026.04.06.716845)
**Reference data:** [doi:10.5281/zenodo.19986368](https://doi.org/10.5281/zenodo.19986368)

---

## Requirements

| Tool | Version | Notes |
|------|---------|-------|
| [Nextflow](https://nextflow.io) | ≥ 23.04 | Requires Java 11+ |
| [Conda](https://conda.io) / [Mamba](https://github.com/mamba-org/mamba) | any | for `local` and `slurm` profiles |
| [Docker](https://www.docker.com) | any | for `local,docker` profile |
| [Singularity](https://sylabs.io) / [Apptainer](https://apptainer.org) | any | for `slurm` profile (default) |
| AWS CLI | ≥ 2 | for `aws_batch` profile and S3 reference storage |
| MAGMA | v1.06b | downloaded separately — see Setup below |

---

## Quick Start

```bash
# 1. Clone
git clone https://github.com/mohammedmsk/GMIP.git
cd GMIP

# 2. Download MAGMA binary (review license at https://ctg.cncr.nl/software/magma)
bash setup_magma.sh

# 3. Download reference files (~5 GB download, ~13 GB on disk)
bash setup_references.sh

# 4. Run
nextflow run gmip_plsr_pipeline/workflows/GMIP_default_wf.nf \
  -c gmip_plsr_pipeline/conf/nextflow.config \
  -profile local \
  --magma_input /path/to/gwas.magma_input.tsv \
  --munged_gwas  /path/to/gwas.sumstats.gz \
  --prefix        MY_GWAS
```

Results are written to `results/` by default. Override with `--outdir /path/to/outdir`.

---

## Setup

### MAGMA binary

MAGMA v1.06b is required but cannot be redistributed here. Download it with:

```bash
bash setup_magma.sh
```

This fetches the Linux x86_64 static binary from the [official MAGMA site](https://ctg.cncr.nl/software/magma) and places it at `bin/magma`. Review the MAGMA license terms before use.

### Reference files

Reference files (~5 GB compressed, ~13 GB uncompressed) are hosted on Zenodo:

```bash
# Default: extract to ./refdir/  (gitignored)
bash setup_references.sh

# Custom local path
bash setup_references.sh --refdir /data/gmip_refdir

# Upload to S3 (for AWS Batch runs)
bash setup_references.sh --refdir s3://my-bucket/gmip_refdir --aws-region us-east-1
```

Reference files include:
- 1000 Genomes Phase 3 EUR LD reference (for MAGMA)
- Full PoPS feature matrices
- S-LDSC baseline model, plink files, frequency files (for Benchmarker)
- Pre-computed gene window LD score files

### Conda environments (optional — only needed for `local` or `slurm,conda` profiles)

Environments are created automatically by Nextflow on first run using the YAMLs in `docker/`. To pre-build them manually:

```bash
conda env create -f docker/environment_gmip.yml
conda env create -f docker/environment_gmip_ldsc.yml
```

---

## Execution Profiles

Select a profile with `-profile <name>`. Profiles can be combined with a comma (e.g., `-profile slurm,conda`).

### Primary profiles

| Profile | Executor | Environment | Typical use |
|---------|----------|-------------|-------------|
| `local` | local | conda | Laptop / workstation |
| `slurm` | SLURM | Singularity | HPC cluster |
| `aws_batch` | AWS Batch | Docker | Cloud |

### Environment overrides

| Profile | Effect |
|---------|--------|
| `conda` | Force conda (combine with `local` or `slurm`) |
| `docker` | Force Docker (combine with `local`) |
| `singularity` | Force Singularity (combine with `slurm`) |

### AWS Batch setup

1. Upload reference files to S3: `bash setup_references.sh --refdir s3://my-bucket/gmip_refdir`
2. Set your queue in `gmip_plsr_pipeline/conf/nextflow.config` (`process.queue`) or pass `--process.queue YOUR_QUEUE`
3. Run:

```bash
nextflow run gmip_plsr_pipeline/workflows/GMIP_default_wf.nf \
  -c gmip_plsr_pipeline/conf/nextflow.config \
  -profile aws_batch \
  -work-dir s3://my-bucket/work \
  --base_refdir  s3://my-bucket/gmip_refdir \
  --magma_input  s3://my-bucket/gwas.magma_input.tsv \
  --munged_gwas  s3://my-bucket/gwas.sumstats.gz \
  --prefix        MY_GWAS
```

---

## Inputs

| Parameter | Description |
|-----------|-------------|
| `--magma_input` | GWAS summary statistics in MAGMA input format (`.tsv`, columns: `SNP`, `CHR`, `BP`, `P`, `N`) |
| `--munged_gwas` | Munged GWAS summary statistics for the Benchmarker (`.sumstats.gz`, LDSC format) |
| `--prefix` | Output file prefix (e.g., trait name) |
| `--base_refdir` | Path to reference files directory (default: `./refdir`) |
| `--outdir` | Output directory (default: `results`) |
| `--gmip_method` | PLSR method string (default: `PLSRegression_nc3`) |
| `--collapse_loci` | Collapse neighbouring genes to loci before PoPS (default: `false`) |

---

## Outputs

```
results/
├── 1_magma/              MAGMA gene-level p-values and raw scores
├── 2_pops/               PoPS per-chromosome predictions
├── 3_gmip/               GMIP-PLSR reprioritized gene scores
├── 4_bm_results/         Benchmarker S-LDSC normalized tau scores
│   ├── part_1/           Per-chromosome LD score annotation files
│   └── part_2/           Partition heritability results
└── report.html           Nextflow execution report
```

The key output is `3_gmip/<prefix>.gmip_*.txt` — a ranked gene list with GMIP-PLSR scores.

---

## Pipeline Overview

```
GWAS summary stats
      │
      ▼
   MAGMA          SNP-to-gene mapping using 1000G EUR LD reference
      │
      ▼
    PoPS           Per-chromosome ridge regression over multi-omics features
      │
      ▼
  GMIP-PLSR        PLSR across chromosomes to handle feature multicollinearity
      │
      ▼
 Benchmarker       S-LDSC partitioned heritability to evaluate reprioritization
```

### Methods

- **MAGMA** ([de Leeuw et al. 2015](https://doi.org/10.1371/journal.pcbi.1004219)) maps SNPs to genes using a window-based approach with LD from 1000 Genomes Phase 3 EUR.
- **PoPS** ([Weeks et al. 2023](https://doi.org/10.1038/s41588-023-01443-6)) scores genes by ridge regression over multi-omics features.
- **GMIP-PLSR** applies Partial Least Squares Regression across PoPS outputs from all chromosomes, exploiting LOCO cross-validation to mitigate multicollinearity.
- **Benchmarker** evaluates reprioritized gene lists via S-LDSC normalized tau scores.

---

## Repository Structure

```
GMIP/
├── gmip_plsr_pipeline/
│   ├── workflows/          Main Nextflow workflow
│   ├── modules/            Process definitions (MAGMA, PoPS, GMIP, Benchmarker)
│   ├── subworkflows/       Benchmarker subworkflow
│   ├── bin/                Pipeline scripts (pops.py, gmip.py, benchmarker/, ldsc/)
│   └── conf/               Nextflow configs
├── gmip_framework/
│   ├── workflows/          Framework workflows (benchmarking variants)
│   ├── modules/            Additional modules
│   ├── bin/                Framework scripts
│   └── conf/               Feature-wise configs (full_pops, ppi_pops, etc.)
├── docker/
│   ├── Dockerfile.gmip         Main pipeline Docker image (Python 3.12)
│   ├── Dockerfile.gmip_ldsc    Benchmarker Docker image (Python 2.7 + LDSC)
│   ├── environment_gmip.yml    Conda environment spec
│   └── environment_gmip_ldsc.yml
├── setup_magma.sh          Download MAGMA binary
├── setup_references.sh     Download / upload reference files
└── LICENSE
```

---

## Citation

If you use GMIP-PLSR in your work, please cite:

> Kanchwala MS, et al. *GMIP-PLSR: A Reproducible Nextflow Pipeline for Post-GWAS Gene Prioritization via Multi-Omics Integration and Partial Least Squares Regression.* bioRxiv (2026). doi:[10.64898/2026.04.06.716845](https://doi.org/10.64898/2026.04.06.716845)

---

## License

MIT License — see [LICENSE](LICENSE).
Copyright © 2026 The University of Texas System.

---

## Contact

Mohammed Shabbir Kanchwala
University of Texas at Dallas
[mohammedmsk@gmail.com](mailto:mohammedmsk@gmail.com)

Issues and pull requests welcome via [GitHub](https://github.com/mohammedmsk/GMIP/issues).
