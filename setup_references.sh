#!/usr/bin/env bash
# Downloads and unpacks the GMIP reference files from Zenodo.
# Supports local destinations and S3 buckets.
#
# Usage:
#   bash setup_references.sh                          # extract to ./refdir (default)
#   bash setup_references.sh --refdir /data/refdir   # custom local path
#   bash setup_references.sh --refdir s3://my-bucket/refdir  # upload to S3
#
# S3 requirements:
#   - AWS CLI installed and configured (aws configure  or  IAM instance role)
#   - Sufficient IAM permissions: s3:PutObject on the target bucket
#   - Set --aws-region if your bucket is not in us-east-1
#
# After setup, run the pipeline with:
#   nextflow run gmip_plsr_pipeline/workflows/GMIP_default_wf.nf \
#     -c gmip_plsr_pipeline/conf/nextflow.config \
#     -profile local \
#     --magma_input /path/to/gwas.magma_input.tsv \
#     --munged_gwas  /path/to/gwas.sumstats.gz \
#     --prefix MY_GWAS
#
# For S3-backed runs also pass:  --base_refdir s3://my-bucket/refdir
# and use:  -work-dir s3://my-bucket/work  (required for aws_batch profile)

set -euo pipefail

# ── Defaults ──────────────────────────────────────────────────────────────────
REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
DEFAULT_REFDIR="${REPO_ROOT}/refdir"
REFDIR="${DEFAULT_REFDIR}"
AWS_REGION="us-east-1"

# Zenodo record — update this URL after uploading the reference archive
ZENODO_URL="https://zenodo.org/records/XXXXXXX/files/gmip_refdir.zip"
ZENODO_MD5=""   # optional: set to expected md5 for verification

# ── Argument parsing ──────────────────────────────────────────────────────────
while [[ $# -gt 0 ]]; do
  case "$1" in
    --refdir)     REFDIR="$2";     shift 2 ;;
    --aws-region) AWS_REGION="$2"; shift 2 ;;
    -h|--help)
      sed -n '2,/^set /p' "${BASH_SOURCE[0]}" | grep '^#' | sed 's/^# \?//'
      exit 0 ;;
    *) echo "Unknown argument: $1" >&2; exit 1 ;;
  esac
done

IS_S3=false
[[ "${REFDIR}" == s3://* ]] && IS_S3=true

# ── Helpers ───────────────────────────────────────────────────────────────────
check_cmd() { command -v "$1" &>/dev/null || { echo "ERROR: '$1' not found. Please install it." >&2; exit 1; }; }

# ── Pre-flight checks ─────────────────────────────────────────────────────────
check_cmd curl
check_cmd unzip
if [[ "${IS_S3}" == true ]]; then
  check_cmd aws
  aws s3 ls "${REFDIR%/}/" &>/dev/null && {
    echo "Reference files already present at ${REFDIR} — skipping download."
    exit 0
  } || true
else
  if [[ -d "${REFDIR}" && -n "$(ls -A "${REFDIR}" 2>/dev/null)" ]]; then
    echo "Reference directory already exists and is non-empty: ${REFDIR}"
    echo "Delete it or choose a different --refdir to re-download."
    exit 0
  fi
fi

# Check Zenodo URL is set
if [[ "${ZENODO_URL}" == *"XXXXXXX"* ]]; then
  echo "ERROR: Zenodo record URL not yet set in this script." >&2
  echo "       Please see: https://zenodo.org/records/XXXXXXX" >&2
  echo "       Or contact the authors for access to the reference files." >&2
  exit 1
fi

# ── Download ──────────────────────────────────────────────────────────────────
TMP_DIR="$(mktemp -d)"
ZIP="${TMP_DIR}/gmip_refdir.zip"

echo "Downloading GMIP reference files (~5 GB) from Zenodo..."
echo "  ${ZENODO_URL}"
curl -fL --retry 3 --progress-bar -o "${ZIP}" "${ZENODO_URL}"

# Optional MD5 check
if [[ -n "${ZENODO_MD5}" ]]; then
  echo "Verifying checksum..."
  echo "${ZENODO_MD5}  ${ZIP}" | md5sum -c -
fi

# ── Extract ───────────────────────────────────────────────────────────────────
echo "Extracting (~13 GB uncompressed)..."
unzip -q "${ZIP}" -d "${TMP_DIR}/extracted"

# The zip contains a top-level 'refdir/' folder — find it
EXTRACTED="$(find "${TMP_DIR}/extracted" -maxdepth 1 -type d | tail -1)"

# ── Install ───────────────────────────────────────────────────────────────────
if [[ "${IS_S3}" == true ]]; then
  echo "Uploading to ${REFDIR} ..."
  aws s3 sync "${EXTRACTED}/" "${REFDIR%/}/" \
    --region "${AWS_REGION}" \
    --no-progress
  echo "Done. Reference files are at: ${REFDIR}"
  echo ""
  echo "Run the pipeline with:"
  echo "  nextflow run gmip_plsr_pipeline/workflows/GMIP_default_wf.nf \\"
  echo "    -c gmip_plsr_pipeline/conf/nextflow.config \\"
  echo "    -profile aws_batch \\"
  echo "    -work-dir s3://YOUR_BUCKET/work \\"
  echo "    --base_refdir ${REFDIR} \\"
  echo "    --magma_input s3://YOUR_BUCKET/gwas.magma_input.tsv \\"
  echo "    --munged_gwas  s3://YOUR_BUCKET/gwas.sumstats.gz \\"
  echo "    --prefix MY_GWAS"
else
  mkdir -p "${REFDIR}"
  mv "${EXTRACTED}"/* "${REFDIR}/"
  echo "Done. Reference files are at: ${REFDIR}"
  echo ""
  echo "Run the pipeline with:"
  echo "  nextflow run gmip_plsr_pipeline/workflows/GMIP_default_wf.nf \\"
  echo "    -c gmip_plsr_pipeline/conf/nextflow.config \\"
  echo "    -profile local \\"
  echo "    --base_refdir ${REFDIR} \\"
  echo "    --magma_input /path/to/gwas.magma_input.tsv \\"
  echo "    --munged_gwas  /path/to/gwas.sumstats.gz \\"
  echo "    --prefix MY_GWAS"
  if [[ "${REFDIR}" == "${DEFAULT_REFDIR}" ]]; then
    echo ""
    echo "  (--base_refdir can be omitted when using the default location)"
  fi
fi

rm -rf "${TMP_DIR}"
