#!/usr/bin/env bash
# Downloads MAGMA v1.06b (Linux x86_64 static binary) from the official source
# and places it at bin/magma, ready for the Nextflow pipeline.
#
# MAGMA is free for non-commercial use. By running this script you agree to
# the MAGMA license terms at https://ctg.cncr.nl/software/magma
#
# Usage (run from the repo root):
#   bash setup_magma.sh

set -euo pipefail

MAGMA_VERSION="1.06b"
MAGMA_URL="https://ctg.cncr.nl/software/MAGMA/prog/magma_v${MAGMA_VERSION}_static.zip"
DEST="bin/magma"
TMP_DIR="$(mktemp -d)"

# Resolve repo root regardless of where the script is called from
REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
DEST="${REPO_ROOT}/bin/magma"

if [[ -x "${DEST}" ]]; then
    echo "MAGMA binary already exists at ${DEST} — skipping download."
    exit 0
fi

echo "Downloading MAGMA v${MAGMA_VERSION} from ${MAGMA_URL} ..."
curl -fsSL --retry 3 -o "${TMP_DIR}/magma.zip" "${MAGMA_URL}"

echo "Extracting ..."
unzip -q "${TMP_DIR}/magma.zip" -d "${TMP_DIR}"

BINARY="$(find "${TMP_DIR}" -type f -name "magma" | head -1)"
if [[ -z "${BINARY}" ]]; then
    echo "ERROR: could not find 'magma' binary in the downloaded zip." >&2
    echo "       Please download manually from https://ctg.cncr.nl/software/magma" >&2
    echo "       and place the binary at: ${DEST}" >&2
    rm -rf "${TMP_DIR}"
    exit 1
fi

mkdir -p "${REPO_ROOT}/bin"
cp "${BINARY}" "${DEST}"
chmod +x "${DEST}"
rm -rf "${TMP_DIR}"

echo "Done. MAGMA v${MAGMA_VERSION} installed at: ${DEST}"
echo ""
echo "By using MAGMA you agree to its license terms:"
echo "  https://ctg.cncr.nl/software/magma"
