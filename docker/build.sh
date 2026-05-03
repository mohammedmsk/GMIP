#!/usr/bin/env bash
# Build and optionally push GMIP Docker images.
# Run from the repo root:  bash docker/build.sh [--push] [--registry <registry>]
#
# Examples:
#   bash docker/build.sh
#   bash docker/build.sh --push --registry ghcr.io/mohammedmsk

set -euo pipefail

REGISTRY=""
PUSH=false
VERSION="${GMIP_VERSION:-latest}"

while [[ $# -gt 0 ]]; do
    case "$1" in
        --push)     PUSH=true; shift ;;
        --registry) REGISTRY="$2"; shift 2 ;;
        --version)  VERSION="$2"; shift 2 ;;
        *) echo "Unknown arg: $1"; exit 1 ;;
    esac
done

tag_and_push() {
    local name="$1"
    local dockerfile="$2"
    local tag="${REGISTRY:+${REGISTRY}/}${name}:${VERSION}"

    echo "==> Building ${tag} from ${dockerfile}"
    docker build -f "${dockerfile}" -t "${tag}" .

    if [[ "${PUSH}" == "true" ]]; then
        echo "==> Pushing ${tag}"
        docker push "${tag}"
    fi
}

tag_and_push "gmip-plsr"   "docker/Dockerfile.gmip"
tag_and_push "gmip-ldsc"   "docker/Dockerfile.gmip_ldsc"

echo ""
echo "Done. To convert to Apptainer .sif files:"
echo "  apptainer build gmip-plsr.sif  docker-daemon://gmip-plsr:${VERSION}"
echo "  apptainer build gmip-ldsc.sif  docker-daemon://gmip-ldsc:${VERSION}"
