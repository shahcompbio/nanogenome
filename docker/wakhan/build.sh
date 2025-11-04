#!/bin/bash

# Build script for Wakhan Docker image (linux/amd64)
# Usage: ./build.sh [image_name] [registry]
# Example: ./build.sh wakhan
# Example: ./build.sh wakhan ghcr.io/shahcompbio

set -e

IMAGE_NAME=${1:-wakhan}
REGISTRY=${2:-}
BUILD_DATE=$(date +%Y%m%d)

# Add registry prefix if provided
if [ -n "$REGISTRY" ]; then
    FULL_IMAGE_NAME="${REGISTRY}/${IMAGE_NAME}"
else
    FULL_IMAGE_NAME="${IMAGE_NAME}"
fi

echo "Building Docker image for linux/amd64..."
docker build --platform linux/amd64 -t ${FULL_IMAGE_NAME}:temp .

echo "Extracting Wakhan commit hash..."
COMMIT_HASH=$(docker run --rm --platform linux/amd64 ${FULL_IMAGE_NAME}:temp cat /opt/wakhan_version.txt)

echo "Wakhan commit: ${COMMIT_HASH}"

# Create version tag
VERSION_TAG="${BUILD_DATE}-${COMMIT_HASH}"

echo "Tagging image as:"
echo "  - ${FULL_IMAGE_NAME}:${VERSION_TAG}"
echo "  - ${FULL_IMAGE_NAME}:${COMMIT_HASH}"
echo "  - ${FULL_IMAGE_NAME}:latest"

docker tag ${FULL_IMAGE_NAME}:temp ${FULL_IMAGE_NAME}:${VERSION_TAG}
docker tag ${FULL_IMAGE_NAME}:temp ${FULL_IMAGE_NAME}:${COMMIT_HASH}
docker tag ${FULL_IMAGE_NAME}:temp ${FULL_IMAGE_NAME}:latest

# Remove temp tag
docker rmi ${FULL_IMAGE_NAME}:temp

echo ""
echo "✓ Build complete!"
echo ""
echo "Available tags:"
docker images ${FULL_IMAGE_NAME} --format "  - {{.Repository}}:{{.Tag}}" | head -3
echo ""
echo "To push to registry:"
echo "  docker push ${FULL_IMAGE_NAME}:${VERSION_TAG}"
echo "  docker push ${FULL_IMAGE_NAME}:${COMMIT_HASH}"
echo "  docker push ${FULL_IMAGE_NAME}:latest"
