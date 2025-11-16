#!/bin/bash

# Build script for Wakhan Docker image (linux/amd64)
# Usage: ./build.sh [image_name] [registry] [--push]
# Example: ./build.sh wakhan
# Example: ./build.sh wakhan ghcr.io/shahcompbio
# Example: ./build.sh wakhan ghcr.io/shahcompbio --push

set -e

IMAGE_NAME=${1:-wakhan}
REGISTRY=${2:-}
PUSH_FLAG=${3:-}
BUILD_DATE=$(date +%Y%m%d)
PLATFORM="linux/amd64"

# Add registry prefix if provided
if [ -n "$REGISTRY" ]; then
    FULL_IMAGE_NAME="${REGISTRY}/${IMAGE_NAME}"
else
    FULL_IMAGE_NAME="${IMAGE_NAME}"
fi

echo "Building Docker image for ${PLATFORM}..."

# Build to extract the commit hash
echo "Building temporary image to extract version info..."
docker build --platform ${PLATFORM} -t ${FULL_IMAGE_NAME}:temp .

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

if [ "$PUSH_FLAG" == "--push" ]; then
    echo ""
    echo "Pushing to registry..."
    docker push ${FULL_IMAGE_NAME}:${VERSION_TAG}
    docker push ${FULL_IMAGE_NAME}:${COMMIT_HASH}
    docker push ${FULL_IMAGE_NAME}:latest
    echo ""
    echo "✓ Push complete!"
else
    echo ""
    echo "To push to registry, run:"
    echo "  docker push ${FULL_IMAGE_NAME}:${VERSION_TAG}"
    echo "  docker push ${FULL_IMAGE_NAME}:${COMMIT_HASH}"
    echo "  docker push ${FULL_IMAGE_NAME}:latest"
    echo ""
    echo "Or rebuild with --push flag:"
    echo "  ./build.sh ${IMAGE_NAME}${REGISTRY:+ $REGISTRY} --push"
fi
