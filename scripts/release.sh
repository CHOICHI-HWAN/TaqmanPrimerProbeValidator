#!/usr/bin/env bash
# Build release archives for every supported platform into dist/.
#   scripts/release.sh v0.1.0
# Pushing a tag (git tag chichipes-v0.1.0 && git push origin chichipes-v0.1.0) makes GitHub Actions run this same build
# and attach the archives to a GitHub Release (.github/workflows/release.yml).
set -euo pipefail
VERSION="${1:?usage: scripts/release.sh vX.Y.Z}"
cd "$(dirname "$0")/.."
rm -rf dist && mkdir -p dist
for target in linux/amd64 linux/arm64 darwin/amd64 darwin/arm64 windows/amd64; do
  os="${target%/*}"; arch="${target#*/}"
  name="chichipes-${VERSION}-${os}-${arch}"
  exe="chichipes"; [ "$os" = windows ] && exe="chichipes.exe"
  mkdir -p "dist/$name"
  CGO_ENABLED=0 GOOS="$os" GOARCH="$arch" go build -trimpath -ldflags "-s -w -X main.Version=${VERSION}" \
    -o "dist/$name/$exe" ./cmd/chichipes
  cp README.md README.ko.md LICENSE THIRD_PARTY_NOTICES.md "dist/$name/"
  cp -r setup.sh setup.ps1 install .gitattributes "dist/$name/"  # the UI installer: ./setup.sh next to the binary
  if [ "$os" = windows ]; then
    (cd dist && zip -qr "$name.zip" "$name")
  else
    tar -C dist -czf "dist/$name.tar.gz" "$name"
  fi
  rm -rf "dist/$name"
done
(cd dist && sha256sum chichipes-* > SHA256SUMS.txt)
ls -la dist
