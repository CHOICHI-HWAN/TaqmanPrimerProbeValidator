#!/usr/bin/env bash
# Pack one built chichipes binary into a release archive in dist/ (the CI runs this once per platform).
#   cargo build --release && scripts/release.sh v0.5.0 linux amd64 target/release/chichipes
# Pushing a tag (git tag chichipes-v0.5.0 && git push origin chichipes-v0.5.0) makes GitHub Actions build every
# platform and attach the archives to a GitHub Release (.github/workflows/release.yml).
set -euo pipefail
VERSION="${1:?usage: scripts/release.sh vX.Y.Z OS ARCH BINARY}"
OS="${2:?os: linux, darwin or windows}"
ARCH="${3:?arch: amd64 or arm64}"
BIN="$(cd "$(dirname "${4:?path to the built binary}")" && pwd)/$(basename "$4")"
cd "$(dirname "$0")/.."
mkdir -p dist
name="chichipes-${VERSION}-${OS}-${ARCH}"
exe="chichipes"; [ "$OS" = windows ] && exe="chichipes.exe"
rm -rf "dist/$name" && mkdir -p "dist/$name"
cp "$BIN" "dist/$name/$exe"
cp README.md README.ko.md LICENSE THIRD_PARTY_NOTICES.md "dist/$name/"
cp -r setup.sh setup.ps1 install .gitattributes "dist/$name/"  # the UI installer: ./setup.sh next to the binary
if [ "$OS" = windows ]; then
  if command -v zip >/dev/null; then (cd dist && zip -qr "$name.zip" "$name"); else (cd dist && 7z a -tzip -bso0 "$name.zip" "$name"); fi
else
  tar -C dist -czf "dist/$name.tar.gz" "$name"
fi
rm -rf "dist/$name"
ls -la dist
