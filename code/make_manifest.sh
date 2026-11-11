#!/bin/bash
# Writes MANIFEST.sha256: the SHA-256 checksums of all the files of the repository except the manifest itself.
# Check it with: sha256sum -c MANIFEST.sha256
set -eu
cd "$(dirname "$0")/.."
git ls-files | grep -v '^MANIFEST.sha256$' | LC_ALL=C sort | tr '\n' '\0' | xargs -0 sha256sum > MANIFEST.sha256
