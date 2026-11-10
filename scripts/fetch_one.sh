#!/bin/bash
# usage: fetch_one.sh <sequences_url> <metadata_url> <out_sequences.fasta> <out_metadata.tsv>
# (URLs come from the virus profile -- see scripts/virus_profile.py source_urls())
set -euo pipefail
seq_url="$1"; meta_url="$2"; out_seq="$3"; out_meta="$4"
tmp_seq_zst="${out_seq}.zst"
tmp_meta_zst="${out_meta}.zst"

curl -fsSL "${seq_url}" -o "$tmp_seq_zst"
zstd -f -d "$tmp_seq_zst" -o "$out_seq"
rm -f "$tmp_seq_zst"

curl -fsSL "${meta_url}" -o "$tmp_meta_zst"
zstd -f -d "$tmp_meta_zst" -o "$out_meta"
rm -f "$tmp_meta_zst"
