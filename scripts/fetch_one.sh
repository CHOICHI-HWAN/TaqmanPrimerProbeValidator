#!/bin/bash
# usage: fetch_one.sh <seg: s|l> <base_url> <out_sequences.fasta> <out_metadata.tsv>
set -euo pipefail
seg="$1"; base_url="$2"; out_seq="$3"; out_meta="$4"
tmp_seq_zst="${out_seq}.zst"
tmp_meta_zst="${out_meta}.zst"

curl -fsSL "${base_url}/${seg}/sequences.fasta.zst" -o "$tmp_seq_zst"
zstd -f -d "$tmp_seq_zst" -o "$out_seq"
rm -f "$tmp_seq_zst"

curl -fsSL "${base_url}/${seg}/metadata.tsv.zst" -o "$tmp_meta_zst"
zstd -f -d "$tmp_meta_zst" -o "$out_meta"
rm -f "$tmp_meta_zst"
