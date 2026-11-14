#!/bin/bash
# Re-run every mfeprimer primer search of the pipeline (found via its .mfe.log) with primer-engine and compare.
# usage: run_all.sh [glob of editions]   -> work/cmp/<edition>_<set>.txt and work/cmp/summary.txt
set -uo pipefail
PE=${CHICHIPES_DIR:-$(cd "$(dirname "$0")/.." && pwd)}
OUT=$PE/work/cmp
mkdir -p "$OUT"
: > "$OUT/summary.txt"
for log in ${EDITIONS:?set EDITIONS to the folder holding the pipeline editions}/${1:-*}/output/primers/*/primer_qc/*_mm*.txt.mfe.log; do
    dir=$(dirname "$log"); ed=$(echo "$log" | sed -E 's#.*/editions/([^/]+)/.*#\1#')
    ref=${log%.mfe.log}; tag="${ed}_$(basename "$ref" .txt)"
    cmd=$(grep -m1 ' Cmd: ' "$log" | sed 's/.* Cmd: //')
    root=$(grep -m1 ' Dir: ' "$log" | sed 's/.* Dir: //')
    args=$(echo "$cmd" | cut -d' ' -f2- | sed -E "s# -o [^ ]+# -o $OUT/$tag.txt#; s# --cpu [0-9]+##")
    ( cd "$root" && $PE/chichipes $args --cpu 24 ) > "$OUT/$tag.log" 2>&1
    echo "### $tag" >> "$OUT/summary.txt"
    python3 -I $PE/harness/compare_spec.py "$ref" "$OUT/$tag.txt" --show 2 >> "$OUT/summary.txt" 2>&1
done
