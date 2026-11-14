#!/bin/bash
# Re-run every probe wide search of the pipeline with `pe search` and compare with mfeprimer's stored output.
set -uo pipefail
PE=${CHICHIPES_DIR:-$(cd "$(dirname "$0")/.." && pwd)}
CHI=${CHICHIPES_BIN:-$PE/target/release/chichipes}  # built with: cargo build --release
OUT=$PE/work/cmp_search
mkdir -p "$OUT"
: > "$OUT/summary.txt"
for ref in ${EDITIONS:?set EDITIONS to the folder holding the pipeline editions}/${1:-*}/output/primers/*/alignment/*_probe_wide_search.txt; do
    ed=$(echo "$ref" | sed -E 's#.*/editions/([^/]+)/.*#\1#'); set=$(basename "$ref" _probe_wide_search.txt)
    root=${EDITIONS:?set EDITIONS to the folder holding the pipeline editions}/$ed
    db=$(grep -m1 ' Cmd: ' "$(dirname "$ref")/../primer_qc/${set}_mm2.txt.mfe.log" | grep -o ' -d [^ ]*' | cut -d' ' -f3)
    mm=$(sed -n 's/# max allowed mismatch: //p' "$ref")
    ( cd "$root" && $CHI search "$db" "${ref%.txt}.fa" "$OUT/${ed}_${set}.txt" -m "$mm" ) 2>/dev/null
    echo "### ${ed}_${set} $(python3 -I $PE/harness/compare_search.py "$ref" "$OUT/${ed}_${set}.txt" --quiet)" >> "$OUT/summary.txt"
done
