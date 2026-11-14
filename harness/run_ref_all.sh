#!/bin/bash
# Re-run the reference-genome searches kept in work/cmp_ref (<ed>_<ref>_<set>_<query>.mfe.txt from mfeprimer search
# -m 6) with `pe search` and compare. usage: run_ref_all.sh  -> work/cmp_ref/summary.txt
set -uo pipefail
PE=${CHICHIPES_DIR:-$(cd "$(dirname "$0")/.." && pwd)}
CHI=${CHICHIPES_BIN:-$PE/target/release/chichipes}  # built with: cargo build --release
ED=${EDITIONS:?set EDITIONS to the folder holding the pipeline editions}
OUT=$PE/work/cmp_ref
: > "$OUT/summary.txt"
for ref in "$OUT"/*.mfe.txt; do
    b=$(basename "$ref" .mfe.txt)
    [[ $b =~ ^([^_]+)_(.+)_(P[0-9]+)_(FR\.fa|probe_wide_search\.fa)$ ]] || { echo "### $b unparsed" >> "$OUT/summary.txt"; continue; }
    ed=${BASH_REMATCH[1]}; db=data/refgenome/${BASH_REMATCH[2]}; set=${BASH_REMATCH[3]}; q=${BASH_REMATCH[4]}
    if [ "$q" = FR.fa ]; then qp=output/primers/$set/primer_qc/${set}_FR.fa; else qp=output/primers/$set/alignment/${set}_probe_wide_search.fa; fi
    ( cd "$ED/$ed" && $CHI search "$db" "$qp" "$OUT/$b.pe.txt" -m 6 ) 2>/dev/null
    echo "### $b $(python3 -I $PE/harness/compare_search.py "$ref" "$OUT/$b.pe.txt" --quiet)" >> "$OUT/summary.txt"
done
