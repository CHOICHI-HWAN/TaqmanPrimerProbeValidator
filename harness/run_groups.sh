#!/bin/bash
# Pool each group's member sets' F/R (ids SET::F / SET::R, as group_cross_binding.py does) into one search, run
# mfeprimer and the engine, compare. usage: run_groups.sh <edition> "<G> <P..>" ...   -> work/groups/<ed>_<G>.*
set -uo pipefail
PE=${CHICHIPES_DIR:-$(cd "$(dirname "$0")/.." && pwd)}
CHI=${CHICHIPES_BIN:-$PE/target/release/chichipes}  # built with: cargo build --release
MFE=${MFEPRIMER:?set MFEPRIMER to the mfeprimer 4.2.4 binary}
ED=${EDITIONS:?set EDITIONS to the folder holding the pipeline editions}/$1; shift
OUT=$PE/work/groups; mkdir -p "$OUT"; chmod 700 "$OUT"
for spec in "$@"; do
    read -r g sets <<< "$spec"
    tag=$(basename "$ED")_$g; fa="$OUT/$tag.FR.fa"; : > "$fa"
    first=${sets%% *}
    db=$(grep -m1 ' Cmd: ' "$ED/output/primers/$first/primer_qc/${first}_mm2.txt.mfe.log" | grep -o ' -d [^ ]*' | cut -d' ' -f3)
    for s in $sets; do
        awk -v s="$s" '/^>/{sub(/^>[^_]*_/, ">" s "::"); print; next} {print}' "$ED/output/primers/$s/primer_qc/${s}_FR.fa" >> "$fa"
    done
    rm -f "${OUT:?}/${tag:?}".mfe.txt* "${OUT:?}/${tag:?}".pe.txt*
    ( cd "$ED" && /usr/bin/time -f "mfe %es" "$MFE" -d "$db" -i "$fa" -o "$OUT/$tag.mfe.txt" --fasta --virus --misMatch 2 --cpu 24 >/dev/null ) 2> "$OUT/$tag.time"
    ( cd "$ED" && /usr/bin/time -f "chichipes %es" "$CHI" -d "$db" -i "$fa" -o "$OUT/$tag.pe.txt" --fasta --virus --misMatch 2 --cpu 24 >/dev/null ) 2>> "$OUT/$tag.time"
    echo "### $tag ($sets) db=$db $(grep -h -o '[a-z]* [0-9.]*s$' "$OUT/$tag.time" | tr '\n' ' ')"
    python3 -I "$PE/harness/compare_spec.py" "$OUT/$tag.mfe.txt" "$OUT/$tag.pe.txt" --show 0
done
