#!/bin/bash
# Builds the Auspice JSON (and the clade assignments) for Nipah locally with the Nextstrain CLI, from the
# official nextstrain/nipah workflow -- nothing is taken from the Auspice web service. The data inputs are the
# workflow's own defaults (data.nextstrain.org/files/workflows/nipah/{sequences.fasta,metadata.tsv}.zst).
#
# usage: scripts/run_nextstrain_nipah.sh [path/to/nextstrain_nipah_clone]
# needs: the Nextstrain env (snakemake, augur, nextclade, iqtree, mafft, snakemake-storage-plugin-http, nextstrain-cli)
#        built by `./setup.sh --nextstrain` (environment_nextstrain.yml); its path is read from config/config.yaml
#        (nextstrain_env_prefix) or $NEXTSTRAIN_ENV. Older machines fall back to /BI/conda/envs/nextstrain*.
set -euo pipefail
REPO="${1:-../nextstrain_nipah}"
[ -d "$REPO/phylogenetic" ] || git clone --depth 1 https://github.com/nextstrain/nipah.git "$REPO"
NS_ENV="${NEXTSTRAIN_ENV:-$(sed -n 's/^nextstrain_env_prefix:[[:space:]]*"\{0,1\}\([^"#]*\)"\{0,1\}.*/\1/p' config/config.yaml 2>/dev/null | head -n 1)}"
if [ -n "$NS_ENV" ] && [ -x "$NS_ENV/bin/nextstrain" ]; then
    export PATH="$NS_ENV/bin:$PATH"
else
    export PATH=/BI/conda/envs/nextstrain_nipah/bin:/BI/conda/envs/nextstrain/bin:$PATH
fi
command -v nextstrain >/dev/null || { echo "[error] nextstrain CLI not found -- run ./setup.sh --nextstrain" >&2; exit 1; }
(cd "$REPO/phylogenetic" && nextstrain build --native . --cores "${CORES:-16}" auspice/nipah_all.json)
mkdir -p data/nextstrain
cp "$REPO/phylogenetic/auspice/nipah_all.json" "$REPO/phylogenetic/results/all/clades.json" \
   "$REPO/phylogenetic/defaults/exclude.txt" "$REPO/phylogenetic/defaults/genemap.gff" "$REPO/phylogenetic/defaults/clades.tsv" data/nextstrain/
echo "nextstrain/nipah $(git -C "$REPO" rev-parse --short HEAD) built $(date +%F)" > data/nextstrain/PROVENANCE.txt
