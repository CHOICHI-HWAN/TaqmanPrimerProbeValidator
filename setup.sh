#!/usr/bin/env bash
# chichipes UI installer: one command that sets up the environment, installs the pathogen dashboards (panels)
# from GitHub, points them at chichipes, downloads their genome data and starts them.
#
#   ./setup.sh                                  # asks which panels; mamba environment; downloads data
#   ./setup.sh --panels lasv,ebov --yes         # no questions
#   ./setup.sh --backend docker                 # everything in a container (Windows, or no mamba wanted)
#
# options:
#   --backend auto|mamba|brew|docker   how to get Python, Snakemake, MAFFT, htslib (default auto = mamba; Windows: docker
#                                      or WSL, see below)
#   --panels LIST     comma list of panel ids, "all" (every panel except the very large ones) or "everything"
#   --dir PATH        install folder (default ~/chichipes-ui)
#   --port-base N     port of the front page; panels follow it (default 8760 -> panels 8765-8783)
#   --cores N         CPUs for the data download/analysis (default: half of this machine)
#   --primers FILE    primer sets exported with tools/export_primer_sets.py (an "edition" column + config/primers.tsv
#                     columns) -- copied into each installed panel
#   --no-data         install and start only; download data later with <dir>/fetch_data.sh
#   --no-start        do not start the web pages
#   --env NAME        conda environment name (default chichipes_ui)
#   --update          update an existing environment / panel checkouts
#   -y, --yes         never ask (defaults: all panels, data downloaded)
#
# Windows: the pipeline tools (Snakemake, MAFFT, htslib, pixy) do not run on Windows itself. Run this script inside
# WSL2 (Ubuntu) -- setup.ps1 prepares WSL -- or use --backend docker with Docker Desktop.
set -euo pipefail
HERE="$(cd "$(dirname "$0")" && pwd -P)"

BACKEND="auto"; PANELS=""; DIR="$HOME/chichipes-ui"; PORT_BASE=8760; CORES=""; PRIMERS=""
NO_DATA=0; NO_START=0; ENV_NAME="chichipes_ui"; UPDATE=0; YES=0
REPO_URL="${CHICHIPES_REPO_URL:-https://github.com/CHOICHI-HWAN/TaqmanPrimerProbeValidator.git}"
while [ $# -gt 0 ]; do
    case "$1" in
        --backend) BACKEND="$2"; shift 2 ;;
        --panels) PANELS="$2"; shift 2 ;;
        --dir) DIR="$2"; shift 2 ;;
        --port-base) PORT_BASE="$2"; shift 2 ;;
        --cores) CORES="$2"; shift 2 ;;
        --primers) PRIMERS="$2"; shift 2 ;;
        --no-data) NO_DATA=1; shift ;;
        --no-start) NO_START=1; shift ;;
        --env) ENV_NAME="$2"; shift 2 ;;
        --update) UPDATE=1; shift ;;
        -y|--yes) YES=1; shift ;;
        -h|--help) sed -n '2,29p' "$0" | sed 's/^# \{0,1\}//'; exit 0 ;;
        *) echo "[error] unknown option: $1 (see --help)" >&2; exit 2 ;;
    esac
done
[ -t 0 ] || YES=1
case "$BACKEND" in auto|mamba|brew|docker|container) ;; *) echo "[error] --backend must be auto, mamba, brew or docker" >&2; exit 2 ;; esac
[ -z "$PRIMERS" ] || PRIMERS="$(cd "$(dirname "$PRIMERS")" && pwd -P)/$(basename "$PRIMERS")"
mkdir -p "$DIR"
DIR="$(cd "$DIR" && pwd -P)"

# shellcheck source=install/common.sh
. "$HERE/install/common.sh"

OS="$(detect_os)"
say "chichipes UI 설치 — $OS, backend: $BACKEND, 설치 위치: $DIR"
case "$OS" in
    ubuntu|debian|linux|wsl) . "$HERE/install/setup_environment_for_ubuntu.sh" ;;
    mac) . "$HERE/install/setup_environment_for_mac.sh" ;;
    windows) . "$HERE/install/setup_environment_for_windows.sh" ;;
    *) die "지원하지 않는 OS입니다: $(uname -s)" ;;
esac

if [ "$BACKEND" = "docker" ]; then
    install_with_docker
    exit 0
fi

prepare_environment          # -> ENV_PREFIX (bin/python, bin/snakemake, mafft, bgzip, ...)
PY="$ENV_PREFIX/bin/python"
install_chichipes            # -> CHI
fetch_entry
choose_panels                # -> SELECTED
for id in $SELECTED; do install_panel "$id"; done
write_registry
[ -z "$PRIMERS" ] || import_primers "$PRIMERS"
write_launchers
[ "$NO_START" = 1 ] || "$DIR/start.sh"
if [ "$NO_DATA" = 1 ]; then
    say "데이터는 나중에 받으세요: $DIR/fetch_data.sh"
else
    start_data_fetch
fi
summary
