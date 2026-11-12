#!/bin/bash
# Starts the front-door page for the LASV / Nipah validation dashboards.
#   ./start.sh [port]      # default 8760 (the dashboards' "전체 파이프라인 목록" link points here)
set -euo pipefail
cd "$(dirname "$(readlink -f "$0")")"
echo "파이프라인 목록: http://localhost:${1:-8760}   (종료: Ctrl+C)"
exec python3 server.py --port "${1:-8760}"
