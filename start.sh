#!/bin/bash
# Starts the LASV primer/probe validation webapp against this project's already-computed
# results (config/, data/, output/) -- no Snakemake/conda env needed for this: server.py
# only reads precomputed TSVs/JSON and uses python3's stdlib (http.server etc.), so a copy
# of this whole project folder onto any machine with a plain python3 is enough to demo it.
# Re-running the pipeline itself (the "데이터 업데이트"/"분석 다시하기" buttons) still needs
# snakemake + mfeprimer, but just browsing/presenting the existing results does not.
#
# Usage:
#   ./start.sh [port]      # default port 8765
set -euo pipefail

PORT="${1:-8765}"

# Always run from this script's own directory, so it works after copying the whole project
# folder anywhere else (USB stick, another machine) regardless of where it's launched from.
cd "$(dirname "$(readlink -f "$0")")"

if ! command -v python3 >/dev/null 2>&1; then
    echo "[error] python3를 찾을 수 없습니다. Python 3.8 이상을 설치한 뒤 다시 실행하세요." >&2
    exit 1
fi

if [ ! -f "config/config.yaml" ]; then
    echo "[error] config/config.yaml이 없습니다 -- 이 스크립트는 프로젝트 폴더 안에서 실행해야 합니다." >&2
    exit 1
fi

if ! python3 - "$PORT" <<'EOF'
import socket
import sys
port = int(sys.argv[1])
s = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
try:
    s.bind(("127.0.0.1", port))
except OSError:
    sys.exit(1)
finally:
    s.close()
EOF
then
    echo "[error] 포트 ${PORT}이 이미 사용 중입니다. 다른 포트로 다시 시도하세요: ./start.sh <port>" >&2
    exit 1
fi

echo "LASV primer/probe validation 대시보드를 시작합니다..."
echo "  브라우저에서 다음 주소를 여세요: http://localhost:${PORT}"
echo "  종료하려면 이 터미널에서 Ctrl+C"
echo ""

# use the conda env that ./setup.sh recorded (snakemake, mafft, pixy ... on PATH for the jobs the web app launches)
ENV_PREFIX="$(sed -n 's/^env_prefix:[[:space:]]*"\{0,1\}\([^"#]*\)"\{0,1\}.*/\1/p' config/config.yaml | head -n 1)"
if [ -n "$ENV_PREFIX" ] && [ -d "$ENV_PREFIX/bin" ]; then
    export PATH="$ENV_PREFIX/bin:$PATH"
else
    echo "[주의] conda 환경이 설정되어 있지 않습니다 (./setup.sh). 이미 PATH에 있는 snakemake/mafft를 사용합니다." >&2
fi

exec python3 webapp/server.py --port "$PORT"
