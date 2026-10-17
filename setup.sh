#!/bin/bash
# First-time setup: builds the conda/mamba environment (snakemake, mafft, pixy, htslib, ...), asks for the
# mfeprimer executable path, and records both in config/config.yaml so the pipeline and the web app use them.
#
#   ./setup.sh                          # interactive
#   ./setup.sh --mfeprimer /path/to/mfeprimer-4.2.4-linux-amd64 --yes
#
# options:
#   --env NAME           conda env name (default: primer_validation)
#   --mfeprimer PATH     mfeprimer executable (otherwise $MFEPRIMER_BIN, the current config value, an
#                        auto-discovered ../*/mfeprimer/ binary, or you are asked)
#   --update             update the env if it already exists (default: leave it alone)
#   --nextstrain         also build the Nextstrain env used by scripts/run_nextstrain_nipah.sh (Nipah clades)
#   --no-nextstrain      skip that question
#   -y, --yes            never ask; take the defaults
#
# Safe to re-run: an existing env is kept, the config lines are simply rewritten. The mfeprimer path can also
# be changed later on the web app's "Analysis Parameters" page.
set -euo pipefail
cd "$(dirname "$(readlink -f "$0")")"

ENV_NAME="primer_validation"; NS_ENV_NAME="${ENV_NAME}_nextstrain"
MFE_ARG=""; UPDATE=0; NS_MODE="ask"; YES=0
while [ $# -gt 0 ]; do
    case "$1" in
        --env) ENV_NAME="$2"; NS_ENV_NAME="${2}_nextstrain"; shift 2 ;;
        --mfeprimer) MFE_ARG="$2"; shift 2 ;;
        --update) UPDATE=1; shift ;;
        --nextstrain) NS_MODE="yes"; shift ;;
        --no-nextstrain) NS_MODE="no"; shift ;;
        -y|--yes) YES=1; shift ;;
        -h|--help) sed -n '2,19p' "$0" | sed 's/^# \{0,1\}//'; exit 0 ;;
        *) echo "[error] 알 수 없는 옵션: $1 (--help 참고)" >&2; exit 2 ;;
    esac
done
[ -t 0 ] || YES=1                         # no terminal to ask on -> behave like --yes

say()  { printf '\n\033[1m== %s\033[0m\n' "$*"; }
warn() { printf '\033[33m[주의] %s\033[0m\n' "$*" >&2; }
die()  { printf '\033[31m[error] %s\033[0m\n' "$*" >&2; exit 1; }
ask()  { local reply; read -r -p "$1" reply || reply=""; printf '%s' "$reply"; }

[ -f config/config.yaml ] || die "config/config.yaml이 없습니다 -- 프로젝트 폴더 안에서 실행하세요."
[ -f environment.yml ]    || die "environment.yml이 없습니다."

# ---- 1. conda / mamba ---------------------------------------------------------------------------
say "1. conda / mamba 찾기"
MGR=""
for c in mamba micromamba conda; do command -v "$c" >/dev/null 2>&1 && { MGR="$c"; break; }; done
[ -n "$MGR" ] || die "mamba / micromamba / conda 중 어느 것도 PATH에 없습니다.
    Miniforge를 먼저 설치하세요: https://github.com/conda-forge/miniforge  (설치 후 새 터미널에서 다시 실행)"
echo "사용: $MGR ($(command -v "$MGR"))"

env_exists() { "$MGR" env list 2>/dev/null | awk 'NF && $1 !~ /^#/ {print $1}' | grep -qx "$1"; }
env_prefix_of() { "$MGR" run -n "$1" python -c 'import sys; print(sys.prefix)' 2>/dev/null | tail -n 1; }

make_env() {   # name, yml
    if env_exists "$1"; then
        if [ "$UPDATE" = 1 ]; then
            echo "환경 '$1'을(를) $2 기준으로 갱신합니다..."
            "$MGR" env update -y -n "$1" -f "$2" --prune
        else
            echo "환경 '$1'이(가) 이미 있어 그대로 씁니다. (갱신하려면 --update)"
        fi
    else
        echo "환경 '$1'을(를) $2 로 만듭니다 (처음엔 수 분 걸립니다)..."
        "$MGR" env create -y -n "$1" -f "$2"
    fi
}

# ---- 2. main environment ------------------------------------------------------------------------
say "2. 분석 환경 ($ENV_NAME)"
make_env "$ENV_NAME" environment.yml
ENV_PREFIX="$(env_prefix_of "$ENV_NAME")"
[ -n "$ENV_PREFIX" ] && [ -d "$ENV_PREFIX/bin" ] || die "환경 '$ENV_NAME'의 경로를 확인하지 못했습니다."
echo "환경 경로: $ENV_PREFIX"

# edits one top-level `key: "value"` line of config/config.yaml (appends it if missing)
set_cfg() {
    "$ENV_PREFIX/bin/python" - "$1" "$2" <<'PY'
import re, sys
key, value = sys.argv[1], sys.argv[2]
path = "config/config.yaml"
text = open(path, encoding="utf-8").read()
line = f'{key}: "{value}"'
if re.search(rf"^{re.escape(key)}:.*$", text, re.M):
    text = re.sub(rf"^{re.escape(key)}:.*$", lambda _m: line, text, count=1, flags=re.M)
else:
    text = text.rstrip("\n") + "\n" + line + "\n"
open(path, "w", encoding="utf-8").write(text)
PY
}
get_cfg() { "$ENV_PREFIX/bin/python" - "$1" <<'PY'
import re, sys
m = re.search(rf'^{re.escape(sys.argv[1])}:\s*"?([^"\n#]*?)"?\s*(?:#.*)?$', open("config/config.yaml", encoding="utf-8").read(), re.M)
print(m.group(1).strip() if m else "")
PY
}
set_cfg env_prefix "$ENV_PREFIX"

# ---- 3. tools in the environment ----------------------------------------------------------------
say "3. 도구 확인"
STATUS_OK=1
check_tool() {   # label, command...
    local label="$1"; shift
    local out
    if out="$(PATH="$ENV_PREFIX/bin:$PATH" "$@" 2>&1 | grep . | tail -n 1)" && [ -n "$out" ]; then
        printf '  ✔ %-10s %s\n' "$label" "$out"
    else
        printf '  ✘ %-10s 실행 실패\n' "$label"; STATUS_OK=0
    fi
}
check_tool snakemake snakemake --version
check_tool mafft mafft --version
check_tool pixy pixy --version
check_tool bgzip bgzip --version
check_tool tabix tabix --version
check_tool zstd zstd --version

# ---- 4. mfeprimer -------------------------------------------------------------------------------
say "4. mfeprimer 실행 파일"
# mfeprimer is a standalone binary (not on conda): download it yourself (https://www.mfeprimer.com/ or its
# GitHub releases) and point to it here.
mfe_ok() {   # path -> 0 if it runs as mfeprimer
    local p="${1/#\~/$HOME}"
    [ -f "$p" ] || { echo "  파일이 없습니다: $p"; return 1; }
    if [ ! -x "$p" ]; then chmod +x "$p" 2>/dev/null && echo "  실행 권한을 추가했습니다." || { echo "  실행 권한이 없습니다: $p"; return 1; }; fi
    local v; v="$("$p" version 2>&1 | head -n 1)" || true
    case "$v" in *mfeprimer*) echo "  ✔ $v"; return 0 ;; *) echo "  mfeprimer가 아닌 것 같습니다 (version 출력: ${v:-없음})"; return 1 ;; esac
}
MFE_DEFAULT="${MFE_ARG:-${MFEPRIMER_BIN:-$(get_cfg mfeprimer_bin)}}"
if [ -z "$MFE_DEFAULT" ] || [ ! -f "${MFE_DEFAULT/#\~/$HOME}" ]; then
    for cand in ../*/mfeprimer/mfeprimer-*/mfeprimer-*-linux-amd64 ../../*/mfeprimer/mfeprimer-*/mfeprimer-*-linux-amd64 "$(command -v mfeprimer 2>/dev/null || true)"; do
        [ -n "$cand" ] && [ -f "$cand" ] && { MFE_DEFAULT="$(readlink -f "$cand")"; break; }
    done
fi
MFE_FINAL=""
if [ "$YES" = 1 ]; then
    if [ -n "$MFE_DEFAULT" ] && mfe_ok "$MFE_DEFAULT"; then MFE_FINAL="${MFE_DEFAULT/#\~/$HOME}"; fi
else
    for _try in 1 2 3; do
        reply="$(ask "mfeprimer 실행 파일 경로${MFE_DEFAULT:+ [$MFE_DEFAULT]} (건너뛰려면 'skip'): ")"
        reply="${reply:-$MFE_DEFAULT}"
        [ "$reply" = "skip" ] && break
        [ -n "$reply" ] || { echo "  경로를 입력하세요."; continue; }
        if mfe_ok "$reply"; then MFE_FINAL="$(readlink -f "${reply/#\~/$HOME}")"; break; fi
    done
fi
if [ -n "$MFE_FINAL" ]; then
    set_cfg mfeprimer_bin "$MFE_FINAL"; echo "config/config.yaml: mfeprimer_bin = $MFE_FINAL"
else
    STATUS_OK=0
    warn "mfeprimer 경로를 설정하지 못했습니다. 웹앱의 'Analysis Parameters' 페이지에서 입력하거나, 'MFEPRIMER_BIN=<경로> ./setup.sh'로 다시 실행하세요."
fi

# ---- 5. Nextstrain environment (only for profiles whose clades come from a local Nextstrain build) ----
PROFILE_PATH="$(get_cfg profile)"
if [ -f environment_nextstrain.yml ] && [ -n "$PROFILE_PATH" ] && grep -q '"auspice_clades"' "$PROFILE_PATH" 2>/dev/null; then
    say "5. Nextstrain 빌드 환경 (클레이드 재계산용, 선택)"
    want="$NS_MODE"
    if [ "$want" = "ask" ]; then
        if [ "$YES" = 1 ]; then want="no"; else
            echo "이 프로필은 클레이드를 로컬 Nextstrain 빌드(scripts/run_nextstrain_nipah.sh)로 얻습니다."
            echo "이미 만들어 둔 data/nextstrain/ 결과를 쓰기만 한다면 필요 없습니다 (용량이 큽니다)."
            [ "$(ask '지금 이 환경도 만들까요? [y/N] ')" = "y" ] && want="yes" || want="no"
        fi
    fi
    if [ "$want" = "yes" ]; then
        make_env "$NS_ENV_NAME" environment_nextstrain.yml
        NS_PREFIX="$(env_prefix_of "$NS_ENV_NAME")"
        [ -n "$NS_PREFIX" ] && { set_cfg nextstrain_env_prefix "$NS_PREFIX"; echo "nextstrain_env_prefix = $NS_PREFIX"; }
    else
        echo "건너뜁니다. 나중에 ./setup.sh --nextstrain 으로 만들 수 있습니다."
    fi
fi

# ---- 6. optional: Chrome for the PDF report ----------------------------------------------------
say "6. 선택 항목"
CHROME=""
for c in google-chrome google-chrome-stable chromium chromium-browser chrome; do command -v "$c" >/dev/null 2>&1 && { CHROME="$c"; break; }; done
if [ -n "$CHROME" ]; then echo "  ✔ PDF 리포트용 브라우저: $CHROME"; else echo "  - Chrome/Chromium이 없어 'PDF 리포트' 버튼만 쓸 수 없습니다 (나머지는 영향 없음)."; fi

say "끝"
if [ "$STATUS_OK" = 1 ]; then
    echo "설정이 완료되었습니다. 다음으로:  ./start.sh   (터미널에 안내되는 주소를 브라우저에서 여세요)"
else
    echo "위의 ✘ / [주의] 항목을 해결한 뒤 ./setup.sh를 다시 실행하세요 (이미 된 단계는 건너뜁니다)."
    exit 1
fi
