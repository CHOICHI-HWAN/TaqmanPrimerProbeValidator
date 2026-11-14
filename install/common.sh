# Shared steps of setup.sh (sourced; bash 3.2 compatible for macOS). Needs the variables setup.sh sets.

say()  { printf '\n\033[1m== %s\033[0m\n' "$*"; }
info() { printf '   %s\n' "$*"; }
warn() { printf '\033[33m[주의] %s\033[0m\n' "$*" >&2; }
die()  { printf '\033[31m[error] %s\033[0m\n' "$*" >&2; exit 1; }
ask_yes() { # ask_yes "question" -> 0 for yes (default yes)
    [ "$YES" = 1 ] && return 0
    printf '%s [Y/n] ' "$1"; read -r a || true
    case "$a" in n|N|no|NO) return 1 ;; *) return 0 ;; esac
}
have() { command -v "$1" >/dev/null 2>&1; }

detect_os() {
    case "$(uname -s)" in
        Darwin) echo mac ;;
        MINGW*|MSYS*|CYGWIN*) echo windows ;;
        Linux)
            if grep -qi microsoft /proc/version 2>/dev/null; then echo wsl
            elif [ -r /etc/os-release ] && grep -qiE '^ID(_LIKE)?=.*(ubuntu|debian)' /etc/os-release; then echo ubuntu
            else echo linux; fi ;;
        *) echo unknown ;;
    esac
}

cpu_count() {
    if have nproc; then nproc
    elif have sysctl; then sysctl -n hw.ncpu
    else echo 2; fi
}
[ -n "$CORES" ] || CORES=$(( $(cpu_count) / 2 )); [ "$CORES" -ge 1 ] || CORES=1

go_os_arch() { # chichipes release naming: <os>-<arch>
    local os arch
    case "$(uname -s)" in Darwin) os=darwin ;; Linux) os=linux ;; *) os=windows ;; esac
    case "$(uname -m)" in x86_64|amd64) arch=amd64 ;; arm64|aarch64) arch=arm64 ;; *) arch="$(uname -m)" ;; esac
    echo "$os-$arch"
}

# ---- conda/mamba (used by the Ubuntu and macOS scripts) --------------------------------------------------------------
MGR=""
find_mamba() {
    local c
    for c in mamba micromamba conda; do have "$c" && { MGR="$(command -v "$c")"; return 0; }; done
    for c in "$HOME/miniforge3/bin/mamba" "$HOME/mambaforge/bin/mamba" "$HOME/miniconda3/bin/conda" "$HOME/anaconda3/bin/conda"; do
        [ -x "$c" ] && { MGR="$c"; return 0; }
    done
    return 1
}

install_miniforge() {
    local os arch url tmp
    case "$(uname -s)" in Darwin) os=MacOSX ;; *) os=Linux ;; esac
    arch="$(uname -m)"; [ "$arch" = "aarch64" ] && [ "$os" = "MacOSX" ] && arch=arm64
    url="https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-$os-$arch.sh"
    ask_yes "mamba가 없습니다. Miniforge(mamba 포함)를 $HOME/miniforge3 에 설치할까요?" || die "mamba / conda가 필요합니다 (또는 --backend docker)."
    tmp="$(mktemp -d)"
    info "내려받는 중: $url"
    curl -fsSL "$url" -o "$tmp/miniforge.sh" || die "Miniforge를 내려받지 못했습니다: $url"
    bash "$tmp/miniforge.sh" -b -p "$HOME/miniforge3" >/dev/null || die "Miniforge 설치 실패"
    rm -rf "$tmp"
    MGR="$HOME/miniforge3/bin/mamba"
    info "설치됨: $HOME/miniforge3 (새 터미널에서 쓰려면: $HOME/miniforge3/bin/conda init)"
}

env_prefix_of() { "$MGR" run -n "$1" python -c 'import sys; print(sys.prefix)' 2>/dev/null | tail -n 1; }

create_env() { # create_env <environment file>
    local file="$1" prefix
    prefix="$(env_prefix_of "$ENV_NAME" || true)"
    if [ -n "$prefix" ] && [ -x "$prefix/bin/snakemake" ] && [ "$UPDATE" = 0 ]; then
        info "기존 환경 사용: $ENV_NAME ($prefix)"
    elif [ -n "$prefix" ] && [ -d "$prefix" ]; then
        say "환경 업데이트: $ENV_NAME"
        "$MGR" env update -n "$ENV_NAME" -f "$file" || die "환경 업데이트 실패 ($file)"
    else
        say "환경 만들기: $ENV_NAME (몇 분 걸립니다)"
        case "$(basename "$MGR")" in
            micromamba) "$MGR" create -y -n "$ENV_NAME" -f "$file" ;;
            *) "$MGR" env create -y -n "$ENV_NAME" -f "$file" ;;
        esac || return 1
    fi
    ENV_PREFIX="$(env_prefix_of "$ENV_NAME")"
    [ -n "$ENV_PREFIX" ] && [ -x "$ENV_PREFIX/bin/snakemake" ] || die "환경 '$ENV_NAME'을 확인하지 못했습니다."
}

# ---- Homebrew backend (macOS, or Linuxbrew): tools from brew, Python packages in a venv -----------------------------
prepare_brew_env() {
    have brew || die "Homebrew가 없습니다: https://brew.sh 에서 설치하거나 --backend mamba를 쓰세요."
    say "Homebrew로 도구 설치 (python@3.12 mafft htslib zstd git)"
    brew install python@3.12 mafft htslib zstd git >/dev/null || die "brew install 실패"
    local bp; bp="$(brew --prefix)"
    ENV_PREFIX="$DIR/env"
    if [ ! -x "$ENV_PREFIX/bin/python" ] || [ "$UPDATE" = 1 ]; then
        "$bp/opt/python@3.12/bin/python3.12" -m venv "$ENV_PREFIX" || die "venv를 만들지 못했습니다."
        "$ENV_PREFIX/bin/pip" install -q --upgrade pip
        "$ENV_PREFIX/bin/pip" install -q "snakemake==9.27.0" pyyaml scikit-allel==1.3.13 numcodecs dask pytest primer3-py==2.3.1 \
            varvamp==1.3.2 "pixy @ git+https://github.com/ksamuk/pixy.git@21c7c33dea315a035d6aff86caba018ab668864a" \
            || die "pip 설치 실패"
    fi
    local t
    for t in mafft bgzip tabix zstd git; do # the panels look for their tools in env_prefix/bin
        [ -e "$ENV_PREFIX/bin/$t" ] || ln -s "$bp/bin/$t" "$ENV_PREFIX/bin/$t" 2>/dev/null || true
    done
    warn "brew 방식은 버전이 고정되지 않아 결과가 기준 서버와 조금 다를 수 있습니다 (권장: mamba)."
}

# ---- chichipes engine -------------------------------------------------------------------------------------------------
install_chichipes() {
    say "chichipes 엔진"
    mkdir -p "$DIR/bin"
    CHI="$DIR/bin/chichipes"
    if [ -x "$HERE/chichipes" ] && "$HERE/chichipes" version >/dev/null 2>&1; then
        cp "$HERE/chichipes" "$CHI"
    elif have go && [ -f "$HERE/go.mod" ]; then
        info "소스에서 빌드 (go)"
        (cd "$HERE" && CGO_ENABLED=0 go build -o "$CHI" ./cmd/chichipes) || die "chichipes 빌드 실패"
    else
        local oa api tag url tmp
        oa="$(go_os_arch)"
        api="https://api.github.com/repos/CHOICHI-HWAN/TaqmanPrimerProbeValidator/releases?per_page=50"
        tag="$("$PY" -c 'import json,sys,urllib.request
rel=json.load(urllib.request.urlopen(sys.argv[1], timeout=30))
print(next(r["tag_name"] for r in rel if r["tag_name"].startswith("chichipes-v") and not r.get("prerelease")))' "$api")" \
            || die "chichipes 릴리즈를 찾지 못했습니다 (인터넷 연결 확인)."
        url="https://github.com/CHOICHI-HWAN/TaqmanPrimerProbeValidator/releases/download/$tag/chichipes-${tag#chichipes-}-$oa.tar.gz"
        info "내려받는 중: $url"
        tmp="$(mktemp -d)"
        curl -fsSL "$url" -o "$tmp/c.tgz" || die "내려받기 실패: $url"
        tar -xzf "$tmp/c.tgz" -C "$tmp" && cp "$tmp"/chichipes-*/chichipes "$CHI" && rm -rf "$tmp"
    fi
    chmod +x "$CHI"
    info "$("$CHI" version)  ->  $CHI"
}

# ---- panels -------------------------------------------------------------------------------------------------------------
fetch_entry() {
    say "패널 목록 (entry)"
    if [ -d "$DIR/entry/.git" ]; then
        git -C "$DIR/entry" checkout -q -- pipelines.json 2>/dev/null || true
        [ "$UPDATE" = 1 ] && { git -C "$DIR/entry" pull -q --ff-only || warn "entry 업데이트 실패 (기존 것 사용)"; }
    else
        git clone -q --depth 1 -b entry "$REPO_URL" "$DIR/entry" || die "entry 브랜치를 받지 못했습니다: $REPO_URL"
    fi
    cp "$DIR/entry/pipelines.json" "$DIR/entry/pipelines.catalog.json"
}

catalog() { # id<TAB>path<TAB>port<TAB>name<TAB>size<TAB>heavy<TAB>note
    "$PY" - "$DIR/entry/pipelines.catalog.json" "$HERE/install/panels.tsv" <<'PY'
import csv, json, sys
extra = {}
for row in csv.reader(open(sys.argv[2], encoding="utf-8"), delimiter="\t"):
    if row and not row[0].startswith("#"):
        extra[row[0]] = (row + ["", "", ""])[1:4]
for p in json.load(open(sys.argv[1], encoding="utf-8"))["pipelines"]:
    size, heavy, note = extra.get(p["id"], ("?", "no", ""))
    print("\t".join([p["id"], p["path"], str(p["port"]), p["name"], size, heavy or "no", note]))
PY
}

choose_panels() {
    local cat ids
    cat="$(catalog)"
    if [ -z "$PANELS" ] && [ "$YES" = 0 ]; then
        say "설치할 패널"
        printf '   %-7s %-7s %-9s %s\n' id 이름 데이터 비고
        printf '%s\n' "$cat" | while IFS="$(printf '\t')" read -r id path port name size heavy note; do
            printf '   %-7s %-7s %-9s %s\n' "$id" "$name" "$size" "$( [ "$heavy" = yes ] && echo "대용량 — all에 포함 안 됨. $note" )"
        done
        printf '\n   쉼표로 고르거나 all (대용량 제외) / everything: [all] '
        read -r PANELS || true
    fi
    [ -n "$PANELS" ] || PANELS=all
    case "$PANELS" in
        all) ids="$(printf '%s\n' "$cat" | awk -F'\t' '$6!="yes"{print $1}')" ;;
        everything) ids="$(printf '%s\n' "$cat" | awk -F'\t' '{print $1}')" ;;
        *) ids="$(printf '%s' "$PANELS" | tr ',' ' ')" ;;
    esac
    SELECTED=""
    for id in $ids; do
        printf '%s\n' "$cat" | awk -F'\t' -v i="$id" '$1==i{f=1} END{exit !f}' || die "알 수 없는 패널: $id"
        SELECTED="$SELECTED $id"
    done
    info "선택: $SELECTED"
}

panel_field() { catalog | awk -F'\t' -v i="$1" -v f="$2" '$1==i{print $f}'; }

cfg_set() { # cfg_set <panel dir> <key> <value>
    "$PY" - "$1/config/config.yaml" "$2" "$3" <<'PY'
import re, sys
path, key, value = sys.argv[1:4]
text = open(path, encoding="utf-8").read()
line = f'{key}: "{value}"'
if re.search(rf"^{re.escape(key)}:.*$", text, re.M):
    text = re.sub(rf"^{re.escape(key)}:.*$", lambda _m: line, text, count=1, flags=re.M)
else:
    text = text.rstrip("\n") + "\n" + line + "\n"
open(path, "w", encoding="utf-8").write(text)
PY
}

install_panel() {
    local id="$1" path dest
    path="$(panel_field "$id" 2)"; dest="$DIR/$path"
    say "패널 $id  ->  $path"
    if [ -d "$dest/.git" ]; then
        [ "$UPDATE" = 1 ] && { git -C "$dest" pull -q --ff-only || warn "$id 업데이트 실패 (기존 것 사용)"; }
    else
        mkdir -p "$(dirname "$dest")"
        git clone -q --depth 1 -b "$id" "$REPO_URL" "$dest" || die "$id 브랜치를 받지 못했습니다."
    fi
    cfg_set "$dest" mfeprimer_bin "$CHI"
    cfg_set "$dest" env_prefix "$ENV_PREFIX"
    if have go && [ -f "$dest/go/scan/main.go" ]; then
        (cd "$dest/go/scan" && CGO_ENABLED=0 go build -o scan . 2>/dev/null) && info "go/scan 빌드됨 (빠른 스캔)" || true
    fi
}

write_registry() {
    "$PY" - "$DIR" "$PORT_BASE" <<'PY'
import json, sys
from pathlib import Path
d, base = Path(sys.argv[1]), int(sys.argv[2])
cat = json.loads((d / "entry" / "pipelines.catalog.json").read_text(encoding="utf-8"))
keep = [p for p in cat["pipelines"] if (d / p["path"] / "Snakefile").exists()]
for p in keep:
    p["port"] = p["port"] - 8760 + base
cat["pipelines"], cat["entry_port"] = keep, base
(d / "entry" / "pipelines.json").write_text(json.dumps(cat, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
print(f"   front page :{base}, panels: " + ", ".join(f"{p['id']}:{p['port']}" for p in keep))
PY
}

import_primers() {
    say "프라이머 세트 가져오기: $1"
    "$PY" - "$DIR" "$1" <<'PY'
import csv, json, sys, time
from pathlib import Path
d, src = Path(sys.argv[1]), sys.argv[2]
rows = list(csv.DictReader(open(src, encoding="utf-8"), delimiter="\t"))
if not rows or "edition" not in rows[0]:
    sys.exit("[error] the file needs an 'edition' column (tools/export_primer_sets.py format)")
fields = [f for f in rows[0] if f != "edition"]
for p in json.loads((d / "entry" / "pipelines.json").read_text(encoding="utf-8"))["pipelines"]:
    mine = [r for r in rows if r["edition"] == p["id"]]
    if not mine:
        continue
    out = d / p["path"] / "config" / "primers.tsv"
    if out.exists():
        out.rename(out.with_name(f"primers.tsv.bak-{time.strftime('%Y%m%d-%H%M%S')}"))
    with open(out, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=fields, delimiter="\t", lineterminator="\n", extrasaction="ignore")
        w.writeheader()
        w.writerows(mine)
    print(f"   {p['id']}: {len(mine)} sets")
PY
}

# ---- launchers ----------------------------------------------------------------------------------------------------------
write_launchers() {
    say "실행 스크립트"
    mkdir -p "$DIR/logs" "$DIR/run"
    cat > "$DIR/start.sh" <<EOF
#!/usr/bin/env bash
# Starts the front page and every installed panel in the background (logs in logs/, stop with ./stop.sh).
DIR="$DIR"; PY="$PY"; HOST="\${CHICHIPES_HOST:-127.0.0.1}"
export PATH="$ENV_PREFIX/bin:\$PATH"
up() { "\$PY" -c 'import socket,sys; s=socket.socket(); s.settimeout(0.3); sys.exit(0 if s.connect_ex(("127.0.0.1", int(sys.argv[1]))) == 0 else 1)' "\$1"; }
start() { # id dir port script
    if up "\$3"; then echo "   \$1: already running on :\$3"; return; fi
    (cd "\$2" && nohup "\$PY" "\$4" --host "\$HOST" --port "\$3" > "\$DIR/logs/web_\$1.log" 2>&1 & echo \$! > "\$DIR/run/\$1.pid")
    echo "   \$1: http://localhost:\$3"
}
"\$PY" - "\$DIR/entry/pipelines.json" <<'PY' | while IFS="\$(printf '\t')" read -r id path port; do start "\$id" "\$DIR/\$path" "\$port" webapp/server.py; done
import json, sys
for p in json.load(open(sys.argv[1], encoding="utf-8"))["pipelines"]:
    print(f"{p['id']}\t{p['path']}\t{p['port']}")
PY
EP=\$("\$PY" -c 'import json,sys; print(json.load(open(sys.argv[1])).get("entry_port", 8760))' "\$DIR/entry/pipelines.json")
start entry "\$DIR/entry" "\$EP" server.py
echo "전체 목록: http://localhost:\$EP"
EOF
    cat > "$DIR/stop.sh" <<EOF
#!/usr/bin/env bash
# Stops what start.sh started.
for f in "$DIR"/run/*.pid; do
    [ -e "\$f" ] || continue
    kill "\$(cat "\$f")" 2>/dev/null && echo "stopped \$(basename "\$f" .pid)"
    rm -f "\$f"
done
EOF
    cat > "$DIR/fetch_data.sh" <<EOF
#!/usr/bin/env bash
# Downloads each installed panel's genome data and runs its analysis, one panel after another.
#   ./fetch_data.sh [panel ...]      (default: every installed panel)  -- progress: logs/data_status.tsv
DIR="$DIR"; PY="$PY"; CORES="$CORES"
export PATH="$ENV_PREFIX/bin:\$PATH"
ids="\$*"
[ -n "\$ids" ] || ids=\$("\$PY" -c 'import json,sys; print(" ".join(p["id"] for p in json.load(open(sys.argv[1]))["pipelines"]))' "\$DIR/entry/pipelines.json")
for id in \$ids; do
    path=\$("\$PY" -c 'import json,sys; print(next(p["path"] for p in json.load(open(sys.argv[1]))["pipelines"] if p["id"]==sys.argv[2]))' "\$DIR/entry/pipelines.json" "\$id")
    printf '%s\trunning\t%s\t\n' "\$id" "\$(date '+%F %T')" >> "\$DIR/logs/data_status.tsv"
    if (cd "\$DIR/\$path" && snakemake --cores "\$CORES" -k > "\$DIR/logs/data_\$id.log" 2>&1); then st=done; else st=failed; fi
    printf '%s\t%s\t\t%s\n' "\$id" "\$st" "\$(date '+%F %T')" >> "\$DIR/logs/data_status.tsv"
    echo "\$id: \$st (log: \$DIR/logs/data_\$id.log)"
done
EOF
    chmod +x "$DIR/start.sh" "$DIR/stop.sh" "$DIR/fetch_data.sh"
    info "$DIR/start.sh · stop.sh · fetch_data.sh"
}

start_data_fetch() {
    say "데이터 받기 (백그라운드, CPU $CORES개)"
    nohup "$DIR/fetch_data.sh" > "$DIR/logs/fetch_data.log" 2>&1 &
    echo $! > "$DIR/run/fetch_data.pid"
    info "진행 상황: $DIR/logs/data_status.tsv · 패널별 로그: $DIR/logs/data_<id>.log"
}

# ---- Docker backend (host side) ----------------------------------------------------------------------------------------
docker_install() {
    have docker || die "docker가 없습니다. Docker Desktop(Windows/macOS) 또는 Docker Engine(Linux)을 설치하세요."
    docker info >/dev/null 2>&1 || die "Docker가 실행 중이 아닙니다 (Docker Desktop을 켜세요)."
    if [ -z "$PANELS" ] && [ "$YES" = 0 ]; then
        say "설치할 패널"
        awk -F'\t' '!/^#/ && NF {printf "   %-7s %-9s %s\n", $1, $2, ($3=="yes" ? "대용량 — all에 포함 안 됨. " $4 : "")}' "$HERE/install/panels.tsv"
        printf '\n   쉼표로 고르거나 all (대용량 제외) / everything: [all] '
        read -r PANELS || true
    fi
    [ -n "$PANELS" ] || PANELS=all
    local args="--panels $PANELS --port-base $PORT_BASE --cores $CORES" last=$((PORT_BASE + 23)) user=""
    [ "$NO_DATA" = 1 ] && args="$args --no-data"
    [ "$UPDATE" = 1 ] && args="$args --update"
    if [ -n "$PRIMERS" ]; then cp "$PRIMERS" "$DIR/primers_import.tsv"; args="$args --primers /work/primers_import.tsv"; fi
    [ "$(uname -s)" = Linux ] && user="--user $(id -u):$(id -g)"   # files in the install folder stay yours
    say "이미지 만들기: chichipes-ui (처음 한 번은 수 분)"
    docker build -t chichipes-ui -f "$HERE/install/docker/Dockerfile" "$HERE" || die "docker build 실패"
    docker rm -f chichipes-ui >/dev/null 2>&1 || true
    say "컨테이너 시작"
    # shellcheck disable=SC2086
    docker run -d --name chichipes-ui --restart unless-stopped $user -v "${DOCKER_DIR:-$DIR}:/work" \
        -p "127.0.0.1:$PORT_BASE-$last:$PORT_BASE-$last" chichipes-ui $args >/dev/null || die "docker run 실패"
    info "설치 진행 상황: docker logs -f chichipes-ui"
    info "전체 목록: http://localhost:$PORT_BASE   (패널·데이터는 $DIR 에 저장)"
    info "멈춤: docker stop chichipes-ui   다시 시작: docker start chichipes-ui"
}

summary() {
    say "완료"
    info "전체 목록: http://localhost:$PORT_BASE"
    info "다시 시작: $DIR/start.sh   멈춤: $DIR/stop.sh   데이터: $DIR/fetch_data.sh"
    info "공개 패널에는 프라이머 세트가 없습니다. 각 패널의 '프라이머 세트' 화면에서 넣거나 --primers 로 가져오세요."
}
