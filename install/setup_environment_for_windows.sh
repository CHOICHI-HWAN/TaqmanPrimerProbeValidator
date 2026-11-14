# Windows (Git Bash / MSYS). Snakemake's tools (MAFFT, htslib, pixy) have no Windows builds, so the panels cannot run
# on Windows itself. Default: Docker Desktop. Alternative: run setup.sh inside WSL2 (Ubuntu).
# (Without Git Bash, use setup.ps1 from PowerShell -- it does the same.)

wsl_help() {
    cat <<'EOF'

  Windows에서 바로는 실행할 수 없습니다 (Snakemake·MAFFT·htslib에 Windows판이 없음). 둘 중 하나를 쓰세요:

  1) Docker Desktop (기본)
     https://www.docker.com/products/docker-desktop/ 에서 설치하고 실행한 뒤:   ./setup.sh

  2) WSL2 (Ubuntu)
     PowerShell(관리자)에서:   wsl --install -d Ubuntu     (처음 한 번, 재부팅 필요)
     그다음 Ubuntu 창에서:      git clone -b chichipes https://github.com/CHOICHI-HWAN/TaqmanPrimerProbeValidator.git chichipes
                                cd chichipes && ./setup.sh
     PowerShell에서는 .\setup.ps1 -Wsl 이 같은 일을 합니다.

  WSL2와 Docker Desktop 모두 BIOS에서 가상화(VT-x/AMD-V)가 켜져 있어야 하고, 처음 설치에는 관리자 권한이 필요합니다.
EOF
}

prepare_environment() {
    wsl_help
    exit 1
}

install_with_docker() {
    export MSYS_NO_PATHCONV=1                # keep /work etc. from being rewritten into C:/Program Files/Git/...
    DOCKER_DIR="$(cygpath -w "$DIR" 2>/dev/null || echo "$DIR")"
    docker_install
}

case "$BACKEND" in
    auto|docker)
        have docker || { wsl_help; exit 1; }
        docker info >/dev/null 2>&1 || die "Docker Desktop을 실행한 뒤 다시 시도하세요."
        BACKEND=docker ;;
    *) wsl_help; exit 1 ;;
esac
