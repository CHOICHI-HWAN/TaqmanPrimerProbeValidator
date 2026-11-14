# Environment for Ubuntu / Debian / WSL2 (and other Linux): mamba (default), Homebrew (Linuxbrew) or Docker.
# Sourced by setup.sh; defines prepare_environment and install_with_docker.

need_system_tools() {
    local missing="" t
    for t in git curl; do have "$t" || missing="$missing $t"; done
    [ -z "$missing" ] && return 0
    if have apt-get; then
        ask_yes "필요한 도구를 설치할까요 (sudo apt-get install$missing)?" || die "설치하세요:$missing"
        sudo apt-get update -q && sudo apt-get install -y -q $missing || die "apt-get 실패"
    else
        die "먼저 설치하세요:$missing"
    fi
}

prepare_environment() {
    case "$BACKEND" in
        container) # inside install/docker's image: the environment is already there
            ENV_PREFIX="$(python -c 'import sys; print(sys.prefix)')"
            ;;
        brew)
            need_system_tools
            prepare_brew_env
            ;;
        *)
            need_system_tools
            find_mamba || install_miniforge
            info "환경 관리자: $MGR"
            create_env "$HERE/install/environment.yml" || {
                warn "고정 버전 환경을 만들지 못했습니다. 버전을 덜 고정한 환경으로 다시 시도합니다."
                create_env "$HERE/install/environment.portable.yml" || die "환경을 만들지 못했습니다."
            }
            ;;
    esac
}

install_with_docker() { docker_install; }
