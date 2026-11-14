# Environment for macOS (Apple Silicon or Intel): mamba (default, Miniforge installed when missing), Homebrew or
# Docker Desktop. Sourced by setup.sh; defines prepare_environment and install_with_docker.

need_mac_tools() {
    if ! xcode-select -p >/dev/null 2>&1; then
        warn "Command Line Tools(git 포함)가 없습니다. 열리는 창에서 설치한 뒤 setup.sh를 다시 실행하세요."
        xcode-select --install 2>/dev/null || true
        exit 1
    fi
    have curl || die "curl이 없습니다."
}

prepare_environment() {
    need_mac_tools
    case "$BACKEND" in
        brew)
            prepare_brew_env
            ;;
        *)
            if [ "$BACKEND" = auto ] && ! find_mamba && have brew; then
                info "mamba가 없고 Homebrew가 있습니다. 권장은 mamba입니다 (--backend brew 로 Homebrew 사용 가능)."
            fi
            find_mamba || install_miniforge
            info "환경 관리자: $MGR"
            # The exact Linux versions of environment.yml are not all built for macOS: use the unpinned file.
            if ! create_env "$HERE/install/environment.portable.yml"; then
                if [ "$(uname -m)" = arm64 ]; then
                    warn "Apple Silicon용 패키지가 없어 Intel(osx-64, Rosetta) 패키지로 다시 시도합니다."
                    have arch && arch -x86_64 /usr/bin/true 2>/dev/null || { softwareupdate --install-rosetta --agree-to-license || die "Rosetta가 필요합니다."; }
                    CONDA_SUBDIR=osx-64 create_env "$HERE/install/environment.portable.yml" || die "환경을 만들지 못했습니다."
                    "$MGR" run -n "$ENV_NAME" conda config --env --set subdir osx-64 >/dev/null 2>&1 || true
                else
                    die "환경을 만들지 못했습니다."
                fi
            fi
            ;;
    esac
}

install_with_docker() { docker_install; }
