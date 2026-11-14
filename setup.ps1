# chichipes UI installer for Windows PowerShell.
#
#   .\setup.ps1 [-Panels all] [-Dir C:\chichipes-ui] [-PortBase 8760] [-NoData]   # Docker Desktop (default)
#   .\setup.ps1 -Wsl --panels lasv,ebov --yes     # WSL2 (Ubuntu) instead: the remaining arguments go to setup.sh
#
# If scripts are blocked:  powershell -ExecutionPolicy Bypass -File .\setup.ps1
param(
    [switch]$Wsl,
    [string]$Panels = "all",
    [string]$Dir = (Join-Path $HOME "chichipes-ui"),
    [int]$PortBase = 8760,
    [switch]$NoData,
    [Parameter(ValueFromRemainingArguments = $true)][string[]]$Rest
)
$ErrorActionPreference = "Stop"
$here = Split-Path -Parent $MyInvocation.MyCommand.Path
function Has($c) { [bool](Get-Command $c -ErrorAction SilentlyContinue) }
function IsAdmin { ([Security.Principal.WindowsPrincipal][Security.Principal.WindowsIdentity]::GetCurrent()).IsInRole([Security.Principal.WindowsBuiltInRole]::Administrator) }
function VirtualizationHint {
    try {
        $v = (Get-CimInstance Win32_Processor | Select-Object -First 1).VirtualizationFirmwareEnabled
        if ($v -eq $false) { Write-Host "[주의] BIOS/UEFI에서 가상화(Intel VT-x / AMD-V)가 꺼져 있습니다. 켠 뒤 다시 실행하세요." -ForegroundColor Yellow }
    } catch {}
}

if (-not $Wsl) {
    if (-not (Has docker)) { Write-Host "Docker Desktop을 설치하세요: https://www.docker.com/products/docker-desktop/"; VirtualizationHint; exit 1 }
    docker info *> $null
    if ($LASTEXITCODE -ne 0) { Write-Host "Docker Desktop을 실행한 뒤 다시 시도하세요."; exit 1 }
    New-Item -ItemType Directory -Force -Path $Dir | Out-Null
    $last = $PortBase + 23
    $cores = [Math]::Max(1, [int]([Environment]::ProcessorCount / 2))
    $runArgs = @("--panels", $Panels, "--port-base", "$PortBase", "--cores", "$cores")
    if ($NoData) { $runArgs += "--no-data" }
    Write-Host "== 이미지 만들기: chichipes-ui (처음 한 번은 수 분)"
    docker build -t chichipes-ui -f (Join-Path $here "install\docker\Dockerfile") $here
    if ($LASTEXITCODE -ne 0) { exit 1 }
    docker rm -f chichipes-ui 2>$null | Out-Null
    docker run -d --name chichipes-ui --restart unless-stopped -v "${Dir}:/work" -p "127.0.0.1:${PortBase}-${last}:${PortBase}-${last}" chichipes-ui @runArgs
    if ($LASTEXITCODE -ne 0) { exit 1 }
    Write-Host "설치 진행 상황: docker logs -f chichipes-ui"
    Write-Host "전체 목록: http://localhost:$PortBase   (패널·데이터: $Dir)"
    exit 0
}

# WSL2 route
$distros = @()
if (Has wsl) { $distros = @((wsl -l -q 2>$null) -replace "`0", "" | Where-Object { $_ -match '\S' }) }
if (-not $distros) {
    Write-Host "WSL2 + Ubuntu가 필요합니다 (Snakemake·MAFFT·htslib은 Windows에서 직접 실행되지 않습니다)."
    VirtualizationHint
    if (IsAdmin) {
        $a = Read-Host "지금 'wsl --install -d Ubuntu'를 실행할까요? 끝나면 재부팅 후 Ubuntu 사용자 이름을 만들고 이 스크립트를 다시 실행하세요 [Y/n]"
        if ($a -notmatch '^[nN]') { wsl --install -d Ubuntu }
    } else {
        Write-Host "관리자 권한 PowerShell에서 실행하세요:  wsl --install -d Ubuntu   (재부팅 후 이 스크립트를 다시 실행)"
        Write-Host "또는 Docker Desktop으로:  .\setup.ps1"
    }
    exit 1
}
$wslPath = (wsl wslpath -a ($here -replace '\\', '/')).Trim()
$quoted = ($Rest | ForEach-Object { "'" + ($_ -replace "'", "'\''") + "'" }) -join " "
Write-Host "== WSL(Ubuntu) 안에서 setup.sh 실행: $wslPath"
wsl bash -lc "cd '$wslPath' && bash ./setup.sh $quoted"
exit $LASTEXITCODE
