# chichipes

[English](README.md)

**chichipes**는 TaqMan 검증 파이프라인이 쓰는 [MFEprimer 4.2.4](https://www.mfeprimer.com) 명령(프라이머쌍 증폭산물
검색과 바이러스 커버리지, 프로브 검색, 다이머·헤어핀 검사)을 그대로 대신하는 실행 파일 하나짜리 엔진이에요.
명령줄도, 출력 파일도, 숫자도 같고, 보통 수십~수백 배 빨라요.

[TaqMan Primer/Probe Validator](https://github.com/CHOICHI-HWAN/TaqmanPrimerProbeValidator) 파이프라인의 검색 엔진이며,
엔진만 필요한 분을 위해 별도 브랜치(`chichipes`)로 배포해요.

- **실행 파일 하나** — Rust로 작성(v0.4까지는 Go). Linux, macOS, Windows용 단일 실행 파일이고 따로 설치할 것이 없어요.
  Linux용은 glibc 2.17 이상(CentOS 7, Ubuntu 14.04 이후)에서 돌아가요.
- **색인 불필요** — DB를 바로 읽어요. `chichipes index`는 다른 도구가 기대하는 `.fai`/`.json`만 만들어요.
- **클린룸 구현** — 공개 열역학 논문과 MFEprimer 출력 비교 실험만으로 규칙을 찾았어요. MFEprimer 코드를 쓰거나
  디컴파일하지 않았어요. 자세한 규칙: [docs/FINDINGS.md](docs/FINDINGS.md).

## 설치

[Releases](https://github.com/CHOICHI-HWAN/TaqmanPrimerProbeValidator/releases)의 `chichipes-v*` 릴리즈에서 운영체제에 맞는 압축 파일을 받아
`chichipes`를 `PATH`에 넣으면 돼요.

```bash
tar xzf chichipes-v0.5.0-linux-amd64.tar.gz
sudo install chichipes-v0.5.0-linux-amd64/chichipes /usr/local/bin/
chichipes version
```

소스에서 빌드 ([Rust](https://rustup.rs) stable과 C 링커 `cc` 필요):

```bash
git clone -b chichipes https://github.com/CHOICHI-HWAN/TaqmanPrimerProbeValidator.git chichipes
cd chichipes
cargo build --release          # -> target/release/chichipes
cargo test --release           # work/에 MFEprimer 기준 실행 결과가 있으면 대조 테스트도 돌아요
```

## 사용법

옵션이 MFEprimer와 같아서, 기존 스크립트는 실행 파일 경로만 바꾸면 돼요.

| 명령 | MFEprimer 대응 | 하는 일 |
|---|---|---|
| `chichipes -d db.fa -i FR.fa -o out.txt --fasta --virus --misMatch 2 [--cpu N]` | `mfeprimer` (spec, virus 모드) | 모든 유전체에서 프라이머쌍이 만드는 증폭산물, 커버리지, 놓친 유전체 |
| `chichipes search db.fa query.fa out.txt -m 8` | `mfeprimer search` | 프로브 등 올리고 결합 위치 (편집거리 미스매치, Tm, ΔG) |
| `chichipes dimer -i oligos.fa -o dimer.txt` | `mfeprimer dimer` | 프라이머·프로브 다이머 |
| `chichipes hairpin -i oligos.fa -o hairpin.txt` | `mfeprimer hairpin` | 헤어핀 |
| `chichipes index -i db.fa` | `mfeprimer index` | `db.fa.fai` / `db.fa.json` 작성 (검색엔 색인 불필요) |
| `chichipes reads -i oligos.fa -r reads.fastq.gz -o out` | — | 시퀀싱 리드에 프라이머·프로브가 있는지 (아래) |
| `chichipes version` | `mfeprimer version` | 버전 |

출력 파일은 MFEprimer와 이름·형식이 같아요: `out.txt`(`Coverage …` 줄 포함 보고서), `out.txt.spec.tsv`,
`out.txt.amp.fasta`, `out.txt.virus.Failed.txt`. 추가로 `out.txt.borderline.tsv`(Tm 기준 ±0.02 °C 결합 목록, 참고용).
형식 상세: [docs/OUTPUT_FORMAT.md](docs/OUTPUT_FORMAT.md).

TaqMan Primer/Probe Validator 파이프라인에서는 `config/config.yaml` 한 줄만 바꾸면 돼요:

```yaml
mfeprimer_bin: "/usr/local/bin/chichipes"
```

## MFEprimer와 얼마나 같은가

17개 바이러스 패널의 실제 프라이머 세트 182개로 MFEprimer 4.2.4 출력과 비교했어요.

| 기능 | 결과 |
|---|---|
| 프라이머쌍 검색 | 서로 다른 증폭산물 2,295,176개: 빠짐 0, 추가 0, F/R Tm·ΔG 동일 |
| 프로브 검색(`-m 8`, 182건)·참조 유전체 검색(`-m 6`, 812건) | 모든 행 동일 |
| 다이머 182건 | 항목·값 동일 (ΔG가 같은 항목끼리 순서만 다를 수 있음 — MFEprimer도 `--cpu`>1이면 실행마다 바뀜) |
| 헤어핀 182건 | 머리말·꼬리말 빼고 바이트 동일 |
| 정밀 비교 (`-j` JSON) | 결합 135,851 · 다이머 15,864 · 헤어핀 29,509건 ΔG·Tm 정확히 일치 |

**일부러 다르게 둔 것** (chichipes가 더 정확):

- **증폭산물 1만 개 초과** — MFEprimer는 앞 1만 개로만 커버리지·Failed 목록을 계산해, 축퇴도가 높은 세트에서 실제로는
  증폭되는 유전체를 놓친 것으로 보고해요. chichipes는 전부 셈.
- **여러 세트를 합친 입력의 중복 올리고** — 어떤 올리고의 펼친 서열이 모두 앞 올리고와 겹치면 MFEprimer는 그 올리고의
  미스매치 변이를 만들지 않아 산물을 놓쳐요. chichipes는 변이를 모두 만들고, 같은 서열의 이름도 모두 보고해요.

알려진 작은 차이: 프로브 Tm이 반올림 경계에서 0.01 °C 다를 수 있음, `ampGC` 열이 일부 행에서 G/C 몇 개 차이,
점수가 같은 헤어핀 후보 선택이 무작위 올리고 약 7천 개 중 1개꼴로 다름. 기본값이 아닌 염 농도·헤어핀
`--maxl`/`--minl`은 아직 검증하지 않았어요.

### 속도

| 작업 | MFEprimer 4.2.4 | chichipes |
|---|---|---|
| 축퇴 프라이머쌍 하나 × 유전체 663개 분절 | 18.2초 | 0.1초 |
| 에볼라바이러스 패널 전체 파이프라인 (유전체 2,900개, 29단계) | 14분 58초 | 16초 |

Rust 엔진(v0.5.0)은 Go 엔진(v0.4.0)과 출력이 바이트 단위로 같고 더 빨라요. 실제 프라이머 검색 182건이 CPU 1개로
326초 → 70초(CPU 24개로 42초 → 20초), 가장 큰 세트(증폭산물 22만 8천 개)가 0.84초 → 0.29초(메모리 140MB → 57MB),
롱리드 4,000개 `--indel` 검사가 13.5초 → 0.5초.

## 이노신 (I)

올리고에 deoxyinosine(`I`)을 쓸 수 있어요. IUPAC 코드는 아니지만 축퇴 프라이머에 자주 쓰여요. chichipes는 `I`를 `N`과
같게 보지 않고, **하나의** 염기가 모든 염기와 고르지 않게 짝짓는 것으로 계산해요.

| | `N` | `I` |
|---|---|---|
| 합성되는 것 | 올리고 4종 혼합 | 올리고 1종 |
| 변이로 펼침 | ×4 | 안 펼침 |
| 미스매치로 셈 | 아니요 | 아니요 (A·C·G·T 모두) |
| 그 자리의 Tm/ΔG | 각 변이가 완전 결합 | I·C > I·A > I·T ≈ I·G, 모두 Watson–Crick 쌍보다 약함 |

이노신 스택 값은 Watkins & SantaLucia (2005) *Nucleic Acids Res.* 33:6258–6267, Table 2의 최근접 이웃 파라미터예요.
논문에 값이 없는 스택(미스매치·`N` 옆의 이노신, dangling end)은 `N`과 같은 규칙으로 염 보정 항만 더해요. 3′ 창에 있는
이노신도 찾아요. 다이머·헤어핀 탐지에서는 I·C와 I·A를 짝으로 셉니다. MFEprimer는 `I`를 받지 않아 비교 대상이 없고,
`I`가 없는 올리고는 이전과 결과가 바이트 단위로 같아요.

## 대시보드(UI) 설치

`setup.sh` 하나로 chichipes를 TaqMan Primer/Probe Validator 전체로 만들 수 있어요. 이 저장소의 병원체 대시보드(패널)를
설치하고, chichipes를 쓰도록 맞추고, 유전체 데이터를 받고, 화면을 띄워요.

```bash
git clone -b chichipes https://github.com/CHOICHI-HWAN/TaqmanPrimerProbeValidator.git chichipes
cd chichipes
./setup.sh                          # 패널을 고르고 설치 → http://localhost:8760
./setup.sh --panels lasv,ebov --yes # 묻지 않고 설치
```

| 컴퓨터 | 명령 | 하는 일 |
|---|---|---|
| Ubuntu / Debian / 기타 Linux | `./setup.sh` | mamba 환경 (없으면 Miniforge 설치) |
| macOS (Apple Silicon·Intel) | `./setup.sh` | mamba 환경, `--backend brew`면 Homebrew |
| Windows | PowerShell에서 `.\setup.ps1` | Docker Desktop (기본), `.\setup.ps1 -Wsl`이면 WSL2 Ubuntu 안에서 `setup.sh` 실행 |
| Docker가 있는 모든 컴퓨터 | `./setup.sh --backend docker` | 컨테이너 하나에서 실행, 패널·데이터는 설치 폴더에 저장 |

Snakemake·MAFFT·htslib은 Windows판이 없어서 Windows에서는 WSL2나 Docker로 실행해요. 둘 다 BIOS에서 가상화가 켜져
있어야 하고, 처음 설치할 때 관리자 권한이 필요해요.

주요 옵션: `--panels 목록|all|everything` (`all`은 대용량인 mpox·인플루엔자 제외), `--dir`(기본 `~/chichipes-ui`),
`--port-base`(목록 화면 포트, 기본 8760), `--no-data`, `--primers 파일`(파이프라인의 `tools/export_primer_sets.py`로
내보낸 세트), `--update`. 설치 폴더에 `start.sh`, `stop.sh`, `fetch_data.sh`가 생기고, 데이터는 백그라운드로 받아요
(`logs/data_status.tsv`). 공개 패널에는 프라이머 세트가 없으니 각 패널의 프라이머 세트 화면에서 넣거나 `--primers`로
가져오세요.

## 시퀀싱 리드 검사 (FASTQ)

`chichipes reads`는 유전체 FASTA를 DB로 쓰듯 FASTQ를 DB로 써서, 어떤 리드에 각 프라이머·프로브가 들어 있는지 알려줘요.
single-end, paired-end, 롱리드 모두 되고, FASTQ/FASTA, gzip 여부는 자동으로 알아봐요.

```bash
chichipes reads -i oligos.fa -r sample.fastq.gz -o sample                      # single-end / 롱리드
chichipes reads -i oligos.fa -1 S_R1.fastq.gz -2 S_R2.fastq.gz -o sample        # paired-end (한 쌍을 한 번으로 셈)
chichipes reads -i oligos.fa -r ont.fastq.gz -o ont -m 3 --indel                # 오류가 많은 롱리드 (ONT, PacBio CLR)
```

올리고 이름을 `SET_F`, `SET_R`, `SET_P`로 주면 세트별 집계도 나와요.

| 파일 | 내용 |
|---|---|
| `<prefix>.summary.tsv` | 올리고별: 들어 있는 리드(쌍) 수, %, 미스매치 0·1·2…개별 수, 가닥 |
| `<prefix>.sets.tsv` | 세트별: F, R, P가 있는 리드, F와 R 모두, F·R·P 모두 |
| `<prefix>.positions.tsv` | 올리고 위치별: 리드가 그 자리에 가진 염기(A/C/G/T/기타/결실)와 일치율 — 변이가 어느 자리에 있는지 보여줌 |
| `<prefix>.hits.tsv` | `--hits`일 때: 리드·올리고마다 한 줄 (위치, 가닥, 미스매치, 리드 염기) |

`-m`은 올리고당 미스매치 수(치환; `--indel`이면 편집 거리)예요. 올리고 전체가 리드 안에 있어야 하고, IUPAC 코드와 이노신은
해당 염기 모두와 맞고, 리드의 `N`은 미스매치예요. 시뮬레이션한 Illumina single-end 2만, paired-end 1만 쌍, ONT형 롱리드
400개(오류 4 %, `--indel`)에서 찾은 리드가 전수 검사 결과와 정확히 같았어요. gzip 150 bp 리드 100만 개에 약 2초(16 스레드).

## 프라이머·프로브 설계

`chichipes design`은 계통 정렬(MAFFT 등으로 만든 aligned FASTA)에서 F·R·프로브를 골라, 통과하는 유전체가 가장 많게
설계해요. 통과 기준은 패널 검증과 같아요.

- 프라이머: 미스매치 `--mm`개 이하이고, 3′ 끝 `--tp`개 염기는 정확히 맞아야 해요.
- 프로브: 완전 일치해야 해요.
- IUPAC 코드는 해당하는 모든 염기와 맞는 것으로 봐요.

```bash
chichipes design cover -a aligned.fa --top 3                       # 커버리지 우선 설계
chichipes design sweep -a aligned.fa --grid mega -o results.jsonl  # 매개변수 조합 수천 개
chichipes design shift -a aligned.fa --sets primers.tsv --shift 50 # 기존 세트 개선
chichipes design tile  -a aligned.fa --step 10                     # 게놈 전체 커버리지 지도
chichipes design labels --sets primers.tsv --groups set_groups.tsv # 형광·퀀처 검사
```

### 명령별 설명

**cover**
- 정렬의 모든 창마다 합의 서열에 IUPAC 코드(2염기)를 최대 `--ambig`개 넣어요. 코드는 통과 유전체가 가장 많이 늘어나는 자리에 들어가요.
- 각 후보가 통과시키는 유전체를 비트로 기록하고, F·R·프로브를 합쳐 통과 수가 최대인 조합을 골라요.
- 그다음 세트는 앞 세트가 놓친 유전체를 대상으로 다시 골라요(탐욕적 집합 덮기).
- 패널의 `scripts/redesign_coverage.py`를 그대로 옮긴 것이에요.
  - LASV S·L의 `hightm` 그리드 240회 실행 모두 설계 결과가 모든 항목(서열, 커버리지, Tm, GC)까지 같았어요.
  - 실행당 2스레드에서 S는 22.7초 → 1.4초, L은 16.6초 → 0.8초예요.
  - Rust 엔진은 이전 Go 엔진과 출력 파일이 같고(LASV에서 cover, shift, tile, 배제, labels), 약 1.5배 빨라요.

**sweep**
- `design_sweep.py`의 그리드(`wide`, `hightm`, `midtm`, `longoligo`, `mega`, `mega2`, `mega3`)를 그대로 돌려요.
- 같은 형식의 `results.jsonl`을 쓰고, 이미 있는 실행은 건너뛰어요.

**shift** (전수 탐색)
- 세트의 각 올리고를 양쪽으로 최대 `--shift` bp까지 `--step` 간격으로 옮겨요.
- 길이 범위의 모든 길이와 IUPAC 코드 0~N개를 모두 만들어, 각 변형을 모든 유전체에 대해 평가해요.
- 올리고별로 상위 `--keep`개를 남긴 뒤 전부 조합해서, 바뀐 올리고 종류(F, R, P, FR, FP, RP, FRP)마다 최고 세트를 보고해요.
- 조건: 각 프라이머의 Tm은 원래 값 ±3 °C 안, 프로브는 프라이머보다 `--probe-above` 이상 높아야 해요.
- LASV S 사용자 세트 25개 결과:
  - ±8 bp: 모든 세트에서 `optimize_user_sets.py` 이상이었어요(45초 → 3.2초).
  - ±50 bp: 4.9초 걸렸고, 여러 세트가 약 0 %에서 55–68 %로 올라갔어요.

**tile** (전수 탐색)
- 모든 보존 창, 모든 길이, IUPAC 코드 0~N개를 전부 평가해요.
- `--step` 간격으로 `[x, x+amp_max)` 구간마다 최고 세트를 찾아요.
- `tile.tsv`는 게놈을 따라 본 커버리지 지도이고, `tile.panel.tsv`는 함께 쓸 때 가장 많이 덮는 상위 `--top` 세트예요.
- LASV S 전체(655개 유전체)에 1.4초 걸렸어요.

### 비표적 배제 (`--exclude other.fa`)

다른 바이러스, 숙주 등을 넣어요.
- `cover`는 비표적 서열에서 프로브까지 포함된 증폭산물이 생기는 세트를 건너뛰고 다음 세트를 골라요.
- 모든 결과에 `off_target` 열이 붙어요: 산물 수, 프로브 포함 산물 수, 올리고별 결합 수.
- 결합 기준은 미스매치 `--mm`개 이하이고 3′ 끝이 맞는 경우예요. 리드 검사와 같은 시드 탐색을 써요.
- 확인 결과: 남미 아레나바이러스 S 129개에서는 0건, LASV 자신을 넣으면 산물을 모두 잡았어요.

### 세트별 검사

- chichipes Tm(SantaLucia & Hicks 2004, spec 검색 완충액)
- 가장 안정한 헤어핀
- 다이머: 프로브의 3′은 퀀처로 막혀 있으므로, 연장될 수 있는 다이머는 프라이머 3′ 끝이 관여한 것만 봐요.
- varVAMP식 페널티(낮을수록 좋음)와 경고

### 형광·퀀처

`--probe-chem`으로 프로브 형식을 정해요.

| 형식 | 예 | 프로브 길이 | 프로브 Tm 보정 |
|---|---|---|---|
| `linear` (기본) | FAM / BHQ1, FAM / TAMRA | 30 nt 이하 (형광과 퀀처가 가까워야 소광됨) | 0 |
| `zen` | FAM / ZEN / IBFQ (내부 퀀처) | 40 nt 이하 | 0 |
| `mgb` | FAM / MGB-NFQ | 13–22 nt | +10 °C (대략값이니, 쓰고 있는 프로브의 제조사 Tm으로 `--probe-tm-offset`을 맞추세요) |

모든 형식에 공통으로 적용되는 규칙:
- 프로브 5′ 끝에 G를 두지 않아요(형광을 소광해요).
- 같은 염기가 4개 이상 연속되지 않게 해요.
- 프로브는 C가 G보다 많은 가닥으로 둬요. 완전 일치 프로브는 어느 가닥이든 통과하는 유전체가 같아서 커버리지는 그대로예요.

`design labels`는 `fluorophore` 열을 읽어서 두 가지를 경고해요.
- 형광에 맞지 않는 퀀처(예: FAM / BHQ2)
- 한 튜브에서 돌린다고 할 때 같은 채널을 쓰는 그룹 내 세트

`6FAM / MGBNFQ`, `SUN`, `CY5` 같은 제조사 표기도 알아봐요.

## 배포 (릴리즈)

```bash
git tag chichipes-v0.5.1 && git push origin chichipes-v0.5.1
```

`.github/workflows/release.yml`은 `chichipes` 브랜치에 푸시할 때마다 테스트하고 빌드해요(Linux amd64/arm64는
`cargo zigbuild`로 glibc 2.17 대상, macOS amd64/arm64, Windows amd64; 압축 파일은 워크플로 아티팩트로 남음).
`chichipes-v*` 태그를 올리면 그 압축 파일을 `SHA256SUMS.txt`와 함께 GitHub Release에 올려요. 로컬에서는
`cargo build --release && scripts/release.sh v0.5.1 linux amd64 target/release/chichipes`로 `dist/`에 만들어요.

## 라이선스와 인용

GPL-3.0-or-later ([LICENSE](LICENSE)), 제3자 고지는 [THIRD_PARTY_NOTICES.md](THIRD_PARTY_NOTICES.md).
chichipes는 MFEprimer와 관계없는 독립 구현이에요. MFEprimer 호환 모드로 얻은 결과를 발표할 때는 MFEprimer도 인용해 주세요:
Wang K. et al., *Nucleic Acids Res.* 2019;47(W1):W610–W613, https://doi.org/10.1093/nar/gkz351.
