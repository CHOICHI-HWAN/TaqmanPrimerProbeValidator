# chichipes

[English](README.md)

**chichipes**는 TaqMan 검증 파이프라인이 쓰는 [MFEprimer 4.2.4](https://www.mfeprimer.com) 명령(프라이머쌍 증폭산물
검색과 바이러스 커버리지, 프로브 검색, 다이머·헤어핀 검사)을 그대로 대신하는 실행 파일 하나짜리 엔진이에요.
명령줄도, 출력 파일도, 숫자도 같고, 보통 수십~수백 배 빨라요.

[TaqMan Primer/Probe Validator](https://github.com/CHOICHI-HWAN/TaqmanPrimerProbeValidator) 파이프라인의 검색 엔진이며,
엔진만 필요한 분을 위해 별도 브랜치(`chichipes`)로 배포해요.

- **의존성 없음** — Go 표준 라이브러리만 사용. Linux, macOS, Windows용 단일 실행 파일.
- **색인 불필요** — DB를 바로 읽어요. `chichipes index`는 다른 도구가 기대하는 `.fai`/`.json`만 만들어요.
- **클린룸 구현** — 공개 열역학 논문과 MFEprimer 출력 비교 실험만으로 규칙을 찾았어요. MFEprimer 코드를 쓰거나
  디컴파일하지 않았어요. 자세한 규칙: [docs/FINDINGS.md](docs/FINDINGS.md).

## 설치

[Releases](https://github.com/CHOICHI-HWAN/TaqmanPrimerProbeValidator/releases)의 `chichipes-v*` 릴리즈에서 운영체제에 맞는 압축 파일을 받아
`chichipes`를 `PATH`에 넣으면 돼요.

```bash
tar xzf chichipes-v0.1.0-linux-amd64.tar.gz
sudo install chichipes-v0.1.0-linux-amd64/chichipes /usr/local/bin/
chichipes version
```

소스에서 빌드 (Go 1.22 이상):

```bash
git clone -b chichipes https://github.com/CHOICHI-HWAN/TaqmanPrimerProbeValidator.git chichipes
cd chichipes
go build -o chichipes ./cmd/chichipes
go test ./...
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

## 배포 (릴리즈)

```bash
git tag chichipes-v0.1.1 && git push origin chichipes-v0.1.1
```

`chichipes-v*` 태그를 올리면 `.github/workflows/release.yml`이 테스트 → `scripts/release.sh`로 Linux/macOS/Windows 압축 파일
빌드 → `SHA256SUMS.txt`와 함께 GitHub Release에 올려요. 로컬에서는 `scripts/release.sh v0.1.1`로 `dist/`에 같은
파일을 만들 수 있어요.

## 라이선스와 인용

GPL-3.0-or-later ([LICENSE](LICENSE)), 제3자 고지는 [THIRD_PARTY_NOTICES.md](THIRD_PARTY_NOTICES.md).
chichipes는 MFEprimer와 관계없는 독립 구현이에요. MFEprimer 호환 모드로 얻은 결과를 발표할 때는 MFEprimer도 인용해 주세요:
Wang K. et al., *Nucleic Acids Res.* 2019;47(W1):W610–W613, https://doi.org/10.1093/nar/gkz351.
