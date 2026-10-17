# TaqMan Primer/Probe Validator — Nipah 에디션

![License](https://img.shields.io/badge/license-MIT-blue) ![Snakemake](https://img.shields.io/badge/snakemake-%E2%89%A59.27-brightgreen) ![Platform](https://img.shields.io/badge/platform-linux--64-lightgrey)

[English README](README.md)  ·  LASV 에디션: [`main`](https://github.com/CHOICHI-HWAN/TaqmanPrimerProbeValidator/tree/main) 브랜치

TaqMan 프라이머/프로브 세트가 **알려진 Nipah virus(NiV) 유전체를 얼마나 잡는지** 계산하고, 결과를 대화형 대시보드로 보여주는 도구입니다.
새 세트를 설계하거나 기존 키트를 점검할 때 *"이 세트가 알려진 유전체의 몇 %를 검출하나, 어떤 계통을 놓치나?"* 를 확인하는 용도입니다.

Nipah 에디션은 비분절 유전체를 다루며, 클레이드는 공식 nextstrain/nipah 워크플로를 로컬에서 빌드해 얻습니다(부분서열은 QC FAIL로 표시). LASV 에디션은 `main` 브랜치에 있습니다.

자세한 문서(파이프라인 단계, 설정, 대시보드 기능, 디렉터리 구조): [docs/PIPELINE.md](docs/PIPELINE.md)

## 특징

1) **정확한 mismatch 계산.** mfeprimer가 보고하는 프라이머 문자열은 5′ 쪽 mismatch를 가리고 시드 단계에서 놓친 위치도 있습니다. 그래서 모든 위치를 **게놈 염기**로 다시 세고,
   mfeprimer와 독립적인 게놈 스캔을 둘 다 돌려 더 나은 결과를 씁니다. 프로브는 amplicon 안에 있을 때만 인정합니다.

2) **단순하고 명확한 판정 규칙.** *검증 통과* = F/R 프라이머 각각 mismatch ≤ 1(3′ 끝 일치) **그리고** 프로브 완전 일치. 결합은 되지만 기준을 조금 넘는 게놈은 따로(“주의”) 보여줍니다.

3) **숙주를 구분한 커버리지.** 기본은 **사람 유래 게놈** 기준이고, 클릭 한 번으로 동물·미상·전체로 바꿉니다. 실험 숙주 게놈은 어디에서도 제외합니다.

4) **클레이드(I / Ia / Ib / II)별 결과.** 클레이드(I / Ia / Ib / II)별 커버리지, 방향을 맞춘 amplicon 다중 정렬(F/Probe/R 위치, 참조 좌표, 뉴클레오타이드 다양성 π 막대).

5) **여러 세트 묶음(패널).** 합집합 커버리지와 **교차 결합**(세트 1의 forward + 세트 2의 reverse)을 확인합니다.

6) **표준물질 확인**(일부만 매치되는 경우 포함; 빨강 = 프로브 또는 프라이머 3′ 끝 mismatch, 노랑 = 감당되는 프라이머 mismatch), dimer·hairpin 점검, 유전자 보존영역, 엑셀 입출력.

7) **자체 점검.** `./test.sh`는 정답을 아는 합성 게놈을 실제 mfeprimer로 돌려 핵심 기능을 확인하고, `./test.sh --full`은 실제 결과까지 감사합니다.

## 개략도

![개략도](docs/schematic/pipeline_schematic.svg)

## 대화형 대시보드

![개요](docs/media/overview.png)

![정렬](docs/media/alignment.png)

<details><summary>움직이는 화면 (숙주 기준 · 세트 상세 · 표준물질)</summary>

![숙주 기준](docs/media/host-scope.gif)

![세트 상세](docs/media/set-detail.gif)

![표준물질](docs/media/standard-match.gif)

</details>

*(업체명, 프라이머·프로브·표준물질 서열과 표준물질 이름은 `docs/demo_blur.js`로 흐리게 처리했습니다.)*

## 설치

### 현재 버전

코드: **v0.19.0 (Nipah fork)**   ·   Ubuntu 24.04, x86_64에서 테스트

### 필요한 것

| | |
|---|---|
| conda / mamba / micromamba | 환경(snakemake, mafft, pixy, htslib …)을 만듭니다 |
| **mfeprimer 4.2.4** | **저장소에 없고 conda로도 설치되지 않습니다. 직접 설치하세요(아래 3–4단계)** |
| 인터넷 | 유전체 내려받기, 선택 항목인 계통수 뷰어 |
| Chrome (선택) | PDF 보고서 내보내기만 |
| Nextstrain CLI (선택) | 클레이드별 커버리지가 필요할 때만: `./setup.sh --nextstrain` 후 `scripts/run_nextstrain_nipah.sh` |

**1)  코드 받기**

```bash
git clone https://github.com/CHOICHI-HWAN/TaqmanPrimerProbeValidator.git
cd TaqmanPrimerProbeValidator
git checkout nipah
```

**2)  환경 만들기** (4단계의 mfeprimer 경로를 물어봅니다. 다시 실행해도 안전합니다)

```bash
./setup.sh                       # = mamba env create -f environment.yml 후 경로를 config/config.yaml에 기록
conda activate primer_validation # snakemake를 직접 돌릴 때만 필요, ./start.sh는 필요 없음
```

**3)  mfeprimer 내려받기** — 단독 실행 파일(MFEprimer v4.2.4, linux-amd64). 공식 사이트 <https://www.mfeprimer.com/> 에서 받으세요.
라이선스에 재배포 허용이 명시되어 있지 않아 이 저장소에는 넣지 않았습니다.

```bash
chmod +x /경로/mfeprimer-4.2.4-linux-amd64
/경로/mfeprimer-4.2.4-linux-amd64 version     # 버전이 나오면 정상
```

**4)  mfeprimer 경로 알려주기** (아래 중 하나):

- `./setup.sh`가 묻는 경로에 입력, 또는 `./setup.sh --mfeprimer /경로/mfeprimer-4.2.4-linux-amd64 --yes`
- `config/config.yaml`에 `mfeprimer_bin: "/경로/mfeprimer-4.2.4-linux-amd64"`
- 대시보드 **"Analysis Parameters"** 페이지에서 입력 (저장할 때 한 번 실행해서 확인)
- 환경변수 `MFEPRIMER_BIN=/경로/mfeprimer-4.2.4-linux-amd64`
- 비워 두면 `../*/mfeprimer/` 또는 `PATH`에서 찾습니다

**5)  설치 점검** (약 10초, 모든 줄이 ✔ 이어야 합니다)

```bash
./test.sh
```

**6)  선택: 계통수 뷰어** ("Nextstrain 계통수" 카드용, MIT 라이선스라 번들하지 않았습니다)

```bash
curl -L -o webapp/static/vendor/peartree.bundle.min.js \
  https://github.com/artic-network/peartree/releases/latest/download/peartree.bundle.min.js
```

# 실행하기

*필수 입력:* **프라이머/프로브 세트** 하나 이상(그리고 선택적으로 **표준물질** 서열). 유전체는 자동으로 내려받습니다.

**대시보드 시작**

```bash
./start.sh            # http://localhost:8767   (127.0.0.1에서만 열림 — 로그인이 없습니다)
```

**프라이머/프로브 세트 입력**: 대시보드 **"📝 프라이머 입력 / 관리"** 페이지(엑셀 다운로드·업로드 가능) 또는 `config/primers.tsv`(탭 구분):

| 열 | 의미 |
|---|---|
| `primer_name` | 표시 이름 (고유 ID `P001, P002 …`는 저장 시 자동 부여) |
| `segment` | `GENOME` (비분절 유전체) |
| `forward_primer`, `reverse_primer`, `probe` | 5′→3′ 서열, IUPAC 축퇴 염기 가능 |
| `fluorophore`, `source`, `remark` | 형광체, 출처(업체·논문), 메모 |
| `proceed` | `True`면 분석에 포함 |

**유전체 내려받기 + 전체 분석**: 사이드바의 **"데이터 업데이트"**, **"분석 다시하기"** 버튼, 또는 터미널에서

```bash
snakemake --cores 4
```

**세트 하나만 다시 분석**

```bash
snakemake --cores 4 output/primers/P001/final/validated_ids.txt
```

**새 세트 하나를 빠르게 분석**(파이프라인 실행 없이): 사이드바 **"+ 새 프라이머 세트 분석"**.

**숙주 기준 선택**: 사이드바 "숙주 기준"(사람 · 동물 · 미상 · 전체). 스크립트에서도 같은 선택을 쓸 수 있습니다.

```bash
curl 'http://localhost:8767/api/sets?host=animal'
```

**신뢰하는 네트워크의 다른 컴퓨터에서 대시보드 열기**

```bash
python3 webapp/server.py --host 0.0.0.0 --port 8767
```

## 결과 파일

| 파일 | 내용 |
|---|---|
| `output/summary/summary.tsv` | 세트별 커버리지와 검증 통과율 (기본 숙주 기준, `host_scope` 열) |
| `output/summary/summary_by_lineage.tsv` | 클레이드(I / Ia / Ib / II)별 검증 통과 게놈 수 |
| `output/summary/group_summary.tsv` | 패널별 합집합 커버리지 |
| `output/primers/<ID>/alignment/mismatch_table.tsv` | 게놈별 F / 프로브 / R mismatch, 3′ 끝 상태, 위치 |
| `output/primers/<ID>/final/validated_ids.txt` | 통과한 accession |
| `output/primers/<ID>/msa/msa.json` | 대시보드에 보이는 amplicon 정렬 |

## 이 저장소에 없는 것

프라이머·프로브·표준물질·게놈(FASTA) **서열과 분석 결과는 커밋하지 않습니다.** `config/`, `data/`, `output/`에 두고(git 제외),
서열이 섞이면 `tests/test_no_sequences_in_repo.py`가 실패합니다. 테스트의 염기서열은 `tests/seqgen.py`가 만들어 씁니다. 공개 저장소의 업체명은 일반 이름(VendorA …)입니다.

# 인용

아직 이 도구의 논문은 없습니다. 이 저장소와 함께 바탕이 된 도구를 인용해 주세요:
**MFEprimer-3.0** (Wang K, et al. *Nucleic Acids Res.* 2019, doi:10.1093/nar/gkz351), MAFFT (Katoh & Standley 2013), pixy (Korunes & Samuk 2021), Snakemake (Mölder et al. 2021), Nextstrain (Hadfield et al. 2018).
서드파티 라이선스: [THIRD_PARTY_LICENSES.md](THIRD_PARTY_LICENSES.md) · [CHANGELOG.md](CHANGELOG.md) · MIT 라이선스, © 2026 Chihwan Choi.
