# Nipah Primer/Probe Validation

**TaqMan 검출 세트(프라이머 2개 + 프로브 1개)가 알려진 Nipah virus(NiV) 유전체를 얼마나 잡는지** 계산하고, 결과를 웹 대시보드로 보여주는 도구입니다.
새 프라이머를 설계하거나 기존 키트를 점검할 때, "이 세트가 현재 알려진 유전체 중 몇 %를 검출하나, 못 잡는 건 어떤 계통인가"를 확인하는 용도입니다.

비분절 유전체를 다루며, 클레이드는 공식 nextstrain/nipah 워크플로를 로컬 Nextstrain CLI로 빌드한 트리에서 가져옵니다(부분서열은 QC FAIL로 표시).

## 하는 일

1. 공개 유전체(서열 + 메타데이터)를 내려받아 품질 필터와 방향 정규화를 합니다.
2. 각 세트를 모든 유전체에 대해 검색합니다(**mfeprimer** + 자체 스캔, mismatch는 게놈 염기로 직접 재계산).
3. **검증 통과** = F/R mismatch ≤ 1(3′ 끝 일치) **그리고** 프로브 완전 일치(프라이머 사이에 위치).
4. 커버리지를 계통별(클레이드 I / Ia / Ib / II), 숙주별(기본 **사람**, 동물·미상·전체로 전환)로 대시보드에 보여줍니다. 표준물질 매치, 정렬, 다양성(π)도 함께 봅니다.

![숙주 기준 전환](docs/media/host-scope.gif)
![세트 상세](docs/media/set-detail.gif)
![표준물질 매치](docs/media/standard-match.gif)

> 업체명·서열·표준물질 이름은 `docs/demo_blur.js`로 흐리게 처리해서 촬영했습니다. 빨강 = 프로브 또는 프라이머 3′ 끝 mismatch, 노랑 = 결합에 지장 없는 프라이머 mismatch.

## 필요한 것

| | |
|---|---|
| OS | Linux x86_64 (Ubuntu 24.04에서 개발·테스트) |
| conda / mamba / micromamba | `setup.sh`가 환경(snakemake, mafft, pixy, htslib …)을 만듭니다 |
| **mfeprimer 4.2.4** | **저장소에 없습니다. 따로 받아서 경로를 지정해야 합니다** → 아래 "mfeprimer 설정" |
| 인터넷 | 데이터 내려받기와 계통수 뷰어(PearTree) 받기 |
| (선택) Chrome | PDF 보고서 내보내기만 |
| (선택) Nextstrain CLI | 클레이드별 커버리지가 필요할 때만: `./setup.sh --nextstrain` 후 `scripts/run_nextstrain_nipah.sh` |

## 빠른 시작

```bash
git clone https://github.com/CHOICHI-HWAN/TaqmanPrimerProbeValidator.git
cd TaqmanPrimerProbeValidator
git checkout nipah         # 이 프로젝트는 nipah 브랜치입니다

./setup.sh                 # 1) 환경 만들기 + mfeprimer 경로 지정 (대화형)
./test.sh                  # 2) 설치 점검 (약 10초): 모든 항목이 ✔ 이어야 합니다
./start.sh                 # 3) 대시보드 → http://localhost:8767
scripts/run_nextstrain_nipah.sh   # (선택) 클레이드별 커버리지가 필요할 때 1회: ./setup.sh --nextstrain 후 실행
```

그다음 대시보드에서:
1. **"📝 프라이머 입력 / 관리"**에서 검증할 세트를 입력합니다(아래 형식). 표준물질이 있으면 **"🧫 표준물질 입력 / 관리"**에도 넣습니다(선택).
2. **"데이터 업데이트"**를 눌러 유전체를 내려받고, **"분석 다시하기"**로 전체를 계산합니다(터미널에서는 `snakemake --cores 4`).
3. **"전체 대시보드"**와 각 세트 화면에서 결과를 봅니다. 새 세트 하나만 빠르게 보려면 **"+ 새 프라이머 세트 분석"**.

대시보드는 기본적으로 **이 컴퓨터에서만**(`127.0.0.1`) 열립니다. 로그인이 없으므로 다른 PC에서 열려면 `--host 0.0.0.0`을 직접 주세요.

### 계통수 뷰어 받기 (선택, 1회)

세트 상세의 "Nextstrain 계통수" 카드에 쓰입니다. 라이선스(MIT) 때문에 저장소에는 넣지 않았습니다.

```bash
curl -L -o webapp/static/vendor/peartree.bundle.min.js \
  https://github.com/artic-network/peartree/releases/latest/download/peartree.bundle.min.js
```

## mfeprimer 설정

mfeprimer는 conda나 pip로 설치되지 않는 **단독 실행 파일**이라 직접 받아야 하고, 재배포가 허용되는지 확인되지 않아 이 저장소에도 넣지 않았습니다.

1. 공식 사이트 <https://www.mfeprimer.com/> 에서 **MFEprimer v4.2.4 (linux-amd64)** 실행 파일을 받습니다.
2. 실행 권한을 주고 동작을 확인합니다.
   ```bash
   chmod +x /경로/mfeprimer-4.2.4-linux-amd64
   /경로/mfeprimer-4.2.4-linux-amd64 version        # 버전이 출력되면 정상
   ```
3. 경로를 알려줍니다. **아래 중 하나**면 됩니다.
   - `./setup.sh` 실행 중 묻는 경로에 입력 (비대화형: `./setup.sh --mfeprimer /경로/mfeprimer-4.2.4-linux-amd64 --yes`)
   - `config/config.yaml`의 `mfeprimer_bin: "/경로/mfeprimer-4.2.4-linux-amd64"`
   - 대시보드의 **"Analysis Parameters"** 페이지에서 입력 (저장할 때 실제로 실행해서 확인합니다)
   - 환경변수 `MFEPRIMER_BIN`
   - 비워 두면 프로젝트 옆 폴더(`../*/mfeprimer/`)나 `PATH`에서 찾습니다.
4. 확인: `./test.sh` 1번 항목에 `mfeprimer (mfeprimer v4.2.4)`가 ✔ 로 나오면 끝입니다.

> 논문에 쓸 때는 MFEprimer-3.0 (Wang K, et al. *Nucleic Acids Res.* 2019, doi:10.1093/nar/gkz351)를 인용하세요.

## 프라이머·프로브 입력 형식

대시보드의 입력 화면(엑셀 다운로드·업로드 가능)이나 `config/primers.tsv`(탭 구분)에 한 줄에 한 세트입니다.

| 열 | 설명 |
|---|---|
| `primer_name` | 표시 이름 (고유한 ID `P001…`는 저장 시 자동 부여) |
| `segment` | `GENOME` (비분절 유전체) |
| `forward_primer`, `reverse_primer`, `probe` | 5′→3′ 서열 (IUPAC 축퇴 염기 가능) |
| `fluorophore`, `source`, `remark` | 형광체, 출처(업체·논문), 메모 |
| `proceed` | `True`면 분석 대상 |

서열과 분석 결과는 `config/`, `data/`, `output/`에만 저장되고 **git에 올라가지 않습니다**(`.gitignore`).

## 더 읽을거리

- [docs/PIPELINE.md](docs/PIPELINE.md): 파이프라인 단계, 설정, 웹 앱 기능, 디렉터리 구조, 숙주 기준 등 자세한 설명
- [CHANGELOG.md](CHANGELOG.md) · [THIRD_PARTY_LICENSES.md](THIRD_PARTY_LICENSES.md) · `profiles/README.md`(다른 바이러스를 추가하는 방법)
- 테스트: `./test.sh` (핵심) / `./test.sh --full` (실제 결과 감사 포함)

> **이 저장소에는 서열이 없습니다.** 프라이머·프로브·표준물질·게놈(FASTA) 서열과 분석 결과는 올리지 않으며(`tests/test_no_sequences_in_repo.py`가 검사),
> 테스트용 염기서열은 `tests/seqgen.py`가 만들어 씁니다. 업체명은 일반 이름(VendorA…)으로 바뀌어 있습니다.

## License

MIT, © 2026 Chihwan Choi. 서드파티(mfeprimer, MAFFT, pixy, Snakemake, Nextstrain …)의 라이선스와 인용은 [THIRD_PARTY_LICENSES.md](THIRD_PARTY_LICENSES.md)를 보세요.
