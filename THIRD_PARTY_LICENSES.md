# Third-party software, data and licenses

이 저장소(GPL-3.0-or-later, `LICENSE`)는 **우리가 작성한 코드만** 담습니다. 아래 도구·라이브러리·데이터는 저장소에 포함되지 않거나(설치/다운로드는 사용자가 수행),
포함된 경우 이 표에 출처와 라이선스를 적었습니다. 검토일 2026-10-14 (설치된 환경의 conda-meta와 각 프로젝트 공식 페이지 기준 — 공개 직전에 한 번 더 확인하세요).

## 1. 저장소에 포함된 서드파티 파일 (`webapp/static/vendor/`)

| 파일 | 버전 | 라이선스 | 저작권 표기 |
|---|---|---|---|
| `d3.min.js` | d3 7.9.0 | ISC | Copyright 2010-2023 Mike Bostock |
| `topojson-client.min.js` | 3.1.0 | ISC | Copyright 2012-2019 Michael Bostock |
| `countries-110m.json` | world-atlas 계열 | ISC (코드) · 지도 원본 Natural Earth | Copyright 2013-2019 Michael Bostock |

ISC 라이선스는 복사본에 저작권 표기와 아래 허가문을 함께 두도록 요구합니다.

```
Permission to use, copy, modify, and/or distribute this software for any purpose with or without fee is hereby
granted, provided that the above copyright notice and this permission notice appear in all copies.

THE SOFTWARE IS PROVIDED "AS IS" AND THE AUTHOR DISCLAIMS ALL WARRANTIES WITH REGARD TO THIS SOFTWARE INCLUDING ALL
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS. IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR ANY SPECIAL, DIRECT,
INDIRECT, OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR PROFITS, WHETHER IN
AN ACTION OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING OUT OF OR IN CONNECTION WITH THE USE OR
PERFORMANCE OF THIS SOFTWARE.
```

`countries-110m.json`의 국경선 데이터 출처인 Natural Earth는 퍼블릭 도메인이라고 알려져 있으나 이 저장소에서 직접 확인하지는 않았습니다
(<https://www.naturalearthdata.com/about/terms-of-use/>). 국경선 표기는 지도 표시용이며 정치적 의미가 없습니다.

### PearTree (저장소에 포함하지 않음)
`webapp/static/vendor/peartree.bundle.min.js`(PearTree, MIT, ARTIC network)는 `.gitignore`로 제외되어 있습니다. 번들 안에는 Bootstrap 5.3(MIT),
Bootswatch(MIT, Thomas Park), Bootstrap Icons 1.11(MIT)이 들어 있습니다. MIT이므로 재배포가 가능하지만, **직접 받도록** 해 두었습니다:

```bash
curl -L -o webapp/static/vendor/peartree.bundle.min.js \
  https://github.com/artic-network/peartree/releases/latest/download/peartree.bundle.min.js
```

주의: 이 번들은 실행 시 `fonts.googleapis.com`에서 Source Sans Pro 글꼴을 불러옵니다(인터넷 필요, 접속 기록이 Google에 남음). 폐쇄망이나
개인정보에 민감한 환경이면 글꼴을 내려받아 self-host 하거나 번들의 `@import`를 제거하세요.

## 2. 실행 시 사용하는 외부 프로그램 (별도 프로세스로 호출, 포함·링크하지 않음)

| 프로그램 | 버전 | 라이선스 | 비고 |
|---|---|---|---|
| **MFEprimer** | 4.2.4 | **명시된 라이선스 파일 없음** | 아래 §3 참고 — **재배포 금지로 간주** |
| Snakemake | 9.27.0 | MIT | |
| MAFFT | 7.526 | BSD-3-Clause | conda-forge 빌드(기본 소스). 확장판(`with-extensions`)의 별도 조건은 해당 없음 |
| pixy | 2.2.3 | MIT (Kieran Samuk) | PyPI 이름이 다른 패키지와 겹쳐 git 커밋으로 고정 설치 (`environment.yml`) |
| scikit-allel / numcodecs / dask | 1.3.13 / 0.16.5 / 2026.8.0 | MIT / MIT / BSD-3-Clause | pixy의 의존성 |
| htslib (bgzip, tabix) | 1.23.1 | MIT/BSD | |
| zstd | 1.5.7 | BSD-3-Clause (GPL-2 이중) | `.zst` 해제 |
| PyYAML | 6.0.3 | MIT | |
| pytest | 9.1.1 | MIT | 테스트용 |
| Python | 3.12 | PSF-2.0 | 우리 코드는 표준 라이브러리만 import |
| Chrome / Chromium | 사용자 설치 | 브라우저별 약관 | PDF 리포트용, 포함하지 않음 |

**Nipah 클레이드 재빌드용 환경**(`environment_nextstrain.yml`, 선택): Augur **AGPL-3.0**, Nextclade MIT, nextstrain-cli MIT, IQ-TREE **GPL-2.0-or-later**,
seqkit MIT, git GPL-2.0-or-later. 이 도구들은 Nextstrain 워크플로가 명령줄로 호출할 뿐 우리 코드에 링크되지 않으므로 우리 코드의 라이선스(MIT)에는
영향이 없습니다. 다만 **Augur/IQ-TREE가 들어 있는 환경이나 이미지를 통째로 배포**하면 그 쪽 소스 제공 의무가 생기므로, 환경은 배포하지 말고
`environment*.yml`만 배포하세요.

주 환경(`environment.yml`) 전이 의존성 146개를 점검한 결과: MIT 64, BSD 34, Apache-2.0 7, 그 밖은 ISC/PSF/Zlib 등 허용형입니다. 카피레프트는
컴파일러 런타임(libgcc 등, GCC 예외 조항 있음), readline/gawk(GPL-3), libiconv/gmp/mpfr(LGPL), 최적화 솔버 COIN-OR(EPL-2.0, Snakemake 의존성)이며
모두 conda가 설치하는 별도 바이너리라 우리 저장소에는 포함되지 않습니다.

## 3. MFEprimer — 가장 주의할 항목

- 저장소에 포함하지 않았고(`mfeprimer_bin`으로 경로만 지정, `setup.sh`가 경로를 물음), 앞으로도 **포함·미러링·재업로드하지 마세요**.
- 공식 사이트(<https://www.mfeprimer.com/license/>)는 "명령줄 도구와 온라인 서버는 비영리·영리 모두 무료"라고만 밝히며, 재배포 조건이나 OSI 오픈소스
  라이선스 전문은 확인되지 않았습니다. 바이너리에도 LICENSE 파일이 없습니다. 재배포·상업적 번들링이 필요하면 저자(quwubin@gmail.com)에게 문의하세요.
- **논문에는 반드시 인용**: Wang K, Li H, Xu Y, Shao Q, Yi J, Wang R, Qu W. *MFEprimer-3.0: quality control for PCR primers.* Nucleic Acids Res. 2019.
  doi:10.1093/nar/gkz351 (프로그램 도움말이 요구하는 인용). 필요하면 2.0(doi:10.1093/nar/gks552), 1.0(doi:10.1093/bioinformatics/btn614)도.
- 과거에 만든 `pipeline_snakemake_demo`(환경 tarball + mfeprimer 바이너리 동봉)와 `packaging/` 꾸러미는 **공개용이 아닙니다**.

## 4. 데이터

| 데이터 | 취급 |
|---|---|
| Nextstrain 데이터 (`data.nextstrain.org`: sequences/metadata, GenBank 유래) | 저장소에 포함하지 않음(`data/`는 gitignore). 실행 때 내려받음. Nextstrain 문서는 다운로드물에 별도 라이선스를 명시하지 않으며, **미발표 데이터는 제공자 허락하에 포함**되어 있어 추가 연구 전 저자(메타데이터의 authors 열)에게 연락하라고 안내합니다. 논문에는 사용한 accession 목록을 보충자료로 싣고 GenBank·원 저자·Nextstrain(Hadfield et al. 2018, *Bioinformatics* 34:4121)을 인용하세요. |
| NCBI GenBank / RefSeq (NC_002728 등) | INSDC 공개 데이터. NCBI는 이용 제한을 두지 않지만 제출자가 권리를 주장할 수 있다는 점을 고지함 → 서열 자체를 재배포하지 말고 accession으로 참조 |
| **nextstrain/nipah 저장소** (`clades.tsv`, `exclude.txt`, `genemap.gff`, 빌드 결과 `nipah_all.json`) | **LICENSE 파일이 없음**(GitHub 라이선스 조회 404, README에도 없음) → 기본 저작권이 유지되는 것으로 보고 `data/nextstrain/`(gitignore)에만 두며 **저장소에 올리지 않음**. 스크립트(`run_nextstrain_nipah.sh`)가 실행 때 GitHub에서 클론함. 인용하고 싶다면 저장소 URL·커밋(`PROVENANCE.txt`)을 적고, 재배포가 필요하면 upstream에 라이선스를 문의하세요. |
| 문헌 프라이머 (대시보드 "📚 문헌 프라이머", 로컬 `config/literature/*.tsv`, 저장소에는 없음) | 논문에서 가져온 짧은 서열과 출처. 항목마다 출처가 필수이며 '원문 대조 완료'로 표시된 것만 검증된 서열입니다. `examples/literature/`에는 재배포가 허용된 논문(Jensen 2018 CC0, Sendow 2013 CC BY, Ranadheera 2018 CC BY 4.0)의 세트만 출처·이용 조건(`SOURCES.md`)과 함께 포함합니다. 그 밖의 서열(저작권이 있는 논문 포함)은 올리지 않습니다. 특허·상용 사용 여부는 확인하지 않았습니다. |
| **우리 프라이머/프로브, KDCA 표준물질, 벤더(VendorD·VendorA 등) 세트** | **비공개**. `config/*`, `data/`, `output/`, `input_primers/`가 gitignore이며 git 이력에도 없음을 확인함. 논문 그림·스크린샷·`output/`·PDF 리포트에 서열이 보이지 않는지 공개 전에 점검하세요. |

## 5. 공개 전 체크리스트

1. `git log --all -p | grep -E "^\+.*[ACGTRYSWKMBDHVN]{22,}"`로 이력에 비공개 서열이 없는지 다시 확인 (2026-10-14 확인 결과: 이력의 긴 서열은 모두 논문에 공개된 프라이머와 그 출처 인용뿐).
2. (처리됨) 예전의 "대외비 — 외부 공유 금지" 배너는 제거했고, 화면 상단과 PDF 표지에는 "in silico 예측 — 실험으로 확인 필요" 안내가 표시됩니다.
3. 논문 Methods/Acknowledgements에 §2–§4의 인용을 넣으세요: MFEprimer, MAFFT(Katoh & Standley 2013), pixy(Korunes & Samuk 2021), Snakemake(Mölder 2021), Nextstrain(Hadfield 2018), Nextclade/Augur(Nipah 클레이드), IQ-TREE(Nipah 클레이드), 사용한 데이터의 GenBank 원 저자.
4. MFEprimer 재배포 여부는 저자 확인 전까지 "사용자가 직접 설치"로 유지.
5. 이 파일의 버전·라이선스 칸은 환경을 갱신할 때마다 다시 확인.
