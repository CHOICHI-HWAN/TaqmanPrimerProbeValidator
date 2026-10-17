# 요구사항 기록 (Nipah) (프로젝트 시작 ~ 현재)

갱신일 2026-10-14 · 총 12건 — 완료 12

> 출처: 프로젝트 시작(2026-09-14)부터의 대화 기록에서 사용자가 요구한 사항을 시간순으로 정리했습니다. 상태는 당시 세션의 처리 결과와 CHANGELOG, 점검 결과 기준이며 오래된 항목을 하나하나 다시 실행해 확인한 것은 아닙니다. 요구 사항 칸의 괄호 번호는 사용자가 직접 붙인 번호입니다.


## 1. Nipah 프로젝트 시작 (10-14)

| # | 날짜 | 구분 | 요구 사항 | 검토 · 수정 내용 | 상태 | 비고 |
|---|---|---|---|---|---|---|
| 1 | 10-14 | 설계 | Nipah를 추가하되 LASV 템플릿은 그대로 두고, 새 폴더에서 새롭게 진행한 뒤 나중에 통합 | `NIPAH/pipeline_nipah`를 LASV 스냅샷(코드만)에서 새 git 저장소로 시작. 변경은 병합을 염두에 두고 프로필(profiles/nipah.json)로 표현 | 완료 | LASV 폴더는 건드리지 않음 |
| 2 | 10-14 | 파이프라인 | 데이터 소스: data.nextstrain.org/files/workflows/nipah/{sequences.fasta,metadata.tsv}.zst, 표준물질·프라이머는 새로 입력 | 프로필 `source`(segment 폴더 없는 flat 구조 지원)로 정의, 수집 규칙이 프로필의 URL 사용. 표준물질/프라이머는 아직 비어 있음(입력 페이지에서 등록) | 완료 | 프라이머·표준물질 서열은 받는 대로 입력 |
| 3 | 10-14 | 파이프라인 | CLADE 정의는 nextstrain 파이프라인을 거쳐야 알 수 있음 → github의 nipah 워크플로를 받아 실제로 돌린 결과를 사용, auspice는 웹에서 받지 말고 로컬 nextstrain CLI로 생성 | 공식 nextstrain/nipah를 clone해 `nextstrain build --native`로 all 빌드 실행(전용 conda 환경 nextstrain_nipah). 만들어진 Auspice JSON에서 tip별 clade_membership(I 3, Ia 73, Ib 20, II 26)을 표로 추출해 메타데이터에 결합. 공식 clades.tsv 규칙 대비 규칙 기반 판정은 122개 중 120개 일치(보조용) | 완료 | 트리 페이지도 로컬 JSON을 사용하도록 서버 수정 |
| 4 | 10-14 | 파이프라인 | 부분 서열은 버리고, 메타데이터에 명시해 QC FAIL로 표기 | 길이 17000bp 미만(306개)을 `QC FAIL: 길이 미달 (부분 서열)`로 기록. 모호염기(N) 과다 28개(허용 1%), Nextstrain 제외 목록 2개도 사유와 함께 QC FAIL. 전체 442개 중 QC PASS 106개 | 완료 | N 허용 비율(1%)은 프로필 qc.max_n_fraction에서 조정 |
| 5 | 10-14 | 웹앱 | auspice 결과를 바탕으로 메타데이터 옆에 CLADE를 추가 | 메타데이터 표에 'QC 사유'와 'Clade' 열을 QC 열 바로 옆에 배치하고, 주석 달린 전체 메타데이터를 `data/filtered/genome_metadata_annotated.tsv`로 생성 | 완료 |  |
| 6 | 10-14 | 웹앱 | Full Genome Alignment에서 consensus가 고정되어 있는데, 제일 중요한 건 프라이머/프로브이니 그것을 기준으로 움직이게 | 정렬 화면의 기준 행을 '기준: 프라이머/프로브'(F/Probe/R 구간은 올리고 자신의 염기, 그 밖은 consensus)로 바꾸고, 점(.)이 그 기준과의 일치를 뜻하게 함. 선택 상자로 'population consensus' 기준으로도 전환 가능(정확한 염기를 `msa.json`의 `seq`에 저장해 브라우저에서 다시 계산). FASTA 내보내기는 점이 아니라 실제 염기를 기록 | 완료 | 기존 `msa.json`은 `msa_alignment` 규칙을 다시 돌려야 전환 상자가 나타남 |
| 7 | 10-14 | 웹앱 | 설명이나 LASV 등 필요한 곳은 변경해서 업데이트 | 프로필에 `ui` 블록(`name_ko`, `web_url`, `lineage_ko`, `lineage_examples`, `lineage_def`, `orientation_note`) 추가. 가이드·투어·용어집·표준물질/다양성 안내·정적 대시보드·PDF 리포트의 문구가 프로필(바이러스 이름, 단일/다중 segment, 길이·N 기준, 참조 accession, Nextstrain 주소)을 따르고, '리니지'↔'클레이드' 표기는 화면 전체에서 프로필 값으로 치환 | 완료 |  |
| 8 | 10-14 | 웹앱 | PRIMER/PROBE 서열 행도 CONSENSUS처럼 고정(LASV에도 적용), 전체화면(최대화)에서 모달창이 안 뜨는 문제 확인 | 정렬 화면에서 눈금·π·PRIMER/PROBE·기준 행이 모두 상단에 고정됨(LASV 동일 적용). 전체화면 중에는 열린 카드만 보이므로 팝업(상세 창·대화상자·이미지 확대)을 <body> 대신 전체화면 요소 안에 붙이도록 수정 | 완료 | 앱 내장 브라우저에서는 실제 전체화면 전환이 되지 않아, 전체화면 요소를 가정한 시뮬레이션으로 확인함 |
| 9 | 10-14 | 문서 | Nipah 파이프라인이 달라졌으니 워크플로 설명 업데이트(메타데이터에 클레이드 추가 등) | 가이드 페이지의 흐름도에 로컬 Nextstrain 빌드 → lineage_table → rename_headers, filter → annotate_metadata 가지를 추가하고, 단계 1b(클레이드 부여)와 단계 2의 QC FAIL/주석 메타데이터 설명, PDF 리포트 워크플로 페이지, README의 Nipah 흐름도를 갱신. 흐름도 상자의 &rarr; 같은 문자가 그대로 보이던 문제와, '리니지→클레이드' 치환이 lineage_table 같은 식별자까지 바꾸던 문제도 수정 | 완료 |  |
| 10 | 10-14 | 설치 | 처음 구동을 위한 conda/mamba 세팅용 setup.sh. mfeprimer는 실행 파일 경로를 넣게 하고 파라미터 설정에도 추가 | `setup.sh`(환경 `primer_validation` 생성·도구 확인·mfeprimer 경로 입력/검증, `config.yaml`의 `env_prefix`/`mfeprimer_bin` 기록)와 `environment_nextstrain.yml`(클레이드 재빌드용, `--nextstrain`) 추가. 'Analysis Parameters' 페이지에 mfeprimer 경로 입력·확인 버튼과 사용 중인 환경/도구 경로 표시, 저장 시 실제 실행으로 검증. `start.sh`와 웹앱이 기록된 환경을 PATH에 올림. LASV에도 동일 적용 | 완료 | 새 환경을 실제로 처음부터 만들어 보지는 않음(해결 가능 여부는 dry-run으로, 스크립트 동작은 기존 환경을 지정해 확인) |
| 11 | 10-14 | 운영 | LASV와 Nipah 분석을 동시에 돌려도 엉키거나 충돌하지 않는지 확인 | 작업 파일·`.snakemake/` 잠금·mfeprimer 인덱스/임시 파일이 모두 프로젝트 폴더 안이고 /tmp 임시 이름은 모두 고유(pid·TemporaryDirectory)하며 포트도 분리(8765/8767)되어 있어 충돌 없음을 확인. 같은 폴더에서 snakemake 두 개는 Snakemake 잠금으로 막힘(웹앱도 작업 1개씩만 허용). 공유되는 것은 읽기 전용(mfeprimer 실행 파일, conda 환경)뿐이며, 동시 실행 시 CPU만 나눠 쓰게 됨 | 완료 | 동시 실행 시 mfeprimer_cpu(64) × snakemake --cores 16이 코어 수를 크게 넘으므로 값을 낮추면 좋음 |
| 12 | 10-14 | 점검 | 전체 프로세스 점검 | 독립 무차별 재계산(`scripts/audit_recount.py`)으로 판정 근거인 mismatch 표를 검증해, F/R 프라이머의 mismatch가 3' 끝 k-mer 밖에서는 0으로 과소 집계되던 결함을 찾아 수정(게놈 염기 기준 재계산). 프로브가 증폭 산물 안에 있는지도 판정에 포함. Nipah 결과 재계산 및 감사 일치 확인 | 완료 | LASV는 그대로 두었고 같은 결함이 있음(감사로 확인) — 결과 재계산 여부는 사용자가 결정 |

## 남은 확인 사항

| 항목 | 내용 |
|---|---|
| 트리에 없는 PASS 게놈 2개(FV537359, FV537360)의 클레이드 | 현재 unclassified. 필요하면 clades.tsv 규칙 기반 보조 판정 추가 |
| 표준물질·프라이머 세트 입력 | 서열을 받는 대로 입력 페이지에서 등록하고 전체 분석 실행 |
| LASV 파이프라인과의 통합 | 병합 전 LASV 회귀 비교 필요 (새 파라미터로 일부 규칙이 재실행 대상이 됨) |
| 가이드의 LASV 수치 예시 | adjustdirection 설명의 55%/20% 같은 LASV 실측 수치는 `ui.orientation_note`가 있는 프로필에만 표시됨. Nipah는 일반 문구 |
