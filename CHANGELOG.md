# Changelog

All notable changes to this project are documented here. Format follows
[Keep a Changelog](https://keepachangelog.com/), versioning is `MAJOR.MINOR` starting at 0.1.

## [Unreleased]

## [0.1.4] - 2026-11-13

### Added
- **커버리지 우선 자동 설계**(`scripts/design_sweep*.py`, `scripts/optimize_user_sets.py`): 여러 격자(wide/hightm/midtm/longoligo 등)로 전체 유전체 대비 프라이머·프로브 후보를 탐색하고, 기존 세트를 부분적으로 개선한 안도 별도 이름(`-OPT-* ★`)으로 등록합니다. 상위 3개 세트에는 🥇🥈🥉 메달을 표시합니다.
- **전체 데이터 백업 · 복원**: 대시보드 우측 상단 아이콘(💾)을 누르면 모달 창에서 이 에디션의 설정·입력 프라이머·프로필·data·output을 하나의 압축 파일(xz, 또는 더 빠른 zstd)로 백업하고 되돌릴 수 있습니다. 진행바가 표시되고, 복원은 기존 폴더를 지우지 않고 `.pre_restore_*/`로 옮겨 둡니다. 코드(스크립트/화면)는 포함하지 않으며, 서열이 들어 있어 저장 위치는 git 밖입니다.
- **세트 상세 페이지에 커버리지 곡선**: 허용하는 mismatch 수가 늘 때 통과하는 유전체 비율을 보여주는 곡선(🧭 프라이머 매핑지도에서 여러 세트를 겹쳐 보던 것과 같은 방식)을 세트 하나만 볼 때도 바로 표시합니다.

### Changed
- **전체 대시보드·사이드바**: 자동 설계(AUTO-DESIGN)로 만들어진 세트를 사용자가 입력한 세트 아래, 별도 "⚙ 자동 설계" 구획으로 분리해서 보여줍니다(이전에는 알파벳 순서상 위로 섞여 보이기도 했습니다).
- **다른 에디션 전환**: 사이드바의 LASV/NiV/... 전환 버튼이 폴더 구조가 다른 에디션(EBOLA/*, SAHF/*)에서도 항상 보이도록 고쳤습니다.
- **같은 창에서 이동**: 전체 목록 페이지로 가거나 다른 에디션 대시보드를 열 때 새 창이 아니라 같은 창에서 열립니다.
- `aggregate_groups`가 세트 그룹이 0개일 때도 실패하지 않도록 수정(기본 헤더에 `host_scope` 포함).

### Removed
- 거의 쓰이지 않는 "+ 새 프라이머 세트 분석" 버튼을 없앴습니다.

### Changed
- **기본 스캔 엔진을 Go로** (`scan_engine: "go"`): 전체 유전체 스캔이 Go 엔진으로 돕니다. 같은 입력을 Python·Go로 각각 전체 분석해 비교한 결과 세트·그룹 판정 표가 모두 같았습니다(NiV, LASV, EBOV, SUDV, BDBV, RESTV, TAFV, MPXV 축소본). 스캔을 쓰는 `alignment_table` 단계는 약 9~60배 빨라졌습니다. Go가 없으면 Python 스캔으로 자동 전환됩니다.

### Fixed
- NCBI 다중 분류군 수신: 아형 부분 필터(`sub=H3`/`N2`)의 오류를 고치고(`subtype_has`, 테스트 추가), 분절의 **참조 계정을 항상 데이터에 포함**합니다(`ref=`; 연도·샘플 조건에 걸려 참조가 빠지면 참조 게놈을 잘라 낼 수 없었음).

### Added
- **한 생물종에 참조 게놈이 여러 개일 때**: 프로파일 분절(`segments.<단위>`)에 `lineage_in: [값…]`을 쓰면 전체를 한 번 받아(`data/raw/_all_*`) 라이니지 열 값으로 분석 단위를 나눕니다(규칙 `fetch_all`·`split_raw`, `scripts/split_raw_by_lineage.py`). 각 단위는 자기 `reference_accession`을 씁니다. 엠폭스 clade I/Ia·Ib·II 참조 3개에 쓰입니다.
- NCBI 출처에서 **여러 분류군·필터**(`taxon` 목록, `since`, `host_filter`, `group_by`, 분절·아형별 요청, 아형 이름이 없는 기록을 HA의 BLAST로 판정하는 `classify_subtype`, 연도·아형당 균주 상한 `cap_per_group`) 지원: 인플루엔자 A형용. 수신 규칙은 URL 인자를 따옴표로 감쌉니다.

### Added
- NCBI Datasets 출처(`source.type = "ncbi"`)가 **분절 바이러스**를 지원합니다(`segment_subdirs: true`, 분절별 `taxon/SEG` 요청). 기록에 `segment`가 있으면 그것으로, 없으면 길이가 한 분절 전장에 가까울 때 길이로 분절을 정합니다(리프트밸리열 바이러스, 분절별 `SEG_LEN`).

## [0.1.3] - 2026-11-10

### Added
- **🦠 감염병 정보** 메뉴: 이 바이러스가 일으키는 감염병의 간단한 정보(병원체·게놈·숙주·분포·잠복기·치명률·백신/치료, 전파·증상)와 WHO 참고 링크. **⏱ 감염 경과 · 바이러스량** 메뉴: 노출 후 발병 시기, 바이러스(RNA/DNA)량이 높아지는 시기(모식도), IgM/IgG 출현 시기, 검체별 검출 시기. 내용은 프로파일 JSON(`disease`, `course`)에 있어 고칠 수 있습니다. 일반 공개 자료의 대략적인 요약이며 임상·방역 판단용이 아닙니다.

### Added
- Pathoplexus 출처에 **세부 계통 열**(`source.sublineage_field`, 예: mpox `lineage`)을 `sub_lineage`로 함께 받습니다. 그룹(리니지) 열은 그대로이고, 메타데이터 표의 Sub-lineage 열과 게놈 카드에 표시만 됩니다. 계통 값이 `unassigned`·`unknown`·`none`이면 `unclassified`로 통일합니다.

## [0.1.2] - 2026-11-10

### Added
- **모든 분석용 게놈에 데이터 출처 표기**: 개요 화면 제목 아래에 "데이터 출처: Pathoplexus (organism …) · 받은 날짜 · 기록별 출처(GenBank/RefSeq 등)" 줄이 생기고(`/api/profile`의 `data_source`), 메타데이터 표에 **데이터 출처** 열이, 게놈 메타데이터 카드(식별)에 출처 항목이 생겼습니다. 기록별 출처는 Nextstrain·Pathoplexus·NCBI Datasets 중 어디서 받았는지와 원본 DB(GenBank / RefSeq / 자체 등록)입니다. (`webapp/data_source.py`, 테스트 `tests/test_data_source.py`)
- 출처 줄에는 **받은 날**과 원본의 **최신 갱신일·최신 등록일**(기록의 `date_updated`·`date_released` 중 가장 최근)도 함께 표시합니다.
- 방향 판별(`scripts/orientation.py`)은 6만 염기가 넘는 게놈(두창류 약 197 kb)에서는 **가운데 2만 염기**만 비교합니다. 전체를 비교하면 1.5만 건에 몇 시간이 걸립니다. 그 미만의 게놈(기존 모든 바이러스)은 결과가 그대로입니다.

## [0.1.1] - 2026-11-10

### Changed
- **모바일(폰) 화면**: 사이드바가 화면 위쪽 절반을 차지하던 것을 **전체 화면 서랍**으로 바꿨습니다. 처음에는 닫혀 있고(왼쪽 위 ☰ 메뉴), 항목을 고르면 자동으로 닫힙니다(폰에서는 열림/닫힘을 기억하지 않음). 뒤로/앞으로 막대는 폰 여백에 맞추고 스크롤해도 위에 붙어 있으며, 가로로 넘치던 표·긴 경로·두 칼럼 배치(문헌·설계 화면)는 칸 안에서 스크롤되거나 한 칸으로 쌓입니다. 데스크톱 화면은 그대로입니다.

## [0.1.0] - 2026-11-10

### Added
- **TAFV 에디션 첫 버전**: NiV 에디션(단일 게놈)을 바탕으로 만든 Orthoebolavirus taiense 검증 대시보드. 게놈과 메타데이터는 **NCBI Datasets API**(`NCBI taxon 186541`, 공개 데이터; Pathoplexus에는 이 종이 없음)에서 받아 오며 (`scripts/fetch_ncbi_virus.py`, 프로파일 `source.type = "ncbi"`), 에볼라에는 클레이드 체계가 없어 **채집 국가별**로 묶습니다(리니지 열 = `country_lineage`). 참조 게놈은 RefSeq .
- 프라이머·표준물질 목록은 비어 있는 상태로 시작합니다(`config/primers.tsv`, `config/standards.tsv`).
