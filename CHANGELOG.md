# Changelog

All notable changes to this project are documented here. Format follows
[Keep a Changelog](https://keepachangelog.com/), versioning is `MAJOR.MINOR` starting at 0.1.

## [Unreleased]

## [0.1.3] - 2026-11-10

### Added
- **🦠 감염병 정보** 메뉴: 이 바이러스가 일으키는 감염병의 간단한 정보(병원체·게놈·숙주·분포·잠복기·치명률·백신/치료, 전파·증상)와 WHO 참고 링크. **⏱ 감염 경과 · 바이러스량** 메뉴: 노출 후 발병 시기, 바이러스(RNA/DNA)량이 높아지는 시기(모식도), IgM/IgG 출현 시기, 검체별 검출 시기. 내용은 프로파일 JSON(`disease`, `course`)에 있어 고칠 수 있습니다. 일반 공개 자료의 대략적인 요약이며 임상·방역 판단용이 아닙니다.

### Added
- Pathoplexus 출처에 **세부 계통 열**(`source.sublineage_field`, 예: mpox `lineage`)을 `sub_lineage`로 함께 받습니다. 그룹(리니지) 열은 그대로이고, 메타데이터 표의 Sub-lineage 열과 게놈 카드에 표시만 됩니다. 계통 값이 `unassigned`·`unknown`·`none`이면 `unclassified`로 통일합니다.

## [0.1.2] - 2026-11-10

### Added
- **모든 분석용 게놈에 데이터 출처 표기**: 개요 화면 제목 아래에 "데이터 출처: Pathoplexus (organism …) · 받은 날짜 · 기록별 출처(GenBank/RefSeq 등)" 줄이 생기고(`/api/profile`의 `data_source`), 메타데이터 표에 **데이터 출처** 열이, 게놈 메타데이터 카드(식별)에 출처 항목이 생겼습니다. 기록별 출처는 Nextstrain·Pathoplexus·NCBI Datasets 중 어디서 받았는지와 원본 DB(GenBank / RefSeq / 자체 등록)입니다. (`webapp/data_source.py`, 테스트 `tests/test_data_source.py`)
- 출처 줄에는 **받은 날**과 원본의 **최신 갱신일·최신 등록일**(기록의 `date_updated`·`date_released` 중 가장 최근)도 함께 표시합니다.
- Pathoplexus에서 기록이 아주 많은 종(예: mpox 약 1.5만 건·3 GB)은 한 번에 받으면 멈추므로 **1000건씩 정렬된 페이지로** 받습니다(`scripts/fetch_pathoplexus.py`). 작은 종은 이전과 같은 한 번 요청입니다.
- 방향 판별(`scripts/orientation.py`)은 6만 염기가 넘는 게놈(두창류 약 197 kb)에서는 **가운데 2만 염기**만 비교합니다. 전체를 비교하면 1.5만 건에 몇 시간이 걸립니다. 그 미만의 게놈(기존 모든 바이러스)은 결과가 그대로입니다.

## [0.1.1] - 2026-11-10

### Changed
- **모바일(폰) 화면**: 사이드바가 화면 위쪽 절반을 차지하던 것을 **전체 화면 서랍**으로 바꿨습니다. 처음에는 닫혀 있고(왼쪽 위 ☰ 메뉴), 항목을 고르면 자동으로 닫힙니다(폰에서는 열림/닫힘을 기억하지 않음). 뒤로/앞으로 막대는 폰 여백에 맞추고 스크롤해도 위에 붙어 있으며, 가로로 넘치던 표·긴 경로·두 칼럼 배치(문헌·설계 화면)는 칸 안에서 스크롤되거나 한 칸으로 쌓입니다. 데스크톱 화면은 그대로입니다.

## [0.1.0] - 2026-11-10

### Added
- **SUDV 에디션 첫 버전**: NiV 에디션(단일 게놈)을 바탕으로 만든 Orthoebolavirus sudanense 검증 대시보드. 게놈과 메타데이터는 **Pathoplexus**(LAPIS, `ebola-sudan`, 공개 데이터)에서 받아 오며 (`scripts/fetch_pathoplexus.py`, 프로파일 `source.type = "pathoplexus"`), 에볼라에는 클레이드 체계가 없어 **채집 국가별**로 묶습니다(리니지 열 = `country_lineage`). 참조 게놈은 RefSeq .
- 프라이머·표준물질 목록은 비어 있는 상태로 시작합니다(`config/primers.tsv`, `config/standards.tsv`).
