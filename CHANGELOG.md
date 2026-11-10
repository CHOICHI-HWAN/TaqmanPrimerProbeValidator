# Changelog

All notable changes to this project are documented here. Format follows
[Keep a Changelog](https://keepachangelog.com/), versioning is `MAJOR.MINOR` starting at 0.1.

## [Unreleased]

## [0.1.1] - 2026-11-10

### Changed
- **모바일(폰) 화면**: 사이드바가 화면 위쪽 절반을 차지하던 것을 **전체 화면 서랍**으로 바꿨습니다. 처음에는 닫혀 있고(왼쪽 위 ☰ 메뉴), 항목을 고르면 자동으로 닫힙니다(폰에서는 열림/닫힘을 기억하지 않음). 뒤로/앞으로 막대는 폰 여백에 맞추고 스크롤해도 위에 붙어 있으며, 가로로 넘치던 표·긴 경로·두 칼럼 배치(문헌·설계 화면)는 칸 안에서 스크롤되거나 한 칸으로 쌓입니다. 데스크톱 화면은 그대로입니다.

## [0.1.0] - 2026-11-10

### Added
- **EBOV 에디션 첫 버전**: NiV 에디션(단일 게놈)을 바탕으로 만든 Orthoebolavirus zairense 검증 대시보드. 게놈과 메타데이터는 **Pathoplexus**(LAPIS, `ebola-zaire`, 공개 데이터)에서 받아 오며 (`scripts/fetch_pathoplexus.py`, 프로파일 `source.type = "pathoplexus"`), 에볼라에는 클레이드 체계가 없어 **채집 국가별**로 묶습니다(리니지 열 = `country_lineage`). 참조 게놈은 RefSeq .
- 프라이머·표준물질 목록은 비어 있는 상태로 시작합니다(`config/primers.tsv`, `config/standards.tsv`).
