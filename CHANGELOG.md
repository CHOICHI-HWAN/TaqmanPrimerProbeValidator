# Changelog

All notable changes to this project are documented here. Format follows
[Keep a Changelog](https://keepachangelog.com/), versioning is `MAJOR.MINOR` starting at 0.1.

## [Unreleased]

### Changed
- **기본 스캔 엔진을 Go로** (`scan_engine: "go"`): 전체 유전체 스캔이 Go 엔진으로 돕니다. 같은 입력을 Python·Go로 각각 전체 분석해 비교한 결과 세트·그룹 판정 표가 모두 같았습니다(NiV, LASV, EBOV, SUDV, BDBV, RESTV, TAFV, MPXV 축소본). 스캔을 쓰는 `alignment_table` 단계는 약 9~60배 빨라졌습니다. Go가 없으면 Python 스캔으로 자동 전환됩니다.

## [0.1.0] - 2026-11-10

### Added
- **FLUA 에디션 첫 버전**: NiV 에디션(단일 게놈)을 바탕으로 만든 Alphainfluenzavirus influenzae 검증 대시보드. 게놈과 메타데이터는 **NCBI Datasets API**(`NCBI taxon 11320`, 공개 데이터; Pathoplexus에는 이 종이 없음)에서 받아 오며 (`scripts/fetch_ncbi_virus.py`, 프로파일 `source.type = "ncbi"`), 에볼라에는 클레이드 체계가 없어 **채집 국가별**로 묶습니다(리니지 열 = `country_lineage`). 참조 게놈은 RefSeq .
- 프라이머·표준물질 목록은 비어 있는 상태로 시작합니다(`config/primers.tsv`, `config/standards.tsv`).
