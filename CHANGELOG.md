# Changelog

All notable changes to this project are documented here. Format follows
[Keep a Changelog](https://keepachangelog.com/), versioning is `MAJOR.MINOR` starting at 0.1.

## [Unreleased]

## [0.1.0] - 2026-11-10

### Added
- **BDBV 에디션 첫 버전**: NiV 에디션(단일 게놈)을 바탕으로 만든 Orthoebolavirus bundibugyoense 검증 대시보드. 게놈과 메타데이터는 **Pathoplexus**(LAPIS, `ebola-bdbv`, 공개 데이터)에서 받아 오며 (`scripts/fetch_pathoplexus.py`, 프로파일 `source.type = "pathoplexus"`), 에볼라에는 클레이드 체계가 없어 **채집 국가별**로 묶습니다(리니지 열 = `country_lineage`). 참조 게놈은 RefSeq .
- 프라이머·표준물질 목록은 비어 있는 상태로 시작합니다(`config/primers.tsv`, `config/standards.tsv`).
