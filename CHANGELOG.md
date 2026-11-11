# Changelog

All notable changes to this project are documented here. Format follows
[Keep a Changelog](https://keepachangelog.com/), versioning is `MAJOR.MINOR` starting at 0.1.

## [Unreleased]

### Added
- **한 생물종에 참조 게놈이 여러 개일 때**: 프로파일 분절(`segments.<단위>`)에 `lineage_in: [값…]`을 쓰면 전체를 한 번 받아(`data/raw/_all_*`) 라이니지 열 값으로 분석 단위를 나눕니다(규칙 `fetch_all`·`split_raw`, `scripts/split_raw_by_lineage.py`). 각 단위는 자기 `reference_accession`을 씁니다. 엠폭스 clade I/Ia·Ib·II 참조 3개에 쓰입니다.
- NCBI 출처에서 **여러 분류군·필터**(`taxon` 목록, `since`, `host_filter`, `group_by`, 분절·아형별 요청, 아형 이름이 없는 기록을 HA의 BLAST로 판정하는 `classify_subtype`, 연도·아형당 균주 상한 `cap_per_group`) 지원: 인플루엔자 A형용. 수신 규칙은 URL 인자를 따옴표로 감쌉니다.

## [0.1.0] - 2026-11-10

### Added
- **RVFV 에디션 첫 버전**: NiV 에디션(단일 게놈)을 바탕으로 만든 Phlebovirus riftense 검증 대시보드. 게놈과 메타데이터는 **NCBI Datasets API**(`NCBI taxon 11588`, 공개 데이터; Pathoplexus에는 이 종이 없음)에서 받아 오며 (`scripts/fetch_ncbi_virus.py`, 프로파일 `source.type = "ncbi"`), 에볼라에는 클레이드 체계가 없어 **채집 국가별**로 묶습니다(리니지 열 = `country_lineage`). 참조 게놈은 RefSeq .
- 프라이머·표준물질 목록은 비어 있는 상태로 시작합니다(`config/primers.tsv`, `config/standards.tsv`).
