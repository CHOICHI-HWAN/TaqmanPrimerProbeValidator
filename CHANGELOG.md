# Changelog

All notable changes to this project are documented here. Format follows
[Keep a Changelog](https://keepachangelog.com/), versioning is `MAJOR.MINOR` starting at 0.1.

## [Unreleased]

## [0.1.0] - 2026-11-10

### Added
- **MPXV 에디션 첫 버전**: NiV 에디션(단일 게놈)을 바탕으로 만든 Orthopoxvirus monkeypox 검증 대시보드. 게놈과 메타데이터는 **Pathoplexus**(LAPIS, `mpox`, 공개 데이터)에서 받아 오며 (`scripts/fetch_pathoplexus.py`, 프로파일 `source.type = "pathoplexus"`), 에볼라에는 클레이드 체계가 없어 **채집 국가별**로 묶습니다(리니지 열 = `country_lineage`). 참조 게놈은 RefSeq .
- 프라이머·표준물질 목록은 비어 있는 상태로 시작합니다(`config/primers.tsv`, `config/standards.tsv`).
- 방향 판별(`scripts/orientation.py`)은 6만 염기가 넘는 게놈(두창류 약 197 kb)에서는 **가운데 2만 염기**만 비교합니다. 전체를 비교하면 1.5만 건에 몇 시간이 걸립니다. 그 미만의 게놈(기존 모든 바이러스)은 결과가 그대로입니다.
- 데이터 출처(이름·받은 날·원본 최신 갱신일·기록별 GenBank/RefSeq)가 개요·메타데이터 표·게놈 카드에 표시됩니다. Pathoplexus mpox(약 1.5만 건)는 1000건씩 페이지로 받습니다. 클레이드(I, Ia, Ib, IIa, IIb)는 Pathoplexus `clade` 열을 씁니다.
- 프라이머 검증은 **subclade(Ia, Ib, IIa, IIb)** 단위로 봅니다. 더 세부적인 **lineage**(IIb의 B.1·C.1…, Ib의 A·A.1…)는 Pathoplexus `lineage` 열을 `sub_lineage`로 받아 메타데이터 표(Sub-lineage 열)와 게놈 카드에 표시만 합니다.
