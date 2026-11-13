# Changelog

All notable changes to this project are documented here. Format follows
[Keep a Changelog](https://keepachangelog.com/), versioning is `MAJOR.MINOR` starting at 0.1.

## [Unreleased]

## [0.1.1] - 2026-11-13

### Added
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

## [0.1.0] - 2026-11-10

### Added
- **FLUA 에디션 첫 버전**: NiV 에디션(단일 게놈)을 바탕으로 만든 Alphainfluenzavirus influenzae 검증 대시보드. 게놈과 메타데이터는 **NCBI Datasets API**(`NCBI taxon 11320`, 공개 데이터; Pathoplexus에는 이 종이 없음)에서 받아 오며 (`scripts/fetch_ncbi_virus.py`, 프로파일 `source.type = "ncbi"`), 에볼라에는 클레이드 체계가 없어 **채집 국가별**로 묶습니다(리니지 열 = `country_lineage`). 참조 게놈은 RefSeq .
- 프라이머·표준물질 목록은 비어 있는 상태로 시작합니다(`config/primers.tsv`, `config/standards.tsv`).
