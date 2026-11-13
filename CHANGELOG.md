# Changelog

All notable changes to this project are documented here. Format follows
[Keep a Changelog](https://keepachangelog.com/), versioning is `MAJOR.MINOR` starting at 0.1.

## [Unreleased]

## [0.2.0] - 2026-11-13

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

### Changed (MPXV 전용)
- **분석 단위를 하나로 통합**: Clade I/Ia(NC_003310), Ib(PP601219), II(NC_063383)로 나뉘어 있던 세 분석 단위를 하나(`genome`)로 합쳐, 같은 세트가 세 번 복제되던 문제(18개 → 6개)를 없앴습니다.
- **참조유전체 여러 개**: 합치면서 빠진 Clade I·Ib 참조유전체를 세트 상세의 "참조유전체 결합위치"에 주 참조(Clade II)와 나란히 보여줍니다(프로파일 `extra_references`).
- **다양성(pixy π) 계산 끄기 옵션**: 유전체가 커서(약 197kb × 게놈 7천여 개) 전체 MSA가 12시간 넘게 걸리던 문제로, 프로파일에서 `diversity.enabled: false`로 끌 수 있게 했습니다.

### Changed
- **기본 스캔 엔진을 Go로** (`scan_engine: "go"`): 전체 유전체 스캔이 Go 엔진으로 돕니다. 같은 입력을 Python·Go로 각각 전체 분석해 비교한 결과 세트·그룹 판정 표가 모두 같았습니다(NiV, LASV, EBOV, SUDV, BDBV, RESTV, TAFV, MPXV 축소본). 스캔을 쓰는 `alignment_table` 단계는 약 9~60배 빨라졌습니다. Go가 없으면 Python 스캔으로 자동 전환됩니다.

### Fixed
- NCBI 다중 분류군 수신: 아형 부분 필터(`sub=H3`/`N2`)의 오류를 고치고(`subtype_has`, 테스트 추가), 분절의 **참조 계정을 항상 데이터에 포함**합니다(`ref=`; 연도·샘플 조건에 걸려 참조가 빠지면 참조 게놈을 잘라 낼 수 없었음).

## [0.1.2] - 2026-11-11

### Changed
- **참조 게놈 3개, 분석 단위 3개**: 엠폭스를 clade 묶음별 세 분석 단위로 나눴습니다 — **Clade I/Ia**(참조 NC_003310), **Clade Ib**(참조 PP601219; C3L(OPG032) 결실), **Clade II**(IIa·IIb; 참조 NC_063383; C3L 없음). 전체를 한 번 받아(`data/raw/_all_*`) 라이니지 열로 나눕니다(`scripts/split_raw_by_lineage.py`, 프로파일 `segments.<단위>.lineage_in`). Ib/IIb 혼재는 Ib로 읽고 클레이드 미지정 기록은 뺍니다. 기존 프라이머 세트 6개는 Clade II에 그대로 두고 Clade I/Ia·Ib에 복사했습니다(P007~P018).

### Added
- NCBI Datasets 출처(`source.type = "ncbi"`)가 **분절 바이러스**를 지원합니다(`segment_subdirs: true`, 분절별 `taxon/SEG` 요청). 기록에 `segment`가 있으면 그것으로, 없으면 길이가 한 분절 전장에 가까울 때 길이로 분절을 정합니다(리프트밸리열 바이러스, 분절별 `SEG_LEN`).

## [0.1.1] - 2026-11-10

### Added
- **🦠 감염병 정보** 메뉴: 이 바이러스가 일으키는 감염병의 간단한 정보(병원체·게놈·숙주·분포·잠복기·치명률·백신/치료, 전파·증상)와 WHO 참고 링크. **⏱ 감염 경과 · 바이러스량** 메뉴: 노출 후 발병 시기, 바이러스(RNA/DNA)량이 높아지는 시기(모식도), IgM/IgG 출현 시기, 검체별 검출 시기. 내용은 프로파일 JSON(`disease`, `course`)에 있어 고칠 수 있습니다. 일반 공개 자료의 대략적인 요약이며 임상·방역 판단용이 아닙니다.

## [0.1.0] - 2026-11-10

### Added
- **MPXV 에디션 첫 버전**: NiV 에디션(단일 게놈)을 바탕으로 만든 Orthopoxvirus monkeypox 검증 대시보드. 게놈과 메타데이터는 **Pathoplexus**(LAPIS, `mpox`, 공개 데이터)에서 받아 오며 (`scripts/fetch_pathoplexus.py`, 프로파일 `source.type = "pathoplexus"`), 에볼라에는 클레이드 체계가 없어 **채집 국가별**로 묶습니다(리니지 열 = `country_lineage`). 참조 게놈은 RefSeq .
- 프라이머·표준물질 목록은 비어 있는 상태로 시작합니다(`config/primers.tsv`, `config/standards.tsv`).
- 방향 판별(`scripts/orientation.py`)은 6만 염기가 넘는 게놈(두창류 약 197 kb)에서는 **가운데 2만 염기**만 비교합니다. 전체를 비교하면 1.5만 건에 몇 시간이 걸립니다. 그 미만의 게놈(기존 모든 바이러스)은 결과가 그대로입니다.
- 데이터 출처(이름·받은 날·원본 최신 갱신일·기록별 GenBank/RefSeq)가 개요·메타데이터 표·게놈 카드에 표시됩니다. Pathoplexus mpox(약 1.5만 건)는 1000건씩 페이지로 받습니다. 클레이드(I, Ia, Ib, IIa, IIb)는 Pathoplexus `clade` 열을 씁니다.
- 프라이머 검증은 **subclade(Ia, Ib, IIa, IIb)** 단위로 봅니다. 더 세부적인 **lineage**(IIb의 B.1·C.1…, Ib의 A·A.1…)는 Pathoplexus `lineage` 열을 `sub_lineage`로 받아 메타데이터 표(Sub-lineage 열)와 게놈 카드에 표시만 합니다.
