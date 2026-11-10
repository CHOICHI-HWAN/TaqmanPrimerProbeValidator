# GenomeCurator(GC) → 프라이머 검증 파이프라인: 필요한 필드와 값 규칙

작성일 2026-10-24. 목적: GC를 이 파이프라인의 **지속적인 단일 데이터 원천**으로 만들기 위해, 파이프라인이 실제로 읽는 필드와 그 값이 지켜야 할 규칙을 한곳에 기록한다.
GC 개발의 기준 문서로 쓴다. 파이프라인 쪽에서 값을 고쳐 쓰는 어댑터는 최소로 두고, 규칙은 **원천(GC)에서** 지킨다.

근거: 파이프라인 코드(`scripts/rename_headers.py`, `filter.py`, `host_scope.py`, `webapp/server.py`)와 GC 코드(`genomecurator/export.py`, `vocab.py`, `datasets.py`, `qc.py`)를 읽고 확인한 것이다.
"GC 현재"의 ✔ = 코드에서 확인, △ = 있으나 형식·동작이 달라 보완 필요, ✘ = 없음, ? = 미확인.

우선순위: **P0** 없으면 연동 불가이거나 결과가 조용히 틀림 · **P1** 재현성·정확도에 필요 · **P2** 있으면 좋음.

## 1. 유전체(accession)별 필드

한 행 = 한 서열 레코드. 분절 바이러스(Lassa)는 분절별로 파일을 나누어 받는다(§3).

| 필드(우리 이름) | 우선 | 용도 | 형식·허용 값 | 어긋나면 | GC 현재 |
| --- | --- | --- | --- | --- | --- |
| `accession` | **P0** | 서열과 메타데이터를 잇는 유일한 키 | **버전 없는** INSDC/RefSeq 계정번호(`KM822128`). 분절 파일 안에서 유일. FASTA 헤더의 첫 토큰과 정확히 같음 | 메타데이터 행이 없는 서열은 **조용히 버려짐**(`dropped_no_metadata`). 헤더가 `KM822128.1`이면 전부 탈락 | △ `accession`·`accession_version` 별도 열 ✔, 그러나 `plain` 내보내기의 FASTA 헤더와 id 열은 버전 포함 |
| `segment` | **P0**(분절 바이러스) | 분절별 파일 분리(`s`, `l`) | 소문자 키 `s` `l`. 단일 게놈 바이러스는 `genome` 또는 빈 값 | 분절이 섞이면 S와 L이 한 모집단으로 분석됨 | △ `segment` 열 ✔, 값 형식 ?, 데이터셋 규칙에 분절 조건 ✘ |
| `lineage` (= clade) | **P0** | 리니지별 커버리지, Venn·곡선, 재설계, 다양성 | 닫힌 값 목록. Lassa `LI` `LII` `LIII` `LIV_V` `LVI` `LVII`, NiV `I` `Ia` `Ib` `II`. 알 수 없으면 **`unclassified`**(빈 값 금지). 분절별로 다를 수 있음(Lassa는 S·L 각각) | 리니지 분석이 `NA`로 뭉침 | ✘ Lassa·Nipah에 부여 로직 없음(Nextclade 데이터셋·clade 사전이 mpox 등만). `clade`·`lineage` 열은 ✔ |
| `is_lab_host` | **P0** | 실험 숙주(세포주) 서열을 모든 범위에서 제외 | `true` 이면 lab. 빈 값은 아님. (받아들이는 값: `1 true t yes y`, 대소문자 무시) | 세포주 서열이 **animal 범위에 섞임** | △ 사전의 `flag = cell line`으로 판정 가능하나 내보내기에 열이 없음 |
| `host` | **P0** | 사람/동물/미상 범위 선택(기본 사람) | 표준 숙주명 문자열(§2) | 빈 값은 `unknown` 범위로 감 | ✔ `host_name`(사전 표준화) → `host`로 이름 변경 필요 |
| `date` | **P0** | 연도별 차트·표 | `YYYY-MM-DD`, 모르는 부분은 `XX`(`2012-XX-XX`, `XXXX-XX-XX`). 연도는 앞 4자리가 숫자일 때만 집계 | 빈 값은 `XXXX-XX-XX`로 처리됨. `2022-07` 같은 표기는 달·일 구분 불가 | △ `nextstrain` 프로파일만 이 형식(`std_date`). `plain`은 원형 `collection_date` |
| `country` | **P0** | 국가별 지도·집계 | 표준 국가명, 미상은 **빈 값**(`Missing` 같은 임의 표기 금지) | 지도에서 빠짐, 집계 키가 갈라짐 | △ 사전 표준화 ✔(`isolation_country`). 우리 지도 이름과의 일치 ? |
| `region` | **P0** | 대륙 단위 집계 | `Africa` `Asia` `Europe` `North America` `South America` `Oceania`. 미상은 빈 값(우리가 `Unknown/Missing`으로 바꿈) | 대륙 집계가 비어 있음 | △ `region_of(country)` ✔이지만 `nextstrain` 프로파일에서만 계산 |
| `strain` 정규화 값 / `virus_id` | **P0**(분절 바이러스) | 같은 바이러스의 S·L 레코드를 잇는 키(파이프라인 결과를 바이러스 단위로 묶어 보는 통합 벤 다이어그램에 사용) | 대소문자·공백을 정리한 strain 이름, 한 바이러스에 레코드가 여럿이면 같은 `virus_id` 공유. 연결 근거와 신뢰도(`exact` / `manual` / `uncertain`) 함께 | 이름만으로 이으면 충돌(`AV`, `CSF` 등 짧은 이름)과 표기 차이(`LASV/H.sapiens-tc/NGA/2013/IRR_005` ↔ `IRR_005`)를 구분 못 함 | △ `strain` ✔, 연결·신뢰도 ✘ |
| `host_taxid` | P1 | 숙주 분류의 우선 근거(`9606`이면 사람) | 정수 NCBI taxid | 자유 텍스트 규칙으로 분류(오분류 여지) | △ 사전에 `taxon_id` 열 ✔, 내보내기에 ✘ |
| `strain` | P1 | 표시 이름 | 문자열 | 표에서 accession으로 대체 | ✔ |
| `accession_version` | P1 | 출처 추적 | `KM822128.1` | — | ✔ |
| `duplicate_of` | P1 | 같은 유전체의 RefSeq/GenBank 중복을 하나만 | 대표 accession 또는 빈 값 | **분모가 부풀어** 커버리지 왜곡. 우리 중복 제거는 같은 계정번호 문자열만 거름 | ? `genbank_accessions`·`refseq_accessions` 열 ✔, 처리 방식 ? |
| `excluded` + `reason` | P1 | 제외 목록(`exclude.txt`) 대체, 이유 기록 | accession과 사유 | 의심 서열이 분석에 남음 | ✔ 제외 기능 있음(사유·되돌리기), 내보내기에 `excluded.tsv` 기록 |
| `length` | P2 | 표시용(우리는 서열에서 다시 잼) | 정수 | — | ✔ `genome_length` |
| `date_released` `date_updated` | P2 | 메타데이터 표 | `YYYY-MM-DD` | — | ✔ `release_date` `update_date` |
| `institution` `authors` | P2 | 메타데이터 표 | 문자열 | — | ✔ `submitter` |

> 분절 바이러스의 **연결(가안)**: 지금 파이프라인은 `webapp/strain_link.py`가 `strain` 이름 일치로 S·L 레코드를 묶는다(분석 코어는 그대로, 결과를 재가공). 실제 데이터에서 이어진 300쌍은 국가·리니지가 모두 일치했지만, 이름 충돌과 표기 차이는 원천에서 해결하는 것이 맞다.

## 2. 값이 지켜야 할 규칙(자세히)

**host → 분류(`scripts/host_scope.py`)**

| 값 | 분류 |
| --- | --- |
| `is_lab_host` 참 | **lab**(모든 범위에서 제외) |
| `host_taxid` = 9606 | human |
| 빈 값, `unknown`, `na`, `n/a`, `missing`, `not collected`, `not provided`, `unspecified`, `-`, `?` | unknown |
| `non-human`·`nonhuman` 포함 | animal |
| `homo sapiens*`, `human`, `patient`, `man`, `woman`, `person` | human |
| 그 밖의 값(`Mastomys natalensis`, `Pteropus lylei` …) | animal |

→ GC는 숙주를 **학명으로 표준화**하고, 세포주·환경 시료는 `host`에 넣지 말고 `is_lab_host = true`로 돌려야 한다. 현재 사전은 세포주를 "숙주 아님"으로 표시만 하고 값은 그대로 둔다(`Vero E6`가 animal로 분류됨).

**date** — 정규식 `^(\d{4}|XXXX)-(\d{2}|XX)-(\d{2}|XX)$`. 범위 표기(`2022-05-01/2022-05-28`)는 한 값으로 줄이거나 부분 `XX` 처리.

**lineage** — 항상 비어 있지 않을 것. 분류 못 하면 `unclassified`. 값 목록은 바이러스별로 닫힌 집합(§1)이고, 새 값이 생기면 GC 사전에 추가한 뒤 내보낸다.

**accession** — 대소문자·공백·`|` 없음. 한 분절 파일 안에서 유일.

**is_lab_host** — 빈 값 또는 `true`만 쓴다(`false`·`0`은 빈 값과 같게 처리되지만 혼동을 피한다).

참고 — 현재 Nextstrain 원천의 실제 값(우리가 쓰던 것): Lassa S 모집단 662개, 리니지 `LII 397 / LIV_V 183 / LIII 64 / LVII 11 / LI 5 / LVI 1 / unclassified 1`, 지역 `Africa 620 / Unknown/Missing 41 / North America 1`. 날짜는 `2012-XX-XX`, `XXXX-XX-XX` 형태. NiV 모집단 `Ia 67 / II 26 / Ib 10 / unclassified 2 / I 1`. 숙주에 `Homo sapiens`, `Mastomys natalensis`, `Pteropus lylei`, 실험 숙주 7건(`is_lab_host=true`).

## 3. 서열(FASTA)과 파일 구성

| 항목 | 우선 | 규칙 | GC 현재 |
| --- | --- | --- | --- |
| 헤더 | **P0** | `>accession`(버전 없음). 공백 뒤 설명은 무시됨 | △ 버전 포함 |
| 분절별 파일 | **P0** | `{seg}_sequences.fasta` + `{seg}_metadata.tsv`. 서열마다 메타 행이 있을 것(행 수 일치) | ✘ 한 파일로 나감 |
| 서열 문자 | P0 | 정렬하지 않은 원서열, IUPAC 문자, 소문자 허용(우리가 대문자로 바꿈), **갭(`-`) 금지** | ✔ 원서열 저장 |
| 방향 | — | 원본 방향 그대로 둔다. 우리가 참조 유전체 방향으로 맞춘다(Lassa S는 662개 중 619개를 뒤집음) | ✔ |
| 기준 유전체 포함 | P1 | 참조 accession(Lassa S `AY628203`, L `AY628202`, NiV `NC_002728`)이 모집단에 있을 것. 없으면 첫 유전체가 방향 기준이 됨 | ? RefSeq 포함 여부 |
| 길이·N 필터 | — | 우리가 처리한다(S ≥ 3000, L ≥ 6000, NiV ≥ 17000, N 비율 한도). GC가 부분 서열을 미리 빼지 않아도 됨 | — |

## 4. 데이터셋·내보내기·재현성

| 항목 | 우선 | 요구 | GC 현재 |
| --- | --- | --- | --- |
| 규약 버전 | **P0** | 내보내기 열 이름·의미를 `schema_version`과 함께 문서로 고정. 바뀌면 버전을 올린다 | ✘ |
| 불변 스냅샷 | **P0** | 같은 해시 = 같은 내용. 결과마다 "어느 데이터셋의 어느 스냅샷으로 계산했는지"를 우리가 기록한다 | △ 매니페스트에 체크섬·규칙·멤버 해시 ✔ |
| 매니페스트 내용 | P1 | 생성 시각, 데이터셋 이름·규칙, 멤버 해시, 파일별 체크섬, 분절별 건수, 제외 건수, 규약 버전, GC 버전 | △ 일부 ✔ |
| 변경 목록(diff) | P1 | 이전 스냅샷 대비 추가·제거·수정된 accession (신선도 알림과 "반영 전후 차이" 화면용) | ✘ |
| 분절 조건 | P1 | 데이터셋 규칙에 `segment` 추가 | ✘ (규칙 항목: clade·lineage·outbreak·region·country·host·completeness·qc_not 등) |
| 비대화형 실행 | P1 | CLI 한 줄(또는 내부 API)로 데이터셋을 내보내고 종료 코드로 성공을 알림 | ✔ `genomecurator export`, `/api/export/...` |
| 공개 데이터만 | — | Pathoplexus 제한 사용 레코드는 나가지 않음 | ✔ |

## 5. 우리 쪽 수용 검사(내보낸 파일을 받을 때 자동으로 확인할 것)

1. 규약 버전이 기대와 같은가.
2. 모든 FASTA id가 메타데이터 `accession`에 있고, 반대도 그러한가(행 수 일치).
3. `accession`이 분절 파일 안에서 유일하고 버전 접미사(`.숫자`)가 없는가.
4. `lineage`가 비어 있지 않고 §1 목록 안에 있는가(밖의 값은 경고로 표시).
5. `date`가 정규식에 맞는가, `region`이 허용 값인가.
6. 기준 유전체가 들어 있는가.
7. 이전 스냅샷과 비교해 건수·리니지 분포가 크게 달라지지 않았는가(크게 달라지면 사람에게 확인).

## 6. GC 개발 우선순위(이 문서 기준)

| 순서 | 내용 | 이유 |
| --- | --- | --- |
| P0-1 | 내보내기 프로파일 `pipeline` 신설: 버전 없는 accession, FASTA 헤더 = accession, 분절별 파일, 표준 date·region·host 열 | 이것만으로 연동이 가능해짐 |
| P0-2 | `is_lab_host` 열(사전의 세포주 판정 노출) | 세포주 서열 혼입 방지 |
| P0-3 | Lassa·Nipah 리니지/클레이드 부여(출처와 근거 기록, 미분류는 `unclassified`) | 리니지 분석의 전제 |
| P0-4 | 규약 문서와 `schema_version` | 원천이 계속 바뀌어도 소비자가 깨지지 않게 함 |
| P1 | 스냅샷 diff, 분절 조건, RefSeq/GenBank 중복 정리, `host_taxid` 노출, 제외 목록 내보내기 | 재현성과 정확도 |
| P2 | completeness·QC 열, 국가 이름을 우리 지도 이름과 맞추는 점검 | 보조 |
