# chichipes (구 primer-engine) 출력 포맷

chichipes은 파이프라인이 지금 쓰는 mfeprimer 4.2.4 명령을 그대로 대신할 수 있도록 **mfeprimer와 같은
파일 이름·열·형식**으로 출력합니다. 아래는 파이프라인 코드(`scripts/*.py`)가 실제로 읽는 부분을 기준으로 정리한 것입니다.
예시의 서열은 모두 가짜(설명용)이며 실제 프라이머가 아닙니다.

## 0. 파이프라인이 mfeprimer를 부르는 곳

| 호출 | 명령 | 읽는 쪽 | 엔진 상태 |
|---|---|---|---|
| 프라이머 쌍 검색 (virus 모드, mm2) | `-d db -i FR.fa -o X.txt --fasta --virus --misMatch 2 --cpu N` | `run_primer.py`, `build_alignment_table.py`, `live_analyze.py`, `group_cross_binding.py`, `audit_pairs.py` | 구현됨, 대조 중 (§1–4) |
| 단일 서열 검색 | `search db.fa oligo.fa out.txt -m M` | `run_refgenome_binding.py`(참조유전체 결합, M=6), `build_alignment_table.py`(프로브 넓은 검색, M=8) | 구현됨, 실제 결과와 100% 일치 (§5) |
| 다이머 | `dimer -i oligos.fa -o dimer.txt` (+옵션) | `run_dimer_hairpin.py` | 구현됨, 실제 182개 보고서 값 100% 일치 (§6) |
| 헤어핀 | `hairpin -i oligos.fa -o hairpin.txt` | `run_dimer_hairpin.py` | 구현됨, 실제 182개 보고서 바이트 단위 일치 (§6) |
| 인덱스 | `index -i db.fa` | `ensure_index` | 필요 없음 — 받기만 하고 무시 (28 GB 인덱스 불필요) |

공통: 좌표의 기준(0/1-based)이 파일마다 다르므로 각 절에 적었습니다. 숫자는 소수 둘째 자리(`%.2f`)입니다.

---

## 1. `<out>.spec.tsv` — 증폭산물 목록 (핵심)

```
#1-based coordinate, should be idential to blat output if align primer sequences to genome
#name	chrom	ampStart	ampEnd	ampGC	ampSize	fpName	fpStart	fpEnd	fpSeq	fpTm	fpGC	fpDg	rpName	rpEnd	rpStart	rpSeq	rpTm	rpGC	rpDg	note
Amp_1	GENOME01	809	897	51.69	89	SET_F.1	809	831	ACGT…AAC	58.85	47.83	-23.12	SET_R.1.8	878	897	CTGA…CTC	61.45	60.00	-23.37	
```

헤더 두 줄은 글자 그대로(오타 `idential` 포함) 씁니다. 행마다 탭 21칸(마지막 `note`는 빈 칸).

| 열 | 의미 |
|---|---|
| `name` | `Amp_<번호>`. mfeprimer는 멀티스레드라 **번호·행 순서가 실행마다 바뀜** → 비교는 번호를 빼고 행 내용으로 함 |
| `chrom` | 유전체 레코드 ID (FASTA 헤더 첫 단어) |
| `ampStart`, `ampEnd` | 증폭산물 1-based, 양끝 포함 |
| `ampGC` | **프라이머 서열을 끼운** 산물(fpSeq + 유전체 가운데 + rpSeq)의 GC% — 유전체 구간 그대로가 아님 |
| `ampSize` | `ampEnd − ampStart + 1` |
| `fpName` | 왼쪽(+가닥에 그대로 붙는) 프라이머의 이름. 오른쪽/왼쪽은 **위치 기준**이라 F+F, R+R, R+F 산물도 나옴 |
| `fpStart`, `fpEnd` | 왼쪽 프라이머 위치, 1-based 양끝 포함 |
| `fpSeq` | 실제로 붙은 **변이 서열**(5'→3'). 3' 9염기 창은 유전체와 정확히 일치하고, 그 밖의 미스매치는 이 글자에 드러나지 않음 |
| `fpTm`, `fpDg` | 그 결합의 Tm(°C), ΔG(kcal/mol, 37 °C) — 표적 문맥(옆 염기, 미스매치) 포함 |
| `fpGC` | fpSeq의 GC% |
| `rpName` | 오른쪽(−가닥에 붙는) 프라이머 이름 |
| `rpEnd`, `rpStart` | **순서 주의**: `rpEnd` = 작은 좌표(프라이머 3' 끝 쪽), `rpStart` = 큰 좌표 = `ampEnd` |
| `rpSeq` | 오른쪽 프라이머를 **+가닥으로 뒤집은**(역상보) 서열 |
| `rpTm`, `rpGC`, `rpDg` | 오른쪽 프라이머 결합 값 (`rpGC`는 프라이머 자체 GC%) |

### 이름 규칙 (`fpName`/`rpName`)
- 축퇴 염기가 없는 프라이머: 원래 이름 `SET_F`, 3' 창 변이는 `SET_F.1` … `SET_F.351`
- 축퇴 프라이머: 펼친 서열 `SET_F.1`, `SET_F.2`, … 각각의 변이 `SET_F.1.1` …
- 펼치는 순서(IUPAC): R→GA, Y→TC, M→AC, K→GT, S→GC, W→TA, B→GTC, D→GTA, H→TAC, V→GAC, N→GATC (뒤 위치가 가장 빨리 바뀜)
- 변이 순서: 미스매치 개수 → 3' 쪽 위치 조합(사전순) → 치환 염기(A→TCG, C→TAG, G→TCA, T→ACG)
- 이미 만든 서열과 같은 펼침/변이는 다시 만들지 않고 번호를 당겨 매김
- 실제 보고서 182개의 프라이머 표가 이름·순서까지 100% 일치 (`internal/primers`)

### 어떤 산물이 행이 되는가
1. 3' 끝 9염기 창의 미스매치 ≤ `--misMatch`
2. 결합 Tm이 `-t`(30) 이상 `-T`(100) 이하 — 양쪽 프라이머 모두
3. 산물 크기 ≤ `-S`(2000), ≥ `-s`(0)
4. 두 프라이머가 겹치면 겹침 < 9염기일 때만
5. 결합 부위의 3' 창에 A/C/G/T 이외 글자(N 등)가 있으면 그 자리는 찾지 않음

## 2. `<out>.amp.fasta` — 증폭산물 서열

```
>Amp_1 SET_F.1+SET_R.1.8==>GENOME01:808-897 data/filtered/genome_filtered.fasta
ACGTACGTAC…TTGCacgtacgtacgt…acgtGATCGATCGATCGATCGATC
```
- 헤더: `>이름 왼쪽+오른쪽==>chrom:시작(0-based)-끝 db경로` — 여기 좌표는 **0-based 시작**(spec.tsv와 1 차이)
- 서열: 왼쪽 프라이머(대문자) + 유전체 가운데(소문자) + 오른쪽 프라이머의 +가닥 서열(대문자)
- 처음 `-M`(10000)개 산물만 씀 → 파이프라인은 서열을 이 파일이 아니라 spec.tsv 좌표로 유전체에서 잘라 씀
- 두 프라이머가 겹치는 산물(전체의 0.0001%): mfeprimer는 서열을 불규칙하게 씀. 엔진은 유전체 구간 대문자로 씀 (알려진 차이, 파이프라인 영향 없음)

## 3. `<out>.virus.Failed.txt` — 산물이 하나도 없는 유전체

```
GENOME07 
GENOME19 
```
- 한 줄에 ID 하나, **뒤에 공백 한 칸**, DB 순서
- 파이프라인은 첫 단어만 읽음 (`parse_primer_coverage`)

## 4. `<out>` (보고서) — 커버리지 줄

파이프라인은 정규식으로 이 한 줄만 읽습니다:
```
Coverage 99.10%, total records in database: 2900, amplicon number: 5812, missed number: 26.
```
`Coverage = (전체 − missed) / 전체 × 100`, `amplicon number` = spec.tsv 행 수.
mfeprimer 보고서의 나머지(프라이머 표, 산물 설명 13만 줄 등)는 파이프라인이 읽지 않아 엔진은 요약만 씁니다.

**산물이 0개일 때** (직접 실행해 확인): mfeprimer는 `spec.tsv`를 **만들지 않고**, `amp.fasta`는 빈 파일,
`Failed.txt`에는 모든 레코드, 보고서에는 `Coverage 0.00%, …, amplicon number: 0, missed number: N.` — 엔진도 똑같이 씀.
(파이프라인은 빠진 spec.tsv를 빈 파일로 채움)

## 4b. `<out>.borderline.tsv` — 경계 사례 (엔진 전용)

mfeprimer에는 없는 파일. 결합 Tm이 `-t` 기준값(30 °C)의 ±0.02 °C 안에 있는 결합을 모두 적어요. 결합 Tm은 이제
mfeprimer와 같은 표로 계산해서 값이 일치해요(FINDINGS §3). 이 파일은 판정이 소수 끝자리에 달린 결합을 알려 주는
참고용이에요.
```
# bindings with Tm within 0.02 C of the cut-off; decision may differ from mfeprimer
#chrom	start	end	strand	name	seq	tm	dg	decision
GENOME01	18082	18107	-	SET_R.339	TGCT…GTCC	30.006	-9.530	kept
```
좌표 1-based 양끝 포함, strand `+` = 프라이머가 −가닥에 붙음(왼쪽), `-` = +가닥에 붙음(오른쪽), Tm·ΔG는 소수 셋째 자리,
`decision`은 엔진의 판정(kept/dropped). 파이프라인은 읽지 않습니다.

## 5. `search` 출력

```
# mfeprimer search
# max allowed mismatch: 8

#Chrom	Start(0-based)	End	QueryName	Size	Strand	MisMatch	DeltaG	Tm	QuerySeq	HitSeq
GENOME01	2859	2879	SET_P.1	20	-	0	-16.97	49.61	ATGA…AATGG	ATGA…AATGG
```
헤더 4줄(빈 줄 포함) 다음 한 줄에 한 hit. 좌표 0-based 반열림 `[Start, End)`, `Size` = End−Start.
- `Strand` `+`: 질의가 +가닥에 그대로 읽힘, `-`: 역상보로 읽힘. `HitSeq`는 질의 방향으로 쓴 유전체 서열(− 가닥이면 IUPAC 상보)
- `MisMatch`: 질의와 HitSeq의 편집거리(삽입·결실 포함)
- `DeltaG`, `Tm`: 길이가 같으면 결합값, 갭이 있으면 `0.00`
- `QueryName`: 축퇴 없으면 이름 그대로, 축퇴면 `NAME.1`, `NAME.2` …
- 찾는 조건: 질의 양끝 9염기가 모두 정확히 일치 + 편집거리 ≤ `-m` (FINDINGS §6)

파이프라인은 유전체마다 MisMatch가 가장 작은 행 하나를 씁니다.

## 6. `dimer` / `hairpin` 출력

두 보고서는 같은 틀이에요. 첫 줄의 시각과 맨 끝 꼬리말(mfeprimer는 인용문·연락처·실행 시간, 엔진은
`chichipes (...)` 한 줄)만 다르고, 그 사이는 mfeprimer와 같은 바이트로 써요.

```
MFEprimer-3.0 Hairpin Reports (2026-11-14 05:48:54)

Primer ID                      Sequence (5'-->3')                    Length     GC      Tm      Dg
                                                                      (bp)      (%)    (°C)  (kcal/mol)

F.1                            <서열>                              18   61.11   60.46    -21.06
...



Hairpin List (1)


Hairpin 1: R

    Score: 7, Tm = 5.63 °C, Delta G = -1.54 kcal/mol
    
    ///////---\\\\\\\----
    <서열>





Parameters
...
```

- 프라이머 표: `%-31s%-36s%8d%8.2f%8.2f%10.2f` (이름, 서열, 길이, GC%, Tm, ΔG). 이 표의 Tm은 기체 상수 1.9872로
  계산해요(`thermo` 명령은 1.987). 축퇴 서열은 펼쳐서 `이름.1`, `이름.2`, … 로 써요.
- 다이머 항목: `Dimer k: A x B`, `Score: s, Delta G = g kcal/mol`, 정렬 그림 세 줄(`5' A 3'` / `:`·`.` 표시 /
  `3' B 5'`). 항목 사이 빈 줄 2개, 마지막 항목 뒤 빈 줄 5개.
- 헤어핀 항목: `Hairpin k: 이름`, `Score: s, Tm = t °C, Delta G = g kcal/mol`, 공백 4칸만 있는 줄, 구조 줄
  (`/` 5′ 쪽 팔, `\` 3′ 쪽 팔, `-` 짝 없음), 서열 줄. 항목 사이 빈 줄 2개, 마지막 항목 뒤 빈 줄 4개.
- 올리고 하나당 헤어핀은 최대 1개. 목록은 ΔG 오름차순.
- 빈 목록: `Dimer List (0)` 다음 줄에 `No dimer found.`, `Hairpin List (0)` 다음 줄에 `No hairpins found.`.
- ΔG가 같은 항목끼리의 순서는 mfeprimer도 `--cpu` > 1에서 실행마다 달라서 재현 대상이 아니에요. 다이머 실제
  보고서 66개의 텍스트 차이가 전부 이 순서 차이예요.
- 파이프라인은 이를 `dimer_hairpin.json`(역할별 최악 ΔG, 그때의 Score·Tm·그림, ΔG ≤ −9 개수)으로 요약해 써요.
- virus 모드 primer 검색도 부산물로 `<out>.dimer.tsv`(`S1Name S1Seq S1GC S2Name S2Seq S2GC Score DeltaG`),
  `<out>.hairpin.tsv`를 쓰지만 파이프라인은 읽지 않아요.

## 7. 비교 기준 (재현 판정)

`harness/compare_spec.py`, `harness/run_all.sh`:
- 행 일치: `(chrom, fpName, fpStart, rpName, rpStart)`가 같은 행 — 번호(`Amp_n`)·순서는 무시
- 값 일치: Tm/ΔG/GC는 ±0.01 이내(반올림 차), 나머지는 글자 그대로
- 실패 목록: 집합과 순서
- 커버리지 줄
