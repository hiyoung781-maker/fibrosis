# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## 프로젝트 성격

MASH 유래 대상성 간경변(F4) 타겟 발굴 → REINVENT4 소분자 lead 설계까지를 3주에 수행하는 학부 심포지엄 연구. 코드베이스가 아니라 **사전등록된 연구 프로토콜의 실행 기록**이다. 소프트웨어 제품이 아니므로 빌드·린트·CI가 없고, 산출물은 `results/` 아래의 TSV/RDS/PNG와 포스터다.

**`plan.md`가 유일한 기준 문서다.** 분석 설계, 게이트 통과 기준, 데이터셋 정의, 대조군 목록, 일정, 변경 기록이 모두 여기 있다. 작업 전에 해당 절을 읽고, 코드가 아니라 plan.md가 정한 숫자를 따른다.

## 사전등록 규칙 (가장 중요)

- plan.md §11이 **사전등록 변경 기록**이다. 설계·기준·대비(contrast) 정의를 바꾸면 **날짜·이유·근거를 §11에 반드시 추가**한다. 기록 없는 변경은 금지.
- **최종 랭킹을 본 뒤에 기준을 바꾸지 않는다.** 특히 기준을 **완화**하는 방향의 변경은 2026-09-17 이후 금지로 선언되어 있다(§11 2차 개정). 기준 미달이면 미달로 보고하고 이유를 분석하는 것이 정답이다 (§9 "가장 조심할 것").
- 게이트 통과/탈락 기준은 사전에 숫자로 고정돼 있다. 결과를 맞추려고 임계값(FDR 0.05, |log2FC| 0.5, B≥10, AUC≥0.75, 백분위 75/25 등)을 조정하지 말 것.
- 사후에 떠오른 아이디어는 §11 "미결 항목"에 두고, 승인 전에는 계획으로 취급하지 않는다.

## 파이프라인 현재 상태

| 게이트 | 내용 | 상태 |
| --- | --- | --- |
| GATE 1 | Open Targets 소분자 tractability 등급 선필터 (킬 게이트) | 구 주분석(F4 vs F0+대조)으로 **통과** (2026-09-16: 등급 B 295, A–C 457). **신 주분석 F4 vs F0로 재실행 필요** |
| GATE 2 | GSE162694 독립 재현 (방향일치 ≥70%, ρ≥0.3) | 구 정의로 계산 시 통과(89.6%, ρ=0.616, 민감도 3으로 공개). **신 주분석으로 재실행 필요** |
| GATE 3 | GSE202379/GSE136103 세포형 귀속 + 세포 구성 판정 + decoupleR TF | 미실행 |
| GATE 4 | Geneformer in silico perturbation 재순위 + §2 대조군 채점 | 미실행 |
| GATE 5 | ChEMBL ≥300 활성 + High-Quality Ligand 보유 확인 | 미실행 |

**주분석 대비는 F4 vs F0다** (2026-09-17 3차 개정, plan.md §11). 대조군(정상 간)은 기준군에 넣지 않는다 — 기준군 안의 정상 비율이 코호트마다 달라(GSE135251 17%, GSE162694 47%) 두 코호트의 기준군이 다른 집단이 되기 때문이다. 구 주분석 F4 vs F0+대조는 민감도 3으로 남아 있다.

| 코호트 | 주분석 F4 vs F0 | 민감도 1 F4 vs F0–F1 | 민감도 2 F3F4 vs F0–F1 | 민감도 3 F4 vs F0+대조 |
| --- | --- | --- | --- | --- |
| GSE135251 | 14 vs 38, DEG 1,868 (완료) | 14 vs 85, 2,057 | 68 vs 85, 874 | 14 vs 46, 1,737 |
| GSE162694 | 12 vs 35, DEG 9,482 | 12 vs 65, 10,274 | 20 vs 65, 4,383 | 12 vs 66, 11,475 |

GSE162694에는 보조 대비 F4 vs 정상(12 vs 31, DEG 12,187, 파일명 `main_F4_vs_ctrl`)도 있다.

추가 대비(F4 vs F3, 단계 추세 스플라인 LRT)는 아직 없다.

## 설계상 반드시 지켜야 할 것

- **채점 범위는 DEG 목록이 아니다.** §2 채점은 소분자 등급 A–C 유전자 전체(≈5,007개) 순위에서 한다. 양성 대조군 5개 중 F4 vs F0 DEG가 0개이기 때문이다(THRB 포함).
- **대조군은 DEG 여부·게이트 통과 여부와 무관하게 GATE 3·4를 끝까지 통과시킨다.** 필터로 떨어뜨리면 채점 자체가 불가능해진다.
- **대조군은 약물 프로그램 단위로 묶는다.** PPARA/D/G = 1단위, CCR2/CCR5 = 1단위, CASP 9개 = 1단위, COL4A3/4/5 = 1단위. 단위 점수는 최상위 구성원. 유전자 수로 세지 말 것.
- **치료 방향을 후보와 대조군에 동일하게 적용한다.** THRB·PPAR은 활성화가 치료 방향, 나머지는 억제. 방향 보정 없으면 AUC가 방향 차이를 측정한다.
- **`sm_evidence`(strong/weak) 열은 순위 계산에 절대 쓰지 않는다.** 용도는 동점 처리와 GATE 5 진입 판단뿐.
- **코호트를 ComBat으로 병합하지 않는다.** 코호트별 독립 DEG 후 순위·방향성만 비교.
- **두 코호트의 대비는 같은 정의로 맞춘다.** 이름이 같아도 그룹 조성(예: 기준군 내 정상 비율)이 다르면 같은 대비가 아니다. 새 대비를 추가할 때 양쪽 코호트의 구성 인원을 표로 나란히 확인한다.
- 단일세포 split은 **기증자 단위**(`GroupKFold(groups=donor_id)`)이며 코드에 `assert`로 강제한다.
- 음성 대조군이 상위에 와도 강등·삭제하지 않는다. 주석(양식/실패 사유/결합 입증)을 달아 결과로 보고한다.
- wet lab은 수행하지 않는다. 실험 자리는 독립 코호트 재현·세포형 귀속·임상 채점·도킹이 대신한다.

## 디렉터리

```
plan.md                     사전등록 프로토콜 (기준 문서)
geo_data/                   GEO 원자료. GSE135251/ GSE162694/ GSE202379/
  GSE135251/GSE135251_samples.tsv   샘플 메타 216행, exclude·main_contrast 열 포함 (main_contrast 열은 구 주분석 = 민감도 3 정의)
  GSE162694/GSE162694_samples.tsv   샘플 메타 143행 (series matrix에서 파싱: gsm·col_id·fibrosis_stage·age·sex). col_id가 raw_counts.csv 열 이름
scripts/                    번호순 실행 스크립트 (01_ → 02_ → …)
results/deseq2/             DEG TSV + dds.rds + sessionInfo.txt + volcano_plots/ (+ plot_volcano.py)
results/gate1/              tractability 전체 표 + 등급 A·B DEG 목록
results/gate2/              재현 지표 CSV·리포트·산점도·RRHO 히트맵 (현재는 구 정의 결과)
opentarget/target/*.parquet Open Targets target 데이터셋 (tractability 필드)
opentarget/patcid*.jsonl    PatCID 특허-분자 매핑 (아직 파이프라인에 미사용)
Geneformer/                 HuggingFace ctheodoris/Geneformer 스냅샷 (vendored, V1-10M·V2-104M·V2-316M 가중치 포함)
REINVENT4/                  MolecularAI/REINVENT4 v4.8.24 클론 (vendored, 유일한 git 저장소)
geo_data.ipynb              탐색용 노트북 (GEO·Open Targets 구조 확인. 재현 경로 아님)
```

`Geneformer/`와 `REINVENT4/`는 외부 도구 스냅샷이다. 우리 코드는 `scripts/`와 `results/`의 플로팅 스크립트에만 넣는다. 외부 도구를 수정해야 하면 이유를 남긴다.

`plan.md`, `geo_data.ipynb`, 결과 주석은 한국어다. 스크립트 주석·출력은 영어다 — 각자의 관례를 유지한다.

## 실행

프로젝트 루트는 `/home/seyoung/academic_symposium/fibrosis`이고 스크립트가 절대경로 `BASE`를 하드코딩한다. conda 환경은 `~/miniconda3/envs/`.

```bash
# DESeq2 (R 4.4.1 + DESeq2 1.46 + apeglm, MulticoreParam(16))
~/miniconda3/envs/geo/bin/Rscript scripts/01_deseq2_GSE135251.R

# GATE 1 tractability 교차 (pandas + pyarrow)
~/miniconda3/envs/geo/bin/python scripts/02_gate1_tractability.py

# GSE162694 독립 DESeq2 (~group)
~/miniconda3/envs/geo/bin/Rscript scripts/03_deseq2_GSE162694.R

# GATE 2 재현 판정 (scipy 필요 → geneformer 환경)
~/miniconda3/envs/geneformer/bin/python scripts/04_gate2_replication.py

# 볼케이노 플롯 재생성 (FILES 목록에 없는 TSV는 SKIP으로 출력)
~/miniconda3/envs/geo/bin/python results/deseq2/plot_volcano.py

# Geneformer (torch 2.14+cu130, GPU 2장 인식)
~/miniconda3/envs/geneformer/bin/python ...
```

`02_gate1_tractability.py`는 stdout으로 게이트 판정(`GATE 1: ... -> PASS/FAIL`), 등급별 교차표, 대조군 유전자 등급·DEG 여부를 출력한다. 판정은 콘솔에만 남으므로 **실행 후 결과를 plan.md §4/§11에 옮겨 적는다.**

REINVENT4 테스트(외부 도구 검증용, 우리 파이프라인 아님):

```bash
cd REINVENT4 && pytest tests/                    # 기본 --device cpu
pytest tests/ -m "not integration"               # 통합 테스트 제외
pytest tests/scoring/test_foo.py::test_bar       # 단일 테스트
```

## 알려진 함정

- **주분석 전환 후 아직 옛 파일명을 읽는 스크립트가 있다.** `02_gate1_tractability.py`와 `04_gate2_replication.py`는 `GSE135251_main_F4_vs_F0ctrl.tsv`(구 주분석)를 입력으로 쓴다. GATE 1·2 재실행 전에 `GSE135251_main_F4_vs_F0.tsv` / `GSE162694_main_F4_vs_F0.tsv`로 바꾼다. 출력 파일명(`GSE135251_main_tractability.tsv` 등)도 구 결과를 덮어쓰지 않게 정한다.
- **`results/deseq2/`에 구 주분석 파일이 남아 있다.** `GSE135251_main_F4_vs_F0ctrl.{tsv,_dds.rds}`는 `GSE135251_sens3_F4_vs_F0ctrl`과 파일 단위로 동일한 중복본이고, `GSE162694_main_F4_vs_F0ctrl.*`와 `volcano_plots/*_main_F4_vs_F0ctrl_volcano.png`도 구 정의 산출물이다. 이름에 `main`이 붙어 있어도 주분석이 아니다.
- **GSE162694의 DEG 수는 기준군을 F0로 바꿔도 크게 줄지 않았다** (F0+정상 11,475 → F0 9,482; GSE135251은 1,868). 기준군 조성만으로는 설명되지 않으므로, 원인이 확인되기 전에는 GSE162694의 DEG 개수를 해석하지 말 것.
- **GSE162694에는 배치 정보가 없다.** GEO series matrix에 배치 필드가 없고 플랫폼도 GPL21290(HiSeq 3000) 하나뿐이라 `~group`으로 돌린다. 논문 보충자료의 배치 정보는 미확인.
- **`geo` 환경에 scipy가 없다.** 통계 검정이 필요한 Python 스크립트는 `geneformer` 환경(pandas·scipy·matplotlib 있음)으로 실행한다.
- **GPU 작업은 이 서버에서 하지 않는다.** 이 서버에는 TITAN Xp 2장뿐이고, Geneformer·REINVENT4 실행은 **A100 서버로 디렉터리를 복제해서** 수행한다. 따라서 GPU 단계 코드는 `BASE` 절대경로 하드코딩에 주의하고(복제 후 경로가 달라진다), 복제 전후 산출물 동기화 방향을 명확히 할 것.
- **`geneformer` 환경에서 `import geneformer`가 실패한다** — transformers 5.17이 설치돼 있고 `SpecialTokensMixin` import가 깨진다. GATE 4 착수 전에 transformers를 Geneformer가 요구하는 버전으로 내려야 한다.
- **REINVENT4는 설치되지 않았고 `priors/` 가중치도 없다.** 어느 conda 환경에도 `reinvent` 모듈이 없다. GATE 5 이후 `REINVENT4/install.py`로 별도 환경을 만들어야 한다(README 참조, Python ≥3.11).
- **프로젝트 루트는 git 저장소가 아니다** (`REINVENT4/.git`만 존재). plan.md Day 17의 "GitHub QR 코드"를 위해서는 먼저 `git init`이 필요하고, 45GB 원자료·모델 가중치는 반드시 커밋에서 제외한다.
- **GATE 3 도구가 아직 하나도 없다** — `geo` 환경에 decoupleR(R/Python 양쪽), LIANA, clusterProfiler, limma가 설치돼 있지 않다(R 라이브러리 84개뿐). GATE 3 착수 전에 설치 계획을 세운다.
