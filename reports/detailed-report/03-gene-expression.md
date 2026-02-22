# 모듈 03: 유전자 발현 예측 (Gene Expression — RNA_SEQ · CAGE · PROCAP)

## 1. 한 줄 요약

> RNA_SEQ(667 tracks), CAGE(546 tracks), PROCAP(12 tracks) 세 가지 유전자 발현 OutputType을 검증했으며, GeneMaskLFCScorer를 통해 유전자 수준 log-fold-change를 정량화한다. PROCAP은 human only(12 tracks, 6 cell lines)이다.

## 2. 왜 이 분석이 필요했나

유전자 발현(gene expression)은 "DNA에서 RNA로 전사되는 과정"을 의미하며, 세포 기능의 핵심 메커니즘이다. AlphaGenome은 세 가지 관점에서 유전자 발현을 예측한다.

1. **RNA_SEQ**: 세포 내 RNA 분자의 양을 측정 (전체 유전자 발현량)
2. **CAGE**: 전사 시작 지점(TSS)의 활성도를 측정 (어디서 전사가 시작되는가)
3. **PROCAP**: RNA Polymerase II의 활성도를 측정 (전사 기계의 위치와 활동)

변이가 유전자 발현에 미치는 영향을 정량화하려면, 이 세 가지 신호를 모두 확인해야 한다. RNA_SEQ만으로는 "왜" 발현이 변했는지 알 수 없고, CAGE와 PROCAP을 함께 보면 TSS 변경인지 Polymerase 활성 변화인지 구분할 수 있다.

## 3. 분석 과정 (The Mechanics)

### Quick Reference

```python
from alphagenome_dna_client import client as dna_client
from alphagenome_dna_client.variant_scorers import GeneMaskLFCScorer

# 1. RNA_SEQ: 유전자 발현 수준 예측
output = dna_model.predict_interval(
    interval=genome.Interval(chromosome="chr19", start=40991281, end=41018398, strand="+"),
    requested_outputs={dna_client.OutputType.RNA_SEQ},
)
output.rna_seq.values.shape  # (27117, 667) — 667 human tracks

# 2. CAGE: 전사 시작 지점 활성도 (strand-specific)
output_cage = dna_model.predict_interval(
    interval=genome.Interval(chromosome="chr1", start=1000000, end=2048576, strand="."),
    requested_outputs={dna_client.OutputType.CAGE},
    ontology_terms=['UBERON:0000955', 'UBERON:0002048'],  # Brain, Lung
)
output_cage.cage.values.shape  # (1048576, 4) — 2 tissues x 2 strands

# 3. PROCAP: RNA Polymerase II 활성도 (human only)
PROCAP_CELL_LINES = {
    'K562': 'EFO:0002067',
    'Caco-2': 'EFO:0001099',
    'A673': 'EFO:0002106',
    'Calu3': 'EFO:0002819',
    'MCF 10A': 'EFO:0001200',
    'HUVEC': 'CL:0002618',
}
output_procap = dna_model.predict_interval(
    interval=interval,
    requested_outputs={dna_client.OutputType.PROCAP},
    ontology_terms=list(PROCAP_CELL_LINES.values()),
)
output_procap.procap.values.shape[-1]  # 12 tracks (6 cell lines x 2 strands)

# 4. 변이 스코어링 (GeneMaskLFCScorer)
scorer = GeneMaskLFCScorer(
    gene_interval=gene_interval,
    output_type=dna_client.OutputType.RNA_SEQ,
)
variant_scores = dna_model.score_variants(
    interval=interval,
    variants=[genome.Variant(chromosome="chr19", position=40991639, ...)],
    scorer=scorer,
)
variant_scores.variant_scores.shape  # (37, 667) — 37 variants x 667 tracks
```

### IPO 요약

| Phase | Input | Process | Output |
|-------|-------|---------|--------|
| **RNA_SEQ** | chr19:40991281-41018398:+ (CYP2B6 locus) | GeneMaskLFCScorer(gene mask × 10 genes) | (37, 667) variant scores, PNG 121 KB |
| **CAGE** | chr1:1000000-2048576:. (Brain, Lung) | 2 tissues × 2 strands (strand-specific) | (1048576, 4) predictions, PNG 60 KB |
| **PROCAP** | chr22:35677410-36725986:. (LARGE locus) | 6 cell lines × 2 strands | 12 tracks, PNG ~100 KB |
| **Batch Result** | chr16:636337:G>A | GeneMaskLFCScorer across 687 ontologies | Neutrophil CAGE +0.3135 (99.5th percentile) |

### RNA_SEQ (667 tracks)

| Item | Value |
|------|-------|
| **Description** | RNA 발현 수준 예측 (Gene expression coverage) |
| **Resolution** | 1 bp |
| **Human Tracks** | 667 |
| **Mouse Tracks** | 173 |
| **Strand** | Strand-specific (+/-) |
| **Scorer** | GeneMaskLFCScorer (gene-level log-fold-change) |
| **Used in** | Quick Start, Viz Tour, CLI, Batch, Analysis Workflow |

#### Quick Start 검증 (CYP2B6 locus)

CYP2B6 유전자(cytochrome P450 효소, 약물 대사 관련)를 대상으로 RNA_SEQ 예측을 수행했다.

```python
interval = genome.Interval(chromosome="chr19", start=40991281, end=41018398, strand="+")
# 29 transcripts, 10 genes in this region

output = dna_model.predict_interval(
    interval=interval,
    requested_outputs={dna_client.OutputType.RNA_SEQ},
)
output.rna_seq.values.shape  # (27117, 667)
```

**GeneMaskLFCScorer 동작 원리**:
1. 유전자 마스크 생성: 10개 유전자 각각에 대해 exon 영역만 True로 마스킹
2. Reference vs Variant 예측: 변이 도입 전후의 RNA_SEQ 신호 예측
3. Log-fold-change 계산: `log2(variant / reference)` (마스크 영역만 평균)

```python
# Quick Start 결과
variant_scores.variant_scores.shape  # (37, 667)
# 37 variants x 667 human tracks
# 각 값은 해당 변이가 특정 조직/세포에서 유전자 발현을 얼마나 변화시키는지 (log2FC)
```

**시각화 결과** (`results/quick_start/cyp2b6_rna_seq.png`):
- X축: 게놈 좌표 (chr19:40991281-41018398)
- Y축: RNA_SEQ 신호 강도
- 667개 트랙 중 대표 조직(간, 폐, 뇌 등) 선택하여 표시
- 변이 위치에 수직선 표시

#### Variant Scoring 예시

Quick Start에서 사용한 37개 변이 중 일부:

```python
variants = [
    genome.Variant(chromosome="chr19", position=40991639, reference_bases="G", alternate_bases="A"),
    genome.Variant(chromosome="chr19", position=40991660, reference_bases="C", alternate_bases="T"),
    # ... 35 more
]
```

**결과 해석**:
- Positive score (+0.5): 변이로 인해 해당 조직에서 RNA 발현량이 1.4배(2^0.5) 증가
- Negative score (-0.5): 변이로 인해 해당 조직에서 RNA 발현량이 0.7배(2^-0.5) 감소
- Zero score (0.0): 발현량 변화 없음

#### Visualization Tour 결과

```python
# results/visualization_tour/01_rna_seq.png
# 파일 크기: 121 KB
# Ontology terms: Sigmoid colon (UBERON:0001159), Transverse colon (UBERON:0001157)
```

### CAGE (546 tracks)

| Item | Value |
|------|-------|
| **Description** | Cap Analysis Gene Expression — TSS activity prediction |
| **Resolution** | 1 bp |
| **Human Tracks** | 546 |
| **Mouse Tracks** | 188 |
| **Strand** | Strand-specific (+/-) |
| **Scorer** | GeneMaskLFCScorer / CenterMaskScorer |
| **Used in** | Quick Start, Viz Tour, Batch |

#### CAGE의 특징

CAGE(Cap Analysis Gene Expression)는 mRNA의 5' cap 구조를 시퀀싱하여 **전사 시작 지점(TSS)**을 1bp 해상도로 식별한다.

**왜 CAGE가 중요한가**:
1. RNA_SEQ는 "얼마나 발현되는가"를 측정
2. CAGE는 "어디서 전사가 시작되는가"를 측정
3. 변이가 TSS를 변경하면 (예: 프로모터 변이) CAGE 신호가 이동하거나 소실됨

**Strand-specific 특성**:
- 각 ontology term마다 +/- 두 개의 트랙 생성
- 예: Brain → Brain(+), Brain(-)
- Quick Start 예시: 2 tissues × 2 strands = 4 tracks

```python
ontology_terms = ['UBERON:0000955', 'UBERON:0002048']  # Brain, Lung
output_cage = dna_model.predict_interval(
    interval=genome.Interval(chromosome="chr1", start=1000000, end=2048576, strand="."),
    requested_outputs={dna_client.OutputType.CAGE},
    ontology_terms=ontology_terms,
)
output_cage.cage.values.shape  # (1048576, 4)
# Brain(+), Brain(-), Lung(+), Lung(-)
```

#### Batch Result 분석

Batch 실험(5 variants × 687 ontologies)에서 CAGE 변이 스코어링 결과:

**chr16:636337:G>A 변이** (가장 강한 positive effect):
- **Neutrophil CAGE**: +0.3135 (99.5th percentile across all scores)
- **Eosinophil CAGE**: +0.3001 (99.2nd percentile)
- **해석**: 이 변이는 호중구와 호산구에서 새로운 TSS를 생성하거나 기존 TSS를 강화함

**CAGE vs RNA_SEQ 비교**:
- RNA_SEQ가 높지만 CAGE가 낮다 → 전사 후 조절(post-transcriptional regulation) 가능성
- CAGE가 높지만 RNA_SEQ가 낮다 → mRNA 불안정성(degradation) 가능성
- 둘 다 높다 → 전사 수준 활성화(transcriptional activation)

#### Visualization Tour 결과

```python
# results/visualization_tour/02_cage.png
# 파일 크기: 60 KB
# Strand-specific TSS peaks 시각화
```

시각화에서 확인할 수 있는 패턴:
1. **Sharp peaks**: TSS 위치 (대부분 유전자 5' 끝)
2. **Bidirectional peaks**: Divergent transcription (양방향 전사)
3. **Multiple peaks**: Alternative TSS 사용

### PROCAP (12 tracks)

| Item | Value |
|------|-------|
| **Description** | RNA Polymerase II activity (Precision Run-On) |
| **Resolution** | 1 bp |
| **Human Tracks** | 12 |
| **Mouse Tracks** | **- (not supported)** |
| **Strand** | Strand-specific (+/-) |
| **Scorer** | GeneMaskLFCScorer |
| **Used in** | PROCAP Visualization |

#### PROCAP의 특징

PROCAP(Precision Run-On and sequencing of capped RNA)은 RNA Polymerase II의 활성 위치를 측정한다.

**CAGE vs PROCAP 차이**:
- **CAGE**: mRNA의 5' cap 구조를 시퀀싱 → "완성된 transcript의 TSS"
- **PROCAP**: 전사 중인 RNA Pol II를 포획 → "현재 전사 중인 위치"
- **의미**: PROCAP이 더 동적(dynamic)이며, pausing/elongation 상태를 반영

**CRITICAL: Mouse NOT Supported**:
```python
# PROCAP은 human genome (GRCh38)에서만 사용 가능
# Mouse genome (GRCm39)에서 PROCAP 요청 시 에러 발생
```

#### 지원 Cell Lines (6개)

| Cell Line | Ontology Term | Description |
|-----------|---------------|-------------|
| **A673** | EFO:0002106 | Ewing sarcoma cell line |
| **Caco-2** | EFO:0001099 | Colorectal adenocarcinoma cell line |
| **K562** | EFO:0002067 | Chronic myelogenous leukemia cell line |
| **Calu3** | EFO:0002819 | Lung adenocarcinoma cell line |
| **MCF 10A** | EFO:0001200 | Breast epithelial cell line |
| **HUVEC** | CL:0002618 | Human umbilical vein endothelial cells |

각 cell line마다 +/- 2개 strand → 6 × 2 = **12 tracks**

#### Script 구현 (`scripts/run_procap_visualization.py`)

```python
PROCAP_CELL_LINES = {
    'A673': 'EFO:0002106',
    'Caco-2': 'EFO:0001099',
    'K562': 'EFO:0002067',
    'Calu3': 'EFO:0002819',
    'MCF 10A': 'EFO:0001200',
    'HUVEC': 'CL:0002618',
}

interval = genome.Interval(chromosome="chr22", start=35677410, end=36725986, strand=".")
# LARGE locus (dystroglycanopathy gene, ~1 Mbp)

output = dna_model.predict_interval(
    interval=interval,
    requested_outputs={dna_client.OutputType.PROCAP},
    ontology_terms=list(PROCAP_CELL_LINES.values()),
)
output.procap.values.shape[-1]  # 12 tracks
```

**결과** (`results/procap_visualization/procap.png`):
- X축: chr22:35677410-36725986 (1.05 Mbp)
- Y축: 12개 트랙 (6 cell lines × 2 strands)
- LARGE 유전자의 exon 위치에서 PROCAP 신호 peak 확인

**해석**:
- Exon-intron junction 근처에서 Pol II pausing 관찰
- K562와 Caco-2에서 강한 신호 → 이 세포주에서 LARGE 유전자 활발히 전사
- HUVEC에서 약한 신호 → 내피세포에서는 발현 낮음

## 4. 왜 이 방법인가

### RNA_SEQ: 전체 발현량의 기준

RNA_SEQ는 세포 내 모든 RNA 분자의 양을 측정하므로, "이 유전자가 얼마나 발현되는가"라는 질문에 답한다. 변이의 기능적 영향을 평가할 때 가장 먼저 확인하는 지표이다.

**장점**:
- 667 human tracks → 다양한 조직/세포에서의 발현 패턴 비교 가능
- 173 mouse tracks → Cross-species conservation 분석 가능
- GeneMaskLFCScorer → 유전자 단위 정량화 (exon만 마스킹)

**한계**:
- "왜" 발현이 변했는지는 알 수 없음 (전사? 안정성? 번역?)
- 대안 스플라이싱(alternative splicing) 정보 부족 → SPLICE_SITES로 보완

### CAGE: TSS 변경 감지

프로모터 변이는 TSS 위치를 변경할 수 있다. CAGE의 strand-specific 특성은 이를 정확히 포착한다.

**왜 strand-specific인가**:
- Sense transcript와 antisense transcript를 구분
- 양방향 프로모터(bidirectional promoter) 감지
- Alternative TSS 사용 패턴 분석

**Batch 결과의 의미**:
- chr16:636337:G>A가 Neutrophil CAGE +0.3135 → 호중구 특이적 TSS 생성
- 이 변이는 염증 반응 관련 유전자의 발현을 조절할 가능성

### PROCAP: 전사 기계의 역학

PROCAP은 "지금 이 순간 어디서 전사가 일어나는가"를 보여준다. CAGE가 "최종 제품"이라면, PROCAP은 "공장 가동 상태"이다.

**왜 6개 cell line만 지원하는가**:
- PROCAP 실험 데이터가 이 6개 세포주에서만 존재 (DeepMind 학습 데이터 제약)
- 향후 추가 cell line 데이터로 확장 가능

**왜 mouse를 지원하지 않는가**:
- Mouse PROCAP 실험 데이터 부족
- Human genome 전용 학습 → GRCm39에서는 예측 불가

## 5. 해석과 한계

### 검증된 내용

1. **RNA_SEQ 667 tracks**: Quick Start(CYP2B6), Viz Tour, Batch, CLI, Analysis Workflow 모두 정상 작동
2. **CAGE 546 tracks**: Strand-specific 동작 확인 (Brain+/-, Lung+/-)
3. **PROCAP 12 tracks**: 6 cell lines × 2 strands, LARGE locus 시각화 성공
4. **GeneMaskLFCScorer**: 10 genes × 37 variants → (37, 667) variant scores
5. **Batch 최고 스코어**: chr16:636337:G>A, Neutrophil CAGE +0.3135 (99.5th percentile)

### 한계

#### 1. PROCAP의 제약
- **Human only**: Mouse genome에서 사용 불가
- **6 cell lines만 지원**: 조직 특이성 분석 제한
- 해결책: RNA_SEQ와 CAGE를 조합하여 간접 추론

#### 2. GeneMaskLFCScorer의 가정
- Exon 영역만 마스킹 → Intronic regulation 무시
- Log-fold-change 평균 → 위치별 차이 손실
- 해결책: CenterMaskScorer(TSS 중심 마스킹) 병행 사용

#### 3. Strand 정보의 복잡성
- CAGE/PROCAP은 strand-specific → 트랙 수 2배
- 일부 ontology는 strand 정보 없음 → "." strand 처리 필요
- 주의: Interval strand와 output strand를 혼동하지 말 것

#### 4. Batch vs CLI 스코어 범위 차이
- Batch: 121,550 rows (687 ontologies × 5 variants)
- CLI: 38,357 rows (19 scorers × 다양한 변이)
- CLI가 더 포괄적이지만, Batch가 조직 특이성 분석에 유리

### 실전 가이드

**변이 우선순위 결정**:
```python
# Step 1: RNA_SEQ로 발현량 변화 확인
rna_scores = variant_scores[variant_scores['scorer_type'] == 'RNA_SEQ']
top_rna = rna_scores[abs(rna_scores['score']) > 0.5]  # |log2FC| > 0.5 → 1.4배 변화

# Step 2: CAGE로 TSS 변경 확인
cage_scores = variant_scores[variant_scores['scorer_type'] == 'CAGE']
top_cage = cage_scores[abs(cage_scores['score']) > 0.3]

# Step 3: 교집합 찾기 (RNA + CAGE 모두 높은 변이)
high_impact = pd.merge(top_rna, top_cage, on='variant_id')
```

**조직 특이성 분석**:
```python
# Batch 결과에서 특정 변이의 조직별 효과 추출
variant_id = "chr16:636337:G>A"
tissue_effects = batch_df[batch_df['variant_id'] == variant_id].sort_values('score', ascending=False)

# Top 5 tissues
print(tissue_effects.head())
# Neutrophil +0.3135
# Eosinophil +0.3001
# ...
```

## 6. 다음으로 (The Bridge)

유전자 발현을 측정하는 3개 OutputType(RNA_SEQ 667 tracks, CAGE 546 tracks, PROCAP 12 tracks)을 검증했다. 이제 핵심 질문이 남았다.

**"유전자가 발현되려면 DNA가 열려있어야 한다."**

Chromatin accessibility(염색질 접근성)가 보장되지 않으면, 전사 기계(RNA Pol II, 전사 인자)가 DNA에 접근할 수 없다. 변이가 chromatin 구조를 변경하면:
- RNA_SEQ/CAGE 신호가 높아도 실제 발현 안 됨 (닫힌 chromatin)
- 반대로 chromatin이 열리면 잠재적 유전자가 활성화됨

**다음 단계**: DNASE(305 tracks)와 ATAC(167 tracks)로 염색질 접근성을 확인한다.

→ [다음: 04-chromatin.md](04-chromatin.md)
