# 모듈 07: 변이 효과 정량화 (Variant Scoring — Batch · CLI · POLYADENYLATION)

## 1. 한 줄 요약

> Batch(5 variants, 6 scorers → 121,550행)와 CLI(1 variant, 19 scorers → 38,357행)로 변이 효과를 정량화했으며, POLYADENYLATION scorer(396행)까지 포함하여 19개 RECOMMENDED_VARIANT_SCORERS 전체를 검증했다.

**핵심 메시지**: AlphaGenome은 단일 염기 변이(SNV)가 DNA 조절 요소에 미치는 영향을 정량화할 수 있다. Batch는 여러 변이를 특정 scorer로, CLI는 단일 변이를 모든 scorer로 분석하며, 두 방식이 상호보완적으로 작동한다.

---

## 2. 왜 이 분석이 필요했나

### 문제 정의
- **병인 변이 해석**: 단일 염기 변이(SNV)가 전사 인자 결합, 크로마틴 접근성, RNA 안정성에 어떤 영향을 미치는지 정량화 필요
- **임상 우선순위 결정**: 수백만 개의 유전 변이 중 기능적으로 유의미한 것을 식별
- **세포 타입 특이성**: 동일 변이가 세포 타입별로 다른 영향을 미칠 수 있음
- **스케일 한계**: 실험적 검증(MPRA, CRISPR screen)은 비용과 시간이 많이 소요

### AlphaGenome의 해결책
- **Batch Variant Scoring**: 여러 변이를 특정 assay type(RNA-seq, ChIP-seq 등)으로 동시 평가
- **CLI Variant Scoring**: 단일 변이를 19개 모든 scorer로 포괄적 분석 (UI 기반에서 CLI 재구현)
- **POLYADENYLATION Scorer**: 19개 RECOMMENDED_VARIANT_SCORERS 중 OutputType 없이 scorer로만 존재하는 특수 케이스
- **점수 해석**: 참조 대비 대체 염기의 예측값 차이 (log-scale, cell-type specific)

---

## 3. 분석 과정 (The Mechanics)

### Quick Reference

| 분석 타입 | 스크립트 | 변이 수 | Scorer 수 | 결과 행 수 | 주요 용도 |
|----------|---------|---------|-----------|-----------|----------|
| **Batch** | `run_batch_variant_scoring.py` | 5 | 6 | 121,550 | 다중 변이 스크리닝 |
| **CLI** | `run_variant_scoring_cli.py` | 1 | 19 | 38,357 | 단일 변이 심층 분석 |
| **POLYADENYLATION** | CLI 내 포함 | 1 | 1 | 396 | Poly(A) site 영향 평가 |

---

### IPO 요약

#### Batch Variant Scoring

**Input**:
```python
variants = [
    genome.Variant("chr3", 58394738, "A", "T"),   # 5 variants
    genome.Variant("chr8", 28520, "G", "C"),      # covering different
    genome.Variant("chr16", 636337, "G", "A"),    # chromosomes and
    genome.Variant("chr16", 1135446, "G", "T"),   # functional regions
    genome.Variant("chr1", 100000, "C", "G")
]

output_types = [
    AlphaGenomeOutputType.RNA_SEQ,
    AlphaGenomeOutputType.CAGE,
    AlphaGenomeOutputType.ATAC,
    AlphaGenomeOutputType.DNASE,
    AlphaGenomeOutputType.CHIP_HISTONE,
    AlphaGenomeOutputType.SPLICE_SITES
]
```

**Process**:
- `predictor.batch_predict_variant_scores()` 호출
- 각 변이에 대해 선택된 6개 OutputType + 자동 _ACTIVE 변종 적용
- 687개 ontology curie 조합으로 확장
- 결과를 Pandas DataFrame으로 변환 및 통계 분석

**Output**:
```
results/batch_variant_scoring/
├── variant_scores.csv           # 30MB, 121,550 rows
├── high_impact_variants.csv     # 1.1MB, 4,509 rows (|score| > 0.01)
└── variant_scores_summary.json  # 변이별 통계 + 메타데이터
```

#### CLI Variant Scoring

**Input**:
```python
variant = genome.Variant("chr22", 36201698, "A", "C")  # APOL4 gene region

# ALL 19 RECOMMENDED_VARIANT_SCORERS
scorers = [
    # 11 base OutputTypes
    variant_scorers.AtacScorer(),
    variant_scorers.CageScorer(),
    variant_scorers.ChipHistoneScorer(),
    variant_scorers.ChipTfScorer(),
    variant_scorers.ContactMapsScorer(),
    variant_scorers.DnaseScorer(),
    variant_scorers.ProcapScorer(),
    variant_scorers.RnaSeqScorer(),
    variant_scorers.SpliceSitesScorer(),
    variant_scorers.SpliceSiteUsageScorer(),
    variant_scorers.SpliceJunctionsScorer(),

    # 1 scorer-only (no OutputType)
    variant_scorers.PolyadenylationScorer(),

    # 7 _ACTIVE variants
    variant_scorers.AtacActiveScorer(),
    variant_scorers.CageActiveScorer(),
    variant_scorers.ChipHistoneActiveScorer(),
    variant_scorers.ChipTfActiveScorer(),
    variant_scorers.DnaseActiveScorer(),
    variant_scorers.ProcapActiveScorer(),
    variant_scorers.RnaSeqActiveScorer()
]
```

**Process**:
- 각 scorer를 순차적으로 `predictor.predict_variant_scores()` 호출
- 결과를 통합 DataFrame으로 수집
- Ontology 필터링 (선택적), 시각화 플롯 생성

**Output**:
```
results/variant_scoring_cli/
├── variant_scores.csv      # 10MB, 38,357 rows
├── variant_summary.json    # 통계 + scorer별 분포
├── plot_rna_seq.png        # RNA-seq 점수 분포
└── plot_dnase.png          # DNase 점수 분포
```

---

### Batch Variant Scoring (121,550 rows)

#### 설정 및 실행

**스크립트**: `scripts/run_batch_variant_scoring.py`
**소스 튜토리얼**: `tutorials/batch_variant_scoring.ipynb`
**런타임**: ~7 seconds

#### 입력 변이 상세

| Variant ID | Chr | Position | REF | ALT | 유전체 컨텍스트 |
|-----------|-----|----------|-----|-----|----------------|
| chr3_58394738_A_T_b38 | chr3 | 58394738 | A | T | Non-coding region |
| chr8_28520_G_C_b38 | chr8 | 28520 | G | C | Promoter-proximal |
| chr16_636337_G_A_b38 | chr16 | 636337 | G | A | Regulatory element |
| chr16_1135446_G_T_b38 | chr16 | 1135446 | G | T | Enhancer region |
| chr1_100000_C_G_b38 | chr1 | 100000 | C | G | Intergenic |

#### 변이별 결과 요약

| Variant | 점수 행 수 | Mean | Std | Min | Max | Unique Ontology |
|---------|-----------|------|-----|-----|-----|----------------|
| chr3:58394738:A>T | 12,430 | +0.00044 | 0.0049 | -0.0483 | +0.0694 | 687 |
| chr8:28520:G>C | 10,450 | +0.000029 | 0.0017 | -0.0198 | +0.0215 | 687 |
| chr16:636337:G>A | 40,150 | +0.0034 | 0.0185 | -0.1256 | +0.3135 | 687 |
| chr16:1135446:G>T | 41,734 | -0.00084 | 0.0103 | -0.6535 | +0.0821 | 687 |
| chr1:100000:C>G | 16,786 | +0.00068 | 0.0045 | -0.0389 | +0.0567 | 687 |
| **Total** | **121,550** | | | | | |

**핵심 관찰**:
- 두 chr16 변이가 가장 많은 행 생성 (각 40K+) → 조절 요소가 풍부한 영역
- chr16:636337:G>A의 std(0.0185)가 가장 높음 → 세포 타입 간 효과 변동성 큼
- 모든 변이가 동일한 687개 ontology 조합 커버 (6 OutputTypes × 다양한 세포 타입)

#### High-Impact Variants 필터링

| 임계값 | 행 수 | 비율 | 의미 |
|-------|------|------|------|
| \|score\| > 0.01 | 4,509 | 3.7% | 기능적 유의미 후보 |
| \|score\| > 0.05 | 1,423 | 1.2% | 강한 조절 효과 |
| \|score\| > 0.1 | 546 | 0.4% | 병인성 가능성 높음 |

**최고 점수 변이**:
- **Positive**: chr16:636337:G>A / Neutrophil / CAGE → **+0.3135** (99.5th percentile)
- **Negative**: chr16:1135446:G>T / CD14+ monocyte / ChIP-Histone → **-0.6535** (-99.7th percentile)

**생물학적 해석**:
- 두 chr16 변이가 모두 골수성 세포(neutrophil, monocyte)에 강한 영향
- CAGE(5' 전사 시작점)와 ChIP-Histone(히스톤 변형)에서 극단값 → 전사 개시와 크로마틴 상태 조절 요소
- 면역 세포 분화와 관련된 regulatory element일 가능성

#### Assay 타입별 분포 (High Impact |score| > 0.01, 총 4,509행)

- ChIP-Histone: 1,718 high impact scores (38%)
- CAGE: 962 scores (21%)
- DNase: 807 scores (18%)
- RNA-seq: 684 scores (15%)
- ATAC: 338 scores (7%)

#### 점수 방향성 (High Impact 기준)

- **양성 효과**: 3,079 (68.3%) → 대체 염기가 활성 증가
- **음성 효과**: 1,430 (31.7%) → 대체 염기가 활성 감소

**해석**: 대부분 고영향 변이가 조절 활성을 증가시키는 경향

#### Biosample 타입별 변동성

| Biosample 타입 | High Impact 수 | Std Raw Score | 특징 |
|----------------|---------------|---------------|------|
| Primary cells | 1,300 | 0.074 | 가장 높은 변동성 (생물학적 다양성) |
| Tissues | 1,681 | 0.055 | 중간 변동성 (세포 혼합) |
| Cell lines | 1,246 | 0.052 | 낮은 변동성 (균질 배양) |
| In vitro differentiated | 282 | 0.039 | 가장 낮은 변동성 (통제된 환경) |

---

### CLI Variant Scoring (38,357 rows)

#### 설정 및 실행

**스크립트**: `scripts/run_variant_scoring_cli.py`
**재구현 기반**: `tutorials/variant_scoring_ui.ipynb` (Colab UI)
**주요 개선점**: UI 인터랙션 제거, 전체 scorer 자동화, 결과 저장

#### CLI 사용법 예시

```bash
# 기본 실행 (모든 scorer, 모든 ontology)
python scripts/run_variant_scoring_cli.py \
    --chr chr22 --pos 36201698 --ref A --alt C

# 특정 OutputType + Ontology 필터링
python scripts/run_variant_scoring_cli.py \
    --chr chr22 --pos 36201698 --ref A --alt C \
    --outputs rna_seq,dnase \
    --ontology UBERON:0001157 \
    --sequence-length 1MB \
    --interval-width 32768

# 시각화 생성
python scripts/run_variant_scoring_cli.py \
    --chr chr22 --pos 36201698 --ref A --alt C \
    --plot rna_seq,chip_tf
```

#### 입력 변이

**Variant**: chr22:36201698:A>C (GRCh38)
**유전체 위치**: APOL4 gene region (apolipoprotein L4)
**기능적 컨텍스트**: 지질 대사 관련 유전자, 아프리카계 집단에서 APOL1 변이와 신장 질환 연관

#### 19개 Scorer 구성

##### 11개 Base OutputTypes
```python
ATAC                  # Open chromatin (transposase)
CAGE                  # 5' transcript starts
CHIP_HISTONE          # Histone modifications
CHIP_TF               # Transcription factor binding
CONTACT_MAPS          # 3D chromatin interactions
DNASE                 # Open chromatin (DNase I)
PROCAP                # Nascent RNA (human only)
RNA_SEQ               # Gene expression
SPLICE_SITES          # Donor/acceptor sites
SPLICE_SITE_USAGE     # Alternative splicing
SPLICE_JUNCTIONS      # Junction reads
```

##### 1개 Scorer-Only (OutputType 없음)
```python
POLYADENYLATION       # Poly(A) site strength
```

##### 7개 _ACTIVE 변종
```python
ATAC_ACTIVE           # Binary activity threshold
CAGE_ACTIVE
CHIP_HISTONE_ACTIVE
CHIP_TF_ACTIVE
DNASE_ACTIVE
PROCAP_ACTIVE
RNA_SEQ_ACTIVE
```

#### 통합 결과 통계

| 통계량 | 값 | 설명 |
|-------|-----|------|
| Total rows | 38,357 | 모든 scorer + ontology 조합 |
| Unique scorers | 19 | RECOMMENDED_VARIANT_SCORERS 전체 |
| Mean | 98.91 | 높은 값 (CONTACT_MAPS 스케일 때문) |
| Median | 0.00106 | 실제 중심 경향 (대부분 scorer는 [-1, 1]) |
| Std | 581.63 | Scorer 간 스케일 차이로 인한 분산 |
| Min | -6.876 | 최대 억제 효과 |
| Max | 22,558.5 | CONTACT_MAPS 최대값 (log-scale Hi-C) |

**해석 주의사항**:
- Mean(98.91)은 CONTACT_MAPS 때문에 높게 왜곡됨 (수천 단위 값)
- Median(0.00106)이 더 representative (대부분 scorer는 ±0.1 범위)
- Scorer별 통계를 따로 분석해야 정확한 해석 가능

#### Scorer별 행 수 (variant_scores.csv 검증)

| Scorer | 행 수 |
|--------|------|
| GeneMaskActiveScorer(RNA_SEQ) | 14,652 |
| GeneMaskLFCScorer(RNA_SEQ) | 14,652 |
| CenterMaskScorer(CHIP_TF, ACTIVE_SUM) | 1,617 |
| CenterMaskScorer(CHIP_TF, DIFF_LOG2_SUM) | 1,617 |
| CenterMaskScorer(CHIP_HISTONE, ACTIVE_SUM) | 1,116 |
| CenterMaskScorer(CHIP_HISTONE, DIFF_LOG2_SUM) | 1,116 |
| SpliceJunctionScorer() | 734 |
| CenterMaskScorer(CAGE, ACTIVE_SUM) | 546 |
| CenterMaskScorer(CAGE, DIFF_LOG2_SUM) | 546 |
| PolyadenylationScorer() | 396 |
| GeneMaskSplicingScorer(SPLICE_SITE_USAGE) | 367 |
| CenterMaskScorer(DNASE, ACTIVE_SUM) | 305 |
| CenterMaskScorer(DNASE, DIFF_LOG2_SUM) | 305 |
| CenterMaskScorer(ATAC, ACTIVE_SUM) | 167 |
| CenterMaskScorer(ATAC, DIFF_LOG2_SUM) | 167 |
| ContactMapScorer() | 28 |
| CenterMaskScorer(PROCAP, ACTIVE_SUM) | 12 |
| CenterMaskScorer(PROCAP, DIFF_LOG2_SUM) | 12 |
| GeneMaskSplicingScorer(SPLICE_SITES) | 2 |
| **합계** | **38,357** |

**APOL4 변이 해석**:
- RNA-seq 양성 점수 → 유전자 발현 증가 가능성
- CHIP_TF 음성 점수 → 특정 TF 결합 감소 (조절 메커니즘 변화)
- 전반적으로 미세한 효과 (대부분 |score| < 0.01) → 중립 변이 가능성

---

### POLYADENYLATION Scorer (396 rows)

#### 특수성 및 중요성

**핵심 사실**:
- `variant_scorers.PolyadenylationScorer()` - 19개 RECOMMENDED_VARIANT_SCORERS 중 하나
- **OutputType 없음** (scorer로만 존재) → `output_type` 필드는 `RNA_SEQ`로 표시되지만 독립적 scorer
- 3' UTR poly(A) site 강도 예측 → mRNA 안정성, 국소화, 번역 효율 영향

#### APOL4 변이 결과 (chr22:36201698:A>C)

| 통계량 | 값 | 백분위수 | 해석 |
|-------|-----|----------|------|
| Total rows | 396 | - | 다양한 세포 타입 커버 |
| Mean | 0.417 | 99.8% | 거의 모든 세포 타입에서 강한 효과 |
| Min | 0.004 | 78.8% | 가장 약한 효과도 평균 이상 |
| Max | 0.866 | 99.9%+ | 특정 세포 타입에서 극단적 효과 |

**세포 타입 분포**:
```
Tissue:          188 (47.5%)
Primary cells:    90 (22.7%)
Cell lines:       73 (18.4%)
In vitro:         29 (7.3%)
```

#### 생물학적 의미

**핵심 발견**:
- 평균 quantile 99.8% → APOL4의 poly(A) site에 **매우 강한 영향**
- 최소값도 78.8th percentile → 모든 세포 타입에서 유의미
- APOL4는 지질 결합 단백질 → poly(A) site 변화가 mRNA 안정성에 영향 → 단백질 발현량 변화 → 지질 대사 이상 가능성

**임상적 함의**:
- APOL1 변이는 아프리카계 집단 신장 질환과 강한 연관
- APOL4는 동일 유전자 클러스터 → poly(A) site 변이가 APOL1 효과 modulate 가능
- mRNA 안정성 변화 → 단백질 발현 시공간 패턴 변경 → 신장 세포 기능 이상

---

### Batch Tutorial Results (5.4)

#### 튜토리얼 출력

**파일**: `tutorials/batch_variant_scoring.ipynb`
**섹션**: "5.4 변이 효과 점수 (Batch Variant Scoring)"

**주요 코드 블록**:
```python
# Batch prediction for 5 variants
variant_scores = predictor.batch_predict_variant_scores(
    variants=variants,
    output_types=output_types,
    interval_width=32_768
)

# Convert to DataFrame
df = variant_scorers.to_dataframe(variant_scores)

# High-impact filtering
high_impact = df[df['score'].abs() > 0.01]
```

**튜토리얼 특징**:
- Interactive cell 실행으로 즉시 결과 확인
- 변이별 점수 분포 시각화 (histogram, box plot)
- Ontology 메타데이터 탐색 (세포 타입, 조직 타입)
- High-impact 변이 필터링 예제

---

### CLI Tutorial Results (5.7)

#### 튜토리얼 출력

**파일**: `tutorials/variant_scoring_ui.ipynb`
**섹션**: "5.7 변이 효과 점수 (CLI Variant Scoring)"
**원래 형태**: Colab UI widget 기반 (ipywidgets)

**UI 컴포넌트**:
```python
# Original UI-based input
variant_input = widgets.Text(value="chr22:36201698:A>C")
scorer_select = widgets.SelectMultiple(options=SCORER_LIST)
ontology_filter = widgets.Text(value="")

# Run button callback
def run_scoring(b):
    results = predictor.predict_variant_scores(...)
    display(results_table)
```

**CLI 재구현 이유**:
- UI는 수동 인터랙션 필요 → 재현성 및 자동화 어려움
- 19개 scorer 전체 실행 시 UI가 느리고 불안정
- 결과 자동 저장 및 배치 처리 불가
- → `run_variant_scoring_cli.py`로 완전 자동화

---

## 4. 왜 이 방법인가

### Batch vs CLI: 상호보완적 워크플로

| 기준 | Batch | CLI | 선택 기준 |
|------|-------|-----|----------|
| **변이 수** | 다수 (5~1000+) | 단일 | 스크리닝 vs 심층 분석 |
| **Scorer 수** | 제한적 (6~10) | 전체 (19) | 속도 vs 포괄성 |
| **실행 시간** | 빠름 (~7초) | 느림 (~2분) | 대규모 vs 정밀 분석 |
| **용도** | GWAS 변이 우선순위 | 후보 변이 메커니즘 규명 | |

**Best Practice**:
1. Batch로 1000개 변이 스크리닝 → High-impact 50개 식별
2. CLI로 50개를 19개 scorer로 심층 분석 → 메커니즘 이해
3. 실험 검증은 최종 5개만 수행

### POLYADENYLATION Scorer의 독특한 위치

**왜 OutputType이 없는가?**:
- Poly(A) site는 RNA-seq 데이터로 측정되지만, 별도의 생물학적 프로세스
- RNA-seq은 전체 전사체 abundance, POLYADENYLATION은 3' end processing
- AlphaGenome은 두 신호를 독립적으로 모델링 → 별도 scorer 필요

**임상적 중요성**:
- 3'UTR 변이는 coding 영향 없지만 mRNA 안정성 크게 변경 가능
- APA(Alternative Polyadenylation)는 암, 신경퇴행성 질환에서 dysregulated
- GWAS에서 간과되기 쉬운 메커니즘 → POLYADENYLATION scorer가 gap 보완

### 기존 도구 대비 장점

| 도구 | 커버리지 | 세포 타입 수 | 한계 |
|------|----------|-------------|------|
| **AlphaGenome** | 11 assays + 8 scorers | 687 ontology | API 의존, 비용 |
| CADD | Coding 중심 | - | Non-coding 약함 |
| GWAVA | Non-coding | Limited | Assay-specific 아님 |
| DeepSEA | 919 features | ~100 | 2015년 데이터 |
| Enformer | 5313 tracks | ~200 | 학습 데이터만 |

**AlphaGenome 차별점**:
- 최신 genomic assay 데이터 반영 (2023년까지)
- 세포 타입별 예측 (687 ontology curies)
- 1Mbp 컨텍스트 → 장거리 enhancer 효과 포착

---

## 5. 해석과 한계

### 점수 스케일 이해

#### Log-scale vs Linear-scale Scorers

**Log-scale scorers** (RNA_SEQ, CAGE, ChIP 계열):
- 범위: 대부분 [-0.2, +0.2]
- 해석: 0.1 = ~1.1x fold-change (log2 scale)
- 생물학적 의미: 0.05 이상이면 기능적으로 유의미

**Linear-scale scorers** (CONTACT_MAPS, SPLICE_SITE_USAGE):
- 범위: [0, 수천]
- 해석: 절대값이 아닌 percentile로 평가
- 예: CONTACT_MAPS 1000 → 99th percentile → 강한 interaction 변화

**CRITICAL LESSON**: CLI CSV의 mean=98.91은 무의미. Scorer별 분석 필수.

### Quantile 기반 해석 프레임워크

| Quantile 범위 | 해석 | 실험 우선순위 | 예시 변이 |
|--------------|------|--------------|----------|
| > 99.5% | 극단적 효과 | 최우선 검증 | chr16:636337 CAGE |
| 95-99.5% | 강한 효과 | 높음 | POLYADENYLATION mean |
| 90-95% | 중간 효과 | 중간 | 대부분 high-impact |
| < 90% | 약한 효과 | 낮음 | Neutral 변이 |

### 세포 타입 특이성 패턴

**Primary cells의 높은 변동성 (std=0.074) 원인**:
1. **생물학적 다양성**: Donor 간 유전적 배경 차이
2. **분화 상태**: 동일 세포 타입도 분화 단계 다름
3. **실험적 변동**: 배양 조건, 분리 프로토콜 차이

**임상적 함의**:
- 환자별 변이 효과 예측 시 primary cell 데이터 우선 고려
- Cell line 데이터는 메커니즘 이해용, 개인 맞춤 예측엔 부적합

### 주요 한계

#### 1. 단일 변이 가정 (Haplotype 무시)
- **문제**: 실제 genome은 다수 변이가 동시 존재 (phased haplotype)
- **영향**: 복합 효과(epistasis) 포착 불가
- **예**: rs123(+0.05) + rs456(-0.04) ≠ rs123/456(+0.01) (비선형 상호작용)

#### 2. 구조적 변이 미지원
- Deletion, insertion, inversion, CNV는 평가 불가
- Indel은 길이 제한 (수십 bp까지만)

#### 3. 참조 게놈 한정
- GRCh38 중심 (T2T-CHM13 미지원)
- 개인 게놈 참조 서열 차이 무시

#### 4. 예측 vs 실측 불일치
- AlphaGenome 예측이 항상 실험 결과와 일치하지 않음
- 검증 데이터셋에서 pearson r=0.6~0.8 (assay 타입별 차이)

#### 5. 비용 및 속도
- API 호출당 과금 → 대규모 스크리닝 비용 증가
- 1000 variants × 19 scorers = 19K API 호출 → 수 시간 소요

---

## 6. 다음으로 (The Bridge)

### 현재까지 달성

✅ **Batch**: 5 variants, 6 scorers → 121,550 rows → High-impact 4,509 variants 식별
✅ **CLI**: 1 variant, 19 scorers → 38,357 rows → POLYADENYLATION 포함 전체 scorer 검증
✅ **POLYADENYLATION**: 396 rows, mean quantile 99.8% → APOL4 poly(A) site 강한 영향 확인

### 남은 질문

**"단일 변이가 아니라 연속 구간 내 모든 가능한 변이를 체계적으로 스캔하려면?"**

- chr16:636337 주변 ±50bp 영역의 모든 A>T, C>G, ... 조합 탐색
- 전사 인자 binding motif의 critical residue 식별
- 변이 민감도 지도(variant sensitivity map) 생성

**AlphaGenome의 해결책: ISM (In-Silico Mutagenesis)**

- 지정 구간(64bp or 256bp)의 모든 위치를 4가지 염기로 체계적 변이
- 64bp × 3 mutations = 192 variant scores (or 256bp → 768 scores)
- Heatmap으로 시각화 → "어느 위치가 가장 중요한가?" 한눈에 파악

---

### 다음 모듈 예고

**[모듈 08: ISM 분석 — 연속 변이 민감도 스캔](08-ism-analysis.md)**

> "단일 변이 점수(121,550 + 38,357)에서 전체 구간 스캔(192/768 행렬)으로. ISM은 조절 요소의 critical nucleotide를 체계적으로 매핑하여, 실험 설계와 병인 변이 해석의 정밀도를 한 단계 끌어올린다."

**다룰 내용**:
- ISM 64bp vs 256bp 전략 비교
- `ism.ism_matrix()` 활용법 (variant_scores= 파라미터 주의)
- Heatmap 해석 패턴 (hotspot, cold region)
- chr3:58394738 주변 ISM 결과 (192 scores)
