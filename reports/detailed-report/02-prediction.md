# 모듈 02: 예측 파이프라인 (Prediction Pipeline)

## 1. 한 줄 요약

> AlphaGenome은 predict_sequence → predict_interval → predict_variant → score_variant의 4단계 파이프라인으로 DNA 서열 예측부터 변이 효과 정량화까지 수행하며, Quick Start에서 37 genes x 667 tracks (14,652 tidy scores)를 생성했다.

---

## 2. 왜 이 분석이 필요했나

- **이전 모듈에서 남은 질문**: 환경 설정과 입력 시스템(Interval, Variant, Ontology)을 이해했다. 이 입력을 실제로 어떻게 예측 파이프라인에 전달하고 결과를 받는가?
- **이 모듈이 답하려는 것**: 4가지 핵심 예측 메서드의 시그니처, 반환값 구조, 실행 흐름을 검증한다.

---

## 3. 분석 과정 (The Mechanics)

### Quick Reference

AlphaGenome의 예측 파이프라인은 추상화 수준에 따라 4단계로 나뉜다.

```
Level 0: predict_sequence(sequence: str)
         → ModelOutput (raw predictions for 1MB string)

Level 1: predict_interval(interval: genome.Interval)
         → ModelOutput (auto-fetch reference genome)

Level 2: predict_variant(interval, variant)
         → VariantOutput (reference + alternate predictions)

Level 3: score_variant(interval, variant, variant_scorers)
         → list[AnnData] (aggregated gene-level scores)
```

**핵심 설계 원칙**:
- Level 0-1: 단일 서열 예측 (변이 비교 없음)
- Level 2: 변이 비교 (사용자가 직접 차이 계산)
- Level 3: 자동 스코어링 (quantile normalization + gene masking)

### IPO 요약

| 단계 | Input | Process | Output | 핵심 수치 |
|------|-------|---------|--------|----------|
| predict_sequence | DNA string (1MB) | 서열 인코딩 → 예측 | ModelOutput | (1048576, N) |
| predict_interval | genome.Interval | 참조 게놈 조회 → 예측 | ModelOutput | (1048576, N) |
| predict_variant | Interval + Variant | REF/ALT 각각 예측 | VariantOutput | REF + ALT 쌍 |
| score_variant | Interval + Variant + Scorers | 유전자 수준 집계 | list[AnnData] | (37, 667) |

---

### predict_sequence()

가장 원시적인 예측 메서드로, 1MB DNA 문자열을 직접 입력받아 예측을 수행한다.

**시그니처:**

```python
output = dna_model.predict_sequence(
    sequence: str,                           # DNA 서열 (정확한 길이 필요)
    requested_outputs: list[OutputType],     # 요청할 출력 타입
    ontology_terms: list[str] = None,        # 조직/세포 ontology 코드
    organism: Organism = HOMO_SAPIENS,       # 생물종
) -> ModelOutput
```

**반환 타입: ModelOutput**
- 각 OutputType에 대응하는 attribute (예: `output.dnase`, `output.cage`)
- 각 attribute는 `TrackData` 객체: `.values` (numpy array), `.metadata` (DataFrame)

**Quick Start 실행 예시:**

```python
# 단일 출력 (DNase for Lung)
output = dna_model.predict_sequence(
    sequence='GATTACA'.center(dna_client.SEQUENCE_LENGTH_1MB, 'N'),
    requested_outputs=[dna_client.OutputType.DNASE],
    ontology_terms=['UBERON:0002048'],  # Lung
)
# output.dnase.values.shape = (1048576, 1)

# 복수 출력 (CAGE + DNase for Lung + Brain)
output = dna_model.predict_sequence(
    sequence='GATTACA'.center(dna_client.SEQUENCE_LENGTH_1MB, 'N'),
    requested_outputs=[dna_client.OutputType.CAGE, dna_client.OutputType.DNASE],
    ontology_terms=['UBERON:0002048', 'UBERON:0000955'],
)
# output.dnase.values.shape = (1048576, 2)  -- Lung + Brain
# output.cage.values.shape  = (1048576, 4)  -- Lung(+/-) + Brain(+/-)
```

CAGE의 shape가 (1048576, 4)인 이유: CAGE는 strand-specific이므로 각 ontology term에 대해 `+`/`-` strand 각각 1개 track이 생성된다. 2 ontology x 2 strands = 4 tracks.

---

### predict_interval()

게놈 좌표(Interval)를 입력으로 받아, 참조 게놈에서 자동으로 서열을 조회한 후 예측을 수행한다.

**시그니처:**

```python
output = dna_model.predict_interval(
    interval: genome.Interval,               # 게놈 구간
    requested_outputs: set[OutputType],      # 요청할 출력 타입
    ontology_terms: list[str] = None,        # 조직/세포 ontology 코드
    organism: Organism = HOMO_SAPIENS,       # 생물종
) -> ModelOutput
```

**Quick Start 실행 예시 (CYP2B6 유전자 RNA-seq):**

```python
from alphagenome.data import gene_annotation, genome

# GENCODE v46 유전자 주석 로드
gtf = pd.read_feather(
    'https://storage.googleapis.com/alphagenome/reference/gencode/'
    'hg38/gencode.v46.annotation.gtf.gz.feather'
)

# CYP2B6 유전자 interval 조회
interval = gene_annotation.get_gene_interval(gtf, gene_symbol='CYP2B6')
# 결과: chr19:40991281-41018398:+

# 1MB로 resize (모델 입력 크기)
interval = interval.resize(dna_client.SEQUENCE_LENGTH_1MB)

# RNA-seq 예측 (간 우엽)
output = dna_model.predict_interval(
    interval=interval,
    requested_outputs=[dna_client.OutputType.RNA_SEQ],
    ontology_terms=['UBERON:0001114'],  # Right liver lobe
)
# output.rna_seq.values.shape = (1048576, 3)
```

**Visualization Tour에서 사용된 구간:**

```python
# chr22 APOL 유전자 클러스터 주변
interval = genome.Interval('chr22', 36_150_498, 36_252_898).resize(
    dna_client.SEQUENCE_LENGTH_1MB
)
# 결과: chr22:35677410-36725986:.
```

---

### predict_variant()

특정 변이의 효과를 예측한다. 참조(reference)와 대체(alternate) 서열 각각에 대해 예측을 수행하고, 두 결과를 `VariantOutput`으로 반환한다.

**시그니처:**

```python
variant_output = dna_model.predict_variant(
    interval: genome.Interval,               # 게놈 구간
    variant: genome.Variant,                 # 변이
    requested_outputs: set[OutputType],      # 요청할 출력 타입
    ontology_terms: list[str] = None,        # 조직/세포 ontology 코드
    organism: Organism = HOMO_SAPIENS,       # 생물종
) -> VariantOutput
```

**반환 타입: VariantOutput**
- `.reference`: ModelOutput (참조 서열 예측 결과)
- `.alternate`: ModelOutput (대체 서열 예측 결과)

**Quick Start 실행 예시:**

```python
variant = genome.Variant(
    chromosome='chr22',
    position=36201698,
    reference_bases='A',
    alternate_bases='C',
)

interval = variant.reference_interval.resize(dna_client.SEQUENCE_LENGTH_1MB)

variant_output = dna_model.predict_variant(
    interval=interval,
    variant=variant,
    requested_outputs=[dna_client.OutputType.RNA_SEQ],
    ontology_terms=['UBERON:0001157'],  # Colon
)

# 결과
variant_output.reference.rna_seq.values.shape  # (1048576, 3)
variant_output.alternate.rna_seq.values.shape  # (1048576, 3)
```

**Analysis Workflow에서 사용 (TAL1 분석):**

```python
tal1_interval = genome.Interval(
    chromosome='chr1', start=47209255, end=47242023, strand='-'
)

output = dna_model.predict_variant(
    interval=tal1_interval.resize(2**20),  # 1MB
    variant=variant,
    requested_outputs={
        dna_client.OutputType.RNA_SEQ,
        dna_client.OutputType.CHIP_HISTONE,
        dna_client.OutputType.DNASE,
    },
    ontology_terms=['CL:0001059'],  # CD34+ common myeloid progenitor
)

# 차이 계산
diff = (output.alternate.rna_seq.filter_to_nonpositive_strand()
        - output.reference.rna_seq.filter_to_nonpositive_strand())
```

---

### score_variant()

변이의 효과를 사전 정의된 scorer를 사용하여 정량적으로 평가한다. 결과는 `AnnData` 객체 리스트로 반환된다.

**시그니처:**

```python
variant_scores = dna_model.score_variant(
    interval: genome.Interval,                # 게놈 구간
    variant: genome.Variant,                  # 변이
    variant_scorers: list[VariantScorer],     # scorer 리스트
    organism: Organism = HOMO_SAPIENS,        # 생물종
) -> list[AnnData]
```

**AnnData 구조:**

| 필드 | 내용 | 예시 |
|------|------|------|
| `.X` | score matrix | (37, 667) -- 37 genes x 667 tracks |
| `.obs` | gene metadata | gene_name, gene_id, strand |
| `.var` | track metadata | ontology_curie, biosample_name |
| `.uns` | variant/scorer metadata | variant, scorer |
| `.layers['quantiles']` | quantile-normalized scores | (37, 667) |

### Scoring Pipeline Overview

Variant scoring은 두 가지 수준에서 이루어진다:

**Level 1: predict_variant (raw predictions)**
- REF/ALT 각각의 raw 예측값 반환
- 사용자가 직접 차이를 계산하고 해석
- 시각화에 적합 (REF vs ALT overlay)

**Level 2: score_variant (aggregated scores)**
- 사전 정의된 scorer가 자동으로 REF-ALT 차이 계산
- Gene-level, center-window 등 다양한 집계 방식
- Quantile normalization으로 게놈 전체 분포 대비 상대적 위치 제공
- AnnData → tidy DataFrame 변환으로 분석/필터링에 적합

**19 RECOMMENDED_VARIANT_SCORERS:**

```python
from alphagenome.models import variant_scorers

all_scorers = variant_scorers.RECOMMENDED_VARIANT_SCORERS
# 19개 키:
#   11개 base: ATAC, CAGE, CHIP_HISTONE, CHIP_TF, CONTACT_MAPS,
#              DNASE, PROCAP, RNA_SEQ, SPLICE_SITES, SPLICE_SITE_USAGE, SPLICE_JUNCTIONS
#   1개 scorer-only: POLYADENYLATION
#   7개 _ACTIVE: ATAC_ACTIVE, CAGE_ACTIVE, CHIP_HISTONE_ACTIVE, CHIP_TF_ACTIVE,
#                DNASE_ACTIVE, PROCAP_ACTIVE, RNA_SEQ_ACTIVE
```

**Scorer 타입별 분류:**

| Scorer 클래스 | 알고리즘 | 사용 OutputType | 주요 파라미터 |
|---------------|---------|----------------|--------------|
| `GeneMaskLFCScorer` | Gene-level log-fold-change | RNA_SEQ, CAGE, PROCAP | `requested_output` |
| `CenterMaskScorer` | Center-window 비교 | DNASE, ATAC, CHIP_HISTONE, CHIP_TF | `requested_output`, `width`, `aggregation_type` |
| `SpliceSitesScorer` | Splice site 확률 변화 | SPLICE_SITES | - |

**tidy_scores() 변환:**

```python
tidy_df = variant_scorers.tidy_scores(
    variant_scores,
    match_gene_strand=True
)
# 결과: (14652, 19) DataFrame
# 컬럼: variant_id, gene_name, ontology_curie, raw_score, quantile_score, ...
```

**Quick Start 결과:**

```python
variant_scorer = variant_scorers.RECOMMENDED_VARIANT_SCORERS['RNA_SEQ']
variant_scores = dna_model.score_variant(
    interval=interval,
    variant=variant,
    variant_scorers=[variant_scorer],
)

variant_scores[0].X.shape  # (37, 667) -- 37 genes x 667 tracks
variant_scores[0].obs['gene_name'].tolist()[:10]
# ['RBFOX2', 'APOL4', 'APOL1', 'MYH9', 'TXN2', 'FOXRED2', 'EIF3D', 'APOL3', 'APOL5', 'APOL2']

tidy = variant_scorers.tidy_scores([variant_scores[0]], match_gene_strand=True)
# Tidy Scores Shape: (14652, 19)
```

---

### Batch Processing Patterns

대량 변이 분석 시 사용되는 배치 처리 패턴이다.

#### Pattern 1: 순차 반복 (run_batch_variant_scoring.py)

```python
results = []
for i, vcf_row in tqdm(vcf.iterrows(), total=len(vcf), desc="Scoring variants"):
    variant = genome.Variant(
        chromosome=str(vcf_row.CHROM),
        position=int(vcf_row.POS),
        reference_bases=vcf_row.REF,
        alternate_bases=vcf_row.ALT,
        name=vcf_row.variant_id,
    )
    interval = variant.reference_interval.resize(sequence_length)

    variant_scores = dna_model.score_variant(
        interval=interval, variant=variant,
        variant_scorers=selected_scorers, organism=organism,
    )
    results.append(variant_scores)

df_scores = variant_scorers.tidy_scores(results)
```

#### Pattern 2: 병렬 처리 (run_analysis_workflow.py)

```python
scores = dna_model.score_variants(
    intervals=eval_df['interval'].to_list(),
    variants=eval_df['variant'].to_list(),
    variant_scorers=[variant_scorers.RECOMMENDED_VARIANT_SCORERS['RNA_SEQ']],
    max_workers=2,
)
```

#### Pattern 3: 전체 Scorer 사용 (CLI)

```python
# 19개 scorer 모두 사용
variant_scores = dna_model.score_variant(
    interval=interval, variant=variant,
    variant_scorers=list(variant_scorers.RECOMMENDED_VARIANT_SCORERS.values()),
    organism=organism,
)
# CLI: 38,357 score rows from 19 scorers
```

**Scorer Coverage 비교:**

| Script | Scorers | Rows | 커버리지 |
|--------|---------|------|---------|
| Batch (`run_batch_variant_scoring.py`) | 6 | 121,550 | RNA_SEQ, CAGE, ATAC, DNASE, CHIP_HISTONE, SPLICE_SITES |
| CLI (`run_variant_scoring_cli.py`) | 19 | 38,357 | **전체** 11 + POLYADENYLATION + 7 _ACTIVE |
| Analysis (`run_analysis_workflow.py`) | 1 | varies | RNA_SEQ only |

---

### Quick Start Tutorial Results (5.1)

**소스**: `tutorials/quick_start.ipynb`
**스크립트**: `results/quick_start/run_quick_start.py`

**9단계 실행 흐름:**

1. **모델 생성**: `dna_client.create(api_key)` -- Success
2. **OutputType 열거**: 11가지 OutputType 확인
3. **predict_sequence (단일)**: DNASE for Lung -- shape (1048576, 1)
4. **predict_sequence (복수)**: CAGE+DNASE for Lung+Brain -- DNASE (1048576, 2), CAGE (1048576, 4)
5. **predict_interval**: CYP2B6 RNA-seq in liver -- chr19:40991281-41018398:+, (1048576, 3)
6. **predict_variant**: chr22:36201698:A>C RNA-seq in colon -- REF (1048576, 3), ALT (1048576, 3)
7. **score_variant**: RNA_SEQ scorer -- (37, 667), Tidy (14652, 19)
8. **ISM**: 64bp window CenterMaskScorer DNASE -- 192 variants
9. **Mouse prediction**: DNase for Lung -- (1048576, 1)

**생성 파일:**

| 파일 | 설명 |
|------|------|
| `results.json` | 전체 실행 결과 메타데이터 |
| `cyp2b6_rna_seq.png` | CYP2B6 RNA-seq 시각화 |
| `variant_effect.png` | REF vs ALT RNA-seq overlay |
| `variant_scores.csv` | Tidy scores (14,652 rows) |
| `run_quick_start.py` | 실행 스크립트 |

---

## 4. 왜 이 방법인가

| 고려한 방법 | 채택 여부 | 이유 |
|------------|----------|------|
| 계층적 API (4 levels) | ✅ 채택 | "간단한 것은 쉽게, 복잡한 것은 가능하게" |
| 단일 통합 메서드 | ❌ | 유연성 부족 |
| AnnData 반환 | ✅ 채택 | scanpy/Seurat 생태계 호환 |
| Custom DataFrame | ❌ | 표준화 부족 |
| tidy_scores 변환 | ✅ 채택 | 필터링/분석에 적합한 long format |

---

## 5. 해석과 한계

**강점:**
- 4단계 파이프라인이 단순 예측부터 유전자 수준 스코어링까지 커버
- AnnData 구조로 scanpy 생태계 연동 가능
- tidy_scores()로 표준적 데이터 분석 워크플로우 지원
- 19개 RECOMMENDED_VARIANT_SCORERS로 모든 OutputType 커버

**한계:**
1. **API 기반**: 모델 내부 접근 불가 (attention weights, embeddings)
2. **Rate Limit**: 대량 분석 시 병렬 처리 제한 (max_workers=2 권장)
3. **Fixed model version**: 재현성 관리에 주의 필요
4. **Quantile normalization**: track 간 직접 비교 시 스케일 차이에 주의

---

## 6. 다음으로 (The Bridge)

> 예측 파이프라인의 6개 핵심 메서드를 검증했고, Quick Start에서 37 genes x 667 tracks의 실제 결과를 확인했다. 이제 **유전자 발현을 측정하는 RNA_SEQ, CAGE, PROCAP의 상세 결과**를 살펴볼 차례이다.

→ [다음: 03-gene-expression.md](03-gene-expression.md)
