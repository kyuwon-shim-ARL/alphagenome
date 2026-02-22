# 모듈 08: In-Silico Mutagenesis (ISM — 64bp · 256bp)

## 1. 한 줄 요약

> ISM은 연속 구간 내 모든 위치에서 3가지 대체 염기를 체계적으로 스캔하며, 64bp(192 variants)와 256bp(768 variants, 9.96초)를 실행하여 DNASE 접근성 변화를 SeqLogo heatmap으로 시각화했다.

---

## 2. 왜 이 분석이 필요했나

- **이전 모듈에서 남은 질문**: Batch(121,550행)와 CLI(38,357행)로 단일 변이의 효과를 정량화했다. 그런데 단일 변이가 아니라 연속 구간 전체를 어떻게 스캔하는가?
- **이 모듈이 답하려는 것**: ISM으로 연속 구간 내 모든 가능한 단일 염기 변이를 체계적으로 스캔하여 position-specific sensitivity를 매핑한다.

**ISM의 핵심 가치:**
- **체계적 스캔**: 연속 구간 내 모든 위치에서 3가지 대체 염기를 시뮬레이션
- **Motif discovery**: 전사인자 결합 부위, 보존된 regulatory element 식별
- **Position-specific sensitivity**: 각 염기 위치의 functional importance 정량화

---

## 3. 분석 과정 (The Mechanics)

### Quick Reference

| 항목 | 64bp (Quick Start) | 256bp (Extended) |
|------|-------------------|-----------------|
| **Sequence interval** | chr20:3753000-3753400 (→16KB) | chr20:3745008-3761392:. |
| **ISM interval** | 64bp (centered) | chr20:3753072-3753328:. |
| **Variants** | 192 (64 × 3) | 768 (256 × 3) |
| **Scorer** | CenterMaskScorer(DNASE, 501, DIFF_MEAN) | 동일 |
| **소요 시간** | ~3초 | 9.96초 |
| **스크립트** | tutorials/quick_start.ipynb | scripts/run_ism_256bp.py |

### IPO 요약

| 단계 | 내용 | 핵심 수치 | 파일 경로 |
|------|------|----------|----------|
| **Input** | interval (16KB) + ism_interval (64/256bp) + scorer | 192 / 768 variants | - |
| **Process** | `score_ism_variants()` → `ism.ism_matrix()` | (4 × width) matrix | - |
| **Output** | SeqLogo heatmap + CSV | mean=-0.0072 (256bp) | `results/ism_256bp/` |

---

### ISM 원리와 API

ISM = In-Silico Mutagenesis: 연속 구간의 모든 위치에서 가능한 대체 염기를 시뮬레이션하여 regulatory impact를 매핑하는 기법.

**Variant 생성 공식:**
```
Variant count = ism_interval.width × 3
```
- 각 위치에서 원본 제외 3가지 대체 염기 (예: A→C/G/T)
- 64bp → 192 variants, 256bp → 768 variants

**API 시그니처:**

```python
ism_scores = dna_model.score_ism_variants(
    interval: genome.Interval,        # 전체 서열 구간 (16KB)
    ism_interval: genome.Interval,    # ISM 스캔 구간 (64/256bp)
    variant_scorers: list[VariantScorer],
)
```

**반환 타입:** ISM 변이 결과 리스트 (각 항목: 변이별 score tuple)

**16KB vs 1MB context:**
- ISM: 16KB 사용 (`SEQUENCE_LENGTH_16KB`) — 다중 variant이므로 효율성 우선
- Single variant scoring: 1MB 사용 — 최대 regulatory context
- Trade-off: ISM은 속도를 위해 약간의 context를 희생

---

### ISM Matrix 변환

```python
from alphagenome.interpretation import ism

ism_matrix_data = ism.ism_matrix(
    variant_scores=[extract_mean_score(score_tuple) for score_tuple in ism_scores],
    variants=[score_tuple[0].uns['variant'] for score_tuple in ism_scores],
)
# 결과: numpy array (4 × width) — A, C, G, T rows
```

**주의사항 (CRITICAL):**
- ✅ `ism.ism_matrix(variant_scores=..., variants=...)` — 올바른 파라미터명
- ❌ `ism.ism_matrix(scores=...)` — 에러 발생
- ❌ `variant_scorers.tidy_scores()` — ISM 결과에 사용 불가

---

### Quick Start ISM (64bp, 192 variants)

| 항목 | 값 |
|------|-----|
| 구간 | chr20:3753000-3753400 (16KB로 resize) |
| ISM window | 64 bp |
| 변이 수 | 192 (64 × 3) |
| Scorer | CenterMaskScorer(DNASE, width=501, DIFF_MEAN) |

```python
sequence_interval = genome.Interval('chr20', 3_753_000, 3_753_400)
sequence_interval = sequence_interval.resize(dna_client.SEQUENCE_LENGTH_16KB)
ism_interval = sequence_interval.resize(64)

dnase_variant_scorer = variant_scorers.CenterMaskScorer(
    requested_output=dna_client.OutputType.DNASE,
    width=501,
    aggregation_type=variant_scorers.AggregationType.DIFF_MEAN,
)

ism_scores = dna_model.score_ism_variants(
    interval=sequence_interval,
    ism_interval=ism_interval,
    variant_scorers=[dnase_variant_scorer],
)
# ISM Variants Scored: 192
```

---

### ISM 256bp Extended Analysis (768 variants)

**스크립트:** `scripts/run_ism_256bp.py`

| 항목 | 값 |
|------|-----|
| Sequence interval | chr20:3745008-3761392:. |
| ISM interval | chr20:3753072-3753328:. |
| ISM window | 256 bp |
| 변이 수 | 768 (256 × 3) |
| 소요 시간 | 9.96초 |
| Scorer | CenterMaskScorer(DNASE, width=501, DIFF_MEAN) |

**구현 상세:**

```python
# 서열 구간 정의 (동일 위치)
sequence_interval = genome.Interval('chr20', 3_753_000, 3_753_400)
sequence_interval = sequence_interval.resize(dna_client.SEQUENCE_LENGTH_16KB)
# 결과: chr20:3745008-3761392:.

# ISM 구간 (64bp → 256bp 확장)
ism_interval = sequence_interval.resize(256)
# 결과: chr20:3753072-3753328:.

# ISM 실행
ism_scores = dna_model.score_ism_variants(
    interval=sequence_interval,
    ism_interval=ism_interval,
    variant_scorers=[dnase_variant_scorer],
)
# 768 variants scored in 9.96s
```

**Score 통계:**

| 통계량 | 값 |
|--------|-----|
| Mean | **-0.0072** |
| Std | 0.038 |
| Min | -0.193 |
| Max | 0.151 |
| Median | -0.00011 |

**해석:**
- **Negative mean (-0.0072)**: 이 구간의 대부분 변이는 DNASE 접근성을 약간 감소시킴
- **Min -0.193**: 특정 위치는 강한 repressive effect — 중요 motif 파괴?
- **Max +0.151**: 일부 위치는 accessibility 증가 — repressor motif 제거?
- **Narrow distribution (std=0.038)**: 대부분 변이는 미약한 영향

---

### ISM Matrix Visualization

**SeqLogo heatmap 생성:**

```python
from alphagenome.interpretation import plot_components

plot_components.plot(
    [plot_components.SeqLogo(
        scores=ism_matrix_data,
        scores_interval=ism_interval,
        ylabel='ISM DNASE (mean)',
    )],
    interval=ism_interval,
    fig_width=35,
)
```

**시각화 특징:**
- X축: Genomic position (256bp)
- Y축: 4개 염기 (A, C, G, T)
- Color intensity: Score magnitude
- Pattern: 특정 motif에서 clustered high/low scores → regulatory element

**결과 파일:**

| 파일 | 설명 |
|------|------|
| `results/ism_256bp/ism_scores.csv` | 768 rows (position, ref, alt, score) |
| `results/ism_256bp/ism_heatmap.png` | SeqLogo heatmap |
| `results/ism_256bp/results.json` | 메타데이터 + 통계 |

---

## 4. 왜 이 방법인가

| 고려한 방법 | 채택 여부 | 이유 |
|------------|----------|------|
| ISM (single-nucleotide scan) | ✅ 채택 | 체계적, motif-agnostic, API batch 효율 |
| Random sampling | ❌ | 불완전한 커버리지 |
| Deep mutational scanning | ❌ | 코돈 수준, 단백질 중심 |
| DNASE scorer for ISM | ✅ 채택 | 305 tracks, 범용적 regulatory region 탐지 |

**ISM window 크기 선택 가이드:**
- 64bp: 단일 TFBS (전사인자 결합 부위)
- 128-256bp: Enhancer/promoter core element (권장)
- >512bp: 비용/시간 급증, context 희석

---

## 5. 해석과 한계

**강점:**
1. **Systematic coverage**: 모든 가능한 single-nucleotide 변화를 평가
2. **Motif-agnostic**: 사전 지식 없이 regulatory element 발견
3. **API 효율성**: 768 variants를 10초 내 처리 (batch API 활용)
4. **SeqLogo 시각화**: Position-specific sensitivity를 직관적으로 표현

**한계:**
1. **Single-nucleotide only**: Indel, structural variant 미지원
2. **Limited context**: 16KB (vs 1MB for single variant scoring)
3. **Haplotype ignorance**: Compound heterozygote 효과 미반영
4. **규모 제한**: 256bp ISM (768 variants)은 관리 가능하지만, 1000bp+ 은 비용/시간 급증
5. **In silico 한정**: 예측은 "potential" impact → wet-lab validation 필수

---

## 6. 다음으로 (The Bridge)

> ISM으로 64bp(192)와 256bp(768) 연속 구간을 체계적으로 스캔하여 DNASE 접근성 변화를 매핑했다. 이론적 도구를 모두 검증했으니, **실제 질병 연구에 어떻게 적용되는가?** T-ALL TAL1 oncogenic variant 분석으로 실전 사례를 확인한다.

→ [다음: 09-tall-workflow.md](09-tall-workflow.md)
