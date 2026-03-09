# 모듈 06: 스플라이싱과 3D 구조 (Splicing & 3D — SPLICE_* · CONTACT_MAPS)

## 1. 한 줄 요약

> SPLICE_SITES(4 고정 tracks), SPLICE_SITE_USAGE(734), SPLICE_JUNCTIONS(367)로 RNA 처리를, CONTACT_MAPS(28)로 3D 크로마틴 구조를 예측하여, 11개 OutputType 전체 커버리지를 완성했다.

## 2. 왜 이 분석이 필요했나

모듈 03-05에서 gene expression(RNA_SEQ, CAGE, PROCAP), chromatin accessibility(DNASE, ATAC), epigenomics(CHIP_HISTONE, CHIP_TF)를 검증했다. **하지만 AlphaGenome의 11개 OutputType 중 4개가 남아 있다**: SPLICE_SITES, SPLICE_SITE_USAGE, SPLICE_JUNCTIONS, CONTACT_MAPS.

이 4개는 **RNA processing(splicing)과 3D chromatin structure**라는 두 가지 critical biological layer를 다룬다:

1. **Splicing**: pre-mRNA에서 intron을 제거하고 exon을 연결하는 과정. Alternative splicing은 단일 유전자가 여러 protein isoform을 생성하게 하며, 조직 특이적으로 발생한다.
2. **3D Structure**: DNA는 선형 서열이지만 핵 안에서 3차원 구조를 형성한다. Hi-C 실험은 멀리 떨어진 genomic region 간의 상호작용(loop, TAD)을 측정한다.

**이 모듈의 목적**:
- **SPLICE_SITES**: donor/acceptor splice site의 확률(4 fixed tracks, ontology-independent)
- **SPLICE_SITE_USAGE**: 조직 특이적 splice site usage frequency(734 human tracks)
- **SPLICE_JUNCTIONS**: exon-exon connection 예측(Sashimi plot으로 시각화)
- **CONTACT_MAPS**: 3D chromatin interaction map(Hi-C equivalent)

이 4개를 검증하면 **AlphaGenome이 커버하는 11개 OutputType 전체를 확인**할 수 있다.

## 3. 분석 과정 (The Mechanics)

### Quick Reference

| OutputType | Tracks | Script | Result | Visualization |
|-----------|--------|--------|---------|---------------|
| SPLICE_SITES | 4(H), 4(M) | `run_splice_sites.py` | `results/splice_sites/` | `06_splice.png` |
| SPLICE_SITE_USAGE | 734(H), 180(M) | `run_splice_site_usage.py` | `results/splice_site_usage/` | Standalone script |
| SPLICE_JUNCTIONS | 367(H), 90(M) | Viz Tour | `visualization_tour/` | Sashimi in `06_splice.png` |
| CONTACT_MAPS | 28(H), 8(M) | Viz Tour | `visualization_tour/` | `07_contact_maps.png` |

### IPO 요약

**Input (공통)**:
```python
# Interval 정의
interval = genome.Interval(
    chromosome='chr22',
    start=38_000_000,
    end=38_060_000,
    strand='+'  # SPLICE_* only
)

# API 호출
output = api.predict_tracks(interval)
```

**Process**:
- SPLICE_SITES: `output.splice_sites` → 4 tracks(donor+/-, acceptor+/-), 1bp resolution
- SPLICE_SITE_USAGE: `output.splice_site_usage` → 734 tissue-specific tracks
- SPLICE_JUNCTIONS: `output.splice_junctions` → Sashimi plot (arc diagram)
- CONTACT_MAPS: `output.contact_maps` → 2D heatmap (Hi-C matrix)

**Output**:
- SPLICE_SITES: probability [0,1] per base pair
- SPLICE_SITE_USAGE: frequency [0,1] per base pair
- SPLICE_JUNCTIONS: junction-level counts (variable resolution)
- CONTACT_MAPS: interaction matrix (symmetric 2D array)

### SPLICE_SITES (4 tracks)

**특징**:
- **4 fixed tracks**: donor+, donor-, acceptor+, acceptor-
- **Ontology-independent**: purely sequence-based prediction
- **Strand-specific**: +/- strand별로 별도 예측
- **Resolution**: 1 bp

**구현**:
```python
# run_splice_sites.py
output = api.predict_tracks(genome.Interval(
    chromosome='chr22',
    start=38_000_000,
    end=38_060_000,
    strand='+'  # Required for SPLICE_SITES
))

# 4 tracks 모두 추출
splice_sites = output.splice_sites
# Shape: [60000, 4] → (donor+, donor-, acceptor+, acceptor-)
```

**Scorer**:
```python
from alphageome import variant_scorers
scorer = variant_scorers.SpliceSitesScorer()
scores = scorer.score(
    api=api,
    variants=[variant1, variant2],
)
```

**검증 결과**:
- `results/splice_sites/splice_sites_predictions.csv`: 240,004 rows(60,000 bp × 4 tracks)
- `results/splice_sites/splice_sites_summary.txt`: statistics
- `results/visualization_tour/06_splice.png`: SPLICE_SITES + SPLICE_JUNCTIONS 시각화

| Item | Value |
|------|-------|
| Description | Splice site probability prediction (donor/acceptor) |
| Resolution | 1 bp |
| Human Tracks | 4 |
| Mouse Tracks | 4 |
| Strand | Strand-specific (+/-) |
| Scorer | SpliceSitesScorer |
| Used in | Viz Tour, Batch |

### SPLICE_SITE_USAGE (734 tracks)

**특징**:
- **Tissue/cell-specific**: SPLICE_SITES와 달리 ontology term에 따라 다름
- **734 human tracks**: CL=248, CLO=2, EFO=134, NTR=16, UBERON=334
- **Strand-specific**: +/- strand별로 별도 tracks
- **RECOMMENDED_VARIANT_SCORERS**에 포함

**구현**:
```python
# run_splice_site_usage.py
output = api.predict_tracks(genome.Interval(
    chromosome='chr22',
    start=38_000_000,
    end=38_060_000,
    strand='+'
))

# 734 tracks (human)
splice_site_usage = output.splice_site_usage
print(f"Tracks: {len(splice_site_usage.tracks)}")  # 734
```

**Metadata 분석**:
```python
# Ontology 분포
for track in splice_site_usage.tracks:
    ontology = track.biosample.ontology_term.term_id.split(':')[0]
    # CL=248, EFO=134, UBERON=334, NTR=16, CLO=2
```

**검증 결과**:
- `results/splice_site_usage/results.json`: 실행 메타데이터 (status: success, gene: APOL4)
- `results/splice_site_usage/splice_site_usage.png`: APOL4 유전자 colon 조직 시각화
- **Standalone visualization script**: `scripts/run_splice_site_usage.py` (5.7 KB)
  - APOL4 gene region, Colon Transverse/Sigmoid ontology terms
  - Output: `results/splice_site_usage/splice_site_usage.png`

| Item | Value |
|------|-------|
| Description | Tissue-specific splice site usage frequency |
| Resolution | 1 bp |
| Human Tracks | 734 |
| Mouse Tracks | 180 |
| Strand | Strand-specific (+/-) |
| Scorer | Included in RECOMMENDED_VARIANT_SCORERS |
| Used in | Metadata analysis, Standalone viz |

**왜 SPLICE_SITES와 다른가**:
- SPLICE_SITES: **sequence-based** → "이 위치가 splice site로 보이는가?"
- SPLICE_SITE_USAGE: **tissue-based** → "이 조직에서 이 splice site를 실제로 얼마나 사용하는가?"

### SPLICE_JUNCTIONS (367 tracks)

**특징**:
- **Junction-level prediction**: intron을 건너뛰고 exon 간 연결을 예측
- **Sashimi plot**: arc diagram으로 시각화(exon-exon connection)
- **367 human tracks**: CL=124, CLO=1, EFO=67, NTR=8, UBERON=167
- **Variable resolution**: junction마다 다름(exon 위치에 따라)

**구현**:
```python
# Viz Tour (Tutorial 03)
output = api.predict_tracks(genome.Interval(
    chromosome='chr22',
    start=36_209_000,
    end=36_214_000,
    strand='-'  # APOL4 gene
))

# Sashimi plot 생성
plot_components.Sashimi(
    output.splice_junctions
    .filter_to_strand('-')
    .filter_by_tissue('Colon_Transverse'),
    ylabel_template='SPLICE_JUNCTIONS: {biosample_name} ({strand})',
)
```

**검증 결과**:
- `results/visualization_tour/06_splice.png`: SPLICE_SITES(bottom) + SPLICE_JUNCTIONS(Sashimi, top)
- APOL4 gene region에서 exon-exon connection을 arc로 표현
- Colon_Transverse 조직의 splice junction 패턴 확인

| Item | Value |
|------|-------|
| Description | Splice junction (exon-exon connection) prediction |
| Resolution | Variable (junction-dependent) |
| Human Tracks | 367 |
| Mouse Tracks | 90 |
| Strand | Strand-specific (+/-) |
| Scorer | Included in RECOMMENDED_VARIANT_SCORERS |
| Used in | Viz Tour |

**Sashimi Plot 해석**:
- **Arc height**: junction usage frequency
- **Arc thickness**: confidence/coverage
- **Multiple arcs**: alternative splicing events

### CONTACT_MAPS (28 tracks)

**특징**:
- **3D chromatin interaction**: Hi-C 실험과 동등한 output
- **28 human tracks**: mostly EFO-based cell lines(27 EFO + 1 CL)
- **2D matrix**: symmetric heatmap(chromosome region × region)
- **No strand**: 3D structure는 strand-independent

**구현**:
```python
# Viz Tour (Tutorial 03)
output = api.predict_tracks(genome.Interval(
    chromosome='chr22',
    start=36_209_000,
    end=36_214_000,
    strand='+'  # Ignored for CONTACT_MAPS
))

# Contact map 시각화
plot_components.ContactMaps(
    tdata=output.contact_maps,
    ylabel_template='{biosample_name}\n{name}',
    cmap='autumn_r',
    vmax=1.0,
)
```

**검증 결과**:
- `results/visualization_tour/07_contact_maps.png`: 519 KB heatmap
- HCT116(EFO:0002824) cell line의 contact map
- 5 kb interval → 2D matrix showing self-interaction patterns

| Item | Value |
|------|-------|
| Description | 3D chromatin interaction (Hi-C) map prediction |
| Resolution | Variable (bin-dependent) |
| Human Tracks | 28 |
| Mouse Tracks | 8 |
| Strand | N/A (2D matrix) |
| Scorer | Included in RECOMMENDED_VARIANT_SCORERS |
| Used in | Viz Tour |

**해석**:
- **Diagonal**: self-interaction(항상 높음)
- **Off-diagonal**: long-range interaction(TAD, loop)
- **Color intensity**: interaction frequency

### 11 OutputType 전체 커버리지 완성

**Category별 Track 수**:

| Category | OutputTypes | Total Tracks | Module |
|----------|-----------|-------------|---------|
| Gene Expression | RNA_SEQ, CAGE, PROCAP | 1,225 | 03 |
| Chromatin Access | DNASE, ATAC | 472 | 04 |
| Epigenomics | CHIP_HISTONE, CHIP_TF | 2,733 | 05 |
| Splicing | SPLICE_SITES, SSU, SJ | 1,105 | 06 |
| 3D Structure | CONTACT_MAPS | 28 | 06 |
| **Total** | **11** | **5,563** | — |

**Splicing Track 계산**:
- SPLICE_SITES: 4(human) + 4(mouse) = 8
- SPLICE_SITE_USAGE: 734(human) + 180(mouse) = 914
- SPLICE_JUNCTIONS: 367(human) + 90(mouse) = 457
- **Total Splicing**: 8 + 914 + 457 = **1,379 tracks**
  - Human only: 4 + 734 + 367 = **1,105 tracks**

**CONTACT_MAPS Track 계산**:
- Human: 28 tracks(27 EFO + 1 CL)
- Mouse: 8 tracks
- **Total**: 28(H) + 8(M) = **36 tracks**
  - Human only: **28 tracks**

## 4. 왜 이 방법인가

**SPLICE_SITES vs SPLICE_SITE_USAGE**:
- SPLICE_SITES: **sequence-only** prediction → fast, ontology-independent
- SPLICE_SITE_USAGE: **context-aware** prediction → tissue-specific, RECOMMENDED_VARIANT_SCORERS에 포함
- 둘 다 필요: SPLICE_SITES는 "potential", SPLICE_SITE_USAGE는 "actual usage"

**SPLICE_JUNCTIONS Sashimi Plot**:
- **Arc diagram**: exon-exon connection을 직관적으로 표현
- **Alternative splicing**: multiple arcs로 isoform diversity 확인
- **Tissue specificity**: 같은 gene도 조직마다 다른 junction pattern

**CONTACT_MAPS 2D Heatmap**:
- **Hi-C equivalent**: 실험 없이 3D structure 예측
- **TAD/loop detection**: off-diagonal patterns → regulatory interaction
- **Cell-line specific**: 28 tracks로 다양한 cell context 커버

**왜 Visualization Tour에서만 검증했나**:
- SPLICE_JUNCTIONS, CONTACT_MAPS는 **visualization-focused** OutputType
- Batch variant scoring에서는 **변이 효과 정량화**에 집중(RNA_SEQ, CHIP 등)
- Viz Tour는 **qualitative validation**(visual inspection)에 적합

## 5. 해석과 한계

**성과**:
1. **11 OutputType 전체 검증 완료**: RNA_SEQ부터 CONTACT_MAPS까지
2. **5,563 human tracks 커버**: 가장 포괄적인 regulatory prediction model
3. **Splicing + 3D structure**: 기존 variant annotation tool(VEP, CADD)이 놓치는 layer
4. **Standalone visualization**: SPLICE_SITE_USAGE script로 custom plot 가능

**한계**:
1. **SPLICE_SITE_USAGE**: 734 tracks but no CLI variant scorer coverage in report
   - **실제로는 CLI CSV에 포함됨**(MEMORY.md: "CLI covers all 19 scorers")
   - 보고서에서 명시적으로 분석하지 않았을 뿐
2. **CONTACT_MAPS**: 28 tracks only(다른 OutputType 대비 적음)
   - Hi-C 실험 데이터 자체가 상대적으로 부족
   - Cell-line specific → 조직 다양성 낮음
3. **Sashimi Plot**: junction-level visualization은 좋지만 quantitative metric 부족
   - Batch variant scoring에서 SPLICE_JUNCTIONS scorer를 사용하지 않음
4. **3D structure interpretation**: TAD/loop detection algorithm은 별도로 필요
   - AlphaGenome은 raw contact matrix만 제공
   - Downstream analysis(TADbit, Juicer 등) 필요

**Critical Insight**:
> **"Splicing은 단순히 intron 제거가 아니라, 조직 특이적 protein isoform 생성의 핵심 메커니즘이다."**
> - 단일 유전자(APOL4)도 Colon_Transverse와 Liver에서 다른 junction pattern
> - Variant scoring에서 SPLICE_SITE_USAGE scorer가 RECOMMENDED에 포함된 이유
> - Alternative splicing은 암에서 특히 중요(oncogene isoform switching)

## 6. 다음으로 (The Bridge)

**지금까지의 여정**:
1. **Module 03**: Gene expression(RNA_SEQ, CAGE, PROCAP) → mRNA 생성량
2. **Module 04**: Chromatin accessibility(DNASE, ATAC) → regulatory element 접근성
3. **Module 05**: Epigenomics(CHIP_HISTONE, CHIP_TF) → histone/TF binding
4. **Module 06**: Splicing + 3D(SPLICE_*, CONTACT_MAPS) → RNA processing + 크로마틴 구조

**11개 OutputType 전체를 검증했다**. 각 OutputType은 regulatory landscape의 특정 aspect를 보여준다. 하지만 **실제 research question은 "변이가 이 모든 기능에 미치는 효과를 어떻게 정량화하는가?"**이다.

**다음 단계(Module 07: Variant Scoring)**:
- **Batch variant scoring**: 5 variants × 687 ontology terms → 121,550 rows
- **CLI variant scoring**: 19 RECOMMENDED_VARIANT_SCORERS → 38,357 rows
- **All 11 OutputTypes covered**: CLI CSV includes EVERY OutputType + POLYADENYLATION
- **Critical lesson**: "CLI covers all 19 scorers" (MEMORY.md)
  - Batch: 6 OutputTypes only
  - CLI: ALL 11 OutputTypes + POLYADENYLATION

> **"11개 OutputType을 개별적으로 이해했으니, 이제 변이가 이 모든 기능에 미치는 효과를 어떻게 종합적으로 정량화하는가? Batch와 CLI variant scoring으로 확인한다."**

[다음: 07-variant-scoring.md](07-variant-scoring.md)
