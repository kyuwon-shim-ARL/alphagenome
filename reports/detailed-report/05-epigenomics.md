# 모듈 05: 에피게놈 (Epigenomics — CHIP_HISTONE · CHIP_TF)

## 1. 한 줄 요약

> CHIP_HISTONE(1,116 tracks)과 CHIP_TF(1,617 tracks)는 가장 많은 human track을 보유한 OutputType으로, CenterMaskScorer(width=2001)로 히스톤 변형과 전사인자 결합을 정량화하며, K562와 HepG2 세포주에서 845 tracks를 분석했다.

## 2. 왜 이 분석이 필요했나

유전자 발현(RNA_SEQ, CAGE)과 크로마틴 접근성(DNASE, ATAC)만으로는 부족했다. **누가 스위치를 조절하는가?**

- **CHIP_HISTONE**: 히스톤 변형은 크로마틴 상태의 **기억 시스템**이다. H3K4ME3는 "이 위치는 프로모터"라는 주석을 달고, H3K27ME3는 "여기는 침묵"이라는 잠금장치를 건다.
- **CHIP_TF**: 전사인자 결합은 **문맥별 실행 명령**이다. 같은 유전자라도 CTCF, POLR2A, EP300가 어디에 붙느냐에 따라 발현 결과가 달라진다.

DNASE/ATAC이 "문이 열렸다"를 알려준다면, CHIP_HISTONE은 "이 방은 언제 열려야 하는지", CHIP_TF는 "누가 열고 있는지"를 알려준다.

### 왜 이 둘이 가장 많은 track을 가지나?

| OutputType | Human Tracks | 이유 |
|-----------|--------------|-----|
| **CHIP_TF** | 1,617 | 수백 종의 TF x 수백 종의 세포주 조합 |
| **CHIP_HISTONE** | 1,116 | 6-10개 히스톤 마커 x 다양한 조직 |
| RNA_SEQ | 667 | 세포주별 1-2개 track |
| DNASE | 305 | 세포주별 1개 track |

CHIP_TF는 **조합 폭발**(combinatorial explosion)의 결과다. 1,447개 track이 EFO 기반 세포주에서 나온다는 것은, ENCODE 프로젝트가 세포주별로 수십 가지 TF ChIP-seq을 체계적으로 수행했음을 의미한다.

## 3. 분석 과정 (The Mechanics)

### Quick Reference

**CHIP_HISTONE**:
```python
from alphagenome.models import dna_client, variant_scorers

output = dna_model.predict_interval(
    interval=interval,
    requested_outputs={dna_client.OutputType.CHIP_HISTONE},
    ontology_terms=['CL:0001059'],  # CD34+ common myeloid progenitor
)
# output.chip_histone.values.shape: (1048576, N)
# N: 해당 ontology의 ChIP-histone track 수

# Variant Scoring
scorer = variant_scorers.CenterMaskScorer(
    requested_output=dna_client.OutputType.CHIP_HISTONE,
    width=2001,  # DNASE의 501보다 넓음 (histone은 넓은 영역에 분포)
    aggregation_type=variant_scorers.AggregationType.DIFF_MEAN,
)
variant_scores = dna_model.score_variant(
    interval=interval,
    variant=variant,
    variant_scorers=[scorer],
)
```

**CHIP_TF**:
```python
# TF 필터링 함수
def filter_to_tfs(chip_tf_data, tf_names):
    """특정 TF만 선택 (대소문자 무관)"""
    metadata = chip_tf_data.metadata
    tf_mask = metadata['transcription_factor'].str.upper().isin(
        [tf.upper() for tf in tf_names]
    )
    filtered_indices = metadata.index[tf_mask]
    return chip_tf_data.select_tracks_by_index(filtered_indices)

# K562 CHIP_TF 예측
output_k562 = dna_model.predict_interval(
    interval=interval,
    requested_outputs={dna_client.OutputType.CHIP_TF},
    ontology_terms=['EFO:0002067'],  # K562 chronic myelogenous leukemia
)

# CTCF만 필터링
k562_ctcf = filter_to_tfs(output_k562.chip_tf, ['CTCF'])
# CTCF + RAD21 co-localization
k562_ctcf_rad21 = filter_to_tfs(output_k562.chip_tf, ['CTCF', 'RAD21'])
```

### IPO 요약

| 단계 | Input | Process | Output |
|-----|-------|---------|--------|
| **CHIP_HISTONE** | interval, ontology | CenterMaskScorer(width=2001) | (1048576, N) array |
| **CHIP_TF** | interval, ontology, TF names | filter_to_tfs() | (1048576, M) filtered array |
| **Variant Scoring** | variant, scorer | score_variant() | AnnData (genes x tracks) |

### CHIP_HISTONE (1,116 tracks)

| 항목 | 내용 |
|------|------|
| **설명** | ChIP-seq histone modification markers |
| **Resolution** | 1 bp |
| **Human Tracks** | 1,116 |
| **Mouse Tracks** | 183 |
| **Strand** | Unstranded (.) |
| **Scorer** | CenterMaskScorer (width=2001) |
| **사용된 분석** | Viz Tour, Batch, Analysis Workflow |

**히스톤 마커별 색상 코딩 (Viz Tour)**:
```python
histone_to_color = {
    'H3K27AC':  '#e41a1c',  # Red - Active enhancer
    'H3K36ME3': '#ff7f00',  # Orange - Gene body (actively transcribed)
    'H3K4ME1':  '#377eb8',  # Blue - Enhancer (poised/active)
    'H3K4ME3':  '#984ea3',  # Purple - Active promoter
    'H3K9AC':   '#4daf4a',  # Green - Active chromatin
    'H3K27ME3': '#ffc0cb',  # Pink - Repressive (Polycomb silencing)
}
```

**생물학적 의미**:
- **H3K4ME3** (Purple): TSS 근처에 집중. "프로모터 활성화" 신호.
- **H3K27AC** (Red): Active enhancer. "이 enhancer가 지금 작동 중".
- **H3K36ME3** (Orange): Gene body에 분포. "전사 진행 중".
- **H3K27ME3** (Pink): Polycomb 단백질이 붙는 위치. "장기 침묵".

**Ontology 분포**:
| Ontology Prefix | Track 수 | 주요 타입 |
|-----------------|---------|----------|
| CL | 323 | Primary cells |
| EFO | 389 | Cell lines (실험적) |
| UBERON | 364 | Tissues |
| CLO | 5 | Cell line ontology |
| NTR | 35 | ENCODE custom terms |

**CenterMaskScorer width=2001 vs DNASE width=501**:
- DNASE는 **point source** (TF binding site는 ~10bp 영역)
- Histone modification은 **broad peak** (nucleosome 단위, ~147bp + 주변 확산)
- 2001bp window는 변이 중심 ±1000bp, 약 6-7개 nucleosome 포함

**Batch Variant Scoring 결과**:
```
chr16:636337:G>A / CD14+ monocyte / ChIP-histone: -0.1265 (99.4th percentile negative)
```
- CD14+ monocyte (단핵구)에서 histone modification이 강하게 감소
- 이는 myeloid 계열 세포에서 regulatory element 비활성화를 시사
- chr16:1135446:G>T의 -0.6535 (극단값)와 동일 myeloid 세포 계열

**결과 파일**:
- `results/visualization_tour/05_chip_histone.png` (566 KB)
- `results/batch_variant_scoring/variant_scores.csv`에 CHIP_HISTONE 포함

### CHIP_TF (1,617 tracks)

| 항목 | 내용 |
|------|------|
| **설명** | ChIP-seq transcription factor binding profiles |
| **Resolution** | 1 bp |
| **Human Tracks** | 1,617 |
| **Mouse Tracks** | 127 |
| **Strand** | Unstranded (.) |
| **Scorer** | CenterMaskScorer |
| **사용된 분석** | ChIP-TF Analysis script (6.3) |

**Why 1,617 tracks?**
- 가장 많은 human track (CHIP_HISTONE의 1.45배)
- 1,447개가 EFO 기반 세포주 (전체의 89%)
- ENCODE 프로젝트: 세포주별 × TF별 조합 ChIP-seq

**Ontology 분포**:
| Ontology Prefix | Track 수 | 설명 |
|-----------------|---------|------|
| EFO | 1,447 | Cell lines (실험용 세포주) |
| CL | 63 | Primary cells |
| UBERON | 100 | Tissues |
| CLO | 4 | Cell line ontology |
| NTR | 3 | ENCODE custom |

**EFO 우위의 이유**:
- TF ChIP-seq은 **실험적 조작이 가능한 세포주**에서 주로 수행
- Primary cell은 배양이 어렵고, tissue는 세포 이질성 문제
- EFO:0002067 (K562), EFO:0001187 (HepG2) 등 ENCODE 핵심 세포주

### ChIP-TF Analysis Script Details (6.3)

**Script**: `scripts/run_chip_tf_analysis.py`

**4개 분석 항목**:
1. **K562 CTCF Binding Profile** (EFO:0002067)
2. **HepG2 CTCF Binding Profile** (EFO:0001187)
3. **CTCF-RAD21 Co-localization** (cohesin complex)
4. **Multi-TF Comparison** (CTCF, RAD21, POLR2A, EP300)

**핵심 함수**: `filter_to_tfs()`
```python
def filter_to_tfs(chip_tf_data, tf_names):
    """
    TrackData에서 특정 TF만 선택.

    Args:
        chip_tf_data: TrackData (CHIP_TF)
        tf_names: list of str (TF names, 대소문자 무관)

    Returns:
        TrackData: Filtered tracks
    """
    metadata = chip_tf_data.metadata
    tf_mask = metadata['transcription_factor'].str.upper().isin(
        [tf.upper() for tf in tf_names]
    )
    filtered_indices = metadata.index[tf_mask]
    return chip_tf_data.select_tracks_by_index(filtered_indices)
```

**분석 대상 세포주**:
| 세포주 | Ontology Code | 설명 | Available TFs |
|--------|---------------|------|---------------|
| **K562** | EFO:0002067 | Chronic myelogenous leukemia | CTCF, RAD21, POLR2A, EP300 |
| **HepG2** | EFO:0001187 | Hepatocellular carcinoma | CTCF |

**CRITICAL: HepG2 Ontology**
- **CORRECT**: `EFO:0001187`
- **WRONG**: `EFO:0002052` (이것은 다른 cell line)
- Script에서 명시적으로 `EFO:0001187` 사용

**분석 1: K562 CTCF Binding**
```python
interval = genome.Interval('chr11', 5_220_000, 5_230_000)
interval = interval.resize(dna_client.SEQUENCE_LENGTH_1MB)

output_k562 = dna_model.predict_interval(
    interval=interval,
    requested_outputs={dna_client.OutputType.CHIP_TF},
    ontology_terms=['EFO:0002067'],  # K562
)

k562_ctcf = filter_to_tfs(output_k562.chip_tf, ['CTCF'])
# k562_ctcf.values.shape: (1048576, N_ctcf)
```
- CTCF: 크로마틴 구조 조직화 (insulator, loop anchor)
- 결과: `results/chip_tf_analysis/chip_tf_k562_ctcf.png` (70 KB)

**분석 2: HepG2 CTCF Binding**
```python
output_hepg2 = dna_model.predict_interval(
    interval=interval,
    requested_outputs={dna_client.OutputType.CHIP_TF},
    ontology_terms=['EFO:0001187'],  # HepG2
)

hepg2_ctcf = filter_to_tfs(output_hepg2.chip_tf, ['CTCF'])
```
- 간암 세포주의 CTCF 결합 패턴
- K562와 비교하여 세포 타입별 CTCF 특이성 확인
- 결과: `results/chip_tf_analysis/chip_tf_hepg2_ctcf.png` (94 KB)

**분석 3: CTCF-RAD21 Co-localization**
```python
k562_ctcf_rad21 = filter_to_tfs(output_k562.chip_tf, ['CTCF', 'RAD21'])
# CTCF + RAD21 동시 필터링
```
- **Cohesin complex**: RAD21은 cohesin의 핵심 구성요소
- CTCF와 cohesin은 TAD (Topologically Associating Domain) 경계 형성
- 두 TF의 결합 위치가 겹치면 → 크로마틴 loop anchor 가능성
- 결과: `results/chip_tf_analysis/chip_tf_ctcf_rad21_coloc.png` (120 KB)

**분석 4: Multi-TF Comparison**
```python
k562_multi_tf = filter_to_tfs(output_k562.chip_tf,
                               ['CTCF', 'RAD21', 'POLR2A', 'EP300'])
```
- **POLR2A** (RNA Polymerase II): 전사 중심
- **EP300**: 히스톤 아세틸전이효소, active enhancer marker
- 4개 TF 동시 비교로 regulatory landscape 전체 파악
- 결과: `results/chip_tf_analysis/chip_tf_multi_tf.png` (182 KB)

**K562 + HepG2 통계**:
- **총 tracks**: 845 across 2 cell types
- K562: 4개 TF 분석 (CTCF, RAD21, POLR2A, EP300)
- HepG2: 1개 TF 분석 (CTCF)
- 4개 plot 모두 성공적으로 생성

**생물학적 통찰**:
1. **CTCF**: K562와 HepG2에서 공통/차이점 비교 가능
2. **CTCF-RAD21**: 크로마틴 3D 구조의 anchor point 추론
3. **POLR2A + EP300**: 활성 전사 중심(POLR2A)과 enhancer(EP300) 관계

## 4. 왜 이 방법인가

### CenterMaskScorer width=2001의 근거

**CHIP_HISTONE에서 width=2001을 사용하는 이유**:

1. **Nucleosome 단위**: 147bp DNA + linker ~50bp = ~200bp/nucleosome
2. **Histone modification 확산**: 한 nucleosome의 변형은 인접 nucleosome에 전파 (chromatin spreading)
3. **2001bp = ±1000bp**: 변이 중심 기준 양쪽 5-7개 nucleosome 포함
4. **Broad peak 특성**: DNASE(narrow peak, width=501)와 달리 histone은 넓은 영역에 분포

**비교**:
| Scorer | Width | 이유 | OutputType |
|--------|-------|------|-----------|
| CenterMaskScorer | 501 | TF binding site (~10bp 영역) | DNASE, ATAC |
| CenterMaskScorer | 2001 | Histone domain (~1kb 영역) | CHIP_HISTONE |

### filter_to_tfs() 함수의 필요성

**문제**: CHIP_TF는 1,617 tracks. 한 번에 모두 시각화하면 해석 불가능.

**해결책**: TF 이름으로 필터링
```python
metadata['transcription_factor'].str.upper().isin([tf.upper() for tf in tf_names])
```
- `.str.upper()`: 대소문자 무관 매칭 (CTCF, ctcf, Ctcf 모두 동일)
- `.isin()`: 여러 TF 동시 선택 (CTCF + RAD21)

**왜 metadata 기반?**
- TrackData.metadata는 pandas DataFrame
- 각 track의 'transcription_factor' 컬럼에 TF 이름 저장
- Boolean indexing으로 원하는 track만 선택

### K562 vs HepG2 선택의 의미

**K562 (EFO:0002067)**:
- Chronic myelogenous leukemia (만성 골수성 백혈병)
- Hematopoietic lineage (조혈 계열)
- ENCODE Tier 1 cell line (가장 많이 연구됨)

**HepG2 (EFO:0001187)**:
- Hepatocellular carcinoma (간암)
- Liver-specific regulatory elements
- Drug metabolism 연구에 중요

**왜 둘을 비교?**
- 서로 다른 lineage (myeloid vs hepatic)
- CTCF 결합 위치의 **공통점** = 구조적 필수 요소 (insulator)
- CTCF 결합 위치의 **차이점** = 세포 타입 특이적 크로마틴 구조

## 5. 해석과 한계

### CHIP_HISTONE 해석

**Batch 변이 chr16:636337:G>A 사례**:
```
CD14+ monocyte ChIP-histone: -0.1265 (99.4th percentile negative)
```

**의미**:
- CD14+ monocyte (단핵구, myeloid 계열)에서 히스톤 변형 감소
- 99.4th percentile negative = 게놈 전체에서 하위 0.6%에 해당하는 강한 감소
- 동일 변이가 CAGE에서 +0.3135 (99.5th percentile positive)

**역설**: CAGE 증가 + ChIP-histone 감소?
- CAGE는 **전사 시작점 활성** (promoter activity)
- ChIP-histone은 **넓은 chromatin domain** 상태
- Possible scenario: Promoter는 활성화되었으나 주변 chromatin은 열리지 않음 → **aberrant activation** (비정상 활성화)

### CHIP_TF 해석

**CTCF-RAD21 Co-localization의 의미**:
- CTCF alone: Insulator, loop anchor candidate
- RAD21 alone: Cohesin 복합체, sister chromatid cohesion
- **CTCF + RAD21 together**: TAD boundary, 크로마틴 loop의 anchor point

**Multi-TF (4개 TF) 패턴**:
1. **POLR2A peak** alone: 전사 중심 (promoter/gene body)
2. **EP300 peak** alone: Enhancer (POLR2A와 떨어진 위치)
3. **EP300 + POLR2A overlap**: Active enhancer-promoter interaction
4. **CTCF + RAD21 + (EP300 or POLR2A)**: Loop anchor가 regulatory element와 접촉

### 한계

**1. Static snapshot**:
- ChIP-seq은 특정 시점의 결합 상태
- 동적 변화 (stimulation 후 변화) 포착 불가
- AlphaGenome 예측도 동일 한계

**2. TF 조합 폭발**:
- 인간 genome: ~1,600개 TF
- CHIP_TF: 1,617 tracks이지만 일부 TF만 커버
- 많은 tissue-specific TF는 데이터 부족

**3. Indirect measurement**:
- ChIP-seq은 **단백질-DNA 결합** 측정
- 기능적 효과 (transcriptional output) 별도 분석 필요
- AlphaGenome은 ChIP signal 예측, 기능 예측은 RNA_SEQ 등과 통합 필요

**4. Histone code의 복잡성**:
- 6개 마커만으로 chromatin state 전체를 표현하기 어려움
- H3K9ME2/3, H4K20ME3 등 다른 마커 미포함
- Bivalent promoter (H3K4ME3 + H3K27ME3) 같은 복합 상태 해석 어려움

**5. K562/HepG2의 대표성**:
- 2개 세포주로 모든 조직 대표 불가
- Primary tissue의 in vivo regulatory landscape와 차이

## 6. 다음으로 (The Bridge)

지금까지:
- **03-gene-expression**: RNA_SEQ(667), CAGE(546), PROCAP(12) → 유전자 발현 예측
- **04-chromatin**: DNASE(305), ATAC(167) → 크로마틴 접근성
- **05-epigenomics**: CHIP_HISTONE(1,116), CHIP_TF(1,617) → 에피게놈 상태

### 현재 위치 (The Status)

11개 OutputType 중 7개 커버:
- ✅ RNA_SEQ, CAGE, PROCAP (gene expression)
- ✅ DNASE, ATAC (chromatin accessibility)
- ✅ CHIP_HISTONE, CHIP_TF (epigenomics)

### 남은 4개 OutputType

| OutputType | Tracks | 내용 |
|-----------|--------|------|
| SPLICE_SITES | 4 | Splice site 확률 (donor/acceptor) |
| SPLICE_SITE_USAGE | 734 | 조직별 splice site 사용 빈도 |
| SPLICE_JUNCTIONS | 367 | Exon-exon junction 예측 |
| CONTACT_MAPS | 28 | 3D chromatin interaction (Hi-C) |

### 다음 단계: RNA Processing & 3D Structure

**왜 splicing이 중요한가?**
- Human genome: 20,000 genes → 100,000+ isoforms
- Alternative splicing: 단백질 다양성의 핵심
- Splicing variant: 질병 원인의 ~15%

**왜 CONTACT_MAPS가 중요한가?**
- CTCF/RAD21 결합 위치를 알았음 (CHIP_TF)
- 하지만 **실제 loop 형성 여부**는 CONTACT_MAPS에서 확인
- 3D 구조 = 1D 서열 + epigenomic marks의 최종 결과

> "히스톤 변형(1,116 tracks)과 전사인자 결합(1,617 tracks)으로 에피게놈 상태를 파악했다. 이로써 유전자 발현(03), 크로마틴 접근성(04), 에피게놈(05) 측면을 모두 검증했다. 남은 것은 **RNA 처리(splicing)와 3D 크로마틴 구조**이다." → [다음: 06-splicing-and-3d.md](06-splicing-and-3d.md)
