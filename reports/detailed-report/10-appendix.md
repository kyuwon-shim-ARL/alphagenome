# 모듈 10: 부록 (Appendix — API Reference · Architecture · Files · Ontology · Timeline)

## 1. 한 줄 요약

> AlphaGenome API는 7개 핵심 메서드, gRPC 기반 client-server 아키텍처, 10개 results 디렉토리, 12개 스크립트로 구성되며, ~156회 API 호출을 ~8분 내에 완료했다.

---

## A. API Reference Summary

### A.1 핵심 메서드 시그니처

#### Model Creation

```python
dna_model = dna_client.create(api_key: str) -> DNAModel
```

#### predict_sequence

```python
output = dna_model.predict_sequence(
    sequence: str,
    requested_outputs: list[dna_client.OutputType],
    ontology_terms: list[str] = None,
    organism: dna_client.Organism = dna_client.Organism.HOMO_SAPIENS,
) -> ModelOutput
```

#### predict_interval

```python
output = dna_model.predict_interval(
    interval: genome.Interval,
    requested_outputs: set[dna_client.OutputType] | list[dna_client.OutputType],
    ontology_terms: list[str] = None,
    organism: dna_client.Organism = dna_client.Organism.HOMO_SAPIENS,
) -> ModelOutput
```

#### predict_variant

```python
variant_output = dna_model.predict_variant(
    interval: genome.Interval,
    variant: genome.Variant,
    requested_outputs: set[dna_client.OutputType] | list[dna_client.OutputType],
    ontology_terms: list[str] = None,
    organism: dna_client.Organism = dna_client.Organism.HOMO_SAPIENS,
) -> VariantOutput
```

#### score_variant

```python
variant_scores = dna_model.score_variant(
    interval: genome.Interval,
    variant: genome.Variant,
    variant_scorers: list[VariantScorer],
    organism: dna_client.Organism = dna_client.Organism.HOMO_SAPIENS,
) -> list[AnnData]
```

#### score_variants (batch parallel)

```python
scores = dna_model.score_variants(
    intervals: list[genome.Interval],
    variants: list[genome.Variant],
    variant_scorers: list[VariantScorer],
    max_workers: int = 2,
) -> list[list[AnnData]]
```

#### score_ism_variants

```python
ism_scores = dna_model.score_ism_variants(
    interval: genome.Interval,
    ism_interval: genome.Interval,
    variant_scorers: list[VariantScorer],
) -> list[tuple[AnnData, ...]]
```

#### output_metadata

```python
metadata = dna_model.output_metadata(
    organism: dna_client.Organism = dna_client.Organism.HOMO_SAPIENS,
) -> OutputMetadata
# OutputMetadata.concatenate() -> pd.DataFrame
```

### A.2 Data Structure 시그니처

#### genome.Interval

```python
genome.Interval(
    chromosome: str,          # 'chr1', 'chr22', etc.
    start: int,              # 0-based, inclusive
    end: int,                # 0-based, exclusive
    strand: str = '.',       # '+', '-', '.'
    name: str = None,        # optional label
)
# Properties: .center(), .width
# Methods: .resize(width), .overlaps(other), .contains(other), .intersect(other)
```

#### genome.Variant

```python
genome.Variant(
    chromosome: str,          # 'chr1', 'chr22', etc.
    position: int,           # 1-based (VCF compatible)
    reference_bases: str,    # 'A', 'AGGGATC', etc.
    alternate_bases: str,    # 'C', 'CGTCAAT', etc.
    name: str = None,        # optional label
)
# Properties: .reference_interval
# Methods: .reference_overlaps(interval), .alternate_overlaps(interval)
```

#### TrackData

```python
track_data.TrackData(
    values: np.ndarray,       # shape (positions, tracks)
    metadata: pd.DataFrame,   # track metadata (name, strand, ...)
    resolution: int = None,   # bp per value
    interval: genome.Interval = None,
)
# Methods:
#   .change_resolution(resolution)
#   .filter_to_positive_strand() / .filter_to_negative_strand() / .filter_to_unstranded()
#   .filter_to_nonpositive_strand()
#   .filter_to_strand(strand)
#   .filter_by_tissue(tissue_name)
#   .resize(width)
#   .slice_by_positions(start, end)
#   .slice_by_interval(interval)
#   .select_tracks_by_name(names)
#   .select_tracks_by_index(idx)
#   .reverse_complement()
# Arithmetic: trackdata1 - trackdata2 (element-wise subtraction)
```

#### Variant Scorers

```python
# GeneMaskLFCScorer: Gene-level log-fold-change
variant_scorers.GeneMaskLFCScorer(
    requested_output: dna_client.OutputType,
)

# CenterMaskScorer: Center-window comparison
variant_scorers.CenterMaskScorer(
    requested_output: dna_client.OutputType,
    width: int,                    # 501 for DNASE, 2001 for CHIP_HISTONE
    aggregation_type: variant_scorers.AggregationType,
)

# AggregationType enum
variant_scorers.AggregationType.DIFF_MEAN
variant_scorers.AggregationType.DIFF_LOG2_SUM

# tidy_scores: AnnData -> DataFrame
variant_scorers.tidy_scores(
    scores: list[list[AnnData]] | list[AnnData],
    match_gene_strand: bool = False,
) -> pd.DataFrame
```

### A.3 Visualization Components

| 컴포넌트 | 용도 | 주요 파라미터 |
|---------|------|-------------|
| `TranscriptAnnotation` | 유전자 구조 표시 | `transcripts` |
| `Tracks` | 예측/점수 track 표시 | `tdata`, `ylabel_template`, `filled`, `track_colors` |
| `OverlaidTracks` | REF vs ALT overlay | `tdata` (dict), `colors` (dict) |
| `Sashimi` | Splice junction arc | `tdata`, `ylabel_template` |
| `ContactMaps` | Hi-C 상호작용 맵 | `tdata`, `cmap`, `vmax` |
| `SeqLogo` | ISM sequence logo | `scores`, `scores_interval`, `ylabel` |
| `VariantAnnotation` | 변이 위치 표시 | `variants`, `alpha` |

**plot() 메인 함수:**

```python
fig = plot_components.plot(
    components: list,                 # 컴포넌트 리스트 (위→아래 순서)
    interval: genome.Interval,       # 표시 구간
    annotations: list = None,        # 추가 주석 (VariantAnnotation 등)
    title: str = None,               # 그림 제목
    fig_width: float = None,         # 그림 너비 (인치)
    despine_keep_bottom: bool = False,
)
```

**ylabel_template 포맷:**

| 패턴 | 용도 |
|------|------|
| `'{biosample_name} ({strand})'` | 일반 track |
| `'{biosample_name}\n{histone_mark}'` | ChIP-histone |
| `'CHIP_TF: {biosample_name}\n{transcription_factor} ({strand})'` | ChIP-TF |
| `'{name} ({strand})'` | SPLICE_SITES |

### A.4 Performance Notes

| 작업 | 소요 시간 | 비고 |
|------|-----------|------|
| Model creation | ~2초 | gRPC 연결 초기화 |
| predict_sequence (1MB) | ~3-5초 | OutputType 수에 따라 변동 |
| predict_interval (1MB) | ~3-5초 | 참조 게놈 조회 포함 |
| predict_variant (1MB) | ~6-10초 | REF + ALT 2회 예측 |
| score_variant (1 scorer) | ~3-5초 | scorer 수에 비례 |
| score_variant (19 scorers) | ~10-15초 | 전체 scorer |
| score_ism_variants (256bp) | ~10초 | 768 variants batch |
| Batch (5 variants, 6 scorers) | ~7초 | ~1.4초/variant |

**Rate Limiting**: 연속 호출 간 1-2초 delay 권장.

### A.5 Error Handling Patterns

```python
# API 키 검증
api_key = os.environ.get('ALPHAGENOME_API_KEY')
if not api_key:
    raise ValueError("ALPHAGENOME_API_KEY environment variable not set")

# Unsupported scorer 필터링
unsupported_scorers = [
    scorer for scorer in selected_scorers
    if organism.value not in variant_scorers.SUPPORTED_ORGANISMS[scorer.base_variant_scorer]
]
for unsupported in unsupported_scorers:
    selected_scorers.remove(unsupported)

# TrackData null 체크
if output.procap is not None and output.procap.values.shape[-1] > 0:
    # 데이터 처리
else:
    print("No data available")
```

### A.6 Known Limitations

1. **서열 길이 제한**: 정확히 지원되는 길이(16KB, 100KB, 500KB, 1MB) 중 하나여야 함
2. **PROCAP 제한**: 마우스 미지원, 6개 세포주만 사용 가능
3. **CONTACT_MAPS 제한**: 28개 human track만 존재 (대부분 EFO 기반)
4. **POLYADENYLATION**: 대응하는 OutputType 없이 scorer로만 존재
5. **Colab 의존성**: variant_scoring_ui.ipynb는 Google Colab widget에 의존 (CLI 대체 필요)
6. **Strand 매칭**: `tidy_scores`에서 `match_gene_strand=True` 시 strand 방향 일치 track만 포함
7. **ISM 규모**: 256bp (768 variants)는 관리 가능하지만, 1000bp+ 은 비용/시간 급증
8. **Ontology term 오류**: 잘못된 ontology code → 빈 TrackData (shape[-1] == 0), 에러 없음
9. **API Rate Limit**: 대량 분석 시 `time.sleep()` 사용 권장
10. **Organism 호환성**: 일부 scorer는 특정 organism만 지원 → `SUPPORTED_ORGANISMS` 확인 필요

### A.7 Best Practices

**1. API 호출 최적화:**

```python
# 단일 호출로 복수 OutputType 요청 (권장)
output = dna_model.predict_interval(
    interval=interval,
    requested_outputs={
        dna_client.OutputType.RNA_SEQ,
        dna_client.OutputType.DNASE,
        dna_client.OutputType.CHIP_HISTONE,
    },
    ontology_terms=ontology_terms,
)
# 개별 호출은 3배 느림
```

**2. Scorer 선택 전략:**

```python
# 탐색: 소수 핵심 scorer
exploratory = ['RNA_SEQ', 'DNASE', 'ATAC']

# 검증: 전체 19 scorer
comprehensive = list(variant_scorers.RECOMMENDED_VARIANT_SCORERS.values())

# 분야별:
splicing = ['SPLICE_SITES', 'SPLICE_SITE_USAGE', 'SPLICE_JUNCTIONS']
chromatin = ['DNASE', 'ATAC', 'CHIP_HISTONE', 'CHIP_TF']
expression = ['RNA_SEQ', 'CAGE', 'PROCAP']
```

**3. 결과 필터링:**

```python
# 고영향 (|raw_score| > 0.01)
high_impact = df[abs(df['raw_score']) > 0.01]

# Quantile 기반 (상위/하위 5%)
extreme = df[abs(df['quantile_score']) > 0.95]

# 조직 특이적
brain = df[df['ontology_curie'].str.startswith('UBERON:0000955')]
```

**4. TrackData Null Safety:**

```python
if output.procap is not None and output.procap.values.shape[-1] > 0:
    fig = plot_components.plot([...])
else:
    print("No data available")
```

---

## B. Reverse-Engineered Architecture

### B.1 gRPC-based Client-Server Pattern

```
┌──────────────────────┐     gRPC      ┌──────────────────────────┐
│   Python Client      │ ──────────>   │  AlphaGenome Server      │
│   (alphagenome pkg)  │ <──────────   │  (Google Cloud)          │
│                      │               │                          │
│  dna_client.create() │               │  Model Inference Engine  │
│  ├─ predict_sequence │               │  ├─ DNA Encoder          │
│  ├─ predict_interval │               │  ├─ Output Decoders (11) │
│  ├─ predict_variant  │               │  ├─ Scorer Algorithms    │
│  ├─ score_variant    │               │  └─ Reference Genome     │
│  ├─ score_variants   │               │      (hg38, mm10)        │
│  └─ score_ism_vars   │               │                          │
└──────────────────────┘               └──────────────────────────┘
```

### B.2 Request/Response Flow

```
1. Client: dna_model.predict_interval(interval, outputs, ontology_terms)
   │
   ├─ 2. interval → 서열 좌표 직렬화
   ├─ 3. gRPC Request → Google Cloud AlphaGenome Server
   ├─ 4. Server: 참조 게놈에서 서열 조회
   ├─ 5. Server: DNA Encoder → 1MB 서열 인코딩
   ├─ 6. Server: OutputType별 Decoder 실행
   │       ├─ RNA_SEQ Decoder → (1048576, N) predictions
   │       ├─ DNASE Decoder   → (1048576, N) predictions
   │       └─ ...
   ├─ 7. Server: ontology_terms 필터 적용
   ├─ 8. gRPC Response: serialized predictions + metadata
   │
   └─ 9. Client: ModelOutput 구성
          ├─ output.rna_seq → TrackData(values, metadata, interval)
          ├─ output.dnase   → TrackData(...)
          └─ ...
```

### B.3 Scorer Algorithm Internals

#### GeneMaskLFCScorer

```
Input: REF predictions, ALT predictions, gene annotation
├─ 1. Gene mask: gene 구간 = 1, 외부 = 0
├─ 2. REF gene signal = sum(REF * gene_mask)
├─ 3. ALT gene signal = sum(ALT * gene_mask)
├─ 4. LFC = log2(ALT_signal / REF_signal)
└─ Output: gene x track matrix (AnnData)
```

#### CenterMaskScorer

```
Input: REF predictions, ALT predictions, center position, width
├─ 1. Center window: [center - width/2, center + width/2]
├─ 2. REF window signal = aggregate(REF[window])
├─ 3. ALT window signal = aggregate(ALT[window])
├─ 4. Aggregation:
│     ├─ DIFF_MEAN: mean(ALT) - mean(REF)
│     └─ DIFF_LOG2_SUM: log2(sum(ALT)) - log2(sum(REF))
└─ Output: 1 x track matrix (AnnData)
```

#### SpliceSitesScorer

```
Input: REF splice_sites, ALT splice_sites
├─ 1. REF splice probabilities (donor+, donor-, acceptor+, acceptor-)
├─ 2. ALT splice probabilities
├─ 3. Difference = ALT_prob - REF_prob
└─ Output: position-level splice site probability changes
```

### B.4 TrackData Internal Structure

```
TrackData
├── values: np.ndarray          # (positions, n_tracks), float32
│   └── positions = interval.width / resolution
├── metadata: pd.DataFrame      # n_tracks rows
│   ├── name, strand, ontology_curie, biosample_name
│   ├── biosample_type, data_source, nonzero_mean
│   └── [assay-specific]: histone_mark, transcription_factor, Assay title
├── resolution: int
└── interval: genome.Interval
```

### B.5 AnnData Structure (score_variant)

```python
adata = variant_scores[0]
adata.X                    # (n_genes, n_tracks) score matrix
adata.obs                  # gene metadata (gene_name, gene_id, strand)
adata.var                  # track metadata (ontology_curie, biosample_name, output_type)
adata.uns                  # variant, scorer metadata
adata.layers['quantiles']  # quantile-normalized scores

# tidy_scores 결과 컬럼 (19개):
# variant_id, scored_interval, gene_name, gene_id, gene_strand,
# ontology_curie, biosample_name, biosample_type, output_type,
# variant_scorer, raw_score, quantile_score, ...
```

---

## C. Complete File Manifest

### C.1 Results Directories

```
results/
├── quick_start/                     (5 files)
│   ├── results.json                 (3.9 KB)   실행 메타데이터
│   ├── run_quick_start.py           (6.8 KB)   실행 스크립트
│   ├── cyp2b6_rna_seq.png          (130 KB)   CYP2B6 RNA-seq
│   ├── variant_effect.png           (147 KB)   REF vs ALT overlay
│   └── variant_scores.csv           (3.8 MB)   14,652 rows
│
├── essential_commands/              (3 files)
│   ├── results.json                 (2.9 KB)
│   ├── run_essential_commands.py     (7.2 KB)
│   └── README.md                    (3.0 KB)
│
├── visualization_tour/              (11 files)
│   ├── 01_rna_seq.png ~ 07_contact_maps.png   7개 시각화
│   ├── results.json, README.md, EXECUTION_SUMMARY.md
│
├── batch_variant_scoring/           (6 files)
│   ├── variant_scores.csv           (30 MB)    121,550 rows
│   ├── high_impact_variants.csv     (1.1 MB)   4,509 rows
│   ├── variant_scores_summary.json  (9.1 KB)
│   ├── README.md, ANALYSIS_SUMMARY.md, EXECUTION_LOG.md
│
├── analysis_workflow/               (Variable)
│   ├── jurkat_variant_effect.png
│   ├── variant_analysis_results.json
│   └── comparison_*.png             (7 그룹)
│
├── tissue_ontology/                 (4 files)
│   ├── track_counts.csv             (352 B)
│   ├── ontology_coverage.csv        (466 B)
│   ├── ontology_terms.json          (30 KB)
│   └── ontology_summary.json        (437 B)
│
├── variant_scoring_cli/             (4 files)
│   ├── variant_scores.csv           (10 MB)    38,357 rows
│   ├── variant_summary.json         (660 B)
│   ├── plot_rna_seq.png, plot_dnase.png
│
├── ism_256bp/                       (3 files)
│   ├── ism_scores.csv               (31 KB)    768 rows
│   ├── ism_heatmap.png              (61 KB)
│   └── results.json                 (988 B)
│
├── procap_visualization/            (2 files)
│   ├── procap.png                   (205 KB)
│   └── results.json                 (932 B)
│
└── chip_tf_analysis/                (5 files)
    ├── chip_tf_k562_ctcf.png        (70 KB)
    ├── chip_tf_hepg2_ctcf.png       (94 KB)
    ├── chip_tf_ctcf_rad21_coloc.png (120 KB)
    ├── chip_tf_multi_tf.png         (182 KB)
    └── results.json                 (1.2 KB)
```

### C.2 Scripts

```
scripts/
├── run_batch_variant_scoring.py     (8.2 KB)
├── run_analysis_workflow.py         (16 KB)
├── run_visualization_tour.py        (13 KB)
├── run_variant_scoring_cli.py       (14 KB)
├── run_ism_256bp.py                 (8.7 KB)
├── run_procap_visualization.py      (6.6 KB)
├── run_chip_tf_analysis.py          (12 KB)
├── run_splice_site_usage.py         (new)
├── tissue_ontology_mapping.py       (9.2 KB)
├── generate_batch_summary.py        (3.0 KB)
├── analyze_batch_results.py         (2.5 KB)
└── verify_install.py                (2.1 KB)
```

### C.3 Tutorials

```
tutorials/
├── quick_start.ipynb                (648 KB)
├── essential_commands.ipynb         (42 KB)
├── visualization_modality_tour.ipynb (3.6 MB)
├── batch_variant_scoring.ipynb      (466 KB)
├── example_analysis_workflow.ipynb   (524 KB)
├── tissue_ontology_mapping.ipynb    (2.2 MB)
└── variant_scoring_ui.ipynb         (26 MB)    Colab 의존적, CLI 대체
```

### C.4 Summary Statistics

| 카테고리 | 수량 |
|---------|------|
| Results 디렉토리 | 10 |
| 총 결과 파일 | ~50+ |
| Scripts | 12 |
| Tutorials | 7 |
| PNG 시각화 | ~20+ |
| CSV 데이터 | 6 |
| JSON 메타데이터 | 10+ |
| 총 결과 데이터 크기 | ~46 MB |

---

## D. Ontology Code Reference

### D.1 자주 사용된 Ontology Codes

**UBERON (Anatomy):**

| Code | Name | 사용 튜토리얼 |
|------|------|-------------|
| UBERON:0002048 | Lung (폐) | Quick Start |
| UBERON:0000955 | Brain (뇌) | Quick Start |
| UBERON:0001114 | Right lobe of liver (간 우엽) | Quick Start |
| UBERON:0001157 | Colon (대장) | Viz Tour, CLI |
| UBERON:0001159 | Sigmoid colon (에스상결장) | Viz Tour |
| UBERON:0001155 | Transverse colon (횡행결장) | Viz Tour |
| UBERON:0000317 | Colonic mucosa (대장 점막) | Viz Tour |

**EFO (Experimental Factors / Cell Lines):**

| Code | Name | 사용 튜토리얼 |
|------|------|-------------|
| EFO:0002067 | K562 (CML 세포주) | ChIP-TF, Batch |
| EFO:0001187 | HepG2 (간암 세포주) | ChIP-TF |
| EFO:0002824 | HCT116 (대장암 세포주) | Viz Tour |
| EFO:0002106 | A673 (Ewing sarcoma) | PROCAP |
| EFO:0001099 | Caco-2 (대장선암) | PROCAP |
| EFO:0002819 | Calu-3 (폐선암) | PROCAP |
| EFO:0001200 | MCF 10A (유방 상피) | PROCAP |
| EFO:0007950 | GM12878 | Batch |

**CL (Cell Ontology):**

| Code | Name | 사용 튜토리얼 |
|------|------|-------------|
| CL:0001059 | CD34+ common myeloid progenitor | Analysis Workflow |
| CL:0002618 | HUVEC (제대 정맥 내피세포) | PROCAP |
| CL:0000775 | Neutrophil | Batch |
| CL:0000771 | Eosinophil | Batch |
| CL:0001054 | CD14+ monocyte | Batch |
| CL:0000863 | Inflammatory macrophage | Batch |

### D.2 Ontology 검색 결과 (Tissue Ontology Tutorial)

| 검색어 | 주요 Ontology Terms |
|--------|-------------------|
| brain | UBERON:0000955, CL:0000540, CL:0000117 |
| liver | UBERON:0002107, UBERON:0002171, EFO:0001187 |
| heart | UBERON:0000948, UBERON:0002082 |
| lung | UBERON:0002048, UBERON:0002049 |
| T cell | CL:0000084, CL:0000624, CL:0000625 |
| neuron | CL:0000540, CL:0000117, CL:0000527 |

---

## E. Execution Timeline

### E.1 튜토리얼 실행 순서

```
2026-02-04  16:32  tutorials/*.ipynb 7개 다운로드
2026-02-04  15:26  scripts/verify_install.py 작성

2026-02-05  09:41  Essential Commands 실행
            09:42  Quick Start 실행
            09:45  Tissue Ontology Mapping 실행
            09:45  Visualization Tour 시작
            09:48  Visualization Tour 완료 (7/7 PNG)
            09:46  Batch Variant Scoring 실행
            09:49  Batch Analysis Summary 작성
            09:50  Analysis Workflow 실행

2026-02-05  10:17  Variant Scoring CLI 스크립트 작성
            10:19  Variant Scoring CLI 실행 완료
            10:21  ChIP-TF Analysis 실행 완료

2026-02-05  16:31  PROCAP Visualization 실행
            16:34  PROCAP 완료
            16:37  ISM 256bp 실행
            16:38  ISM 256bp 완료
```

### E.2 API 호출 통계

| 스크립트 | API 호출 수 | 소요 시간 |
|---------|-----------|----------|
| Quick Start | ~6 | ~30초 |
| Essential Commands | 0 (로컬) | ~1초 |
| Visualization Tour | 7 | ~40초 |
| Batch Variant Scoring | 5 | ~7초 |
| Analysis Workflow | ~130 | ~5분 |
| Tissue Ontology | 2 | ~10초 |
| Variant Scoring CLI | 2 | ~20초 |
| ISM 256bp | 1 (batch) | ~10초 |
| PROCAP Visualization | 1 | ~5초 |
| ChIP-TF Analysis | 2 | ~15초 |
| **총계** | **~156** | **~8분** |

---

> **이전 모듈**: [09-tall-workflow.md](09-tall-workflow.md) — T-ALL 실전 워크플로우
>
> **전체 목차**: [00-INDEX.md](00-INDEX.md)
