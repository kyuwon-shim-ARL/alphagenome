# AlphaGenome 튜토리얼 실행 상세 보고서

**작성일**: 2026년 2월 6일
**실행 완료**: 7개 핵심 튜토리얼 + 4개 추가 스크립트
**OutputType 커버리지**: 11가지 전체
**핵심 수치**: 121,550 변이 점수, 38,357 CLI 점수, 768 ISM 변이, 7 시각화 PNG, 128 분석 변이

---

## 1. Executive Summary

### 1.1 AlphaGenome 소개

AlphaGenome은 Google DeepMind가 개발한 DNA regulatory sequence prediction 모델이다. 최대 1,048,576 bp (1 Mbp) 길이의 게놈 서열을 입력으로 받아, 11가지 서로 다른 생물학적 출력(OutputType)을 예측한다. 이 모델은 인간(Homo sapiens)과 마우스(Mus musculus) 두 종을 지원하며, 5,559개의 고유한 ontology term에 걸쳐 5,563개 이상의 human track을 제공한다.

AlphaGenome의 핵심 기능은 다음과 같다:

- **DNA 서열 기반 예측**: 원시 DNA 서열 또는 참조 게놈 좌표(Interval)로부터 RNA 발현, 크로마틴 접근성, 전사인자 결합, 히스톤 변형, 3D 크로마틴 구조 등을 예측
- **변이 효과 분석**: 단일 변이(SNV, insertion, deletion)가 모든 출력에 미치는 영향을 정량적으로 평가
- **In Silico Mutagenesis (ISM)**: 연속 구간 내 모든 가능한 단일 염기 변이를 체계적으로 스캔
- **조직/세포 특이적 예측**: 687개 이상의 ontology term을 통해 세포/조직 특이적 예측 수행

### 1.2 실행 범위

총 **7개 튜토리얼**과 **4개 추가 분석 스크립트**를 실행했다:

| # | 항목 | 유형 | 소스 | 상태 |
|---|------|------|------|------|
| 1 | Quick Start | Tutorial | `tutorials/quick_start.ipynb` | 완료 |
| 2 | Essential Commands | Tutorial | `tutorials/essential_commands.ipynb` | 완료 |
| 3 | Visualization Modality Tour | Tutorial | `tutorials/visualization_modality_tour.ipynb` | 완료 |
| 4 | Batch Variant Scoring | Tutorial | `tutorials/batch_variant_scoring.ipynb` | 완료 |
| 5 | Example Analysis Workflow | Tutorial | `tutorials/example_analysis_workflow.ipynb` | 완료 |
| 6 | Tissue Ontology Mapping | Tutorial | `tutorials/tissue_ontology_mapping.ipynb` | 완료 |
| 7 | Variant Scoring UI (CLI 대체) | Tutorial/Script | `scripts/run_variant_scoring_cli.py` | 완료 |
| 8 | ISM 256bp Analysis | Script | `scripts/run_ism_256bp.py` | 완료 |
| 9 | PROCAP Visualization | Script | `scripts/run_procap_visualization.py` | 완료 |
| 10 | ChIP-TF Analysis | Script | `scripts/run_chip_tf_analysis.py` | 완료 |
| 11 | Batch Results Analysis | Script | `scripts/analyze_batch_results.py` | 완료 |

### 1.3 핵심 수치 요약

| 항목 | 수치 |
|------|------|
| 총 변이 점수 행 (Batch) | 121,550 |
| 고영향 변이 행 (Batch) | 4,509 |
| 총 변이 점수 행 (CLI) | 38,357 |
| ISM 변이 (Quick Start) | 192 (64bp x 3) |
| ISM 변이 (256bp) | 768 (256bp x 3) |
| 분석 워크플로우 변이 | 128 (32 oncogenic + 96 background) |
| 시각화 PNG (Viz Tour) | 7 |
| ChIP-TF Plot | 4 |
| PROCAP Plot | 1 |
| OutputType 커버리지 | 11/11 |
| Human track 수 | 5,563 |
| Mouse track 수 | 1,038 |
| Human ontology term 수 | 5,559 |

### 1.4 실행 환경

```
프로젝트 경로: /home/kyuwon/projects/alphagenome
Python: 3.11+ (uv-managed virtual environment)
참조 게놈: hg38 (GENCODE v46)
API: AlphaGenome gRPC Client (Google DeepMind)
주요 라이브러리:
  - alphagenome (Google DeepMind API Client)
  - pandas, numpy (데이터 처리)
  - matplotlib, plotnine (시각화)
  - anndata (AnnData 구조)
  - tqdm (진행률 표시)
```

---

## 2. Phase I: Input Architecture

AlphaGenome의 입력 시스템은 4가지 핵심 데이터 구조로 구성된다: `genome.Interval`, `genome.Variant`, raw DNA sequence string, 그리고 Ontology system.

### 2.1 genome.Interval

`genome.Interval`은 게놈 좌표 구간을 나타내는 dataclass이다.

**필드 구조:**

| 필드 | 타입 | 설명 | 필수 |
|------|------|------|------|
| `chromosome` | `str` | 염색체 이름 (예: `'chr1'`) | Yes |
| `start` | `int` | 시작 위치 (0-based, inclusive) | Yes |
| `end` | `int` | 종료 위치 (0-based, exclusive) | Yes |
| `strand` | `str` | 가닥 방향 (`'+'`, `'-'`, `'.'`) | No (기본값: `'.'`) |
| `name` | `str` | 구간 이름 | No |

**좌표 체계**: 0-based half-open 구간 (start 포함, end 미포함). 이는 BED 포맷과 동일하다.

**주요 메서드:**

| 메서드 | 설명 | 반환 타입 |
|--------|------|-----------|
| `resize(width)` | 중심을 기준으로 구간 크기 변경 | `Interval` |
| `overlaps(other)` | 두 구간 겹침 여부 확인 | `bool` |
| `contains(other)` | 다른 구간 완전 포함 여부 | `bool` |
| `intersect(other)` | 겹치는 구간 반환 | `Interval` |
| `center` / `center()` | 구간 중심 좌표 | `int` |
| `width` | 구간 너비 (end - start) | `int` |

**코드 예시:**

```python
from alphagenome.data import genome

# Interval 생성 (keyword arguments)
interval = genome.Interval(chromosome='chr1', start=1_000, end=1_010)
# 결과: chr1:1000-1010:.

# 속성 확인
interval.center()  # 1005
interval.width     # 10

# resize: 중심 기준 크기 변경
resized = interval.resize(100)
# 결과: chr1:955-1055:.

# 겹침 검사
second = genome.Interval(chromosome='chr1', start=1_005, end=1_015)
interval.overlaps(second)   # True
interval.contains(second)   # False
interval.intersect(second)  # chr1:1005-1010:.

# Positional arguments도 가능
interval = genome.Interval('chr22', 36_150_498, 36_252_898)

# Strand 지정
interval = genome.Interval(
    chromosome='chr19', start=40991281, end=41018398, strand='+'
)
# 결과: chr19:40991281-41018398:+
```

**Essential Commands 결과에서 검증된 값:**
- 생성: `chr1:1000-1010:.`
- center: `1005`
- width: `10`
- resize(100): `chr1:955-1055:.`
- overlaps: `True`
- contains: `False`
- intersect: `chr1:1005-1010:.`

### 2.2 genome.Variant

`genome.Variant`는 게놈 변이를 나타내는 dataclass이다.

**필드 구조:**

| 필드 | 타입 | 설명 | 필수 |
|------|------|------|------|
| `chromosome` | `str` | 염색체 이름 | Yes |
| `position` | `int` | 변이 위치 (1-based, VCF 호환) | Yes |
| `reference_bases` | `str` | 참조 대립유전자 | Yes |
| `alternate_bases` | `str` | 대체 대립유전자 | Yes |
| `name` | `str` | 변이 이름 | No |

**위치 체계**: 1-based (VCF 포맷과 동일). Interval의 0-based와 다르므로 주의가 필요하다.

**변이 유형:**

| 유형 | 조건 | 예시 |
|------|------|------|
| SNV | `len(ref) == 1 and len(alt) == 1` | `A>C` |
| Insertion | `len(ref) < len(alt)` | `T>CGTCAAT` |
| Deletion | `len(ref) > len(alt)` | `AGGGATC>C` |

**주요 속성:**

| 속성/메서드 | 설명 |
|-------------|------|
| `reference_interval` | 참조 대립유전자 구간 (Interval) |
| `reference_overlaps(interval)` | 참조 대립유전자와 구간 겹침 |
| `alternate_overlaps(interval)` | 대체 대립유전자와 구간 겹침 |

**코드 예시:**

```python
from alphagenome.data import genome

# SNV 생성
variant = genome.Variant(
    chromosome='chr3',
    position=10000,
    reference_bases='A',
    alternate_bases='C',
)
# 결과: chr3:10000:A>C

# Insertion 생성
insertion = genome.Variant(
    chromosome='chr3',
    position=10000,
    reference_bases='T',
    alternate_bases='CGTCAAT',
)
# 결과: chr3:10000:T>CGTCAAT

# Deletion 생성
deletion = genome.Variant(
    chromosome='chr3',
    position=10000,
    reference_bases='AGGGATC',
    alternate_bases='C',
)
# 결과: chr3:10000:AGGGATC>C

# Quick Start에서 사용된 실제 변이
variant = genome.Variant(
    chromosome='chr22',
    position=36201698,
    reference_bases='A',
    alternate_bases='C',
)

# reference_interval
# chr3:9999-10000:. (1-based position -> 0-based interval)

# Overlap 검사
interval = genome.Interval(chromosome='chr3', start=10_005, end=10_010)
variant.reference_overlaps(interval)   # False
variant.alternate_overlaps(interval)   # True (insertion이 해당 구간에 겹침)
```

**Essential Commands 결과에서 검증된 값:**
- SNV: `chr3:10000:A>C`
- Insertion: `chr3:10000:T>CGTCAAT`
- Deletion: `chr3:10000:AGGGATC>C`
- reference_interval: `chr3:9999-10000:.`
- reference_overlaps: `False`
- alternate_overlaps: `True`

### 2.3 DNA Sequence Input

AlphaGenome은 raw DNA string을 직접 입력으로 받을 수 있다. `predict_sequence()` 메서드에 사용된다.

**지원 서열 길이:**

| 상수명 | 길이 (bp) | 용도 |
|--------|-----------|------|
| `SEQUENCE_LENGTH_16KB` | 16,384 | ISM 분석, 빠른 테스트 |
| `SEQUENCE_LENGTH_100KB` | 100,000 | 중간 규모 분석 |
| `SEQUENCE_LENGTH_500KB` | 500,000 | 대규모 분석 |
| `SEQUENCE_LENGTH_1MB` | 1,048,576 | 전체 맥락 분석 (권장) |

서열은 정확히 해당 길이와 일치해야 한다. 짧은 서열은 'N'으로 패딩한다:

```python
from alphagenome.models import dna_client

# 서열 패딩 (Quick Start에서 사용된 방법)
sequence = 'GATTACA'.center(dna_client.SEQUENCE_LENGTH_1MB, 'N')
# 1,048,576 bp 길이의 서열 (중앙에 'GATTACA', 나머지 'N')

# SUPPORTED_SEQUENCE_LENGTHS 딕셔너리로 접근 가능
seq_len = dna_client.SUPPORTED_SEQUENCE_LENGTHS['SEQUENCE_LENGTH_1MB']
# 1048576
```

**주의사항:**
- 서열 길이가 정확히 일치하지 않으면 API 오류 발생
- 'N' 패딩은 해당 위치의 예측에 영향을 미치지 않음
- ISM 분석에는 16KB 길이가 적합 (연산 효율)
- 변이 분석에는 1MB 길이가 권장 (충분한 맥락)

**길이 선택 가이드:**

| 길이 | 맥락 범위 | 적합한 용도 | 비고 |
|------|-----------|------------|------|
| 16KB | ~8KB 양쪽 | ISM, 빠른 프로토타이핑 | score_ism_variants 권장 |
| 100KB | ~50KB 양쪽 | 중규모 유전자 분석 | 일반적 유전자 크기 커버 |
| 500KB | ~250KB 양쪽 | Large locus 분석 | enhancer-promoter 상호작용 포함 |
| 1MB | ~500KB 양쪽 | 전체 맥락 분석 | 최대 규제 맥락 제공, 권장 |

대부분의 튜토리얼과 스크립트에서 1MB (`SEQUENCE_LENGTH_1MB` = 1,048,576)를 사용한다. 이는 장거리 regulatory interaction을 포착하기 위한 것이다. ISM은 예외적으로 16KB를 사용하는데, ISM은 변이 자체가 다수이므로 더 좁은 맥락에서도 충분한 정보를 얻을 수 있다.

### 2.4 Ontology System

AlphaGenome은 5가지 ontology 체계를 사용하여 조직/세포 타입을 지정한다.

**5가지 Ontology 시스템:**

| 약칭 | 정식명 | 설명 | Human Term 수 |
|------|--------|------|--------------|
| CL | Cell Ontology | 세포 타입 | 1,394 |
| CLO | Cell Line Ontology | 세포주 | 65 |
| EFO | Experimental Factor Ontology | 실험 조건/세포주 | 2,402 |
| NTR | New Term Request | 신규 용어 | 93 |
| UBERON | Uberon Anatomy Ontology | 해부학적 구조 | 1,605 |
| **합계** | | | **5,559** |

**Ontology Curie 포맷**: `PREFIX:숫자코드`

```python
# 튜토리얼에서 사용된 주요 ontology code
ontology_examples = {
    'UBERON:0002048': 'lung (폐)',
    'UBERON:0000955': 'brain (뇌)',
    'UBERON:0001114': 'right lobe of liver (간 우엽)',
    'UBERON:0001157': 'colon (대장)',
    'UBERON:0001159': 'sigmoid colon (에스상결장)',
    'UBERON:0001155': 'transverse colon (횡행결장)',
    'EFO:0002067':   'K562 (만성 골수성 백혈병 세포주)',
    'EFO:0001187':   'HepG2 (간세포암 세포주)',
    'EFO:0002824':   'HCT116 (대장암 세포주)',
    'EFO:0002106':   'A673 (Ewing sarcoma 세포주)',
    'EFO:0001099':   'Caco-2 (대장선암 세포주)',
    'EFO:0002819':   'Calu-3 (폐선암 세포주)',
    'EFO:0001200':   'MCF 10A (유방 상피 세포주)',
    'CL:0002618':    'HUVEC (제대 정맥 내피세포)',
    'CL:0001059':    'CD34+ common myeloid progenitor',
}
```

**Ontology별 OutputType 커버리지 (Human):**

| OutputType | CL | CLO | EFO | NTR | UBERON | 합계 |
|------------|-----|-----|------|------|--------|------|
| ATAC | 15 | 40 | 65 | 3 | 44 | 167 |
| CAGE | 286 | 0 | 38 | 0 | 222 | 546 |
| DNASE | 94 | 12 | 95 | 13 | 91 | 305 |
| RNA_SEQ | 238 | 1 | 130 | 15 | 283 | 667 |
| CHIP_HISTONE | 323 | 5 | 389 | 35 | 364 | 1,116 |
| CHIP_TF | 63 | 4 | 1,447 | 3 | 100 | 1,617 |
| SPLICE_SITE_USAGE | 248 | 2 | 134 | 16 | 334 | 734 |
| SPLICE_JUNCTIONS | 124 | 1 | 67 | 8 | 167 | 367 |
| CONTACT_MAPS | 1 | 0 | 27 | 0 | 0 | 28 |
| PROCAP | 2 | 0 | 10 | 0 | 0 | 12 |
| **합계** | **1,394** | **65** | **2,402** | **93** | **1,605** | **5,559** |

> 참고: SPLICE_SITES는 총 4개 track이지만 ontology 기반이 아닌 고정 track이므로 위 표에서 제외됨.

### 2.5 Output Metadata

`output_metadata()` 메서드는 사용 가능한 모든 track에 대한 메타데이터를 반환한다.

```python
from alphagenome.models import dna_client

dna_model = dna_client.create(api_key)

# Human metadata 조회
human_metadata = dna_model.output_metadata(
    organism=dna_client.Organism.HOMO_SAPIENS
)

# 모든 track을 하나의 DataFrame으로 합침
metadata_df = human_metadata.concatenate()
# 컬럼: output_type, ontology_curie, biosample_name, biosample_type,
#        strand, data_source, nonzero_mean, ...
```

**Human Track 수 (OutputType별):**

| OutputType | Human Track 수 | Mouse Track 수 |
|------------|----------------|----------------|
| ATAC | 167 | 18 |
| CAGE | 546 | 188 |
| DNASE | 305 | 67 |
| RNA_SEQ | 667 | 173 |
| CHIP_HISTONE | 1,116 | 183 |
| CHIP_TF | 1,617 | 127 |
| SPLICE_SITES | 4 | 4 |
| SPLICE_SITE_USAGE | 734 | 180 |
| SPLICE_JUNCTIONS | 367 | 90 |
| CONTACT_MAPS | 28 | 8 |
| PROCAP | 12 | - |
| **합계** | **5,563** | **1,038** |

> 참고: Track 수(5,563)와 ontology term 수(5,559)는 일부 track이 같은 ontology term을 공유하므로 다를 수 있다.

**메타데이터 필드 예시 (DNase Lung):**

```python
{
    'name': 'UBERON:0002048 DNase-seq',
    'strand': '.',
    'Assay title': 'DNase-seq',
    'ontology_curie': 'UBERON:0002048',
    'biosample_name': 'lung',
    'biosample_type': 'tissue',
    'biosample_life_stage': 'embryonic',
    'data_source': 'encode',
    'endedness': 'paired',
    'genetically_modified': False,
    'nonzero_mean': 0.4275
}
```

---

## 3. Phase II: Processing Pipeline

### 3.1 Model Creation

AlphaGenome 모델 생성은 `dna_client.create()` 함수를 통해 이루어진다.

```python
from alphagenome.models import dna_client

# 모델 생성
api_key = os.environ.get('ALPHAGENOME_API_KEY')
dna_model = dna_client.create(api_key)
```

`dna_client.create()`는 내부적으로 gRPC 클라이언트를 초기화하고 Google Cloud 기반의 AlphaGenome 서버에 연결한다. 반환되는 `DNAModel` 객체는 모든 예측 메서드의 진입점이다.

**Organism enum:**

```python
dna_client.Organism.HOMO_SAPIENS     # 인간
dna_client.Organism.MUS_MUSCULUS     # 마우스
```

Quick Start에서 마우스 예측 검증:

```python
output = dna_model.predict_sequence(
    sequence='GATTACA'.center(dna_client.SEQUENCE_LENGTH_1MB, 'N'),
    organism=dna_client.Organism.MUS_MUSCULUS,
    requested_outputs=[dna_client.OutputType.DNASE],
    ontology_terms=['UBERON:0002048'],  # Lung
)
# Mouse DNase Shape: (1048576, 1)
```

**OutputType enum (11가지):**

```python
dna_client.OutputType.ATAC
dna_client.OutputType.CAGE
dna_client.OutputType.CHIP_HISTONE
dna_client.OutputType.CHIP_TF
dna_client.OutputType.CONTACT_MAPS
dna_client.OutputType.DNASE
dna_client.OutputType.PROCAP
dna_client.OutputType.RNA_SEQ
dna_client.OutputType.SPLICE_JUNCTIONS
dna_client.OutputType.SPLICE_SITES
dna_client.OutputType.SPLICE_SITE_USAGE
```

### 3.2 predict_sequence() Workflow

`predict_sequence()`는 raw DNA string을 입력으로 받아 예측을 수행한다.

**시그니처:**

```python
output = dna_model.predict_sequence(
    sequence: str,                           # DNA 서열 (정확한 길이 필요)
    requested_outputs: list[OutputType],     # 요청할 출력 타입
    ontology_terms: list[str] = None,        # 조직/세포 ontology 코드
    organism: Organism = HOMO_SAPIENS,       # 생물종
)
```

**반환 타입**: `ModelOutput`
- 각 OutputType에 대응하는 attribute를 가짐 (예: `output.dnase`, `output.cage`)
- 각 attribute는 `TrackData` 객체: `.values` (numpy array), `.metadata` (DataFrame)

**실행 예시 (Quick Start):**

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

CAGE의 shape가 (1048576, 4)인 이유: CAGE는 strand-specific이므로 각 ontology term에 대해 `+` strand와 `-` strand 각각 1개 track씩 생성된다. 2개 ontology term x 2 strand = 4 track.

### 3.3 predict_interval() Workflow

`predict_interval()`은 게놈 좌표(Interval)를 입력으로 받아, 참조 게놈에서 자동으로 서열을 조회한 후 예측을 수행한다.

**시그니처:**

```python
output = dna_model.predict_interval(
    interval: genome.Interval,               # 게놈 구간
    requested_outputs: set[OutputType],      # 요청할 출력 타입 (set 또는 list)
    ontology_terms: list[str] = None,        # 조직/세포 ontology 코드
    organism: Organism = HOMO_SAPIENS,       # 생물종
)
```

**반환 타입**: `ModelOutput` (predict_sequence와 동일)

**실행 예시 (Quick Start - CYP2B6 유전자 RNA-seq):**

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

### 3.4 predict_variant() Workflow

`predict_variant()`는 특정 변이의 효과를 예측한다. 참조(reference) 서열과 대체(alternate) 서열 각각에 대해 예측을 수행하고, 두 결과를 `VariantOutput`으로 반환한다.

**시그니처:**

```python
variant_output = dna_model.predict_variant(
    interval: genome.Interval,               # 게놈 구간
    variant: genome.Variant,                 # 변이
    requested_outputs: set[OutputType],      # 요청할 출력 타입
    ontology_terms: list[str] = None,        # 조직/세포 ontology 코드
    organism: Organism = HOMO_SAPIENS,       # 생물종
)
```

**반환 타입**: `VariantOutput`
- `.reference`: `ModelOutput` (참조 서열 예측 결과)
- `.alternate`: `ModelOutput` (대체 서열 예측 결과)

**실행 예시 (Quick Start):**

```python
# 변이 정의
variant = genome.Variant(
    chromosome='chr22',
    position=36201698,
    reference_bases='A',
    alternate_bases='C',
)

# 변이 중심 1MB interval
interval = variant.reference_interval.resize(dna_client.SEQUENCE_LENGTH_1MB)

# 예측 수행
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
# TAL1 locus interval
tal1_interval = genome.Interval(
    chromosome='chr1', start=47209255, end=47242023, strand='-'
)

# 복수 출력 타입 요청
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

### 3.5 score_variant() Workflow

`score_variant()`는 변이의 효과를 사전 정의된 scorer를 사용하여 정량적으로 평가한다. 결과는 `AnnData` 객체 리스트로 반환된다.

**시그니처:**

```python
variant_scores = dna_model.score_variant(
    interval: genome.Interval,                # 게놈 구간
    variant: genome.Variant,                  # 변이
    variant_scorers: list[VariantScorer],     # scorer 리스트
    organism: Organism = HOMO_SAPIENS,        # 생물종
)
```

**반환 타입**: `list[AnnData]`

각 AnnData 객체의 구조:
- `.X`: score matrix (genes x tracks)
- `.obs`: gene metadata (gene_name, gene_id, strand 등)
- `.var`: track metadata (ontology_curie, biosample_name 등)
- `.uns`: variant metadata, scorer info
- `.layers['quantiles']`: quantile-normalized scores

**RECOMMENDED_VARIANT_SCORERS (19개 키):**

```python
from alphagenome.models import variant_scorers

all_scorers = variant_scorers.RECOMMENDED_VARIANT_SCORERS
# 19개 키:
#   11개 base scorer:
#     ATAC, CAGE, CHIP_HISTONE, CHIP_TF, CONTACT_MAPS,
#     DNASE, PROCAP, RNA_SEQ, SPLICE_SITES, SPLICE_SITE_USAGE, SPLICE_JUNCTIONS
#   1개 scorer-only (대응하는 OutputType 없음):
#     POLYADENYLATION
#   7개 _ACTIVE variant:
#     ATAC_ACTIVE, CAGE_ACTIVE, CHIP_HISTONE_ACTIVE, CHIP_TF_ACTIVE,
#     DNASE_ACTIVE, PROCAP_ACTIVE, RNA_SEQ_ACTIVE
```

**Scorer 타입별 분류:**

| Scorer 클래스 | 알고리즘 | 사용 OutputType | 주요 파라미터 |
|---------------|---------|----------------|--------------|
| `GeneMaskLFCScorer` | Gene-level log-fold-change | RNA_SEQ, CAGE, PROCAP | `requested_output` |
| `CenterMaskScorer` | Center-window 비교 | DNASE, ATAC, CHIP_HISTONE, CHIP_TF | `requested_output`, `width`, `aggregation_type` |
| `SpliceSitesScorer` | Splice site 확률 변화 | SPLICE_SITES | - |

**tidy_scores() 변환:**

```python
# AnnData -> DataFrame 변환
tidy_df = variant_scorers.tidy_scores(
    variant_scores,
    match_gene_strand=True  # gene strand와 track strand 매칭
)
# 결과: (rows, 19 columns) DataFrame
# 컬럼: variant_id, gene_name, ontology_curie, raw_score, quantile_score, ...
```

**실행 예시 (Quick Start):**

```python
variant_scorer = variant_scorers.RECOMMENDED_VARIANT_SCORERS['RNA_SEQ']
variant_scores = dna_model.score_variant(
    interval=interval,
    variant=variant,
    variant_scorers=[variant_scorer],
)

# 결과
variant_scores[0].X.shape  # (37, 667) -- 37 genes x 667 tracks
variant_scores[0].obs['gene_name'].tolist()[:10]
# ['RBFOX2', 'APOL4', 'APOL1', 'MYH9', 'TXN2', 'FOXRED2', 'EIF3D', 'APOL3', 'APOL5', 'APOL2']

tidy = variant_scorers.tidy_scores([variant_scores[0]], match_gene_strand=True)
# Tidy Scores Shape: (14652, 19)
```

### 3.6 ISM (In-Silico Mutagenesis)

ISM은 연속 구간 내 모든 위치에서 3가지 대체 염기(원래 염기 제외)를 체계적으로 생성하여 변이 효과를 스캔한다.

**시그니처:**

```python
ism_scores = dna_model.score_ism_variants(
    interval: genome.Interval,               # 전체 서열 구간
    ism_interval: genome.Interval,           # ISM 스캔 구간
    variant_scorers: list[VariantScorer],    # scorer 리스트
)
```

**반환 타입**: ISM 변이 결과 리스트 (각 항목: 변이별 score tuple)

**변이 수 계산**: `ism_interval.width * 3`
- 64bp window: 64 x 3 = 192 변이
- 256bp window: 256 x 3 = 768 변이

**ISM Matrix 변환:**

```python
from alphagenome.interpretation import ism

# score와 variant 추출
ism_matrix_data = ism.ism_matrix(
    variant_scores=[extract_mean_score(score_tuple) for score_tuple in ism_scores],
    variants=[score_tuple[0].uns['variant'] for score_tuple in ism_scores],
)
# 결과: numpy array (4 x width) -- A, C, G, T rows
```

**Quick Start ISM (64bp):**

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

**ISM 256bp 확장 분석:**

```python
# 동일 구간, 256bp ISM window
ism_interval = sequence_interval.resize(256)
# ism_interval: chr20:3753072-3753328:.

ism_scores = dna_model.score_ism_variants(
    interval=sequence_interval,
    ism_interval=ism_interval,
    variant_scorers=[dnase_variant_scorer],
)
# Variants scored: 768
# Elapsed time: 9.96s
```

**ISM 256bp 통계:**

| 통계량 | 값 |
|--------|-----|
| Mean raw_score | -0.0072 |
| Std raw_score | 0.038 |
| Min raw_score | -0.193 |
| Max raw_score | 0.151 |
| Median raw_score | -0.00011 |

### 3.7 Batch Processing Patterns

대량 변이 분석 시 사용되는 배치 처리 패턴이다.

**기본 패턴 (Batch Variant Scoring):**

```python
import pandas as pd
from tqdm import tqdm
from alphagenome.data import genome
from alphagenome.models import dna_client, variant_scorers

# VCF DataFrame에서 변이 반복
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
        interval=interval,
        variant=variant,
        variant_scorers=selected_scorers,
        organism=organism,
    )
    results.append(variant_scores)

# 전체 결과 tidy DataFrame으로 변환
df_scores = variant_scorers.tidy_scores(results)
```

**병렬 처리 패턴 (Analysis Workflow):**

```python
# score_variants: 복수 변이 동시 처리
scores = dna_model.score_variants(
    intervals=eval_df['interval'].to_list(),
    variants=eval_df['variant'].to_list(),
    variant_scorers=[variant_scorers.RECOMMENDED_VARIANT_SCORERS['RNA_SEQ']],
    max_workers=2,
)
```

**Scorer 선택 패턴 (Batch):**

```python
# 대소문자 무관 키 매칭으로 scorer 선택
scorer_selections = {
    'rna_seq': True,
    'cage': True,
    'atac': True,
    'dnase': True,
    'chip_histone': True,
    'polyadenylation': False,
    'splice_sites': True,
    'splice_site_usage': False,
    'splice_junctions': False,
}

all_scorers = variant_scorers.RECOMMENDED_VARIANT_SCORERS
selected_scorers = [
    all_scorers[key]
    for key in all_scorers
    if scorer_selections.get(key.lower(), False)
]
# Batch script: 6개 scorer 사용 (RNA_SEQ, CAGE, ATAC, DNASE, CHIP_HISTONE, SPLICE_SITES)
# + 자동 선택된 _ACTIVE variants
```

**CLI에서 전체 scorer 사용:**

```python
# 19개 scorer 모두 사용
variant_scores = dna_model.score_variant(
    interval=interval,
    variant=variant,
    variant_scorers=list(variant_scorers.RECOMMENDED_VARIANT_SCORERS.values()),
    organism=organism,
)
# CLI: 38,357 score rows from 19 scorers
```

---

## 4. Phase III: Output Analysis

### 4.1 11 OutputTypes

#### 4.1.1 RNA_SEQ

| 항목 | 내용 |
|------|------|
| **설명** | RNA 발현 수준 예측 (Gene expression coverage) |
| **Resolution** | 1 bp |
| **Human Track 수** | 667 |
| **Mouse Track 수** | 173 |
| **Strand** | Strand-specific (+/-) |
| **Scorer 타입** | GeneMaskLFCScorer (gene-level log-fold-change) |
| **사용된 튜토리얼** | Quick Start, Viz Tour, CLI, Batch, Analysis Workflow |

RNA_SEQ는 가장 널리 사용되는 OutputType으로, 유전자 발현 수준을 1bp 해상도로 예측한다. 667개 human track은 다양한 조직과 세포 타입에 걸쳐 있으며, GeneMaskLFCScorer는 유전자 구간 내 예측값의 log-fold-change를 계산한다.

Quick Start 결과:
- CYP2B6 유전자 interval `chr19:40991281-41018398:+`에서 RNA-seq 예측
- 29개 transcript, 10개 gene scored
- Variant Scores shape: (37, 667)

시각화 결과: `results/quick_start/cyp2b6_rna_seq.png`, `results/visualization_tour/01_rna_seq.png`

#### 4.1.2 CAGE

| 항목 | 내용 |
|------|------|
| **설명** | Cap Analysis Gene Expression -- 전사 시작점(TSS) 활성 예측 |
| **Resolution** | 1 bp |
| **Human Track 수** | 546 |
| **Mouse Track 수** | 188 |
| **Strand** | Strand-specific (+/-) |
| **Scorer 타입** | GeneMaskLFCScorer (gene-level) / CenterMaskScorer (center-window) |
| **사용된 튜토리얼** | Quick Start (predict_sequence), Viz Tour, Batch |

CAGE는 FANTOM5 데이터 기반의 TSS 활성을 예측한다. Strand-specific이므로 각 ontology term에 대해 `+` strand와 `-` strand 각각 별도의 track이 생성된다.

Quick Start 결과:
- Brain(+/-) + Lung(+/-) = 4 tracks
- shape: (1048576, 4)

Batch 결과에서 chr16:636337:G>A 변이의 CAGE score가 가장 높은 양성 효과를 보임:
- Neutrophil CAGE: +0.3135 (99.5th percentile)
- Eosinophil CAGE: +0.3001 (99.2nd percentile)

시각화 결과: `results/visualization_tour/02_cage.png`

#### 4.1.3 DNASE

| 항목 | 내용 |
|------|------|
| **설명** | DNase-seq 크로마틴 접근성 예측 |
| **Resolution** | 1 bp |
| **Human Track 수** | 305 |
| **Mouse Track 수** | 67 |
| **Strand** | Unstranded (.) |
| **Scorer 타입** | CenterMaskScorer (center-window comparison) |
| **사용된 튜토리얼** | Quick Start, Viz Tour, Batch, CLI, ISM, Analysis Workflow |

DNASE는 크로마틴 접근성을 예측하며, CenterMaskScorer를 사용하여 변이 중심 주변의 신호 변화를 정량화한다. ISM 분석의 기본 scorer로도 사용된다.

Quick Start 결과:
- Lung tissue: (1048576, 1) shape
- ISM에서 CenterMaskScorer(DNASE, width=501, DIFF_MEAN) 사용

Batch 결과에서 chr16:1135446:G>T의 DNASE score:
- Inflammatory macrophage: -0.5527 (DNase accessibility 감소)
- Common myeloid progenitor: -0.5501

시각화 결과: `results/visualization_tour/03_dnase.png`, `results/variant_scoring_cli/plot_dnase.png`

#### 4.1.4 ATAC

| 항목 | 내용 |
|------|------|
| **설명** | ATAC-seq 크로마틴 접근성 예측 |
| **Resolution** | 1 bp |
| **Human Track 수** | 167 |
| **Mouse Track 수** | 18 |
| **Strand** | Unstranded (.) |
| **Scorer 타입** | CenterMaskScorer |
| **사용된 튜토리얼** | Viz Tour, Batch |

ATAC은 DNASE와 유사하게 크로마틴 접근성을 측정하지만, 다른 실험 프로토콜(Tn5 transposase)을 사용한다. Track 수가 DNASE보다 적지만 (167 vs 305), 특히 EFO 기반 세포주에서 높은 커버리지를 보인다 (40 CLO tracks).

Batch 결과: chr3:58394738:A>T에서 ATAC score +0.047 (EFO:0007950, 76.8th percentile)

시각화 결과: `results/visualization_tour/04_atac.png`

#### 4.1.5 CHIP_HISTONE

| 항목 | 내용 |
|------|------|
| **설명** | ChIP-seq 히스톤 변형 마커 예측 |
| **Resolution** | 1 bp |
| **Human Track 수** | 1,116 |
| **Mouse Track 수** | 183 |
| **Strand** | Unstranded (.) |
| **Scorer 타입** | CenterMaskScorer (width=2001) |
| **사용된 튜토리얼** | Viz Tour, Batch, Analysis Workflow |

CHIP_HISTONE은 가장 많은 track을 가진 OutputType 중 하나이다 (1,116). H3K27AC, H3K36ME3, H3K4ME1, H3K4ME3, H3K9AC, H3K27ME3 등의 히스톤 마커를 예측한다.

Viz Tour에서 색상 코딩:
```python
histone_to_color = {
    'H3K27AC':  '#e41a1c',  # Red - Active enhancer
    'H3K36ME3': '#ff7f00',  # Orange - Gene body
    'H3K4ME1':  '#377eb8',  # Blue - Enhancer
    'H3K4ME3':  '#984ea3',  # Purple - Active promoter
    'H3K9AC':   '#4daf4a',  # Green - Active chromatin
    'H3K27ME3': '#ffc0cb',  # Pink - Repressive
}
```

Batch 결과: chr16:636337:G>A에서 CD14+ monocyte ChIP-histone -0.1265 (99.4th percentile negative)
- CenterMaskScorer width=2001 (DNASE의 501보다 넓음)

시각화 결과: `results/visualization_tour/05_chip_histone.png`

#### 4.1.6 CHIP_TF

| 항목 | 내용 |
|------|------|
| **설명** | ChIP-seq 전사인자(TF) 결합 프로파일 예측 |
| **Resolution** | 1 bp |
| **Human Track 수** | 1,617 |
| **Mouse Track 수** | 127 |
| **Strand** | Unstranded (.) |
| **Scorer 타입** | CenterMaskScorer |
| **사용된 튜토리얼** | ChIP-TF Analysis |

CHIP_TF는 가장 많은 human track (1,617)을 가진 OutputType이다. 대부분(1,447)이 EFO 기반 세포주에서 측정된 데이터이다. CTCF, RAD21, POLR2A, EP300 등 다양한 전사인자의 결합 프로파일을 예측한다.

ChIP-TF Analysis에서 분석된 세포주와 TF:
- K562 (`EFO:0002067`): CTCF, RAD21, POLR2A, EP300
- HepG2 (`EFO:0001187`): CTCF
- K562 + HepG2 합산: 845 tracks across 2 cell types

핵심 분석:
1. K562 CTCF binding profile
2. HepG2 CTCF binding profile
3. CTCF-RAD21 co-localization (cohesin complex 구성요소 분석)
4. Multi-TF comparison (4개 TF 동시 비교)

시각화 결과:
- `results/chip_tf_analysis/chip_tf_k562_ctcf.png`
- `results/chip_tf_analysis/chip_tf_hepg2_ctcf.png`
- `results/chip_tf_analysis/chip_tf_ctcf_rad21_coloc.png`
- `results/chip_tf_analysis/chip_tf_multi_tf.png`

#### 4.1.7 SPLICE_SITES

| 항목 | 내용 |
|------|------|
| **설명** | Splice site 확률 예측 (donor/acceptor) |
| **Resolution** | 1 bp |
| **Human Track 수** | 4 |
| **Mouse Track 수** | 4 |
| **Strand** | Strand-specific (+/-) |
| **Scorer 타입** | SpliceSitesScorer |
| **사용된 튜토리얼** | Viz Tour, Batch |

SPLICE_SITES는 4개의 고정 track (donor+, donor-, acceptor+, acceptor-)으로 구성된다. Ontology term에 독립적이며, 순수하게 서열 기반으로 splice site 확률을 예측한다.

시각화 결과: `results/visualization_tour/06_splice.png` (SPLICE_JUNCTIONS와 함께)

#### 4.1.8 SPLICE_SITE_USAGE

| 항목 | 내용 |
|------|------|
| **설명** | 조직별 splice site 사용 빈도 예측 |
| **Resolution** | 1 bp |
| **Human Track 수** | 734 |
| **Mouse Track 수** | 180 |
| **Strand** | Strand-specific (+/-) |
| **Scorer 타입** | RECOMMENDED_VARIANT_SCORERS에 포함 |
| **사용된 튜토리얼** | (메타데이터 분석) |

SPLICE_SITE_USAGE는 SPLICE_SITES와 달리 조직/세포 특이적인 splice site 사용 패턴을 예측한다. 734개의 human track이 다양한 ontology term에 걸쳐 있다.

#### 4.1.9 SPLICE_JUNCTIONS

| 항목 | 내용 |
|------|------|
| **설명** | Splice junction (exon-exon connection) 예측 |
| **Resolution** | Variable (junction-dependent) |
| **Human Track 수** | 367 |
| **Mouse Track 수** | 90 |
| **Strand** | Strand-specific (+/-) |
| **Scorer 타입** | RECOMMENDED_VARIANT_SCORERS에 포함 |
| **사용된 튜토리얼** | Viz Tour |

SPLICE_JUNCTIONS는 Sashimi plot으로 시각화되며, exon 간 연결 패턴을 보여준다.

Viz Tour에서 APOL4 유전자 주변 splice junction 분석:
```python
plot_components.Sashimi(
    output.splice_junctions
    .filter_to_strand('-')
    .filter_by_tissue('Colon_Transverse'),
    ylabel_template='SPLICE_JUNCTIONS: {biosample_name} ({strand})',
)
```

시각화 결과: `results/visualization_tour/06_splice.png`

#### 4.1.10 CONTACT_MAPS

| 항목 | 내용 |
|------|------|
| **설명** | 3D 크로마틴 상호작용 (Hi-C) 맵 예측 |
| **Resolution** | Variable (bin-dependent) |
| **Human Track 수** | 28 |
| **Mouse Track 수** | 8 |
| **Strand** | N/A (2D matrix) |
| **Scorer 타입** | RECOMMENDED_VARIANT_SCORERS에 포함 |
| **사용된 튜토리얼** | Viz Tour |

CONTACT_MAPS는 Hi-C 데이터를 기반으로 한 3D 크로마틴 구조를 예측한다. 28개의 human track은 주로 EFO 기반 세포주(27개)에서 측정된 데이터이다.

Viz Tour에서 HCT116 세포주 contact map 시각화:
```python
plot_components.ContactMaps(
    tdata=output.contact_maps,
    ylabel_template='{biosample_name}\n{name}',
    cmap='autumn_r',
    vmax=1.0,
)
```

시각화 결과: `results/visualization_tour/07_contact_maps.png`

#### 4.1.11 PROCAP

| 항목 | 내용 |
|------|------|
| **설명** | RNA Polymerase II 활성 (Precision Run-On) 예측 |
| **Resolution** | 1 bp |
| **Human Track 수** | 12 |
| **Mouse Track 수** | - (미지원) |
| **Strand** | Strand-specific (+/-) |
| **Scorer 타입** | GeneMaskLFCScorer |
| **사용된 튜토리얼** | PROCAP Visualization |

PROCAP은 12개 human track (6개 세포주 x 2 strand)으로 구성된 가장 작은 OutputType이다. RNA Polymerase II의 TSS 근처 활성을 base-pair 해상도로 측정한다.

**지원 세포주 (6개):**

| 세포주 | Ontology Code | 설명 |
|--------|---------------|------|
| A673 | EFO:0002106 | Ewing sarcoma |
| Caco-2 | EFO:0001099 | 대장선암 |
| K562 | EFO:0002067 | 만성 골수성 백혈병 |
| Calu3 | EFO:0002819 | 폐선암 |
| MCF 10A | EFO:0001200 | 유방 상피 |
| HUVEC | CL:0002618 | 제대 정맥 내피세포 |

PROCAP Visualization 결과:
- 12 tracks 성공적으로 조회
- 6개 세포주 모두 데이터 확인
- 시각화: `results/procap_visualization/procap.png`

### 4.2 Variant Scoring Results

#### 4.2.0 Scoring Pipeline Overview

Variant scoring은 두 가지 수준에서 이루어진다:

**Level 1: predict_variant (raw predictions)**
- REF 서열과 ALT 서열 각각에 대해 모든 track의 raw 예측값 반환
- 사용자가 직접 차이를 계산하고 해석해야 함
- 시각화에 적합 (REF vs ALT overlay)

**Level 2: score_variant (aggregated scores)**
- 사전 정의된 scorer가 자동으로 REF-ALT 차이를 계산
- Gene-level, center-window 등 다양한 집계 방식
- Quantile normalization으로 게놈 전체 분포 대비 상대적 위치 제공
- AnnData -> tidy DataFrame 변환으로 분석/필터링에 적합

```python
# Level 1: Raw prediction
variant_output = dna_model.predict_variant(
    interval=interval, variant=variant,
    requested_outputs=[dna_client.OutputType.RNA_SEQ],
    ontology_terms=['UBERON:0001157'],
)
# variant_output.reference.rna_seq.values  -> (1048576, N) numpy array
# variant_output.alternate.rna_seq.values  -> (1048576, N) numpy array

# Level 2: Scored prediction
variant_scores = dna_model.score_variant(
    interval=interval, variant=variant,
    variant_scorers=[variant_scorers.RECOMMENDED_VARIANT_SCORERS['RNA_SEQ']],
)
# variant_scores[0].X  -> (n_genes, n_tracks) score matrix
# variant_scores[0].layers['quantiles']  -> quantile-normalized scores
```

#### 4.2.1 Batch Variant Scoring

5개 변이를 6개 scorer 유형으로 분석하여 121,550개 score row를 생성했다.

**변이별 결과:**

| 변이 | Score 행 수 | Mean Raw Score | Std | Unique Ontology |
|------|------------|----------------|-----|-----------------|
| chr3:58394738:A>T | 12,430 | +0.00044 | 0.0049 | 687 |
| chr8:28520:G>C | 10,450 | +0.000029 | 0.0017 | 687 |
| chr16:636337:G>A | 40,150 | +0.0034 | 0.0185 | 687 |
| chr16:1135446:G>T | 41,734 | -0.00084 | 0.0103 | 687 |
| chr1:100000:C>G | 16,786 | +0.00068 | 0.0045 | 687 |
| **합계** | **121,550** | | | |

**고영향 변이 필터링 (|raw_score| > 0.01):**

| 카테고리 | 임계값 | 수 | 비율 |
|---------|--------|-----|------|
| High impact | |score| > 0.01 | 4,509 | 3.7% |
| Very high impact | |score| > 0.05 | 1,423 | 1.2% |
| Extreme impact | |score| > 0.1 | 546 | 0.4% |

**Assay 타입별 분포:**
- ChIP-Histone: 1,718 high impact scores (38%)
- CAGE: 962 scores (21%)
- DNase: 807 scores (18%)
- RNA-seq: 684 scores (15%)
- ATAC: 338 scores (7%)

**방향성 분석:**
- 양성 효과: 3,079 (68.3%)
- 음성 효과: 1,430 (31.7%)

**가장 강한 양성 효과:**
- chr16:636337:G>A / Neutrophil / CAGE: +0.3135 (99.5th percentile)

**가장 강한 음성 효과:**
- chr16:1135446:G>T / CD14+ monocyte / ChIP-Histone: -0.6535 (-99.7th percentile)

**생물학적 해석:**

두 개의 chr16 변이가 가장 강한 효과를 보인다는 것은 주목할 만하다:

1. **chr16:636337:G>A**: 면역 세포(특히 myeloid lineage)에서 전사 시작점 활성을 강하게 증가시킨다. CAGE score가 neutrophil에서 +0.3135로, 99.5th percentile에 해당한다. 이는 이 변이가 면역 관련 유전자의 프로모터 활성을 강화할 가능성을 시사한다.

2. **chr16:1135446:G>T**: 동일한 myeloid 세포 계열에서 chromatin accessibility와 histone modification을 감소시킨다. CD14+ monocyte에서 ChIP-histone score -0.6535는 극단적인 값으로, 규제 요소(regulatory element)의 비활성화를 의미한다.

두 변이 모두 myeloid 세포에 특이적이라는 점은, chr16 해당 구간에 면역 세포 분화/기능에 중요한 regulatory element가 집중되어 있음을 시사한다.

**Biosample 타입별 변동성:**

| Biosample 타입 | High Impact 수 | Std Raw Score |
|----------------|---------------|---------------|
| Primary cells | 1,300 | 0.074 |
| Tissues | 1,681 | 0.055 |
| Cell lines | 1,246 | 0.052 |
| In vitro differentiated | 282 | 0.039 |

Primary cell이 가장 높은 변동성(std=0.074)을 보이는 것은, 세포 상태에 따른 규제 효과의 context-dependency를 반영한다.

결과 파일:
- `results/batch_variant_scoring/variant_scores.csv` (121,550 rows, ~30 MB)
- `results/batch_variant_scoring/high_impact_variants.csv` (4,509 rows, ~1.1 MB)
- `results/batch_variant_scoring/variant_scores_summary.json`

#### 4.2.2 CLI Variant Scoring

단일 변이 `chr22:36201698:A>C`를 19개 전체 scorer로 분석했다.

| 항목 | 값 |
|------|-----|
| Variant | chr22:36201698:A>C |
| Total scores | 38,357 |
| Unique scorers | 19 |
| Mean raw score | 98.91 |
| Median raw score | 0.00106 |
| Std raw score | 581.63 |
| Min raw score | -6.876 |
| Max raw score | 22,558.5 |

> 참고: Mean raw score가 98.91로 높은 이유는 POLYADENYLATION, _ACTIVE scorer 등 다양한 스케일의 scorer가 포함되었기 때문이다. Median (0.00106)이 더 대표적인 중심 경향이다.

**POLYADENYLATION Scorer 결과:**

`PolyadenylationScorer()`는 19개 scorer 중 하나로 실행되었으며, APOL4 유전자(variant chr22:36201698:A>C)에 대해 396개의 score row를 생성했다.

| 항목 | 값 |
|------|-----|
| Total rows | 396 |
| Gene | APOL4 (ENSG00000100336) |
| Unique ontology curies | 285 |
| Mean raw score | 0.417 |
| Min raw score | 0.004 |
| Max raw score | 0.866 |
| Mean quantile score | 0.998 |
| Min quantile score | 0.788 |
| Max quantile score | 1.000 |

**Biosample type 분포:**
- Tissue: 188 rows
- Primary cells: 90 rows
- Cell lines: 73 rows
- In vitro differentiated: 29 rows

높은 quantile score(평균 0.998, 최소 0.788)는 이 변이가 거의 모든 세포 타입에서 polyadenylation에 유의미한 영향을 미친다는 것을 나타낸다. 특히 tissue 샘플에서 가장 많은 track이 분석되었으며, 다양한 조직 및 세포 계열에 걸쳐 광범위한 영향력을 시사한다.

> 참고: PolyadenylationScorer는 RECOMMENDED_VARIANT_SCORERS에 포함되어 있으며, output_type 필드는 RNA_SEQ로 표시되지만 별도의 OutputType enum은 존재하지 않는다(scorer-only).

시각화 결과:
- `results/variant_scoring_cli/plot_rna_seq.png`
- `results/variant_scoring_cli/plot_dnase.png`

### 4.3 ISM Analysis Results

#### 4.3.1 Quick Start ISM (64bp)

| 항목 | 값 |
|------|-----|
| 구간 | chr20:3753000-3753400 (16KB로 resize) |
| ISM window | 64 bp |
| 변이 수 | 192 (64 x 3) |
| Scorer | CenterMaskScorer(DNASE, width=501, DIFF_MEAN) |

#### 4.3.2 ISM 256bp 확장 분석

| 항목 | 값 |
|------|-----|
| Sequence interval | chr20:3745008-3761392:. |
| ISM interval | chr20:3753072-3753328:. |
| ISM window | 256 bp |
| 변이 수 | 768 (256 x 3) |
| 소요 시간 | 9.96초 |
| Scorer | CenterMaskScorer(DNASE, width=501, DIFF_MEAN) |

**Score 통계:**

| 통계 | 값 |
|------|-----|
| Mean | -0.0072 |
| Std | 0.038 |
| Min | -0.193 |
| Max | 0.151 |
| Median | -0.00011 |

ISM heatmap (`results/ism_256bp/ism_heatmap.png`)은 SeqLogo 컴포넌트를 사용하여 생성되었으며, 각 위치에서의 DNASE 접근성 변화를 시각화한다. 음수 값(mean=-0.0072)은 해당 구간의 대부분의 변이가 DNASE 접근성을 약간 감소시키는 경향을 보임을 의미한다.

결과 파일:
- `results/ism_256bp/ism_scores.csv` (768 rows)
- `results/ism_256bp/ism_heatmap.png`
- `results/ism_256bp/results.json`

### 4.4 Tissue Ontology Coverage

Tissue Ontology Mapping 튜토리얼에서 생성된 커버리지 분석 결과이다.

**Human Ontology 시스템별 고유 term 수:**

| Ontology | 고유 term 수 | 설명 |
|----------|-------------|------|
| CL | 227 | Cell types |
| CLO | 52 | Cell lines |
| EFO | 198 | Experimental factors |
| NTR | 18 | New term requests |
| UBERON | 209 | Anatomical structures |

> 참고: 위 수치는 `ontology_summary.json`의 `human_ontology_counts` (metadata에서 추출한 고유 ontology curie 수)이며, ontology_coverage.csv의 OutputType별 총합(5,559)과는 다른 관점의 집계이다.

**Mouse Ontology:**

| Ontology | 고유 term 수 |
|----------|-------------|
| CL | 57 |
| EFO | 24 |
| NTR | 12 |
| UBERON | 85 |

결과 파일:
- `results/tissue_ontology/ontology_coverage.csv`
- `results/tissue_ontology/track_counts.csv`
- `results/tissue_ontology/ontology_terms.json`
- `results/tissue_ontology/ontology_summary.json`

---

## 5. Tutorial-by-Tutorial Execution Results

### 5.1 Quick Start (`results/quick_start/`)

**소스**: `tutorials/quick_start.ipynb`
**스크립트**: `results/quick_start/run_quick_start.py`

Quick Start 튜토리얼은 AlphaGenome API의 전체 workflow를 순차적으로 시연한다.

**실행 단계:**

1. **모델 생성**: `dna_client.create(api_key)` -- Success
2. **OutputType 열거**: 11가지 OutputType 확인
3. **predict_sequence (단일)**: DNASE for Lung -- shape (1048576, 1)
4. **predict_sequence (복수)**: CAGE + DNASE for Lung + Brain
   - DNASE shape: (1048576, 2)
   - CAGE shape: (1048576, 4)
5. **predict_interval**: CYP2B6 RNA-seq in liver
   - CYP2B6 Interval: `chr19:40991281-41018398:+`
   - RNA-seq shape: (1048576, 3)
   - 29 transcripts in interval
6. **predict_variant**: chr22:36201698:A>C RNA-seq in colon
   - REF shape: (1048576, 3)
   - ALT shape: (1048576, 3)
7. **score_variant**: RNA_SEQ scorer
   - Variant Scores shape: (37, 667) -- 37 genes x 667 tracks
   - 10 genes scored: RBFOX2, APOL4, APOL1, MYH9, TXN2, FOXRED2, EIF3D, APOL3, APOL5, APOL2
   - Tidy Scores shape: (14652, 19)
8. **ISM**: 64bp window, CenterMaskScorer DNASE
   - 192 variants scored
9. **Mouse prediction**: DNase for Lung
   - shape: (1048576, 1)

**생성 파일:**

| 파일 | 설명 |
|------|------|
| `results.json` | 전체 실행 결과 메타데이터 |
| `cyp2b6_rna_seq.png` | CYP2B6 RNA-seq 시각화 |
| `variant_effect.png` | REF vs ALT RNA-seq overlay |
| `variant_scores.csv` | Tidy scores (14,652 rows) |
| `run_quick_start.py` | 실행 스크립트 |

### 5.2 Essential Commands (`results/essential_commands/`)

**소스**: `tutorials/essential_commands.ipynb`
**스크립트**: `results/essential_commands/run_essential_commands.py`

Essential Commands 튜토리얼은 AlphaGenome의 핵심 데이터 구조를 테스트한다. API 호출 없이 로컬에서 실행되는 순수 데이터 구조 조작이다.

**3가지 테스트 영역:**

1. **Interval Operations**
   - 생성, center, width, resize
   - overlaps, contains, intersect

2. **Variant Operations**
   - SNV, Insertion, Deletion 생성
   - reference_interval
   - reference_overlaps, alternate_overlaps

3. **TrackData Operations**
   - TrackData 생성 (values + metadata)
   - Resolution 변경 (downsample/upsample)
   - Strand filtering (positive, negative, unstranded)
   - Resize (crop/pad)
   - Slicing (by position, by interval)
   - Track selection (by name, by index)
   - Reverse complement

**생성 파일:**

| 파일 | 설명 |
|------|------|
| `results.json` | 모든 연산 결과 (JSON) |
| `run_essential_commands.py` | 실행 스크립트 |
| `README.md` | 튜토리얼 설명 |

### 5.3 Visualization Modality Tour (`results/visualization_tour/`)

**소스**: `tutorials/visualization_modality_tour.ipynb`
**스크립트**: `scripts/run_visualization_tour.py`

7가지 OutputType에 대한 시각화를 생성했다. 모든 시각화는 chr22 APOL 유전자 클러스터 주변 interval `chr22:35677410-36725986:.`에서 수행되었다.

**시각화 결과 (7개 모두 성공):**

| # | OutputType | 파일 | Ontology | 크기 |
|---|-----------|------|----------|------|
| 1 | RNA_SEQ | `01_rna_seq.png` | Sigmoid colon, Transverse colon | 121 KB |
| 2 | CAGE | `02_cage.png` | Sigmoid colon, Transverse colon | 60 KB |
| 3 | DNASE | `03_dnase.png` | Transverse colon, Sigmoid colon | 61 KB |
| 4 | ATAC | `04_atac.png` | Transverse colon, Sigmoid colon | 59 KB |
| 5 | CHIP_HISTONE | `05_chip_histone.png` | Colon (3 UBERON terms) | 566 KB |
| 6 | SPLICE | `06_splice.png` | Colon (2 UBERON terms) | 81 KB |
| 7 | CONTACT_MAPS | `07_contact_maps.png` | HCT116 (EFO:0002824) | 519 KB |

**주요 시각화 기법:**

```python
# TranscriptAnnotation: 유전자 구조 표시
plot_components.TranscriptAnnotation(transcripts)

# Tracks: 예측 결과 track 표시
plot_components.Tracks(tdata=output.rna_seq, ylabel_template='...')

# OverlaidTracks: REF vs ALT overlay
plot_components.OverlaidTracks(
    tdata={'REF': ref_data, 'ALT': alt_data},
    colors={'REF': 'dimgrey', 'ALT': 'red'},
)

# Sashimi: Splice junction arc diagram
plot_components.Sashimi(output.splice_junctions)

# ContactMaps: Hi-C 상호작용 맵
plot_components.ContactMaps(tdata=output.contact_maps, cmap='autumn_r')

# SeqLogo: ISM heatmap
plot_components.SeqLogo(scores=ism_matrix_data, scores_interval=ism_interval)

# VariantAnnotation: 변이 위치 표시
plot_components.VariantAnnotation([variant])
```

**생성 파일:**

| 파일 | 설명 |
|------|------|
| `01_rna_seq.png` ~ `07_contact_maps.png` | 7개 시각화 |
| `results.json` | 실행 결과 메타데이터 |
| `README.md` | 시각화 가이드 |
| `EXECUTION_SUMMARY.md` | 실행 요약 |

### 5.4 Batch Variant Scoring (`results/batch_variant_scoring/`)

**소스**: `tutorials/batch_variant_scoring.ipynb`
**스크립트**: `scripts/run_batch_variant_scoring.py`

5개 변이를 6개 scorer 유형 (RNA_SEQ, CAGE, ATAC, DNASE, CHIP_HISTONE, SPLICE_SITES + 자동 선택된 _ACTIVE variants)으로 분석했다.

**입력 변이:**

| Variant ID | Chromosome | Position | REF | ALT |
|-----------|-----------|----------|-----|-----|
| chr3_58394738_A_T_b38 | chr3 | 58394738 | A | T |
| chr8_28520_G_C_b38 | chr8 | 28520 | G | C |
| chr16_636337_G_A_b38 | chr16 | 636337 | G | A |
| chr16_1135446_G_T_b38 | chr16 | 1135446 | G | T |
| chr1_100000_C_G_b38 | chr1 | 100000 | C | G |

**핵심 결과:**
- 총 121,550 score rows
- 4,509 high impact rows (|raw_score| > 0.01, 3.7%)
- 5 variants x 687 unique ontology curies x 5 scorer types
- Runtime: ~7초

**가장 영향력 있는 변이**: chr16:636337:G>A
- 2,424 high impact scores (전체의 53.7%)
- 면역 세포(neutrophil, eosinophil)에서 CAGE score 강한 양성 효과
- 전사 시작점 활성 증가 시사

**가장 억제적인 변이**: chr16:1135446:G>T
- 1,133 high impact scores
- Myeloid 세포에서 chromatin accessibility 감소
- CD14+ monocyte ChIP-histone: -0.6535 (극단적 음성 효과)

**생성 파일:**

| 파일 | 크기 | 설명 |
|------|------|------|
| `variant_scores.csv` | ~30 MB | 전체 결과 (121,550 rows) |
| `high_impact_variants.csv` | ~1.1 MB | 고영향 필터 (4,509 rows) |
| `variant_scores_summary.json` | ~9 KB | 변이별 요약 통계 |
| `README.md` | ~5 KB | 방법론 문서 |
| `ANALYSIS_SUMMARY.md` | ~6 KB | 생물학적 해석 |
| `EXECUTION_LOG.md` | ~6 KB | 실행 로그 |

### 5.5 Example Analysis Workflow (`results/analysis_workflow/`)

**소스**: `tutorials/example_analysis_workflow.ipynb`
**스크립트**: `scripts/run_analysis_workflow.py`

TAL1 유전자 주변의 T-ALL (T-cell Acute Lymphoblastic Leukemia) 관련 비코딩 변이를 분석하는 실전 워크플로우이다.

**TAL1 Locus:**
- 구간: `chr1:47212072-47239296` (TAL1 주변)
- TAL1 Interval: `chr1:47209255-47242023:-` (음성 가닥)
- Ontology: `CL:0001059` (CD34+ common myeloid progenitor)

**변이 데이터:**
- 32개 oncogenic variants (실제 T-ALL 환자 데이터, 4개 연구)
  - Mansour_2014: 8 variants (Jurkat, MOLT-3, patients)
  - Liu_2020: 3 variants
  - Liu_2017: 18 variants
  - Smith_2023: 3 variants (new 3' enhancer)
- 96개 background variants (oncogenic variant당 3개 랜덤 생성)
- **총 128개 variants**

**분석 결과:**

| 메트릭 | 값 |
|--------|-----|
| Total variants | 128 |
| Oncogenic variants | 32 |
| Background variants | 96 |
| Mean oncogenic TAL1 effect | +0.334 |
| Mean background TAL1 effect | -0.020 |

Oncogenic variant의 평균 TAL1 효과(+0.334)가 background variant의 평균(-0.020)보다 현저히 높다는 것은, T-ALL oncogenic 변이가 TAL1 발현을 유의미하게 증가시킴을 보여준다.

**Oncogenic Variant 상세:**

변이들은 크게 4개 위치 그룹으로 나뉜다:

| 그룹 | 위치 | 변이 유형 | 변이 수 | 설명 |
|------|------|---------|--------|------|
| MUTE site | chr1:47239290-47239296 | Insertion (1-18bp) | 26 | 주요 TAL1 enhancer insertion site |
| Intergenic SNV | chr1:47230639 | SNV (C>T) | 2 | 유전자간 SNV |
| 3' enhancer 1 | chr1:47212072 | Insertion (21bp) | 1 | 새로운 3' enhancer |
| 3' enhancer 2 | chr1:47212074 | Insertion (6bp) | 1 | 새로운 3' enhancer |

MUTE site (chr1:47239296)에서의 insertion이 가장 많으며(26개), 이는 TAL1 promoter 근처에 반복적으로 insertion이 발생하여 enhancer를 생성하는 것으로 해석된다.

**Background Variant 생성 전략:**

```python
def generate_background_variants(variant, max_number=100):
    """Oncogenic variant의 alternate_bases 길이와 동일한 랜덤 서열 생성"""
    nucleotides = np.array(list('ACGT'), dtype='<U1')
    # alternate_bases 길이와 동일한 길이의 랜덤 서열 생성
    # oncogenic alternate와 동일한 서열은 제외
    ...
```

이 방법은 동일 위치, 동일 길이의 insertion이지만 서열이 다른 "neutral" 변이를 생성하여 비교 대조군으로 사용한다. Oncogenic effect (+0.334)과 background effect (-0.020)의 명확한 차이는 서열 특이적 효과를 입증한다.

**시각화 (plotnine/ggplot2 스타일):**

분석 워크플로우는 plotnine (Python ggplot2)을 사용하여 각 변이 그룹별 rain plot + density plot을 생성한다:

```python
plt_ = (
    gg.ggplot(subplot_df)
    + gg.aes(x='tal1_diff_in_cd34')
    + gg.geom_col(...)     # Rain plot (individual variant scores)
    + gg.geom_density(...)  # Background distribution
    + gg.facet_wrap('~output + plot_group', nrow=1, scales='free_x')
    + gg.scale_fill_manual({True: '#FAA41A', False: 'gray'})
    + gg.coord_flip()
)
```

Orange (#FAA41A)는 oncogenic variant, gray는 background variant를 나타낸다.

**생성 파일:**
- `jurkat_variant_effect.png`: Jurkat 세포주 변이의 RNA-seq, DNase, ChIP-histone 효과
- `variant_analysis_results.json`: 정량적 분석 결과
- `comparison_*.png`: Oncogenic vs background variant density 비교 차트

### 5.6 Tissue Ontology Mapping (`results/tissue_ontology/`)

**소스**: `tutorials/tissue_ontology_mapping.ipynb`
**스크립트**: `scripts/tissue_ontology_mapping.py`

AlphaGenome에서 사용되는 조직/세포 ontology 체계를 탐색하고 매핑하는 튜토리얼이다.

**실행 단계:**

1. **Output Metadata 탐색**: Human 5,563 tracks, Mouse 1,038 tracks
2. **Ontology Term 추출**: 5가지 ontology 시스템 (CL, CLO, EFO, NTR, UBERON)
3. **Tissue 검색**: brain, liver, heart, lung, T cell, neuron 검색
4. **Ontology Coverage 분석**: OutputType x Ontology 교차표 생성
5. **요약 저장**: Track count 테이블, ontology coverage matrix

**핵심 수치:**
- Human 총 ontology term: 5,559 (CL=1,394, CLO=65, EFO=2,402, NTR=93, UBERON=1,605)
- Mouse 총 ontology types: CL, EFO, NTR, UBERON (CLO 없음)

**생성 파일:**

| 파일 | 설명 |
|------|------|
| `track_counts.csv` | OutputType별 Human/Mouse track 수 |
| `ontology_coverage.csv` | OutputType x Ontology 교차표 |
| `ontology_terms.json` | 전체 ontology term 목록 (~30 KB) |
| `ontology_summary.json` | 요약 통계 |

### 5.7 Variant Scoring CLI (`results/variant_scoring_cli/`)

**소스**: `tutorials/variant_scoring_ui.ipynb` (Colab 의존적 UI를 CLI로 대체)
**스크립트**: `scripts/run_variant_scoring_cli.py`

`variant_scoring_ui.ipynb`는 Google Colab의 interactive widget에 의존하므로, CLI 기반으로 동일한 기능을 재구현했다.

**CLI 사용법:**

```bash
python scripts/run_variant_scoring_cli.py \
    --chr chr22 --pos 36201698 --ref A --alt C \
    --outputs rna_seq,dnase \
    --ontology UBERON:0001157 \
    --sequence-length 1MB \
    --interval-width 32768
```

**옵션:**

| 옵션 | 기본값 | 설명 |
|------|--------|------|
| `--chr` | (필수) | 염색체 |
| `--pos` | (필수) | 위치 (1-based) |
| `--ref` | (필수) | 참조 대립유전자 |
| `--alt` | (필수) | 대체 대립유전자 |
| `--outputs` | `rna_seq,dnase` | 출력 타입 (쉼표 구분) |
| `--ontology` | None | Ontology terms |
| `--sequence-length` | `1MB` | 서열 길이 |
| `--interval-width` | `32768` | 시각화 구간 너비 |
| `--no-visualize` | False | 시각화 건너뛰기 |
| `--filter-strand` | None | 특정 strand 필터 |

**실행 결과:**
- Variant: chr22:36201698:A>C
- Interval: chr22:35677410-36725986:.
- 19 scorers 전체 사용 (RECOMMENDED_VARIANT_SCORERS)
  - 11개 base scorer: ATAC, CAGE, CHIP_HISTONE, CHIP_TF, CONTACT_MAPS, DNASE, PROCAP, RNA_SEQ, SPLICE_SITES, SPLICE_SITE_USAGE, SPLICE_JUNCTIONS
  - 1개 scorer-only: POLYADENYLATION (396 rows)
  - 7개 _ACTIVE variant: ATAC_ACTIVE, CAGE_ACTIVE, CHIP_HISTONE_ACTIVE, CHIP_TF_ACTIVE, DNASE_ACTIVE, PROCAP_ACTIVE, RNA_SEQ_ACTIVE
- 38,357 score rows 생성
- 2개 시각화 (RNA_SEQ, DNASE)

**생성 파일:**

| 파일 | 설명 |
|------|------|
| `variant_scores.csv` | 전체 scores (~10 MB, 38,357 rows) |
| `variant_summary.json` | 요약 통계 |
| `plot_rna_seq.png` | RNA-seq REF vs ALT overlay |
| `plot_dnase.png` | DNase REF vs ALT overlay |

---

## 6. Additional Scripts

### 6.1 ISM 256bp Analysis (`scripts/run_ism_256bp.py`)

Quick Start의 64bp ISM을 256bp로 확장한 분석이다.

**구현 상세:**

```python
# 서열 구간 정의 (동일 위치)
sequence_interval = genome.Interval('chr20', 3_753_000, 3_753_400)
sequence_interval = sequence_interval.resize(dna_client.SEQUENCE_LENGTH_16KB)
# 결과: chr20:3745008-3761392:.

# ISM 구간 (64bp -> 256bp 확장)
ism_interval = sequence_interval.resize(256)
# 결과: chr20:3753072-3753328:.

# CenterMaskScorer 설정
dnase_variant_scorer = variant_scorers.CenterMaskScorer(
    requested_output=dna_client.OutputType.DNASE,
    width=501,
    aggregation_type=variant_scorers.AggregationType.DIFF_MEAN,
)

# ISM 실행
ism_scores = dna_model.score_ism_variants(
    interval=sequence_interval,
    ism_interval=ism_interval,
    variant_scorers=[dnase_variant_scorer],
)
# 768 variants scored in 9.96s
```

**ISM Matrix 시각화:**

```python
from alphagenome.interpretation import ism

# Score 및 variant 추출
ism_matrix_data = ism.ism_matrix(
    variant_scores=[extract_mean_score(s) for s in ism_scores],
    variants=[s[0].uns['variant'] for s in ism_scores],
)

# SeqLogo 시각화
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

**결과 요약:**
- 768 variants (256bp x 3 alternates)
- Score 분포: mean=-0.0072, std=0.038, range=[-0.193, 0.151]
- 결과 디렉토리: `results/ism_256bp/`

### 6.2 PROCAP Visualization (`scripts/run_procap_visualization.py`)

PROCAP OutputType을 6개 세포주에 걸쳐 시각화하는 스크립트이다.

**구현 상세:**

```python
# PROCAP 세포주 정의
PROCAP_CELL_LINES = {
    'A673': 'EFO:0002106',
    'Caco-2': 'EFO:0001099',
    'K562': 'EFO:0002067',
    'Calu3': 'EFO:0002819',
    'MCF 10A': 'EFO:0001200',
    'HUVEC': 'CL:0002618',
}

# 예측 수행
output = dna_model.predict_interval(
    interval=interval,   # chr22:35677410-36725986:.
    requested_outputs={dna_client.OutputType.PROCAP},
    ontology_terms=list(PROCAP_CELL_LINES.values()),
)

# Track 정보
output.procap.values.shape[-1]  # 12 tracks (6 cell lines x 2 strands)
```

**결과:**
- 12 tracks 성공 조회
- 6개 세포주 모두 데이터 확인 (A673, Caco-2, K562, Calu3, MCF 10A, HUVEC)
- 시각화: `results/procap_visualization/procap.png`

### 6.3 ChIP-TF Analysis (`scripts/run_chip_tf_analysis.py`)

CHIP_TF OutputType을 활용한 전사인자 결합 프로파일 분석이다.

**4가지 분석:**

1. **K562 CTCF Binding Profile**: K562 세포주에서 CTCF 전사인자 결합 위치 시각화
2. **HepG2 CTCF Binding Profile**: HepG2 세포주에서 동일 분석
3. **CTCF-RAD21 Co-localization**: K562에서 CTCF와 RAD21(cohesin subunit) 공동 위치 분석
4. **Multi-TF Comparison**: CTCF, RAD21, POLR2A, EP300 4개 TF 동시 비교

**구현 핵심:**

```python
# TF 필터링 함수
def filter_to_tfs(chip_tf_data, tf_names):
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
    ontology_terms=['EFO:0002067'],  # K562
)

# CTCF 필터링
k562_ctcf = filter_to_tfs(output_k562.chip_tf, ['CTCF'])

# CTCF-RAD21 co-localization
k562_ctcf_rad21 = filter_to_tfs(output_k562.chip_tf, ['CTCF', 'RAD21'])
```

**결과:**
- 4개 plot 모두 성공
- K562 available TFs: CTCF, EP300, POLR2A, RAD21
- 845 tracks across K562 + HepG2
- 결과 디렉토리: `results/chip_tf_analysis/`

### 6.4 Variant Scoring CLI (`scripts/run_variant_scoring_cli.py`)

Section 5.7에서 상세 기술. `tutorials/variant_scoring_ui.ipynb`의 Colab 의존 UI를 완전한 CLI 도구로 재구현했다.

**핵심 차별점:**
- argparse 기반 CLI 인터페이스
- 19개 전체 scorer 사용 (Batch의 6개 대비, POLYADENYLATION 포함)
- 시각화 옵션 (--no-visualize)
- Strand 필터 옵션
- Custom output directory 지원

---

## 7. API Reference Summary

### 7.1 핵심 메서드 시그니처

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

### 7.2 Data Structure 시그니처

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
# Mutation: .resize_inplace(width)
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
    width: int,                    # 비교 window 크기 (501 for DNASE, 2001 for CHIP_HISTONE)
    aggregation_type: variant_scorers.AggregationType,
)

# AggregationType enum
variant_scorers.AggregationType.DIFF_MEAN
variant_scorers.AggregationType.DIFF_LOG2_SUM

# tidy_scores: AnnData -> DataFrame 변환
variant_scorers.tidy_scores(
    scores: list[list[AnnData]] | list[AnnData],
    match_gene_strand: bool = False,
) -> pd.DataFrame
```

### 7.3 Performance Notes

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

**Rate Limiting**: API에 rate limit이 존재하며, 연속 호출 간 1-2초 delay를 권장한다.

### 7.4 Error Handling Patterns

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

# try-except 패턴 (CLI 스크립트)
try:
    # API 호출
    ...
except Exception as e:
    print(f"Error: {e}", file=sys.stderr)
    traceback.print_exc()
    # 부분 결과 저장
    with open(results_file, 'w') as f:
        json.dump({'error': str(e), 'partial_results': results}, f)
    return 1
```

### 7.5 Visualization Components Reference

AlphaGenome의 `plot_components` 모듈은 다양한 시각화 컴포넌트를 제공한다.

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
    components: list,                 # 컴포넌트 리스트 (위에서 아래 순서)
    interval: genome.Interval,       # 표시 구간
    annotations: list = None,        # 추가 주석 (VariantAnnotation 등)
    title: str = None,               # 그림 제목
    fig_width: float = None,         # 그림 너비 (인치)
    despine_keep_bottom: bool = False,  # 하단 축만 유지
)
```

**ylabel_template 포맷 문자열:**

metadata의 컬럼명을 중괄호로 참조한다:
- `'{biosample_name} ({strand})'` -- 일반적인 track
- `'{biosample_name}\n{histone_mark}'` -- ChIP-histone
- `'CHIP_TF: {biosample_name}\n{transcription_factor} ({strand})'` -- ChIP-TF
- `'{name} ({strand})'` -- SPLICE_SITES

### 7.6 Known Limitations

1. **서열 길이 제한**: 정확히 지원되는 길이(16KB, 100KB, 500KB, 1MB) 중 하나여야 함
2. **PROCAP 제한**: 마우스 미지원, 6개 세포주만 사용 가능
3. **CONTACT_MAPS 제한**: 28개 human track만 존재 (대부분 EFO 기반 세포주)
4. **POLYADENYLATION**: 대응하는 OutputType 없이 scorer로만 존재
5. **Colab 의존성**: variant_scoring_ui.ipynb는 Google Colab widget에 의존 (CLI 대체 필요)
6. **Strand 매칭**: tidy_scores에서 match_gene_strand=True 사용 시 strand 방향이 일치하는 track만 포함
7. **ISM 규모**: 256bp ISM은 768 variants로 관리 가능하지만, 1000bp 이상은 비용/시간이 급증
8. **Ontology term 의존성**: 잘못된 ontology code를 지정하면 빈 TrackData가 반환될 수 있으며, 에러 메시지 없이 shape[-1] == 0인 결과가 반환됨
9. **API Rate Limit**: 연속 호출 시 rate limiting이 적용되므로, 대량 분석 시 time.sleep() 사용 권장
10. **Organism 호환성**: 일부 scorer는 특정 organism만 지원하므로, SUPPORTED_ORGANISMS 확인 필요

### 7.7 Best Practices

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

# 개별 호출 (비권장 - 3배 느림)
# output_rna = dna_model.predict_interval(..., requested_outputs={RNA_SEQ})
# output_dna = dna_model.predict_interval(..., requested_outputs={DNASE})
# output_chip = dna_model.predict_interval(..., requested_outputs={CHIP_HISTONE})
```

**2. Scorer 선택 전략:**

```python
# 탐색 단계: 소수의 핵심 scorer 사용
exploratory_scorers = ['RNA_SEQ', 'DNASE', 'ATAC']

# 검증 단계: 전체 scorer 사용
comprehensive_scorers = list(variant_scorers.RECOMMENDED_VARIANT_SCORERS.values())

# 특정 분야 집중: 관련 scorer만 선택
splicing_scorers = ['SPLICE_SITES', 'SPLICE_SITE_USAGE', 'SPLICE_JUNCTIONS']
chromatin_scorers = ['DNASE', 'ATAC', 'CHIP_HISTONE', 'CHIP_TF']
expression_scorers = ['RNA_SEQ', 'CAGE', 'PROCAP']
```

**3. 결과 필터링 전략:**

```python
# 고영향 변이 필터링 (|raw_score| > 0.01)
high_impact = df_scores[abs(df_scores['raw_score']) > 0.01]

# Quantile 기반 필터링 (상위/하위 5%)
extreme = df_scores[abs(df_scores['quantile_score']) > 0.95]

# 조직 특이적 필터링
brain_scores = df_scores[
    df_scores['ontology_curie'].str.startswith('UBERON:0000955')
]

# Scorer 타입별 필터링
dnase_scores = df_scores[
    df_scores['variant_scorer'].str.contains('DNASE')
]
```

**4. TrackData Null Safety:**

```python
# 항상 None 체크 수행
if output.procap is not None and output.procap.values.shape[-1] > 0:
    # 정상 처리
    fig = plot_components.plot([...])
else:
    print("No data available for requested outputs/ontology_terms")
```

---

## Appendix A: Reverse-Engineered Architecture

### A.1 gRPC-based Client-Server Pattern

AlphaGenome은 gRPC 기반의 client-server 아키텍처를 사용한다.

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

### A.2 Request/Response Flow

```
1. Client: dna_model.predict_interval(interval, outputs, ontology_terms)
   │
   ├─ 2. Client: interval -> 서열 좌표 직렬화
   │
   ├─ 3. gRPC Request: (coordinates, outputs, ontology_terms, organism)
   │                    -> Google Cloud AlphaGenome Server
   │
   ├─ 4. Server: 참조 게놈에서 서열 조회
   │
   ├─ 5. Server: DNA Encoder -> 1MB 서열 인코딩
   │
   ├─ 6. Server: 요청된 OutputType별 Decoder 실행
   │       ├─ RNA_SEQ Decoder -> (1048576, N_tracks) predictions
   │       ├─ DNASE Decoder   -> (1048576, N_tracks) predictions
   │       └─ ...
   │
   ├─ 7. Server: ontology_terms 필터 적용 (track 선택)
   │
   ├─ 8. gRPC Response: serialized predictions + metadata
   │
   └─ 9. Client: ModelOutput 객체 구성
          ├─ output.rna_seq -> TrackData(values, metadata, interval)
          ├─ output.dnase   -> TrackData(values, metadata, interval)
          └─ ...
```

### A.3 Scorer Algorithm Internals

#### GeneMaskLFCScorer

```
Input: REF predictions, ALT predictions, gene annotation
│
├─ 1. Gene mask 생성: gene 구간 내 positions = 1, 외부 = 0
├─ 2. REF gene signal = sum(REF * gene_mask)
├─ 3. ALT gene signal = sum(ALT * gene_mask)
├─ 4. LFC = log2(ALT_signal / REF_signal)
│
└─ Output: gene x track matrix (AnnData)
```

#### CenterMaskScorer

```
Input: REF predictions, ALT predictions, center position, width
│
├─ 1. Center window 정의: [center - width/2, center + width/2]
├─ 2. REF window signal = aggregate(REF[window])
├─ 3. ALT window signal = aggregate(ALT[window])
├─ 4. Aggregation:
│     ├─ DIFF_MEAN: mean(ALT[window]) - mean(REF[window])
│     └─ DIFF_LOG2_SUM: log2(sum(ALT[window])) - log2(sum(REF[window]))
│
└─ Output: 1 x track matrix (AnnData)
```

#### SpliceSitesScorer

```
Input: REF splice_sites predictions, ALT splice_sites predictions
│
├─ 1. REF splice probabilities (donor+, donor-, acceptor+, acceptor-)
├─ 2. ALT splice probabilities
├─ 3. Difference = ALT_prob - REF_prob
│
└─ Output: position-level splice site probability changes
```

### A.4 TrackData Internal Structure

TrackData는 AlphaGenome의 핵심 데이터 컨테이너이다.

```
TrackData
├── values: np.ndarray          # shape: (positions, n_tracks)
│   ├── dtype: float32
│   └── positions = interval.width / resolution
│
├── metadata: pd.DataFrame      # n_tracks rows
│   ├── name: str              # track 이름
│   ├── strand: str            # '+', '-', '.'
│   ├── ontology_curie: str    # 'UBERON:0002048', 'EFO:0002067', ...
│   ├── biosample_name: str    # 'lung', 'K562', ...
│   ├── biosample_type: str    # 'tissue', 'cell line', 'primary cell', ...
│   ├── data_source: str       # 'encode', 'fantom', 'roadmap', ...
│   ├── nonzero_mean: float    # 비영점 평균 (data quality 지표)
│   └── [assay-specific columns]
│       ├── histone_mark: str           # CHIP_HISTONE only
│       ├── transcription_factor: str   # CHIP_TF only
│       └── Assay title: str            # 실험 방법
│
├── resolution: int             # bp per value (None if unknown)
│
└── interval: genome.Interval   # 대응하는 게놈 구간
```

**Resolution과 Position 관계:**

```python
# resolution=1 (1bp per value): 대부분의 OutputType
positions = interval.width  # e.g., 1048576 for 1MB

# resolution=128 (128bp per value): CONTACT_MAPS 등
positions = interval.width // 128  # e.g., 8192 for 1MB

# Resolution 변환
tdata_lowres = tdata.change_resolution(resolution=2)  # Downsample
tdata_hires = tdata_lowres.change_resolution(resolution=1)  # Upsample
```

**TrackData 산술 연산:**

```python
# Element-wise subtraction (변이 효과 계산)
diff = output.alternate.rna_seq - output.reference.rna_seq
# 결과: TrackData with diff.values = alt.values - ref.values
```

### A.5 AnnData Structure (score_variant 반환값)

```python
# AnnData 구조 (단일 scorer)
adata = variant_scores[0]

adata.X                    # numpy array: (n_genes, n_tracks) score matrix
adata.obs                  # pd.DataFrame: gene metadata
  # columns: gene_name, gene_id, strand, ...
adata.var                  # pd.DataFrame: track metadata
  # columns: ontology_curie, biosample_name, biosample_type, output_type, ...
adata.uns                  # dict: scorer metadata
  # keys: variant, scorer, ...
adata.layers['quantiles']  # numpy array: quantile-normalized scores

# tidy_scores 변환 결과 컬럼 (19개):
# variant_id, scored_interval, gene_name, gene_id, gene_strand,
# ontology_curie, biosample_name, biosample_type, output_type,
# variant_scorer, raw_score, quantile_score, ...
```

---

## Appendix B: Complete File Manifest

### B.1 Results Directories

```
results/
├── quick_start/                     (5 files)
│   ├── results.json                 (3.9 KB)   실행 메타데이터
│   ├── run_quick_start.py           (6.8 KB)   실행 스크립트
│   ├── cyp2b6_rna_seq.png          (130 KB)   CYP2B6 RNA-seq 시각화
│   ├── variant_effect.png           (147 KB)   REF vs ALT overlay
│   └── variant_scores.csv           (3.8 MB)   Tidy scores (14,652 rows)
│
├── essential_commands/              (3 files)
│   ├── results.json                 (2.9 KB)   연산 결과
│   ├── run_essential_commands.py     (7.2 KB)   실행 스크립트
│   └── README.md                    (3.0 KB)   설명서
│
├── visualization_tour/              (11 files)
│   ├── 01_rna_seq.png              (121 KB)   RNA-seq 시각화
│   ├── 02_cage.png                  (60 KB)    CAGE 시각화
│   ├── 03_dnase.png                 (61 KB)    DNase 시각화
│   ├── 04_atac.png                  (59 KB)    ATAC 시각화
│   ├── 05_chip_histone.png          (566 KB)   ChIP-histone 시각화
│   ├── 06_splice.png                (81 KB)    Splice sites/junctions
│   ├── 07_contact_maps.png          (519 KB)   Contact maps
│   ├── results.json                 (1.2 KB)   실행 결과
│   ├── README.md                    (5.0 KB)   가이드
│   └── EXECUTION_SUMMARY.md         (4.8 KB)   실행 요약
│
├── batch_variant_scoring/           (6 files)
│   ├── variant_scores.csv           (30 MB)    전체 결과 (121,550 rows)
│   ├── high_impact_variants.csv     (1.1 MB)   고영향 필터 (4,509 rows)
│   ├── variant_scores_summary.json  (9.1 KB)   변이별 요약
│   ├── README.md                    (5.0 KB)   방법론
│   ├── ANALYSIS_SUMMARY.md          (5.8 KB)   생물학적 해석
│   └── EXECUTION_LOG.md             (6.0 KB)   실행 로그
│
├── analysis_workflow/               (Variable)
│   ├── jurkat_variant_effect.png                Jurkat 변이 효과
│   ├── variant_analysis_results.json            분석 결과
│   └── comparison_*.png                         비교 차트 (7개 그룹)
│
├── tissue_ontology/                 (4 files)
│   ├── track_counts.csv             (352 B)    OutputType별 track 수
│   ├── ontology_coverage.csv        (466 B)    OutputType x Ontology 교차표
│   ├── ontology_terms.json          (30 KB)    전체 ontology term
│   └── ontology_summary.json        (437 B)    요약 통계
│
├── variant_scoring_cli/             (4 files)
│   ├── variant_scores.csv           (10 MB)    전체 scores (38,357 rows)
│   ├── variant_summary.json         (660 B)    요약 통계
│   ├── plot_rna_seq.png             (156 KB)   RNA-seq 시각화
│   └── plot_dnase.png               (80 KB)    DNase 시각화
│
├── ism_256bp/                       (3 files)
│   ├── ism_scores.csv               (31 KB)    ISM scores (768 rows)
│   ├── ism_heatmap.png              (61 KB)    SeqLogo heatmap
│   └── results.json                 (988 B)    실행 결과
│
├── procap_visualization/            (2 files)
│   ├── procap.png                   (205 KB)   PROCAP 시각화
│   └── results.json                 (932 B)    실행 결과
│
└── chip_tf_analysis/                (5 files)
    ├── chip_tf_k562_ctcf.png        (70 KB)    K562 CTCF
    ├── chip_tf_hepg2_ctcf.png       (94 KB)    HepG2 CTCF
    ├── chip_tf_ctcf_rad21_coloc.png (120 KB)   CTCF-RAD21 공위치
    ├── chip_tf_multi_tf.png         (182 KB)   Multi-TF 비교
    └── results.json                 (1.2 KB)   실행 결과
```

### B.2 Scripts

```
scripts/
├── run_batch_variant_scoring.py     (8.2 KB)   Batch variant scoring 파이프라인
├── run_analysis_workflow.py         (16 KB)    TAL1 분석 워크플로우
├── run_visualization_tour.py        (13 KB)    시각화 투어
├── run_variant_scoring_cli.py       (14 KB)    CLI variant scoring
├── run_ism_256bp.py                 (8.7 KB)   ISM 256bp 확장
├── run_procap_visualization.py      (6.6 KB)   PROCAP 시각화
├── run_chip_tf_analysis.py          (12 KB)    ChIP-TF 분석
├── tissue_ontology_mapping.py       (9.2 KB)   Ontology 매핑
├── generate_batch_summary.py        (3.0 KB)   Batch 요약 생성
├── analyze_batch_results.py         (2.5 KB)   Batch 결과 분석
└── verify_install.py                (2.1 KB)   설치 검증
```

### B.3 Tutorials

```
tutorials/
├── quick_start.ipynb                (648 KB)
├── essential_commands.ipynb         (42 KB)
├── visualization_modality_tour.ipynb (3.6 MB)
├── batch_variant_scoring.ipynb      (466 KB)
├── example_analysis_workflow.ipynb   (524 KB)
├── tissue_ontology_mapping.ipynb    (2.2 MB)
└── variant_scoring_ui.ipynb         (26 MB)    (Colab 의존적, CLI 대체)
```

### B.4 Summary Statistics

| 카테고리 | 수량 |
|---------|------|
| Results 디렉토리 | 10 |
| 총 결과 파일 | ~50+ |
| Scripts | 11 |
| Tutorials | 7 |
| PNG 시각화 | ~20+ |
| CSV 데이터 | 6 |
| JSON 메타데이터 | 10+ |
| 총 결과 데이터 크기 | ~46 MB |

---

## Appendix C: Ontology Code Reference

### C.1 자주 사용된 Ontology Codes

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
| UBERON:0002170 | Adipose tissue (지방 조직) | Batch (결과) |
| UBERON:0002171 | Liver (간) | Batch (결과) |
| UBERON:0008953 | Spleen (비장) | Batch (결과) |

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
| EFO:0002059 | HeLa-S3 | Batch (결과) |
| EFO:0002791 | OCI-LY7 | Batch (결과) |
| EFO:0007950 | GM12878 | Batch (결과) |

**CL (Cell Ontology):**

| Code | Name | 사용 튜토리얼 |
|------|------|-------------|
| CL:0001059 | CD34+ common myeloid progenitor | Analysis Workflow |
| CL:0002618 | HUVEC (제대 정맥 내피세포) | PROCAP |
| CL:0000775 | Neutrophil | Batch (결과) |
| CL:0000771 | Eosinophil | Batch (결과) |
| CL:0001054 | CD14+ monocyte | Batch (결과) |
| CL:0000863 | Inflammatory macrophage | Batch (결과) |
| CL:0000134 | Mesenchymal stem cell | Batch (결과) |
| CL:0002551 | Endothelial cell | Batch (결과) |

### C.2 Ontology 검색 결과 (Tissue Ontology Tutorial)

| 검색어 | 매칭 수 | 주요 Ontology Terms |
|--------|---------|-------------------|
| brain | 다수 | UBERON:0000955, CL:0000540, CL:0000117 |
| liver | 다수 | UBERON:0002107, UBERON:0002171, EFO:0001187 |
| heart | 다수 | UBERON:0000948, UBERON:0002082 |
| lung | 다수 | UBERON:0002048, UBERON:0002049 |
| T cell | 다수 | CL:0000084, CL:0000624, CL:0000625 |
| neuron | 다수 | CL:0000540, CL:0000117, CL:0000527 |

---

## Appendix D: Execution Timeline

### D.1 튜토리얼 실행 순서

```
2026-02-04  16:32  tutorials/quick_start.ipynb 다운로드
            16:33  tutorials/essential_commands.ipynb 다운로드
            16:33  tutorials/batch_variant_scoring.ipynb 다운로드
            16:33  tutorials/example_analysis_workflow.ipynb 다운로드
            16:33  tutorials/tissue_ontology_mapping.ipynb 다운로드
            16:33  tutorials/visualization_modality_tour.ipynb 다운로드
            16:33  tutorials/variant_scoring_ui.ipynb 다운로드

2026-02-04  15:26  scripts/verify_install.py 작성

2026-02-05  09:41  Essential Commands 실행
            09:42  Quick Start 실행
            09:45  Tissue Ontology Mapping 실행
            09:45  Visualization Tour 실행 시작
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

### D.2 API 호출 통계

| 스크립트 | API 호출 수 (추정) | 소요 시간 |
|---------|-------------------|----------|
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

*이 보고서는 AlphaGenome API (Google DeepMind, 2025)를 사용하여 생성된 결과를 바탕으로 작성되었습니다. 모든 코드 예시는 실제 실행된 스크립트에서 추출되었으며, 수치는 결과 파일에서 검증되었습니다.*
