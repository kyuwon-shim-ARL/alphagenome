# 모듈 01: 환경 설정과 입력 시스템 (Setup & Input Architecture)

## 1. 한 줄 요약

> AlphaGenome은 genome.Interval(0-based)과 genome.Variant(1-based VCF), 5,559 ontology term, 11 OutputType으로 구성된 입력 시스템을 갖추며, 5,563 human track과 1,038 mouse track을 지원한다.

## 2. 왜 이 분석이 필요했나

게놈 DNA 서열이 어떻게 RNA 발현, 크로마틴 구조, 전사인자 결합 등 다양한 regulatory 현상을 만들어내는지 이해하는 것은 현대 생물학의 핵심 과제다. 하지만 이러한 예측을 위해서는 복잡한 데이터 구조와 입력 형식을 정확히 이해해야 한다.

AlphaGenome은 Google DeepMind가 개발한 DNA regulatory sequence prediction 모델로, 최대 1,048,576 bp (1 Mbp) 길이의 게놈 서열을 입력받아 11가지 생물학적 출력(OutputType)을 예측한다. 이 모델을 효과적으로 활용하기 위해서는 다음 질문들에 답할 수 있어야 한다:

- **게놈 좌표를 어떻게 지정하는가?** (Interval vs Variant, 0-based vs 1-based)
- **조직/세포 특이성을 어떻게 표현하는가?** (Ontology system)
- **어떤 DNA 서열 길이를 선택해야 하는가?** (16KB ~ 1MB)
- **어떤 종류의 예측이 가능한가?** (11 OutputType, 5,563 track)
- **예측 결과의 메타데이터는 어떤 정보를 포함하는가?** (output_metadata)

이 모듈은 AlphaGenome의 입력 시스템 전체를 체계적으로 분석하여, 이후 예측 파이프라인과 결과 해석을 위한 토대를 제공한다.

## 3. 분석 과정 (The Mechanics)

### 3.1 Quick Reference Table

AlphaGenome 입력 시스템의 핵심 컴포넌트를 한눈에 파악하기 위한 참조표다.

| 컴포넌트 | 설명 | 주요 필드 | 좌표 체계 |
|---------|------|-----------|----------|
| `genome.Interval` | 게놈 구간 | chromosome, start, end, strand | 0-based half-open |
| `genome.Variant` | 게놈 변이 | chromosome, position, ref, alt | 1-based (VCF) |
| DNA Sequence | 원시 DNA 문자열 | 16KB / 100KB / 500KB / 1MB | N/A |
| Ontology Term | 조직/세포 식별자 | CL / CLO / EFO / NTR / UBERON | N/A |
| OutputType | 예측 유형 | RNA_SEQ, DNASE, CAGE, ... (11종) | N/A |

### 3.2 IPO 요약 (Input-Process-Output)

AlphaGenome 입력 시스템의 전체 데이터 흐름을 IPO 관점에서 정리한다.

| 단계 | 항목 | 내용 |
|------|------|------|
| **Input** | Genome Coordinate | `genome.Interval` (0-based) 또는 `genome.Variant` (1-based) |
|  | DNA Sequence | Raw DNA string (16KB ~ 1MB) 또는 참조 게놈 자동 조회 |
|  | Organism | HOMO_SAPIENS (human) / MUS_MUSCULUS (mouse) |
|  | Ontology Terms | 조직/세포 식별자 리스트 (예: `['UBERON:0002048']`) |
|  | OutputType | 예측 유형 리스트 (예: `[OutputType.RNA_SEQ]`) |
| **Process** | Model Creation | `dna_client.create(api_key)` -- gRPC 클라이언트 초기화 |
|  | Metadata Query | `output_metadata()` -- 사용 가능한 track 정보 조회 |
| **Output** | Track Metadata | DataFrame (5,563 human / 1,038 mouse tracks) |
|  | ModelOutput | 예측 결과 (다음 모듈에서 설명) |

### 3.3 genome.Interval: 게놈 구간 지정

`genome.Interval`은 게놈 상의 연속 구간을 나타내는 dataclass다.

#### 필드 구조

| 필드 | 타입 | 설명 | 필수 | 기본값 |
|------|------|------|------|--------|
| `chromosome` | `str` | 염색체 이름 (예: `'chr1'`) | Yes | - |
| `start` | `int` | 시작 위치 (0-based, inclusive) | Yes | - |
| `end` | `int` | 종료 위치 (0-based, exclusive) | Yes | - |
| `strand` | `str` | 가닥 방향 (`'+'`, `'-'`, `'.'`) | No | `'.'` |
| `name` | `str` | 구간 이름 | No | `None` |

**좌표 체계**: 0-based half-open 구간 (start 포함, end 미포함). 이는 BED 포맷, Python slicing과 동일하다.

#### 주요 메서드

| 메서드 | 설명 | 반환 타입 | 예시 |
|--------|------|-----------|------|
| `resize(width)` | 중심을 기준으로 구간 크기 변경 | `Interval` | `interval.resize(100)` |
| `overlaps(other)` | 두 구간 겹침 여부 확인 | `bool` | `interval.overlaps(second)` |
| `contains(other)` | 다른 구간 완전 포함 여부 | `bool` | `interval.contains(second)` |
| `intersect(other)` | 겹치는 구간 반환 | `Interval` | `interval.intersect(second)` |
| `center()` | 구간 중심 좌표 | `int` | `interval.center()` |
| `width` | 구간 너비 (end - start) | `int` | `interval.width` |

#### 코드 예시

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
# 계산: center=1005, new_start=1005-50=955, new_end=1005+50=1055

# 겹침 검사
second = genome.Interval(chromosome='chr1', start=1_005, end=1_015)
interval.overlaps(second)   # True (1005-1010 겹침)
interval.contains(second)   # False (second가 1010을 넘어감)
interval.intersect(second)  # chr1:1005-1010:.

# Positional arguments도 가능
interval = genome.Interval('chr22', 36_150_498, 36_252_898)

# Strand 지정
interval = genome.Interval(
    chromosome='chr19', start=40991281, end=41018398, strand='+'
)
# 결과: chr19:40991281-41018398:+
```

#### Essential Commands 결과에서 검증된 값

| 연산 | 입력 | 결과 |
|------|------|------|
| 생성 | `Interval('chr1', 1000, 1010)` | `chr1:1000-1010:.` |
| center | `interval.center()` | `1005` |
| width | `interval.width` | `10` |
| resize(100) | `interval.resize(100)` | `chr1:955-1055:.` |
| overlaps | `interval.overlaps(Interval('chr1', 1005, 1015))` | `True` |
| contains | `interval.contains(Interval('chr1', 1005, 1015))` | `False` |
| intersect | `interval.intersect(Interval('chr1', 1005, 1015))` | `chr1:1005-1010:.` |

### 3.4 genome.Variant: 게놈 변이 지정

`genome.Variant`는 게놈 변이를 나타내는 dataclass다.

#### 필드 구조

| 필드 | 타입 | 설명 | 필수 | 기본값 |
|------|------|------|------|--------|
| `chromosome` | `str` | 염색체 이름 | Yes | - |
| `position` | `int` | 변이 위치 (1-based, VCF 호환) | Yes | - |
| `reference_bases` | `str` | 참조 대립유전자 | Yes | - |
| `alternate_bases` | `str` | 대체 대립유전자 | Yes | - |
| `name` | `str` | 변이 이름 | No | `None` |

**위치 체계**: 1-based (VCF 포맷과 동일). **Interval의 0-based와 다르므로 주의가 필요하다.**

#### 변이 유형

| 유형 | 조건 | 예시 | 설명 |
|------|------|------|------|
| SNV | `len(ref) == 1 and len(alt) == 1` | `A>C` | Single Nucleotide Variant |
| Insertion | `len(ref) < len(alt)` | `T>CGTCAAT` | DNA 삽입 |
| Deletion | `len(ref) > len(alt)` | `AGGGATC>C` | DNA 결실 |

#### 주요 속성

| 속성/메서드 | 설명 | 반환 타입 |
|-------------|------|-----------|
| `reference_interval` | 참조 대립유전자 구간 (Interval, 0-based) | `Interval` |
| `reference_overlaps(interval)` | 참조 대립유전자와 구간 겹침 | `bool` |
| `alternate_overlaps(interval)` | 대체 대립유전자와 구간 겹침 | `bool` |

#### 코드 예시

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

# reference_interval (1-based position -> 0-based interval)
# position=10000 (1-based) -> start=9999 (0-based)
# chr3:9999-10000:. (0-based half-open)

# Overlap 검사
interval = genome.Interval(chromosome='chr3', start=10_005, end=10_010)
variant.reference_overlaps(interval)   # False
variant.alternate_overlaps(interval)   # True (insertion이 해당 구간에 겹침)
```

#### Essential Commands 결과에서 검증된 값

| 연산 | 입력 | 결과 |
|------|------|------|
| SNV | `Variant('chr3', 10000, 'A', 'C')` | `chr3:10000:A>C` |
| Insertion | `Variant('chr3', 10000, 'T', 'CGTCAAT')` | `chr3:10000:T>CGTCAAT` |
| Deletion | `Variant('chr3', 10000, 'AGGGATC', 'C')` | `chr3:10000:AGGGATC>C` |
| reference_interval | `variant.reference_interval` | `chr3:9999-10000:.` |
| reference_overlaps | `variant.reference_overlaps(Interval('chr3', 10005, 10010))` | `False` |
| alternate_overlaps | `variant.alternate_overlaps(Interval('chr3', 10005, 10010))` | `True` |

### 3.5 DNA Sequence Input: 서열 길이와 패딩

AlphaGenome은 raw DNA string을 직접 입력으로 받을 수 있다. `predict_sequence()` 메서드에 사용된다.

#### 지원 서열 길이

| 상수명 | 길이 (bp) | 용도 | 추천 사용 케이스 |
|--------|-----------|------|-----------------|
| `SEQUENCE_LENGTH_16KB` | 16,384 | ISM 분석, 빠른 테스트 | ISM, 단기 regulatory element |
| `SEQUENCE_LENGTH_100KB` | 100,000 | 중간 규모 분석 | 일반적 유전자 크기 커버 |
| `SEQUENCE_LENGTH_500KB` | 500,000 | 대규모 분석 | Large locus, enhancer-promoter 상호작용 |
| `SEQUENCE_LENGTH_1MB` | 1,048,576 | 전체 맥락 분석 (권장) | 최대 규제 맥락, 장거리 상호작용 |

서열은 정확히 해당 길이와 일치해야 한다. 짧은 서열은 'N'으로 패딩한다.

#### 코드 예시

```python
from alphagenome.models import dna_client

# 서열 패딩 (Quick Start에서 사용된 방법)
sequence = 'GATTACA'.center(dna_client.SEQUENCE_LENGTH_1MB, 'N')
# 1,048,576 bp 길이의 서열 (중앙에 'GATTACA', 나머지 'N')

# SUPPORTED_SEQUENCE_LENGTHS 딕셔너리로 접근 가능
seq_len = dna_client.SUPPORTED_SEQUENCE_LENGTHS['SEQUENCE_LENGTH_1MB']
# 1048576
```

#### 주의사항

- 서열 길이가 정확히 일치하지 않으면 API 오류 발생
- 'N' 패딩은 해당 위치의 예측에 영향을 미치지 않음
- ISM 분석에는 16KB 길이가 적합 (연산 효율)
- 변이 분석에는 1MB 길이가 권장 (충분한 맥락)

#### 길이 선택 가이드

| 길이 | 맥락 범위 | 적합한 용도 | 비고 |
|------|-----------|------------|------|
| 16KB | ~8KB 양쪽 | ISM, 빠른 프로토타이핑 | score_ism_variants 권장 |
| 100KB | ~50KB 양쪽 | 중규모 유전자 분석 | 일반적 유전자 크기 커버 |
| 500KB | ~250KB 양쪽 | Large locus 분석 | enhancer-promoter 상호작용 포함 |
| 1MB | ~500KB 양쪽 | 전체 맥락 분석 | 최대 규제 맥락 제공, 권장 |

대부분의 튜토리얼과 스크립트에서 1MB (`SEQUENCE_LENGTH_1MB` = 1,048,576)를 사용한다. 이는 장거리 regulatory interaction을 포착하기 위한 것이다. ISM은 예외적으로 16KB를 사용하는데, ISM은 변이 자체가 다수이므로 더 좁은 맥락에서도 충분한 정보를 얻을 수 있다.

### 3.6 Ontology System: 조직/세포 특이성

AlphaGenome은 5가지 ontology 체계를 사용하여 조직/세포 타입을 지정한다.

#### 5가지 Ontology 시스템

| 약칭 | 정식명 | 설명 | Human Term 수 |
|------|--------|------|--------------|
| CL | Cell Ontology | 세포 타입 | 1,394 |
| CLO | Cell Line Ontology | 세포주 | 65 |
| EFO | Experimental Factor Ontology | 실험 조건/세포주 | 2,402 |
| NTR | New Term Request | 신규 용어 | 93 |
| UBERON | Uberon Anatomy Ontology | 해부학적 구조 | 1,605 |
| **합계** | | | **5,559** |

**Ontology Curie 포맷**: `PREFIX:숫자코드`

#### 주요 Ontology 코드 예시

```python
# 튜토리얼에서 사용된 주요 ontology code
ontology_examples = {
    # UBERON (조직/해부학적 구조)
    'UBERON:0002048': 'lung (폐)',
    'UBERON:0000955': 'brain (뇌)',
    'UBERON:0001114': 'right lobe of liver (간 우엽)',
    'UBERON:0001157': 'colon (대장)',
    'UBERON:0001159': 'sigmoid colon (에스상결장)',
    'UBERON:0001155': 'transverse colon (횡행결장)',

    # EFO (세포주)
    'EFO:0002067':   'K562 (만성 골수성 백혈병 세포주)',
    'EFO:0001187':   'HepG2 (간세포암 세포주)',
    'EFO:0002824':   'HCT116 (대장암 세포주)',
    'EFO:0002106':   'A673 (Ewing sarcoma 세포주)',
    'EFO:0001099':   'Caco-2 (대장선암 세포주)',
    'EFO:0002819':   'Calu-3 (폐선암 세포주)',
    'EFO:0001200':   'MCF 10A (유방 상피 세포주)',

    # CL (세포 타입)
    'CL:0002618':    'HUVEC (제대 정맥 내피세포)',
    'CL:0001059':    'CD34+ common myeloid progenitor',
}
```

#### Ontology별 OutputType 커버리지 (Human)

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

#### 주요 패턴

- **CHIP_TF**: EFO 기반 세포주가 압도적 (1,447/1,617 = 89.5%)
- **CHIP_HISTONE**: 가장 많은 track (1,116), 모든 ontology에 고르게 분포
- **CAGE/RNA_SEQ**: UBERON(조직) 기반 track이 많음 (CAGE: 222/546, RNA_SEQ: 283/667)
- **CONTACT_MAPS/PROCAP**: 소규모 track (각 28, 12), 주로 EFO 기반

### 3.7 Output Metadata: Track 정보 조회

`output_metadata()` 메서드는 사용 가능한 모든 track에 대한 메타데이터를 반환한다.

#### 메타데이터 조회 예시

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

#### Human Track 수 (OutputType별)

| OutputType | Human Track 수 | Mouse Track 수 | 비고 |
|------------|----------------|----------------|------|
| ATAC | 167 | 18 | |
| CAGE | 546 | 188 | Strand-specific (x2) |
| DNASE | 305 | 67 | |
| RNA_SEQ | 667 | 173 | Strand-specific (x2) |
| CHIP_HISTONE | 1,116 | 183 | 가장 많은 track |
| CHIP_TF | 1,617 | 127 | 가장 많은 track |
| SPLICE_SITES | 4 | 4 | 고정 track (ontology 독립) |
| SPLICE_SITE_USAGE | 734 | 180 | Strand-specific (x2) |
| SPLICE_JUNCTIONS | 367 | 90 | Strand-specific (x2) |
| CONTACT_MAPS | 28 | 8 | |
| PROCAP | 12 | - | Human only, 6 cell lines x 2 strand |
| **합계** | **5,563** | **1,038** | |

> 참고: Track 수(5,563)와 ontology term 수(5,559)는 일부 track이 같은 ontology term을 공유하므로 다를 수 있다.

#### 메타데이터 필드 예시 (DNase Lung)

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

#### SPLICE_SITES: Ontology 독립 고정 Track

SPLICE_SITES는 4개의 고정 track (donor+, donor-, acceptor+, acceptor-)으로 구성된다. Ontology term에 독립적이며, 순수하게 서열 기반으로 splice site 확률을 예측한다.

```python
# SPLICE_SITES는 ontology_terms 파라미터를 무시하고 항상 4개 track 반환
output = dna_model.predict_interval(
    interval=interval,
    requested_outputs=[dna_client.OutputType.SPLICE_SITES],
    # ontology_terms 불필요
)
# output.splice_sites.values.shape = (1048576, 4)
# 4 tracks: donor+, donor-, acceptor+, acceptor-
```

#### PROCAP: Human Only 제한

PROCAP은 12개 human track (6개 세포주 x 2 strand)으로 구성된 가장 작은 OutputType이다. **마우스는 지원하지 않는다.**

**지원 세포주 (6개):**

| 세포주 | Ontology Code | 설명 |
|--------|---------------|------|
| A673 | EFO:0002106 | Ewing sarcoma |
| Caco-2 | EFO:0001099 | 대장선암 |
| K562 | EFO:0002067 | 만성 골수성 백혈병 |
| Calu3 | EFO:0002819 | 폐선암 |
| MCF 10A | EFO:0001200 | 유방 상피 |
| HUVEC | CL:0002618 | 제대 정맥 내피세포 |

### 3.8 Model Creation: 예측 파이프라인 진입점

AlphaGenome 모델 생성은 `dna_client.create()` 함수를 통해 이루어진다.

#### 모델 생성 코드

```python
from alphagenome.models import dna_client
import os

# 모델 생성
api_key = os.environ.get('ALPHAGENOME_API_KEY')
dna_model = dna_client.create(api_key)
```

`dna_client.create()`는 내부적으로 gRPC 클라이언트를 초기화하고 Google Cloud 기반의 AlphaGenome 서버에 연결한다. 반환되는 `DNAModel` 객체는 모든 예측 메서드의 진입점이다.

#### Organism enum

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

#### OutputType enum (11가지)

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

### 3.9 Essential Commands 튜토리얼 결과

**소스**: `tutorials/essential_commands.ipynb`
**스크립트**: `results/essential_commands/run_essential_commands.py`
**특징**: API 호출 없이 로컬에서 실행되는 순수 데이터 구조 조작

Essential Commands 튜토리얼은 AlphaGenome의 핵심 데이터 구조를 테스트한다.

#### 3가지 테스트 영역

**1. Interval Operations**
   - 생성, center, width, resize
   - overlaps, contains, intersect

**2. Variant Operations**
   - SNV, Insertion, Deletion 생성
   - reference_interval
   - reference_overlaps, alternate_overlaps

**3. TrackData Operations**
   - TrackData 생성 (values + metadata)
   - Resolution 변경 (downsample/upsample)
   - Strand filtering (positive, negative, unstranded)
   - Resize (crop/pad)
   - Slicing (by position, by interval)
   - Track selection (by name, by index)
   - Reverse complement

#### 생성 파일

| 파일 | 설명 |
|------|------|
| `results.json` | 모든 연산 결과 (JSON) |
| `run_essential_commands.py` | 실행 스크립트 |
| `README.md` | 튜토리얼 설명 |

### 3.10 Tissue Ontology Mapping 튜토리얼 결과

**소스**: `tutorials/tissue_ontology_mapping.ipynb`
**스크립트**: `scripts/tissue_ontology_mapping.py`
**목적**: 사용 가능한 ontology term의 전체 커버리지 분석

#### Human Ontology 시스템별 고유 term 수

| Ontology | 고유 term 수 | 설명 |
|----------|-------------|------|
| CL | 227 | Cell types |
| CLO | 52 | Cell lines |
| EFO | 198 | Experimental factors |
| NTR | 18 | New term requests |
| UBERON | 209 | Anatomical structures |

> 참고: 위 수치는 `ontology_summary.json`의 `human_ontology_counts` (metadata에서 추출한 고유 ontology curie 수)이며, ontology_coverage.csv의 OutputType별 총합(5,559)과는 다른 관점의 집계이다.

#### Mouse Ontology

| Ontology | 고유 term 수 |
|----------|-------------|
| CL | 57 |
| EFO | 24 |
| NTR | 12 |
| UBERON | 85 |

#### 생성 파일

| 파일 | 설명 |
|------|------|
| `ontology_coverage.csv` | OutputType별 ontology 커버리지 |
| `track_counts.csv` | OutputType별 track 수 |
| `ontology_terms.json` | 고유 ontology term 리스트 |
| `ontology_summary.json` | 전체 요약 통계 |

## 4. 왜 이 방법인가

### 4.1 0-based vs 1-based: 좌표 체계 선택

AlphaGenome이 Interval에는 0-based, Variant에는 1-based를 사용하는 이유는 생물정보학 표준 포맷과의 호환성 때문이다:

- **genome.Interval (0-based)**: BED 포맷, Python slicing과 동일 → 구간 연산이 직관적
- **genome.Variant (1-based)**: VCF 포맷과 동일 → 기존 변이 데이터베이스와 호환

이 설계는 사용자가 이미 익숙한 파일 포맷을 그대로 사용할 수 있게 하면서, 프로그래밍 언어(Python)의 관습도 존중한다.

### 4.2 5가지 Ontology 시스템: 생물학적 복잡성 반영

단일 ontology 대신 5가지 시스템을 사용하는 이유:

1. **CL (Cell Ontology)**: 기본 세포 타입 (T cell, neuron 등)
2. **UBERON (Anatomy)**: 조직/장기 구조 (lung, liver 등)
3. **EFO (Experimental Factor)**: 세포주와 실험 조건 (K562, HepG2 등)
4. **CLO (Cell Line Ontology)**: 특정 세포주 메타데이터
5. **NTR (New Term Request)**: 신규 개념

이 구조는 primary tissue, immortalized cell line, in vitro differentiated cell 등 다양한 생물학적 맥락을 정확히 표현할 수 있게 한다.

### 4.3 서열 길이 선택 (16KB ~ 1MB): 맥락과 효율의 균형

1MB가 권장되는 이유:

- **장거리 enhancer**: 100KB 이상 떨어진 enhancer-promoter 상호작용도 포착 가능
- **TAD (Topologically Associating Domain)**: 일반적 TAD 크기 (~1MB)와 일치
- **3D 크로마틴 구조**: CONTACT_MAPS 예측에 충분한 맥락

ISM이 16KB를 사용하는 이유:

- **변이 수**: 256bp ISM window = 768 변이 → 1MB 맥락에서는 연산 비용 과다
- **로컬 효과**: 단일 염기 변이의 영향은 대부분 ±10KB 내에서 발생
- **효율**: 16KB 서열로 충분한 정보 제공 (Quick Start 검증)

### 4.4 11 OutputType: 다층 regulatory landscape

AlphaGenome이 11가지 OutputType을 지원하는 이유:

| 레이어 | OutputType | 측정 대상 |
|-------|-----------|----------|
| **Transcription** | RNA_SEQ, CAGE, PROCAP | 유전자 발현, TSS 활성 |
| **Chromatin Access** | DNASE, ATAC | 개방 크로마틴 구조 |
| **Histone Marks** | CHIP_HISTONE | 규제 상태 (active/repressive) |
| **TF Binding** | CHIP_TF | 전사인자 결합 프로파일 |
| **RNA Splicing** | SPLICE_SITES, SPLICE_SITE_USAGE, SPLICE_JUNCTIONS | Alternative splicing 패턴 |
| **3D Structure** | CONTACT_MAPS | 크로마틴 상호작용 |

이 다층 구조는 단일 DNA 변이가 어떤 regulatory 레벨에 영향을 미치는지 체계적으로 분석할 수 있게 한다.

## 5. 해석과 한계

### 5.1 주요 발견

1. **좌표 체계 혼용의 복잡성**: Interval (0-based)과 Variant (1-based)를 혼용하면서 `variant.reference_interval`로 자동 변환 제공 → 사용자 편의성 증대하지만 초보자에게는 혼란 가능성

2. **Ontology 커버리지 불균형**: CHIP_TF는 EFO 중심(89.5%), RNA_SEQ는 UBERON 중심(42.4%) → 분석 설계 시 사용 가능한 track 사전 확인 필수

3. **PROCAP Human-only 제한**: 12개 track, 6개 세포주만 지원 → 마우스 연구나 다른 조직 분석에는 대안 필요

4. **서열 길이 제약**: 1MB 고정 크기 → 초대형 유전자 locus (>1MB)는 여러 조각으로 나눠 분석해야 함

### 5.2 한계

#### 기술적 한계

- **API 기반 아키텍처**: 로컬 GPU 불필요하지만 인터넷 연결 필수, 대량 분석 시 API quota 제약
- **고정 서열 길이**: 16KB/100KB/500KB/1MB 중 선택 → 가변 길이 입력 불가
- **좌표 체계 변환 부담**: 외부 데이터 통합 시 0-based ↔ 1-based 변환 오류 가능성

#### 생물학적 한계

- **Ontology 의존성**: 5,559 term 외의 희귀 조직/세포는 가장 유사한 term으로 근사해야 함
- **세포주 편향**: CHIP_TF는 EFO 기반 세포주 중심 → primary tissue 분석에 제약
- **종 제한**: Human, Mouse만 지원 → 다른 모델 생물(zebrafish, C. elegans 등) 불가

### 5.3 실무 권장사항

1. **Interval vs Variant 선택**:
   - 구간 분석 (gene expression across locus) → `genome.Interval` + `predict_interval()`
   - 변이 효과 분석 → `genome.Variant` + `predict_variant()` 또는 `score_variant()`

2. **서열 길이 선택**:
   - 기본값: 1MB (최대 맥락)
   - ISM 분석: 16KB (효율)
   - 빠른 프로토타이핑: 100KB (중간 타협)

3. **Ontology 선택**:
   - Primary tissue 분석 → UBERON 우선 확인
   - 세포주 실험 → EFO 또는 CLO 확인
   - 특정 세포 타입 → CL 확인
   - `output_metadata()`로 사용 가능한 track 사전 확인 필수

4. **메타데이터 활용**:
   - `biosample_type` 필터링으로 tissue vs cell line vs primary cell 구분
   - `nonzero_mean` 필드로 해당 ontology term에서의 baseline 활성도 확인
   - `data_source` (ENCODE, GTEx 등) 확인으로 데이터 품질 평가

## 6. 다음으로 (The Bridge)

> "환경 설정과 입력 시스템을 이해했다. genome.Interval과 genome.Variant로 게놈 좌표와 변이를 정의하고, 5,559 ontology term으로 세포/조직을 지정할 수 있다. 11가지 OutputType과 5,563개 human track의 구조도 파악했다. 그런데 이 입력을 **실제로 어떻게 예측 파이프라인에 전달하고 결과를 받는가?** → [다음: 02-prediction.md](02-prediction.md)"

지금까지 AlphaGenome의 입력 시스템 전체를 살펴봤다:
- genome.Interval (0-based, BED 호환)
- genome.Variant (1-based, VCF 호환)
- DNA sequence (16KB ~ 1MB)
- Ontology system (5,559 term)
- 11 OutputType (5,563 human track)
- Model creation (`dna_client.create()`)

다음 모듈에서는 이 입력들이 어떻게 4가지 예측 메서드 (`predict_sequence()`, `predict_interval()`, `predict_variant()`, `score_variant()`)를 통과하여 실제 regulatory 예측 결과로 변환되는지, 그리고 ISM(In-Silico Mutagenesis)과 배치 처리 패턴을 어떻게 활용하는지 분석한다.
