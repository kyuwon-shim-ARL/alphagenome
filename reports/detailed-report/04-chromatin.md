# 모듈 04: 크로마틴 접근성 (Chromatin Accessibility — DNASE · ATAC)

## 1. 한 줄 요약

> DNASE(305 tracks)과 ATAC(167 tracks)은 동일한 크로마틴 접근성을 서로 다른 실험 프로토콜로 측정하며, CenterMaskScorer(width=501, DIFF_MEAN)로 변이 중심 주변의 신호 변화를 정량화한다.

---

## 2. 왜 이 분석이 필요했나

- **이전 모듈에서 남은 질문**: 유전자 발현(RNA_SEQ, CAGE, PROCAP)을 확인했다. 유전자가 발현되려면 DNA가 열려있어야 한다.
- **이 모듈이 답하려는 것**: DNA 접근성(open chromatin)은 어떻게 측정하며, 변이가 접근성을 어떻게 바꾸는가?

크로마틴이 열려있는 영역(open chromatin)은 전사인자와 조절 단백질의 결합이 가능한 위치이다. 변이가 크로마틴 접근성을 바꾸면 유전자 조절 기능이 손상되거나 활성화된다.

AlphaGenome은 두 가지 표준 프로토콜로 크로마틴 접근성을 예측한다:
- **DNASE-seq**: DNase I 효소로 열린 크로마틴을 분해하여 측정 (305 human tracks)
- **ATAC-seq**: Tn5 transposase로 열린 크로마틴에 어댑터를 삽입하여 측정 (167 human tracks)

---

## 3. 분석 과정 (The Mechanics)

### Quick Reference

| 항목 | DNASE | ATAC |
|------|-------|------|
| **스크립트** | Quick Start, Viz Tour, Batch, CLI, ISM, Analysis | Viz Tour, Batch |
| **핵심 함수** | CenterMaskScorer(DNASE, 501, DIFF_MEAN) | CenterMaskScorer(ATAC, 501, DIFF_MEAN) |
| **실행 명령** | `dna_model.predict_interval(...)` | 동일 |

### IPO 요약

| 단계 | Input | Process | Output |
|------|-------|---------|--------|
| 예측 | genome.Interval + ontology_terms | predict_interval/predict_variant | TrackData (1048576, N) |
| 스코어링 | Interval + Variant | CenterMaskScorer (width=501) | AnnData (score per track) |
| 시각화 | TrackData | plot_components.Tracks | PNG |

---

### DNASE (305 tracks)

| 항목 | 내용 |
|------|------|
| **설명** | DNase-seq 크로마틴 접근성 예측 |
| **Resolution** | 1 bp |
| **Human Track 수** | 305 |
| **Mouse Track 수** | 67 |
| **Strand** | Unstranded (.) |
| **Scorer 타입** | CenterMaskScorer (center-window comparison, width=501) |
| **사용된 튜토리얼** | Quick Start, Viz Tour, Batch, CLI, ISM, Analysis Workflow (6/11) |

DNASE는 AlphaGenome 튜토리얼에서 **가장 널리 사용되는 OutputType**으로, 11개 튜토리얼/스크립트 중 6개에서 사용된다.

**Quick Start 결과:**
- Lung tissue: (1048576, 1) shape
- ISM에서 CenterMaskScorer(DNASE, width=501, DIFF_MEAN) 기본 scorer로 사용

**Batch 결과 (chr16:1135446:G>T):**

| Ontology | Cell/Tissue Type | Score | 해석 |
|----------|------------------|-------|------|
| CL:0000863 | Inflammatory macrophage | -0.5527 | 강한 접근성 감소 |
| CL:0000049 | Common myeloid progenitor | -0.5501 | 강한 접근성 감소 |

- Negative score: 변이가 크로마틴 접근성을 **감소**시킴 → 전사 억제 가능
- Myeloid 세포에서 특이적으로 강한 효과 → chr16 해당 구간의 면역 regulatory element 존재 시사

**시각화:**
- `results/visualization_tour/03_dnase.png` (61 KB)
- `results/variant_scoring_cli/plot_dnase.png`

**Ontology 분포:**

| Ontology | Track 수 |
|----------|---------|
| CL | 94 |
| CLO | 12 |
| EFO | 95 |
| NTR | 13 |
| UBERON | 91 |
| **합계** | **305** |

---

### ATAC (167 tracks)

| 항목 | 내용 |
|------|------|
| **설명** | ATAC-seq 크로마틴 접근성 예측 |
| **Resolution** | 1 bp |
| **Human Track 수** | 167 |
| **Mouse Track 수** | 18 |
| **Strand** | Unstranded (.) |
| **Scorer 타입** | CenterMaskScorer |
| **사용된 튜토리얼** | Viz Tour, Batch (2/11) |

ATAC은 DNASE와 유사하게 크로마틴 접근성을 측정하지만, 다른 실험 프로토콜(Tn5 transposase)을 사용한다. Track 수가 DNASE보다 적지만(167 vs 305), 특히 CLO 기반 세포주에서 높은 커버리지를 보인다(40 CLO tracks).

**Batch 결과 (chr3:58394738:A>T):**
- ATAC +0.047 (EFO:0007950 GM12878, 76.8th percentile)
- Positive score: 크로마틴 접근성 증가

**시각화:**
- `results/visualization_tour/04_atac.png` (59 KB)

**Ontology 분포:**

| Ontology | Track 수 |
|----------|---------|
| CL | 15 |
| CLO | 40 |
| EFO | 65 |
| NTR | 3 |
| UBERON | 44 |
| **합계** | **167** |

---

### DNASE vs ATAC 비교

| 항목 | DNASE | ATAC |
|------|-------|------|
| **프로토콜** | DNase I 효소 (크로마틴 분해) | Tn5 transposase (어댑터 삽입) |
| **Human tracks** | 305 | 167 |
| **Mouse tracks** | 67 | 18 |
| **CenterMask width** | 501 | 501 |
| **Tutorial coverage** | 6/11 (55%) | 2/11 (18%) |
| **Ontology: CL** | 94 | 15 |
| **Ontology: CLO** | 12 | **40** |
| **Ontology: EFO** | 95 | 65 |
| **Ontology: NTR** | 13 | 3 |
| **Ontology: UBERON** | 91 | 44 |

**공통점:**
- 동일한 생물학적 현상 측정 (open chromatin regions)
- 동일한 Scorer 전략 (CenterMaskScorer, width=501)
- 1bp resolution, unstranded

**차이점:**
- DNASE: 더 많은 트랙 (305), 더 높은 primary tissue 커버리지 (UBERON 91 vs 44)
- ATAC: 더 높은 cell line 커버리지 (CLO 40 vs 12), 최신 프로토콜
- DNASE가 AlphaGenome의 기본 크로마틴 접근성 분석 도구

---

## 4. 왜 이 방법인가

| 고려한 방법 | 채택 여부 | 이유 |
|------------|----------|------|
| DNASE + ATAC 둘 다 | ✅ 채택 | 상보적 — DNASE는 primary tissue, ATAC은 cell line 커버리지 |
| DNASE만 | ❌ | CLO 커버리지 부족 (12 tracks) |
| CenterMaskScorer width=501 | ✅ 채택 | 변이 주변 국소 크로마틴 상태 포착, 노이즈 감소 |
| ISM에서 DNASE 기본 scorer | ✅ 채택 | 305개 track으로 범용적 regulatory region 탐지 |

---

## 5. 해석과 한계

**실전 해석 가이드:**

| Score 범위 | 해석 |
|-----------|------|
| < -0.5 | Strong loss of accessibility → 조절 손상 가능성 높음 |
| -0.5 ~ -0.1 | Moderate loss → 추가 증거 필요 |
| -0.1 ~ +0.1 | Neutral → 크로마틴 접근성에 영향 없음 |
| +0.1 ~ +0.5 | Moderate gain → 조절 활성화 가능성 |
| > +0.5 | Strong gain → ectopic activation 위험 |

**한계:**
1. **실험 조건 의존성**: 크로마틴 접근성은 세포 유형/발달 단계에 따라 변함. 305+167 tracks로도 모든 조직/시점 커버 불가
2. **CenterMaskScorer 한계**: 501bp 윈도우 = 국소 영역만 고려. 장거리 enhancer-promoter loop은 CONTACT_MAPS 필요
3. **DNASE-ATAC 불일치 가능성**: 동일 세포에서도 프로토콜 차이로 신호 강도 다를 수 있음

---

## 6. 다음으로 (The Bridge)

> 크로마틴이 열려있는 위치를 DNASE(305)와 ATAC(167)으로 확인했다. 열린 크로마틴 위에 **히스톤 변형 마커와 전사인자가 어떻게 결합**하는가? CHIP_HISTONE(1,116)과 CHIP_TF(1,617)로 살펴볼 차례이다.

→ [다음: 05-epigenomics.md](05-epigenomics.md)
