# AlphaGenome Tutorial Execution & Analysis

Google DeepMind AlphaGenome API를 활용한 DNA regulatory sequence prediction 튜토리얼 실행, 변이 스코어링 분석, 종합 보고서 생성 프로젝트.

## Setup

```bash
# Install dependencies
uv sync

# Set API key
echo 'ALPHAGENOME_API_KEY=your_key_here' > .env
```

Requires Python >= 3.11 and [uv](https://docs.astral.sh/uv/) package manager.

## Directory Structure

```
tutorials/          7 official AlphaGenome Jupyter notebooks
scripts/            12 standalone Python scripts (extracted from tutorials)
results/            Execution results per script (PNGs, JSONs)
reports/
  detailed_report.md           Monolithic detailed report (2,572 lines)
  detailed-report/             Modular 11-chapter version
  learnlm_report.html         Interactive LearnLM HTML report
examples/           Basic prediction example
```

## Scripts

| Script | Purpose | OutputType |
|--------|---------|------------|
| `run_visualization_tour.py` | 7 OutputType visualizations | RNA_SEQ, CAGE, DNASE, ATAC, CHIP_HISTONE, SPLICE, CONTACT_MAPS |
| `run_batch_variant_scoring.py` | 5-variant batch scoring (121,550 rows) | RNA_SEQ, CAGE, DNASE, ATAC, CHIP_HISTONE |
| `run_variant_scoring_cli.py` | All 19 RECOMMENDED_VARIANT_SCORERS (38,357 rows) | All 11 + POLYADENYLATION |
| `run_chip_tf_analysis.py` | Transcription factor binding (CTCF, RAD21) | CHIP_TF |
| `run_analysis_workflow.py` | TAL1 locus T-ALL disease analysis (128 variants) | Multi |
| `run_ism_256bp.py` | In-silico mutagenesis (768 variants) | DNASE |
| `run_procap_visualization.py` | RNA Pol II activity (6 cell lines) | PROCAP |
| `run_splice_site_usage.py` | Splice site usage for APOL4 gene | SPLICE_SITE_USAGE |
| `tissue_ontology_mapping.py` | Ontology term mapping (5,559 terms) | -- |
| `verify_install.py` | Installation verification | -- |

## Running Scripts

```bash
# Run any script
uv run python scripts/<script_name>.py

# Results are saved to results/<script_name>/
```

All scripts include built-in API rate limiting (`API_DELAY`) and save results to `results/`.

## Results

Large CSV files (variant scores) are gitignored. Regenerate via:

```bash
uv run python scripts/run_batch_variant_scoring.py    # ~30 min
uv run python scripts/run_variant_scoring_cli.py      # ~20 min
```

## Reports

- **Detailed Report**: `reports/detailed_report.md` or modular version in `reports/detailed-report/`
- **LearnLM Report**: Open `reports/learnlm_report.html` in a browser for interactive experience
