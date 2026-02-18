# NCBI to Galaxy: Rice Variant Calling Automation

[![CI](https://github.com/jdetras/ncbl2galaxy/actions/workflows/ci.yml/badge.svg)](https://github.com/jdetras/ncbl2galaxy/actions/workflows/ci.yml)
[![Security](https://github.com/jdetras/ncbl2galaxy/actions/workflows/security.yml/badge.svg)](https://github.com/jdetras/ncbl2galaxy/actions/workflows/security.yml)

Automates discovery of rice-related sequencing runs from NCBI literature and enqueues Galaxy variant-calling workflows with resumable execution and security-focused engineering defaults.

## Highlights
- NCBI PubMed -> SRA -> ENA run resolution pipeline
- Single-end and paired-end workflow routing
- Sample-level grouping (`sample_accession`) with optional per-sample histories
- Resume mode using state file (`.ncbi_to_galaxy_state.json`)
- Retry/backoff for external APIs (NCBI, ENA, Galaxy)
- CI and security checks (CodeQL, Bandit, pip-audit, Gitleaks)

## Repository layout
- `ncbi_to_galaxy.py`: main automation script
- `workflows/rice_variant_calling.ga`: single-end workflow template
- `workflows/rice_variant_calling_paired.ga`: paired-end workflow template
- `docs/ARCHITECTURE.md`: architecture and compatibility matrix
- `SECURITY.md`: vulnerability reporting policy

## Prerequisites
- Python 3.10+
- Galaxy account and API key
- Imported Galaxy workflows from this repo
- Stable reference FASTA choice (recommended: pin one build and URL, e.g. IRGSP-1.0)

## Setup
```bash
python3 -m venv .venv
source .venv/bin/activate
pip install -e .[dev]
pre-commit install
```

## Import Galaxy workflows
1. In Galaxy, go to `Workflows` -> `Import`.
2. Upload:
- `workflows/rice_variant_calling.ga`
- `workflows/rice_variant_calling_paired.ga`
3. Keep names unchanged unless you pass explicit workflow IDs.
4. If tools are missing, remap local equivalents in the workflow editor.

## Run
Dry run (discovery only):
```bash
python3 ncbi_to_galaxy.py --retmax 100 --max-runs 20 --dry-run
```

Full run (single + paired, sample-aware histories, resumable):
```bash
python3 ncbi_to_galaxy.py \
  --query '(rice[Title/Abstract] OR "Oryza sativa"[MeSH Terms])' \
  --retmax 100 \
  --max-runs 50 \
  --galaxy-url 'https://usegalaxy.org' \
  --galaxy-api-key 'YOUR_GALAXY_API_KEY' \
  --single-workflow-name 'Rice Variant Calling (BWA-MEM2 + FreeBayes)' \
  --paired-workflow-name 'Rice Variant Calling Paired (BWA-MEM2 + FreeBayes)' \
  --reference-url 'https://example.org/IRGSP-1.0.fasta' \
  --history-per-sample
```

## Key runtime options
- `--state-file`: path for resume-state tracking
- `--reset-state`: ignore previous state file
- `--history-per-sample`: isolate each biological sample in its own history
- `--reference-dataset-id`: use existing Galaxy dataset instead of `--reference-url`
- `--single-workflow-id` / `--paired-workflow-id`: explicit workflow IDs

## Quality and security checks
Local:
```bash
ruff check .
python -m unittest discover -s tests -p "test_*.py" -v
bandit -r . -x tests,.venv
pip-audit
```

GitHub:
- CI workflow: `.github/workflows/ci.yml`
- Security workflow: `.github/workflows/security.yml`
- Dependabot config: `.github/dependabot.yml`

## Governance
- License: `LICENSE`
- Contributing guide: `CONTRIBUTING.md`
- Security policy: `SECURITY.md`
