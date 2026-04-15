# SmartVariant: Automated Genomic Mutation Analyzer

## Overview
SmartVariant analyzes DNA sequence pairs (reference + query), detects variants (SNP/INS/DEL), predicts functional impact, and provides both web UI and REST API. This step includes ClinVar loader, protein effect heuristics, batch API, and Docker deployment.

## Quickstart

1. install dependencies:
   ```bash
   pip install -r requirements.txt
   ```
2. run app:
   ```bash
   python app.py
   ```
3. open: http://127.0.0.1:5000/

## API

- `POST /api/analyze` JSON:
  `{ "reference_seq": "...", "query_seq": "..." }`
- `POST /api/batch` JSON:
  `{ "pairs": [ {"reference_seq":"...","query_seq":"..."}, ... ]}`
- `POST /api/clinvar` form upload field `clinvar_file` (TSV style)

## Docker

```bash
docker build -t smartvariant .
docker run -p 5000:5000 smartvariant
```

## Step 5 (Enterprise Enhancements)

- Added localized ClinVar fetch + VCF ingestion: `utils/clinvar.py` including `download_clinvar_vcf`, `parse_clinvar_vcf`, `variants_to_training_records`
- Enhanced model pipeline for robust feature set: `model/impact_predictor.py` with `is_indel`, `xx > feature`, and `xgboost` optional mode
- Added new API endpoint:
  - `POST /api/clinvar/download` { url, algorithm } to ingest external VCF and retrain
- Added CI pipeline and regression tests for ClinVar pipeline data
- Added full project requirements for reproducibility and venv install:
  - `xgboost`, `requests`, `joblib`, `pytest`

## How to run the new ClinVar endpoint

1. Put your ClinVar VCF URL in the request:
   ```bash
   curl -X POST http://127.0.0.1:5000/api/clinvar/download \
     -H 'Content-Type: application/json' \
     -d '{"url":"https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz","algorithm":"rf"}'
   ```
2. The server downloads and retrains model, returns `training` stats.

