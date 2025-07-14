# ADC Ontology Enrichment Pipeline


## Key Features

- **Field Standardization:** Normalizes drug names, target antigens, cancer indications, payloads, linkers, mechanisms of action, company, trial design, and biomarker strategy.
- **Ontology Integration:** Maps terms to HGNC, ChEMBL, Disease Ontology (DOID), TACA, and other resources.
- **Automated Fuzzy Matching:** Applies string similarity and alias matching to maximize coverage and reduce manual review.
- **Comprehensive Dictionaries:** Builds and maintains versioned dictionaries for all key fields, supporting updates and transparency.
- **Automated Analysis & Reporting:** Quantifies match rates, field coverage, and unknown patterns; exports markdown reports with actionable recommendations.
- **Quality Monitoring:** Supports regular, automated quality assessment to guide ongoing improvements.


## Codebase Overview

### Enrichment Modules
- `antigen.py`: Maps antigens to HGNC/TACA, with fuzzy and alias matching.
- `drug.py`: Standardizes drug names and mechanisms using ChEMBL.
- `disease_enhanced.py`: Normalizes cancer indications via DOID/NCIT, with acronym expansion and fuzzy matching.
- `company.py`: Cleans company names using rule-based logic and drug name extraction.
- `trial_design.py`: Categorizes trial designs into comprehensive research-based groups.
- `biomarker_strategy.py`: Classifies biomarker strategies using ADC-specific categories and keyword matching.
- `payload_linker.py`: Standardizes payload and linker fields using ChEMBL and curated dictionaries.

### Pipeline Orchestration
- `unified_enrichment.py`: Runs all enrichment modules, merges outputs, and produces a unified, fully enriched JSON.
- `run_pipeline.py`: Main entry point for executing the full pipeline.

### Analysis & Reporting
- `quick_ontology_analysis.py`: Analyzes final output, calculates match rates, field coverage, and unknown patterns; exports a markdown report (`ontology_analysis_report.md`).
- `export_ontology.py`: Exports selected fields and ontology data for downstream analysis.
- `test_pipeline.py`: Validates pipeline output and provides quality metrics.


## Outputs

- **Enriched JSON:** `aacrArticle_fully_enriched.json` (all fields, all ontology results)
- **Ontology-Only Export:** `aacrArticle_ontology_only.json` (core fields + ontology)
- **Individual Enriched Files:** Per-module minimal outputs for debugging
- **Dictionaries:** Versioned JSONs for each ontology/field
- **Automated Report:** `ontology_analysis_report.md` (key stats, match rates, recommendations)


## Workflow

1. **Prepare Input:** Place raw ADC data in `aacrArticle.json`.
2. **Run Pipeline:**
   ```bash
   python run_pipeline.py
   ```


## References

This project integrates data from:
- [HGNC](https://www.genenames.org/)
- [ChEMBL](https://www.ebi.ac.uk/chembl/)
- [BioPortal](https://bioportal.bioontology.org/) — DOID / NCIT
- [Disease Ontology GitHub](https://github.com/DiseaseOntology/HumanDiseaseOntology)
- TACA literature curation (tumor glycan antigens)


