# 🧬 Unified ADC Data Enrichment Pipeline

This unified pipeline orchestrates all individual enrichment scripts to create a comprehensive, standardized dataset with full ontology mappings for ADC (Antibody-Drug Conjugate) research data.

## 🚀 Quick Start

```bash
# Run the complete pipeline
python run_pipeline.py

# Or run the unified script directly
python unified_enrichment.py
```

## 📁 What You Get

### 1. Organized Dictionary Files
All generated dictionaries are organized in the `dictionaries/` folder:

```
dictionaries/
├── antigen/
│   ├── aacrArticle_hgnc.json
│   └── aacrArticle_hgnc_unknowns.json
├── disease/
│   └── doid_cancer_leaf_paths.json
├── drug/
│   └── chembl_drug_dictionary.json
└── payload_linker/
    ├── chembl_payload_dictionary.json
    └── chembl_linker_dictionary.json
```

### 2. Comprehensive Enriched JSON
The final output `aacrArticle_fully_enriched.json` contains:

- **Original data** preserved exactly as input
- **Ontology mappings** for all biomedical entities
- **Standardized structure** for easy analysis
- **Match status** for quality assessment

## 📊 Data Structure

Each drug entry in the enriched JSON now has this structure:

```json
{
  "id": "entry_id",
  "createdAt": "timestamp",
  "title": "paper_title",
  "abstract": "paper_abstract",
  "url": "paper_url",
  "extractedDrugs": [
    {
      // Original fields preserved
      "extractedAt": "timestamp",
      "drugName": "original_drug_name",
      "drugNameConfidence": 1.0,
      "company": "company_name",
      "cancerIndication": ["NSCLC"],
      "targetAntigen": ["EGFR", "cMET"],
      "payload": ["MMAE"],
      "linker": ["valine-citrulline"],
      "phase": "Phase 1",
      
      // NEW: Comprehensive ontology mappings
      "ontology": {
        "drug": {
          "chembl_id": "CHEMBL123456",
          "preferred_name": "Standardized Drug Name",
          "max_phase": 3,
          "mechanism_of_action": ["ADC", "Antibody-drug conjugate"],
          "targets": ["EGFR", "cMET"],
          "indications": ["NSCLC"],
          "match_status": "chembl_match"
        },
        "antigen": {
          "hgnc_symbol": "EGFR",
          "hgnc_id": "HGNC:3236",
          "ensembl_gene_id": "ENSG00000146648",
          "synonyms": ["ERBB", "HER1"],
          "locus_type": "protein-coding gene",
          "gene_group": ["Receptor tyrosine kinases"],
          "taca_subtype": null,
          "taca_family": null,
          "match_status": "hgnc_match"
        },
        "disease": {
          "doid_id": "DOID:3908",
          "doid_label": "non-small cell lung carcinoma",
          "ncit_id": "C2926",
          "ncit_label": "Non-Small Cell Lung Carcinoma",
          "synonyms": ["NSCLC", "non-small cell lung cancer"],
          "hierarchy_path": ["disease", "cancer", "lung cancer"],
          "match_status": "doid_match"
        },
        "payload": {
          "chembl_id": "CHEMBL123456",
          "preferred_name": "Monomethyl auristatin E",
          "max_phase": 4,
          "molecule_type": "Small molecule",
          "match_status": "chembl_match"
        },
        "linker": {
          "chembl_id": "CHEMBL123456",
          "preferred_name": "Valine-citrulline",
          "max_phase": null,
          "molecule_type": "Small molecule",
          "match_status": "chembl_match"
        }
      }
    }
  ]
}
```

## 🔧 Configuration

Edit `pipeline_config.json` to customize:

- **Input/output file paths**
- **Script execution order**
- **Data structure fields**
- **Pipeline settings**

```json
{
  "input_files": {
    "main_json": "aacrArticle.json",
    "hgnc_tsv": "hgnc_complete_set.tsv",
    "taca_json": "taca.json"
  },
  "output_files": {
    "final_enriched": "aacrArticle_fully_enriched.json"
  },
  "scripts": {
    "antigen": {"enabled": true},
    "disease": {"enabled": true},
    "drug": {"enabled": true},
    "payload_linker": {"enabled": true}
  }
}
```

## 📋 Match Status Values

Each ontology field includes a `match_status` indicating quality:

| Status | Meaning | Quality |
|--------|---------|---------|
| `chembl_match` | Exact ChEMBL match | High |
| `hgnc_match` | Exact HGNC match | High |
| `doid_match` | Exact DOID match | High |
| `fuzzy_match` | Fuzzy string match | Medium |
| `alias_match` | Alias/synonym match | Medium |
| `taca_match` | TACA antigen match | Medium |
| `unknown` | No match found | Low |

## 🛠️ Individual Scripts

The pipeline runs these scripts in order:

### 1. `antigen.py`
- Maps antigen names to HGNC database
- Falls back to TACA for tumor antigens
- Outputs: `aacrArticle_hgnc.json`

### 2. `disease.py`
- Processes Disease Ontology (DOID)
- Extracts cancer hierarchy
- Outputs: `doid_cancer_leaf_paths.json`

### 3. `drug.py`
- Queries ChEMBL for ADC drugs
- Adds mechanism of action
- Outputs: `chembl_drug_dictionary.json`

### 4. `payload_linker.py`
- Maps payloads and linkers to ChEMBL
- Filters for small molecules
- Outputs: `chembl_payload_dictionary.json`, `chembl_linker_dictionary.json`

## 📈 Quality Metrics

The pipeline provides several quality indicators:

1. **Match Coverage**: Percentage of terms successfully mapped
2. **Unknown Terms**: Terms that couldn't be mapped
3. **Fuzzy Matches**: Terms matched with <100% confidence
4. **Ontology Distribution**: Which ontologies were most useful

## 🔍 Troubleshooting

### Common Issues

1. **Missing Required Files**
   ```
   ❌ Missing required files: ['hgnc_complete_set.tsv']
   ```
   **Solution**: Download HGNC data from https://www.genenames.org/download/

2. **Script Execution Errors**
   ```
   ❌ antigen.py failed: ModuleNotFoundError
   ```
   **Solution**: Install required packages: `pip install pandas tqdm owlready2 chembl_webresource_client`

3. **API Rate Limits**
   ```
   ❌ drug.py failed: HTTP 429 Too Many Requests
   ```
   **Solution**: Add delays between API calls or use cached results

### Log Files

Check `enrichment_pipeline.log` for detailed error information:

```bash
tail -f enrichment_pipeline.log
```

## 📊 Output Analysis

### Dictionary Statistics
```python
import json

# Load final enriched data
with open('aacrArticle_fully_enriched.json', 'r') as f:
    data = json.load(f)

# Count ontology matches
total_drugs = 0
drug_matches = 0
antigen_matches = 0

for entry in data:
    for drug in entry['extractedDrugs']:
        total_drugs += 1
        ontology = drug.get('ontology', {})
        
        if ontology.get('drug', {}).get('match_status') != 'unknown':
            drug_matches += 1
        if ontology.get('antigen', {}).get('match_status') != 'unknown':
            antigen_matches += 1

print(f"Drug match rate: {drug_matches/total_drugs:.1%}")
print(f"Antigen match rate: {antigen_matches/total_drugs:.1%}")
```

## 🎯 Use Cases

### 1. Drug Discovery Analysis
- Find all ADCs targeting specific antigens
- Analyze payload-linker combinations
- Track clinical development phases

### 2. Target Validation
- Identify most common antigen targets
- Map to gene families and pathways
- Find novel therapeutic targets

### 3. Clinical Trial Design
- Analyze disease indications
- Identify biomarker strategies
- Compare trial designs

### 4. Literature Mining
- Extract drug-entity relationships
- Map to standardized ontologies
- Enable semantic search

## 🤝 Contributing

To add new enrichment scripts:

1. Create your script (e.g., `biomarker.py`)
2. Add to `pipeline_config.json`
3. Update `unified_enrichment.py` to handle your output
4. Test with `python run_pipeline.py`

## 📚 References

- **HGNC**: https://www.genenames.org/
- **ChEMBL**: https://www.ebi.ac.uk/chembl/
- **Disease Ontology**: https://disease-ontology.org/
- **TACA**: Tumor-associated carbohydrate antigens

---

For questions or issues, check the log file or create an issue in the repository. 