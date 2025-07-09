# ADC Drug Name Cleaning & Ontology Normalization

This project standardizes key biomedical fields in a dataset of Antibody-Drug Conjugates (ADCs). We integrate scraped clinical trial data with public biomedical ontologies (HGNC, ChEMBL, BioPortal, Disease Ontology) to construct clean, canonical dictionaries and enrich datasets for downstream analysis.


## 📌 Project Goals

- Normalize messy fields (e.g. drug name, antigen, disease indication)
- Canonicalize to authoritative ontologies (HGNC, BioPortal, ChEMBL)
- Generate modular dictionaries and apply them to real trial data
- Preserve ontological hierarchies (esp. DOID) for structured reference


## ✅ Progress Tracker by Task

### **STEP 1: Identify Fields for Standardization**

Focus columns:
- Drug Name / Drug Alias
- Target Antigen
- Cancer Indication (DO / NCIT)
- Payload
- Linker
- Mechanism of Action
- Company, Trial Design, Biomarker Strategy, etc.

✅ **In Progress** — Columns identified; drug, antigen, and disease mapping underway

---

### **STEP 2: Enumerate Raw Terms**

✅ `drug.py`: expands and deduplicates drug aliases  
✅ `disease.py`: extracts unique DO terms and queries BioPortal  
⬜ Payloads, linkers, biomarkers, trial features pending


### **STEP 3: Build Term Dictionaries**

✅ **Antigen Dictionary** via HGNC + TACA fallback  
✅ **Drug Dictionary** via ChEMBL with all aliases (uppercased + deduplicated)  
✅ **Disease Dictionary** via BioPortal (DOID/NCIT) with canonical labels and synonyms  
✅ **DOID Hierarchy Tree** extracted from `doid.owl` using `owlready2`  
⬜ Fuzzy matching & enrichment for payload/linker to be added


### **STEP 4: Apply Dictionaries to Dataset**

✅ `drug.py`: adds `drugNameChembl` + `mechanismOfActionChembl` to input JSON  
✅ `disease.py`: adds `cancerIndicationLabel`, `Ontology`, `MatchID`, `Synonyms`  
✅ `antigen.py`: enriches each drug with HGNC/TACA ontology and match status  
⬜ Remaining fields will be added (e.g., payload, biomarker, combo)


### **STEP 5: Flag Unknown or Ambiguous Terms**

✅ `antigen.py`: status = `"unknown"` or `"taca_match"`  
✅ `disease.py`: status = `"unknown"` for unmatched DO terms  
✅ `drug.py`: filters for ADCs by molecule_type  
⬜ General-purpose rule-based edge case tagging TBD



### **STEP 6: Fuzzy Matching**

✅ `antigen.py`: fuzzy match via `difflib` on HGNC aliases  
✅ `drug.py`: alias fallback + LRU cache to avoid redundant ChEMBL calls  
⬜ Fuzzy matcher for other text fields to be added



### **STEP 7: Documentation**

✅ This README (updated with all pipelines + deliverables)  
✅ Modular, documented codebase across all scripts  
⬜ `NOTES.md` to describe edge cases, fallback logic, and update guides



## 🧩 Codebase Overview

### `antigen.py`
- Maps antigen names to HGNC or TACA references
- Classifies each hit as canonical, alias match, fuzzy, or unknown
- Output: `aacrArticle_hgnc.json`, plus minimal + unknown JSONs

### `disease.py`
- Loads `doid.owl` with `owlready2`
- Builds cancer-specific ancestor trees for all **leaf-level** disease terms
- Traces complete `paths_to_root` (CURIE + label form)
- Output: `doid_cancer_leaf_paths.json`

### `drug.py`
- Maps drug aliases to ChEMBL
- Extracts:
  - Preferred Name
  - ChEMBL ID
  - Mechanism of Action
  - All Aliases (uppercased, deduplicated)
- Filters for `Antibody drug conjugate` types only
- Enriches each drug entry with `drugNameChembl` + `mechanismOfActionChembl`
- Output: `chembl_drug_dictionary.json`, enriched JSON

## 🗂 Deliverables (In Progress)

| Deliverable               | Format        | Status     |
|---------------------------|---------------|------------|
| Raw term lists            | Excel/CSV     | ✅ Partial |
| Cleaned input JSON        | JSON          | ✅ Enriched |
| Term dictionaries         | JSON / CSV    | ✅ Drug, Antigen, Disease |
| DOID tree/hierarchy       | JSON          | ✅ Leaf paths exported |
| Unknowns tagged           | JSON/CSV cols | ✅ Drug / Disease / Antigen |
| Fuzzy match scripts       | Python        | ✅ Antigen, Drug |
| Final documentation       | Markdown      | ✅ This README |


## 🚀 Next Steps

1. Complete payload/linker/biomarker term dictionaries
2. Apply fuzzy + rule-based logic for those fields
3. Merge all cleaned values into master trial dataset
4. Create filtering/ontology-based browser or dashboard
5. Publish `NOTES.md` with edge cases, fallback rules, update pipeline


## 🧠 Acknowledgments

This project integrates data from:
- 🧬 [HGNC](https://www.genenames.org/)
- 💊 [ChEMBL](https://www.ebi.ac.uk/chembl/)
- 🧠 [BioPortal](https://bioportal.bioontology.org/) — DOID / NCIT
- 🍬 TACA (Tumor-Associated Carbohydrate Antigen) glycan references
- 🧾 [Disease Ontology (doid.owl)](https://github.com/DiseaseOntology/HumanDiseaseOntology)

