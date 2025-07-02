# 🧪 ADC Drug Ontology Normalization Pipeline

This project implements a comprehensive pipeline for normalizing ADC (Antibody-Drug Conjugate) drug data using multiple biomedical APIs to standardize drug names, target antigens, and disease indications.

## 📊 Data Source
- **Input**: `ADC drugs and trials_Apr2025.xlsx - OG List.csv`
- **Format**: CSV with header row (skipped during import)

---

## 🧪 FIELD 1: Drug Name / Drug Alias

**Goal:** Canonicalize drug names and retrieve comprehensive metadata including mechanism of action, indications, and target information.

### ✅ Implemented API:

**[ChEMBL WebResource Client](https://chembl.gitbook.io/chembl-interface-documentation/web-services)** – Comprehensive drug database with full metadata.

### 🔧 Key Features:

* **Comprehensive Drug Information**: ChEMBL ID, preferred names, synonyms, development phase
* **Mechanism of Action**: Detailed MOA data with target information
* **Drug Indications**: EFO and MeSH mappings for approved uses
* **Target Metadata**: Full target information including UniProt and HGNC cross-references
* **Drug Safety**: Withdrawal status, black box warnings, ATC classifications

### 🛠️ Implementation:

```python
from chembl_webresource_client.new_client import new_client

def fetch_full_chembl_data(drug_name):
    # Search molecule
    results = molecule_client.search(drug_name)
    if not results:
        return {'Found': False, ...}
    
    # Extract comprehensive metadata:
    # - Molecule Info (ChEMBL ID, preferred name, max phase, etc.)
    # - Mechanism of Action (MOA, action type, target info)
    # - Drug Indications (EFO/MeSH mappings)
    # - Target Metadata (UniProt, HGNC cross-references)
```

**Data Retrieved:**
- **Molecule Info**: ChEMBL ID, preferred name, max phase, approval date, drug type
- **Mechanism of Action**: MOA description, action type, target ChEMBL ID, disease efficacy
- **Drug Indications**: EFO ID/term, MeSH ID/heading, max phase for indication
- **Target Metadata**: Target components, UniProt accessions, HGNC symbols

---

## 🧬 FIELD 2: Target Antigen

**Goal:** Normalize antigen targets (e.g., "TROP2", "HER2", "CD19") to canonical gene symbols with HGNC identifiers.

### ✅ Implemented API:

**[HGNC REST API](https://rest.genenames.org/)** – Official gene nomenclature database.

### 🔧 Key Features:

* **Gene Symbol Standardization**: Maps aliases to approved HGNC symbols
* **Cross-references**: Ensembl gene IDs, UniProt mappings
* **Synonym Management**: Comprehensive alias symbol lists
* **Status Tracking**: Canonical vs unknown status for quality control

### 🛠️ Implementation:

```python
HGNC_API = "https://rest.genenames.org/search/{}"

def query_hgnc(symbol):
    """Query HGNC API for symbol or alias"""
    url = HGNC_API.format(symbol)
    response = requests.get(url, headers={"Accept": "application/json"})
    
    # Returns: approved_symbol, hgnc_id, ensembl_gene_id, synonyms, status
```

**Data Retrieved:**
- **Approved Symbol**: Official HGNC gene symbol
- **HGNC ID**: Unique identifier in HGNC database
- **Ensembl Gene ID**: Cross-reference to Ensembl
- **Synonyms**: All known alias symbols
- **Status**: "canonical" or "unknown" for quality tracking

---

## 🧾 FIELD 3: Disease Ontology (DO) / Indication

**Goal:** Normalize indication terms to standardized disease ontologies (DOID, NCIt) with BioPortal integration.

### ✅ Implemented API:

**[BioPortal API](https://data.bioontology.org/)** – Multi-ontology search including DOID and NCIt.

### 🔧 Key Features:

* **Multi-ontology Search**: Simultaneous search across DOID and NCIt
* **Fuzzy Matching**: Flexible matching for variant disease names
* **Synonym Support**: Comprehensive disease term synonyms
* **Ontology Tracking**: Source ontology identification

### 🛠️ Implementation:

```python
BIOPORTAL_SEARCH_URL = "https://data.bioontology.org/search"

def query_doid_bioportal(disease_term):
    """Search BioPortal DOID/NCIt for disease"""
    params = {
        "q": disease_term,
        "ontologies": "DOID,NCIT",
        "require_exact_match": "false"
    }
    
    # Returns: label, ontology, match_id, synonyms, status
```

**Data Retrieved:**
- **Label**: Preferred disease term from ontology
- **Ontology**: Source ontology (DOID or NCIt)
- **Match ID**: Unique identifier in the ontology
- **Synonyms**: Alternative disease terms
- **Status**: "canonical" or "unknown" for quality tracking

---

## 🔄 Integration Plan Overview

| Field           | API(s)                    | Output Entities                | Implementation Status |
| --------------- | ------------------------- | ------------------------------ | -------------------- |
| Drug Name/Alias | ChEMBL WebResource Client | ChEMBL ID, MOA, Indications, Targets | ✅ **Implemented** |
| Target Antigen  | HGNC REST API             | HGNC ID, Gene Symbol, Synonyms | ✅ **Implemented** |
| Indication (DO) | BioPortal API             | DOID, NCIt, UMLS CUI           | ✅ **Implemented** |

---

## 🧰 Data Pipeline Design

### 1. **Input Processing**
```python
df = pd.read_csv("ADC drugs and trials_Apr2025.xlsx - OG List.csv", header=1)
```

### 2. **Multi-API Integration**
- **ChEMBL**: Comprehensive drug metadata extraction
- **HGNC**: Gene symbol standardization
- **BioPortal**: Disease ontology mapping

### 3. **Data Quality Control**
- **Status Tracking**: Each field tagged as "canonical" or "unknown"
- **Error Handling**: Graceful handling of API failures
- **Comprehensive Metadata**: Rich context for each normalized entity

### 4. **Output Format**
- **Structured Data**: JSON-like dictionaries with full metadata
- **Cross-references**: Multiple identifier systems (ChEMBL, HGNC, DOID, NCIt)
- **Quality Indicators**: Status flags for data validation

---

## 🚀 Usage Examples

### Drug Name Processing
```python
drug_names = ["trastuzumab", "sacituzumab govitecan", "brentuximab vedotin"]
chembl_data = [fetch_full_chembl_data(name) for name in drug_names]
```

### Target Antigen Processing
```python
antigens = ["TROP2", "HER2", "CD19", "unknown antigen"]
results = [query_hgnc(clean_target_antigen_name(ag)) for ag in antigens]
```

### Disease Indication Processing
```python
diseases = ["Triple Negative Breast Cancer", "NSCLC", "Bladder Ca"]
results = [query_doid_bioportal(clean_disease_name(d)) for d in diseases]
```

---

## 📋 Dependencies

```python
# Required packages
pandas
chembl_webresource_client
requests
urllib
json
pprint
```

---

## 🔧 Configuration

- **BioPortal API Key**: Set `BIOPORTAL_API_KEY` for enhanced BioPortal access
- **ChEMBL Client**: No API key required for basic usage
- **HGNC API**: No authentication required

---

Would you like to extend this pipeline with additional features or integrate it with your specific ADC dataset?
