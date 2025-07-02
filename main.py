import pandas as pd

df = pd.read_csv("ADC drugs and trials_Apr2025.xlsx - OG List.csv", header=1)
df.T

# Step 1: Drug Names
from chembl_webresource_client.new_client import new_client
import pandas as pd

molecule_client = new_client.molecule
mechanism_client = new_client.mechanism
indication_client = new_client.drug_indication
target_client = new_client.target

def fetch_full_chembl_data(drug_name):
    # Search molecule
    results = molecule_client.search(drug_name)
    if not results:
        return {
            'Input Name': drug_name,
            'Found': False,
            'Molecule Info': None,
            'Mechanism of Action': [],
            'Drug Indications': [],
            'Target Metadata': []
        }

    mol_info = results[0]
    chembl_id = mol_info['molecule_chembl_id']
    molecule_details = molecule_client.get(chembl_id)

    # --- Molecule Info ---
    molecule_info = {
        'ChEMBL ID': chembl_id,
        'Preferred Name': mol_info.get('pref_name'),
        'Max Phase': molecule_details.get('max_phase'),
        'First Approval': molecule_details.get('first_approval'),
        'Drug Type': molecule_details.get('drug_type'),
        'Molecule Type': molecule_details.get('molecule_type'),
        'Withdrawn': molecule_details.get('withdrawn_flag'),
        'Black Box Warning': molecule_details.get('black_box_warning'),
        'ATC Codes': [x['level5'] for x in molecule_details.get('atc_classifications', []) if 'level5' in x],
        'USAN Stem': molecule_details.get('usan_stem'),
        'Indication Class': molecule_details.get('indication_class')
    }

    # --- Mechanism of Action ---
    moas_raw = mechanism_client.filter(molecule_chembl_id=chembl_id)
    moas = []
    targets_seen = set()
    for m in moas_raw:
        moas.append({
            'Mechanism of Action': m.get('mechanism_of_action'),
            'Action Type': m.get('action_type'),
            'Target ChEMBL ID': m.get('target_chembl_id'),
            'Target Pref Name': m.get('target_pref_name'),
            'Target Type': m.get('target_type'),
            'Disease Efficacy': m.get('disease_efficacy')
        })
        if m.get('target_chembl_id'):
            targets_seen.add(m['target_chembl_id'])

    # --- Drug Indications ---
    indications_raw = indication_client.filter(molecule_chembl_id=chembl_id)
    indications = [{
        'EFO ID': ind.get('efo_id'),
        'EFO Term': ind.get('efo_term'),
        'MeSH ID': ind.get('mesh_id'),
        'MeSH Term': ind.get('mesh_heading'),
        'Max Phase for Indication': ind.get('max_phase_for_ind')
    } for ind in indications_raw]

    # --- Target Metadata ---
    targets = []
    for tid in targets_seen:
        try:
            t = target_client.get(tid)
            targets.append({
                'Target ChEMBL ID': t.get('target_chembl_id'),
                'Pref Name': t.get('pref_name'),
                'Target Type': t.get('target_type'),
                'Organism': t.get('organism'),
                'Target Components': [
                    {
                        'Accession': c.get('accession'),
                        'Component Type': c.get('component_type'),
                        'Organism': c.get('organism'),
                        'UniProt': c.get('target_component_xrefs', [{}])[0].get('xref_id'),
                        'HGNC': c.get('target_component_xrefs', [{}])[0].get('xref_name')
                    }
                    for c in t.get('target_components', [])
                ]
            })
        except:
            continue

    return {
        'Input Name': drug_name,
        'Found': True,
        'Molecule Info': molecule_info,
        'Mechanism of Action': moas,
        'Drug Indications': indications,
        'Target Metadata': targets
    }

# ✅ Example Usage
drug_names = ["trastuzumab", "sacituzumab govitecan", "brentuximab vedotin"]
chembl_data = [fetch_full_chembl_data(name) for name in drug_names]

# Optional: look at one entry
import pprint
pprint.pprint(chembl_data[0])



# Step 2: Target Name
import requests
import pandas as pd

HGNC_API = "https://rest.genenames.org/search/{}"
HEADERS = {"Accept": "application/json"}

def clean_target_antigen_name(name):
    """Normalize input (e.g., strip, uppercase)"""
    return name.strip().upper()

def query_hgnc(symbol):
    """Query HGNC API for symbol or alias"""
    url = HGNC_API.format(symbol)
    response = requests.get(url, headers=HEADERS)
    if response.status_code == 200:
        hits = response.json().get("response", {}).get("docs", [])
        if hits:
            match = hits[0]  # take the first match
            return {
                "input": symbol,
                "approved_symbol": match.get("symbol"),
                "hgnc_id": match.get("hgnc_id"),
                "ensembl_gene_id": match.get("ensembl_gene_id"),
                "synonyms": match.get("alias_symbol", []),
                "status": "canonical"
            }
    return {
        "input": symbol,
        "approved_symbol": None,
        "hgnc_id": None,
        "ensembl_gene_id": None,
        "synonyms": [],
        "status": "unknown"
    }

# EXAMPLE USAGE
antigens = ["TROP2", "HER2", "CD19", "unknown antigen"]
results = [query_hgnc(clean_target_antigen_name(ag)) for ag in antigens]
df_antigens = pd.DataFrame(results)
print(df_antigens)



# Step 3
import requests
import urllib.parse
import json

BIOPORTAL_API_KEY = ""
BIOPORTAL_SEARCH_URL = "https://data.bioontology.org/search"


def clean_disease_name(name):
    return name.strip().lower()

def query_doid_bioportal(disease_term):
    """Search BioPortal DOID/NCIt for disease"""
    params = {
        "q": disease_term,
        "ontologies": "DOID,NCIT",
        "require_exact_match": "false",
        "apikey": BIOPORTAL_API_KEY
    }
    response = requests.get(BIOPORTAL_SEARCH_URL, params=params)
    if response.status_code == 200:
        data = response.json()
        if data["collection"]:
            match = data["collection"][0]
            return {
                "input": disease_term,
                "label": match.get("prefLabel"),
                "ontology": match.get("links", {}).get("ontology"),
                "match_id": match.get("@id"),
                "synonyms": match.get("synonym", []),
                "status": "canonical"
            }
    return {
        "input": disease_term,
        "label": None,
        "ontology": None,
        "match_id": None,
        "synonyms": [],
        "status": "unknown"
    }

# Example usage
diseases = ["Triple Negative Breast Cancer", "NSCLC", "Bladder Ca", "cancer of unknown primary"]
results = [query_doid_bioportal(clean_disease_name(d)) for d in diseases]

# Output as pretty JSON
print(json.dumps(results, indent=2))
