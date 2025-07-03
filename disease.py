import requests
import pandas as pd
import json
from tqdm import tqdm
import os

# Load and clean dataframe
df = pd.read_csv("ADC drugs and trials_Apr2025.xlsx - OG List.csv", header=1)
df = df.iloc[:, 1:]

# --- BioPortal Config ---

BIOPORTAL_API_KEY = os.environ.get("BIOPORTAL_API_KEY")
BIOPORTAL_SEARCH_URL = "https://data.bioontology.org/search"

# --- Clean disease names ---
def clean_disease_name(name):
    return str(name).strip().lower()

# --- Query BioPortal for disease term ---
def query_doid_bioportal(disease_term):
    params = {
        "q": disease_term,
        "ontologies": "DOID,NCIT",
        "require_exact_match": "false",
        "apikey": BIOPORTAL_API_KEY
    }
    try:
        response = requests.get(BIOPORTAL_SEARCH_URL, params=params, timeout=10)
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
    except Exception:
        pass

    return {
        "input": disease_term,
        "label": None,
        "ontology": None,
        "match_id": None,
        "synonyms": [],
        "status": "unknown"
    }

# --- Apply to unique DO entries ---
unique_diseases = df['DO'].dropna().unique()
doid_results = {
    d: query_doid_bioportal(clean_disease_name(d))
    for d in tqdm(unique_diseases, desc="Querying BioPortal DOID/NCIT")
}

# --- Map back to dataframe ---
df['doid_label'] = df['DO'].map(lambda x: doid_results.get(x, {}).get('label'))
df['doid_match_id'] = df['DO'].map(lambda x: doid_results.get(x, {}).get('match_id'))
df['doid_ontology'] = df['DO'].map(lambda x: doid_results.get(x, {}).get('ontology'))
df['doid_synonyms'] = df['DO'].map(lambda x: doid_results.get(x, {}).get('synonyms'))

df[['DO','doid_synonyms']].dropna().to_csv("do.csv")