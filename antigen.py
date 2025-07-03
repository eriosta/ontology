import requests
import pandas as pd
from tqdm import tqdm

# Load and clean dataframe
df = pd.read_csv("ADC drugs and trials_Apr2025.xlsx - OG List.csv", header=1)
df = df.iloc[:, 1:]

# --- HGNC Config ---
HGNC_API = "https://rest.genenames.org/search/{}"
HEADERS = {"Accept": "application/json"}

# --- Helper: Normalize input symbols ---
def clean_target_antigen_name(name):
    return str(name).strip().upper()

# --- Helper: Query HGNC API ---
def query_hgnc(symbol):
    url = HGNC_API.format(symbol)
    try:
        response = requests.get(url, headers=HEADERS, timeout=10)
        if response.status_code == 200:
            hits = response.json().get("response", {}).get("docs", [])
            if hits:
                match = hits[0]
                return {
                    "input": symbol,
                    "hgnc_symbol": match.get("symbol"),
                    "hgnc_id": match.get("hgnc_id"),
                    "ensembl_gene_id": match.get("ensembl_gene_id"),
                    "status": "canonical"
                }
    except Exception:
        pass
    return {
        "input": symbol,
        "hgnc_symbol": None,
        "hgnc_id": None,
        "ensembl_gene_id": None,
        "status": "unknown"
    }

# --- Apply to unique antigens ---
unique_antigens = df['Target Antigen'].dropna().unique()

hgnc_results = {
    ag: query_hgnc(clean_target_antigen_name(ag))
    for ag in tqdm(unique_antigens, desc="Querying HGNC")
}

# --- Map back to dataframe ---
df['hgnc_symbol'] = df['Target Antigen'].map(
    lambda x: hgnc_results.get(x, {}).get('hgnc_symbol')
)
df['hgnc_id'] = df['Target Antigen'].map(
    lambda x: hgnc_results.get(x, {}).get('hgnc_id')
)
df['ensembl_gene_id'] = df['Target Antigen'].map(
    lambda x: hgnc_results.get(x, {}).get('ensembl_gene_id')
)

df[['Target Antigen','hgnc_symbol']].dropna()
