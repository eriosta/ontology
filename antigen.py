# TO-DO: Use GlyTouCan, UniCarbKB / GlyGen, or MeSH (Medical Subject Headings) to standardize the glycan antigens
# HGNC is for genes encoding protein products

import json
import pandas as pd
import re
import difflib
from tqdm import tqdm

# ------------------------------------
# CONFIGURATION
# ------------------------------------
CONFIG = {
    "HGNC_TSV": "hgnc_complete_set.tsv",
    "JSON_INPUT": "aacrArticle.json",
    "JSON_OUTPUT": "aacrArticle_hgnc.json",
    "FUZZY_CUTOFF": 0.85
}


# ------------------------------------
# UTILITIES
# ------------------------------------
def normalize_symbol(s):
    return re.sub(r'[^A-Z0-9]', '', str(s).upper())


def load_json(path):
    with open(path, "r") as f:
        return json.load(f)


def save_json(data, path, indent=2):
    with open(path, "w") as f:
        json.dump(data, f, indent=indent)


# ------------------------------------
# HGNC PROCESSING
# ------------------------------------
def load_hgnc_data(tsv_path):
    df = pd.read_csv(tsv_path, sep="\t", low_memory=False)
    return df
    # return df[df["locus_type"] == "gene with protein product"]


def build_hgnc_maps(hgnc_df):
    symbol_map = {normalize_symbol(row["symbol"]): row for _, row in hgnc_df.iterrows()}
    alias_map = {}
    for _, row in hgnc_df.iterrows():
        aliases = re.split(r"[|,]", str(row.get("alias_symbol", "")))
        for alias in aliases:
            alias_clean = normalize_symbol(alias)
            if alias_clean:
                alias_map[alias_clean] = row
    return symbol_map, alias_map


def query_hgnc(symbol, symbol_map, alias_map, cutoff=0.85):
    cleaned = normalize_symbol(symbol)

    if cleaned in symbol_map:
        row = symbol_map[cleaned]
        status = "canonical"
    elif cleaned in alias_map:
        row = alias_map[cleaned]
        status = "alias_match"
    else:
        all_keys = list(set(symbol_map.keys()) | set(alias_map.keys()))
        match = difflib.get_close_matches(cleaned, all_keys, n=1, cutoff=cutoff)
        if match:
            key = match[0]
            row = symbol_map.get(key)
            if row is None:
                row = alias_map.get(key)
            status = "fuzzy_match"
        else:
            return {
                "input": symbol,
                "hgnc_symbol": None,
                "hgnc_id": None,
                "ensembl_gene_id": None,
                "synonyms": [],
                "status": "unknown"
            }

    return {
        "input": symbol,
        "hgnc_symbol": row["symbol"],
        "hgnc_id": row["hgnc_id"],
        "ensembl_gene_id": row.get("ensembl_gene_id"),
        "synonyms": re.split(r"[|,]", str(row.get("alias_symbol", ""))),
        "status": status
    }


# ------------------------------------
# JSON PROCESSING
# ------------------------------------
def extract_unique_antigens(data):
    antigens = []
    for entry in data:
        for drug in entry.get("extractedDrugs", []):
            ag = drug.get("targetAntigenCanonicalized")
            if isinstance(ag, str):
                antigens.append(ag)
            elif isinstance(ag, list):
                antigens.extend([a for a in ag if isinstance(a, str)])
    return sorted(set(antigens))


def enrich_json(data, hgnc_lookup_results):
    for entry in data:
        for drug in entry.get("extractedDrugs", []):
            ag = drug.get("targetAntigenCanonicalized")
            if isinstance(ag, str):
                drug["HGNC"] = [hgnc_lookup_results.get(ag)]
            elif isinstance(ag, list):
                drug["HGNC"] = [hgnc_lookup_results.get(a) for a in ag if isinstance(a, str)]
    return data


def reduce_json(data):
    minimal_data = []
    for entry in data:
        reduced_entry = {
            "id": entry.get("id"),
            "extractedDrugs": []
        }
        for drug in entry.get("extractedDrugs", []):
            reduced_entry["extractedDrugs"].append({
                "drugName": drug.get("drugName"),
                "cancerIndication": drug.get("cancerIndication"),
                "payload": drug.get("payload"),
                "linker": drug.get("linker"),
                "phase": drug.get("phase"),
                "targetAntigenCanonicalized": drug.get("targetAntigenCanonicalized"),
                "HGNC": drug.get("HGNC")
            })
        minimal_data.append(reduced_entry)
    return minimal_data

def export_unknowns(minimal_data, output_path):
    unknowns = []
    for entry in minimal_data:
        for drug in entry.get("extractedDrugs", []):
            for hgnc_result in (drug.get("HGNC") or []):
                if hgnc_result and hgnc_result.get("status") == "unknown":
                    unknowns.append({
                        "entry_id": entry.get("id"),
                        "drugName": drug.get("drugName"),
                        "cancerIndication": drug.get("cancerIndication"),
                        "payload": drug.get("payload"),
                        "linker": drug.get("linker"),
                        "phase": drug.get("phase"),
                        "targetAntigenCanonicalized": drug.get("targetAntigenCanonicalized"),
                        "HGNC": hgnc_result
                    })

    unknown_output_path = output_path.replace(".json", "_unknowns.json")
    print(f"🛠️  Writing {len(unknowns)} unknown entries to {unknown_output_path}")
    save_json(unknowns, unknown_output_path)
    
    
# ------------------------------------
# MAIN EXECUTION
# ------------------------------------
def main():
    print("🔄 Loading HGNC data...")
    hgnc_df = load_hgnc_data(CONFIG["HGNC_TSV"])
    symbol_map, alias_map = build_hgnc_maps(hgnc_df)

    print("📥 Loading input JSON...")
    data = load_json(CONFIG["JSON_INPUT"])

    print("🧬 Extracting unique antigens...")
    unique_antigens = extract_unique_antigens(data)
    print(f"✅ Found {len(unique_antigens)} unique antigens")

    print("🔍 Querying HGNC...")
    hgnc_results = {
        ag: query_hgnc(ag, symbol_map, alias_map, CONFIG["FUZZY_CUTOFF"])
        for ag in tqdm(unique_antigens, desc="HGNC Lookup")
    }

    print("🧩 Enriching JSON with HGNC results...")
    enriched_data = enrich_json(data, hgnc_results)

    print("📝 Reducing JSON for output...")
    minimal_data = reduce_json(enriched_data)
    export_unknowns(minimal_data, CONFIG["JSON_OUTPUT"])

    print(f"💾 Writing output to {CONFIG['JSON_OUTPUT']}")
    save_json(minimal_data, CONFIG["JSON_OUTPUT"])

    print("🎉 Done!")


# ------------------------------------
# RUN SCRIPT
# ------------------------------------
if __name__ == "__main__":
    main()


