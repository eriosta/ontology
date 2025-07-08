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
    # this comes from: https://www.genenames.org/download/
    "HGNC_TSV": "hgnc_complete_set.tsv",
    # TACA 
    "TACA_JSON": "taca.json",
    # input/output json
    "JSON_INPUT": "aacrArticle.json",
    "JSON_OUTPUT": "aacrArticle_hgnc.json",
    # fuzzy matching cutoff
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


def expand_antigen_name(name):
    """
    Splits strings like "Folate Receptor Alpha (FRα)" into ["Folate Receptor Alpha", "FRα"]
    Returns a list with just the name if no parentheses.
    """
    match = re.match(r"^(.*?)\s*\((.*?)\)$", name)
    if match:
        outside, inside = match.groups()
        return [outside.strip(), inside.strip()]
    return [name.strip()]


# ------------------------------------
# HGNC PROCESSING
# ------------------------------------
def load_hgnc_data(tsv_path):
    df = pd.read_csv(tsv_path, sep="\t", low_memory=False)
    excluded_loci = [
        "pseudogene", "RNA, long non-coding", "RNA, micro", "RNA, transfer", "RNA, small nucleolar",
        "immunoglobulin pseudogene", "T cell receptor pseudogene", "RNA, ribosomal", "RNA, small nuclear",
        "RNA, miscellaneous", "RNA, Y", "RNA, vault"
    ]
    return df[~df["locus_type"].isin(excluded_loci)]


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


def query_hgnc(symbol, symbol_map, alias_map, cutoff=0.85, original=None):
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
                "locus_type": None,
                "gene_group": [],
                "status": "unknown"
            }

    gene_group_raw = row.get("gene_group")
    gene_group_list = [g.strip() for g in str(gene_group_raw).split("|") if g.strip()] if pd.notna(gene_group_raw) else []

    return {
        "input": symbol,
        "hgnc_symbol": row["symbol"],
        "hgnc_id": row["hgnc_id"],
        "ensembl_gene_id": row.get("ensembl_gene_id"),
        "synonyms": re.split(r"[|,]", str(row.get("alias_symbol", ""))),
        "locus_type": row.get("locus_type"),
        "gene_group": gene_group_list,
        "status": status
    }


# ------------------------------------
# TACA PROCESSING
# ------------------------------------
def query_taca(symbol, taca_db, cutoff=0.85):
    candidates = []
    for family in taca_db:
        for subtype in family["subtypes"]:
            ratio = difflib.SequenceMatcher(None, symbol.lower(), subtype.lower()).ratio()
            if ratio >= cutoff:
                candidates.append({
                    "input": symbol,
                    "matched_subtype": subtype,
                    "matched_family": family["family"],
                    "references": family.get("references", []),
                    "score": ratio
                })
    if candidates:
        # Return the best match
        best = sorted(candidates, key=lambda x: -x["score"])[0]
        return {
            "input": best["input"],
            "hgnc_symbol": None,
            "hgnc_id": None,
            "ensembl_gene_id": None,
            "synonyms": [],
            "locus_type": "TACA",
            "gene_group": [best["matched_family"]],
            "status": "taca_match",
            "taca_subtype": best["matched_subtype"],
            "references": best["references"]
        }
    else:
        return None

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
            match_results = []

            if isinstance(ag, str):
                results = [hgnc_lookup_results.get(ag)]
            elif isinstance(ag, list):
                results = [hgnc_lookup_results.get(a) for a in ag if isinstance(a, str)]
            else:
                results = []

            for result in results:
                if not result:
                    continue

                status = result.get("status")
                is_match = status != "unknown"
                match_entry = {
                    "input": result.get("input"),
                    "is_match": is_match,
                    "match_type": None,
                    "HGNC": None,
                    "TACA": None
                }

                if status == "taca_match":
                    match_entry["match_type"] = "TACA"
                    match_entry["TACA"] = {
                        "subtype": result.get("taca_subtype"),
                        "family": ", ".join(result.get("gene_group", [])) if result.get("gene_group") else None
                    }
                elif is_match:
                    match_entry["match_type"] = "HGNC"
                    match_entry["HGNC"] = {
                        "symbol": result.get("hgnc_symbol"),
                        "hgnc_id": result.get("hgnc_id"),
                        "ensembl_gene_id": result.get("ensembl_gene_id"),
                        "synonyms": result.get("synonyms"),
                        "locus": result.get("locus_type"),
                        "family": ", ".join(result.get("gene_group", [])) if result.get("gene_group") else None
                    }

                match_results.append(match_entry)

            drug["targetOntology"] = match_results

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
                "targetOntology": drug.get("targetOntology")
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
                        "HGNC": hgnc_result,
                        "locus_type": hgnc_result.get("locus_type"),
                        "gene_group": hgnc_result.get("gene_group")
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
    
    print("🔄 Loading TACA data...")
    taca_db = load_json(CONFIG["TACA_JSON"])["TACA_classifications"]

    print("📥 Loading input JSON...")
    data = load_json(CONFIG["JSON_INPUT"])

    print("🧬 Extracting unique antigens...")
    unique_antigens = extract_unique_antigens(data)
    print(f"✅ Found {len(unique_antigens)} unique antigens")

    print("🔍 Querying HGNC...")
    hgnc_results = {}
    for ag in tqdm(unique_antigens, desc="HGNC Lookup"):
        candidates = expand_antigen_name(ag)
        if not candidates:
            candidates = [ag]
            
        best_result = None
        for candidate in candidates:
            result = query_hgnc(candidate, symbol_map, alias_map, CONFIG["FUZZY_CUTOFF"], original=ag)
            if result["status"] != "unknown":
                best_result = result
                break

        if not best_result:
            best_result = query_hgnc(ag, symbol_map, alias_map, CONFIG["FUZZY_CUTOFF"], original=ag)
            if best_result["status"] == "unknown":
                taca_result = query_taca(ag, taca_db, 0.50)
                if taca_result:
                    best_result = taca_result

        hgnc_results[ag] = best_result

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