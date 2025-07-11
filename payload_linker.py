import json
import re
from time import sleep
from tqdm import tqdm
from functools import lru_cache
from chembl_webresource_client.new_client import new_client

# ------------------
# CONFIGURATION
# ------------------
INPUT_JSON = "aacrArticle.json"
OUTPUT_JSON = "dictionaries/payload_linker/aacrArticle_chembl_payload_linker_enriched.json"
PAYLOAD_DICT_OUTPUT = "chembl_payload_dictionary.json"
LINKER_DICT_OUTPUT = "chembl_linker_dictionary.json"

# ------------------
# LOAD DATA
# ------------------
with open(INPUT_JSON, "r") as f:
    data = json.load(f)

# ------------------
# EXTRACT UNIQUE PAYLOADS AND LINKERS
# ------------------
payload_set = set()
linker_set = set()

for entry in data:
    for drug in entry.get("extractedDrugs", []):
        # Extract payloads
        payload = drug.get("payload")
        if isinstance(payload, str):
            payload_set.add(payload.strip())
        elif isinstance(payload, list):
            for p in payload:
                payload_set.add(str(p).strip())
        
        # Extract linkers
        linker = drug.get("linker")
        if isinstance(linker, str):
            linker_set.add(linker.strip())
        elif isinstance(linker, list):
            for l in linker:
                linker_set.add(str(l).strip())

payload_list = sorted(payload_set)
linker_list = sorted(linker_set)

print(f"Found {len(payload_list)} unique payloads and {len(linker_list)} unique linkers")

# ------------------
# UTILS
# ------------------
def expand_component_name(name):
    match = re.match(r"^(.*?)\s*\((.*?)\)$", name)
    if match:
        outside, inside = match.groups()
        return [outside.strip(), inside.strip()]
    return [name.strip()]

# ------------------
# QUERY ChEMBL
# ------------------
molecule_client = new_client.molecule
mechanism_client = new_client.mechanism
indication_client = new_client.drug_indication
target_client = new_client.target

@lru_cache(maxsize=None)
def fetch_full_chembl_data(component_name):
    try:
        results = molecule_client.search(component_name)
    except Exception as e:
        # Handle API errors gracefully
        print(f"Warning: API error for '{component_name}': {str(e)}")
        return None
    if not results:
        return None

    mol_info = results[0]
    chembl_id = mol_info['molecule_chembl_id']
    molecule_details = molecule_client.get(chembl_id)

    # Filter for small molecule drugs only
    if molecule_details.get('molecule_type') != 'Small molecule':
        return None

    if not mol_info.get('pref_name') and molecule_details.get('max_phase') is None:
        return None

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

    return {
        'ChEMBL ID': chembl_id,
        'Preferred Name': mol_info.get('pref_name'),
        'All Aliases': [],
        'Molecule Info': molecule_info
    }

def best_component_match(raw_name):
    variants = expand_component_name(raw_name)
    for name in variants:
        # Skip very short search terms that will cause API errors
        if len(name.strip()) < 3:
            continue
        result = fetch_full_chembl_data(name)
        if result and result.get("ChEMBL ID"):
            return result
    return None

# ------------------
# BUILD PAYLOAD LOOKUP DICT
# ------------------
print("\n🔍 Processing Payloads...")
payload_dict = {}
for payload in tqdm(payload_list, desc="Querying ChEMBL for payloads"):
    # Skip very short payload names
    if len(payload.strip()) < 3:
        print(f"Skipping short payload name: '{payload}'")
        continue
        
    result = best_component_match(payload)
    if result and result.get("ChEMBL ID"):
        chembl_id = result["ChEMBL ID"]
        if chembl_id not in payload_dict:
            result['All Aliases'] = [payload.upper()]
            payload_dict[chembl_id] = result
        else:
            if payload.upper() not in payload_dict[chembl_id]['All Aliases']:
                payload_dict[chembl_id]['All Aliases'].append(payload.upper())
                payload_dict[chembl_id]['All Aliases'] = list(set(payload_dict[chembl_id]['All Aliases']))

# ------------------
# BUILD LINKER LOOKUP DICT
# ------------------
print("\n🔗 Processing Linkers...")
linker_dict = {}
for linker in tqdm(linker_list, desc="Querying ChEMBL for linkers"):
    # Skip very short linker names
    if len(linker.strip()) < 3:
        print(f"Skipping short linker name: '{linker}'")
        continue
        
    result = best_component_match(linker)
    if result and result.get("ChEMBL ID"):
        chembl_id = result["ChEMBL ID"]
        if chembl_id not in linker_dict:
            result['All Aliases'] = [linker.upper()]
            linker_dict[chembl_id] = result
        else:
            if linker.upper() not in linker_dict[chembl_id]['All Aliases']:
                linker_dict[chembl_id]['All Aliases'].append(linker.upper())
                linker_dict[chembl_id]['All Aliases'] = list(set(linker_dict[chembl_id]['All Aliases']))

# ------------------
# SAVE DICTIONARIES
# ------------------
import os
os.makedirs("dictionaries/payload_linker", exist_ok=True)

with open("dictionaries/payload_linker/chembl_payload_dictionary.json", "w") as f:
    json.dump(payload_dict, f, indent=2)

with open("dictionaries/payload_linker/chembl_linker_dictionary.json", "w") as f:
    json.dump(linker_dict, f, indent=2)

# ------------------
# ENRICH JSON INPUT
# ------------------
for entry in data:
    for drug in entry.get("extractedDrugs", []):
        # Process payloads
        payload = drug.get("payload")
        if payload:
            # Handle both string and list payloads
            payload_list = []
            if isinstance(payload, str):
                payload_list = [payload]
            elif isinstance(payload, list):
                payload_list = payload

            # Find matching ChEMBL records for payloads
            payload_matches = []
            for p in payload_list:
                for chembl_id, result in payload_dict.items():
                    if p.upper() in result['All Aliases']:
                        payload_matches.append(result)
                        break

            if payload_matches:
                # Add payload fields
                drug['payloadNameChembl'] = payload_matches[0].get('Preferred Name')
                
                # Also keep the original ontology fields for backward compatibility
                drug["payloadOntology"] = payload_matches
                drug["payloadCanonicalized"] = [r.get("pref_name") for r in payload_matches if r and r.get("pref_name")]

        # Process linkers
        linker = drug.get("linker")
        if linker:
            # Handle both string and list linkers
            linker_list = []
            if isinstance(linker, str):
                linker_list = [linker]
            elif isinstance(linker, list):
                linker_list = linker

            # Find matching ChEMBL records for linkers
            linker_matches = []
            for l in linker_list:
                for chembl_id, result in linker_dict.items():
                    if l.upper() in result['All Aliases']:
                        linker_matches.append(result)
                        break

            if linker_matches:
                # Add linker fields
                drug['linkerNameChembl'] = linker_matches[0].get('Preferred Name')
                
                # Also keep the original ontology fields for backward compatibility
                drug["linkerOntology"] = linker_matches
                drug["linkerCanonicalized"] = [r.get("pref_name") for r in linker_matches if r and r.get("pref_name")]

# ------------------
# SAVE ENRICHED JSON
# ------------------
with open(OUTPUT_JSON, "w") as f:
    json.dump(data, f, indent=2)

print(f"\n✅ Saved {len(payload_dict)} unique payload mappings to dictionaries/payload_linker/chembl_payload_dictionary.json")
print(f"✅ Saved {len(linker_dict)} unique linker mappings to dictionaries/payload_linker/chembl_linker_dictionary.json")
print(f"✅ Saved enriched JSON to {OUTPUT_JSON}")