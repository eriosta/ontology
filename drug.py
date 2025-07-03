import pandas as pd
from chembl_webresource_client.new_client import new_client
from tqdm import tqdm
from functools import lru_cache
import re

# --- Load data ---
df = pd.read_csv("ADC drugs and trials_Apr2025.xlsx - OG List.csv", header=1)
df = df.iloc[:, 1:]

# Step 1: Make sure Drug Alias is a list
df['Drug Alias'] = df['Drug Alias'].apply(
    lambda x: [x] if isinstance(x, str) else ([] if pd.isna(x) else x)
)

# Step 2: Append Drug Name to Drug Alias list
df['Drug Alias'] = df.apply(
    lambda row: list(set(row['Drug Alias'] + [row['Drug Name']])) if pd.notnull(row['Drug Name']) else row['Drug Alias'],
    axis=1
)

# Step 3: Drop the original Drug Name column
df.drop(columns=['Drug Name'], inplace=True)

# Optional: Preview rows with aliases
print(df[df['Drug Alias'].apply(len) > 0][['Drug Alias']])

# --- Initialize ChEMBL clients ---
molecule_client = new_client.molecule
mechanism_client = new_client.mechanism
indication_client = new_client.drug_indication
target_client = new_client.target

# --- Fetch full ChEMBL drug record ---
def fetch_full_chembl_data(drug_name):
    try:
        results = molecule_client.search(drug_name)
    except Exception:
        return None
    if not results:
        return None

    mol_info = results[0]
    chembl_id = mol_info['molecule_chembl_id']
    molecule_details = molecule_client.get(chembl_id)

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

    # --- Mechanism of Action ---
    moas_raw = list(mechanism_client.filter(molecule_chembl_id=chembl_id))
    moas = []
    targets_seen = set()
    for m in moas_raw:
        target_id = m.get('target_chembl_id')
        moas.append({
            'Mechanism of Action': m.get('mechanism_of_action'),
            'Action Type': m.get('action_type'),
            'Target ChEMBL ID': target_id,
            'Disease Efficacy': m.get('disease_efficacy')
        })
        if target_id:
            targets_seen.add(target_id)

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
    target_pref_name_map = {}
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
            if t.get('target_chembl_id') and t.get('pref_name'):
                target_pref_name_map[t['target_chembl_id']] = t['pref_name']
        except:
            continue

    # --- Choose Primary Target Name ---
    primary_target_name = None
    for m in moas_raw:
        tid = m.get('target_chembl_id')
        if tid and tid in target_pref_name_map:
            primary_target_name = target_pref_name_map[tid]
            break  # take the first found

    return {
        'Input Name': drug_name,
        'Found': True,
        'Primary Target Name': primary_target_name,
        'Molecule Info': molecule_info,
        'Mechanism of Action': moas,
        'Drug Indications': indications,
        'Target Metadata': targets
    }

# --- LRU Cache to avoid re-fetching the same name ---
@lru_cache(maxsize=None)
def fetch_full_chembl_data_cached(name):
    return fetch_full_chembl_data(name)

# --- Create lookup columns for best-matching ChEMBL info ---
best_matches = []
for aliases in tqdm(df['Drug Alias'], desc="Matching aliases to ChEMBL"):
    match = None
    matched_alias = None
    for alias in aliases:
        result = fetch_full_chembl_data_cached(alias)
        if result and result.get("Found"):
            match = result
            matched_alias = alias
            break
    best_matches.append((matched_alias, match))

# --- Unpack into new columns ---
df['chembl_matched_alias'] = [alias for alias, _ in best_matches]
df['chembl_id'] = [res['Molecule Info']['ChEMBL ID'] if res else "NA" for _, res in best_matches]
df['chembl_name'] = [res['Molecule Info']['Preferred Name'] if res else "NA" for _, res in best_matches]
df['chembl_target'] = [res.get('Primary Target Name', 'NA') if res else "NA" for _, res in best_matches]

# --- Save result ---
df.to_csv("chembl_data.csv", index=False)