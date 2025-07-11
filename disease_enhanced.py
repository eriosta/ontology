#!/usr/bin/env python3
"""
Enhanced Disease Ontology Enrichment Script

This script enriches disease/cancer indication data with:
- Acronym expansion (NSCLC -> non-small cell lung cancer)
- Fuzzy matching for similar terms
- Multiple indications per drug
- DOID and NCIT ontology mappings
- Comprehensive JSON enrichment
"""

import json
import re
import difflib
from collections import defaultdict
from typing import Dict, List, Any, Optional, Tuple
import pandas as pd
from tqdm import tqdm

# Configuration
CONFIG = {
    "INPUT_JSON": "aacrArticle.json",
    "OUTPUT_JSON": "dictionaries/disease/aacrArticle_disease_enriched.json",
    "DOID_HIERARCHY": "dictionaries/disease/doid_cancer_leaf_paths.json",
    "FUZZY_CUTOFF": 0.70,
    "ACRONYM_CUTOFF": 0.60
}

# Common cancer acronyms and their expansions
CANCER_ACRONYMS = {
    "NSCLC": "non-small cell lung cancer",
    "SCLC": "small cell lung cancer", 
    "AML": "acute myeloid leukemia",
    "ALL": "acute lymphoblastic leukemia",
    "DLBCL": "diffuse large B-cell lymphoma",
    "CRC": "colorectal cancer",
    "CRPC": "castrate-resistant prostate cancer",
    "ESCC": "esophageal squamous cell carcinoma",
    "EWS": "Ewing sarcoma",
    "B-ALL": "B-cell acute lymphoblastic leukemia",
    "Ph-like": "Philadelphia chromosome-like",
    "TKI": "tyrosine kinase inhibitor",
    "EGFR-TKI": "EGFR tyrosine kinase inhibitor",
    "HER2": "human epidermal growth factor receptor 2",
    "ER+": "estrogen receptor positive",
    "HER2-": "human epidermal growth factor receptor 2 negative",
    "CD30+": "CD30 positive",
    "CDH17+": "CDH17 positive",
    "CEACAM5+": "CEACAM5 positive",
    "CRLF2+": "CRLF2 positive",
    "ETV6-NTRK3+": "ETV6-NTRK3 fusion positive",
    "B7-H3+": "B7-H3 positive"
}

class DiseaseEnricher:
    """Enhanced disease ontology enricher"""
    
    def __init__(self):
        self.doid_hierarchy = {}
        self.doid_labels = {}
        self.doid_synonyms = {}
        self.expanded_terms = []
        self.fuzzy_matches = {}
        
    def load_doid_hierarchy(self):
        """Load DOID hierarchy from disease.py output"""
        print("🔄 Loading DOID hierarchy...")
        
        try:
            with open(CONFIG["DOID_HIERARCHY"], 'r') as f:
                self.doid_hierarchy = json.load(f)
            
            # Create label and synonym mappings
            for doid, data in self.doid_hierarchy.items():
                label = data.get("label", "")
                self.doid_labels[doid] = label
                
                # Create synonyms from label variations
                synonyms = self.generate_synonyms(label)
                self.doid_synonyms[doid] = synonyms
                
            print(f"✅ Loaded {len(self.doid_hierarchy)} DOID terms")
            
        except FileNotFoundError:
            print("⚠️  DOID hierarchy file not found. Please check the path in CONFIG['DOID_HIERARCHY'].")
            return False
        
        return True
    
    def generate_synonyms(self, label: str) -> List[str]:
        """Generate synonyms for a disease label"""
        synonyms = [label.lower()]
        
        # Common variations
        variations = [
            label.replace("cancer", "carcinoma"),
            label.replace("carcinoma", "cancer"),
            label.replace("tumor", "cancer"),
            label.replace("cancer", "tumor"),
            label.replace("malignant", "").strip(),
            label.replace("malignancy", "").strip(),
            # Add general terms
            label.replace("estrogen-receptor positive ", ""),
            label.replace("estrogen-receptor negative ", ""),
            label.replace("progesterone-receptor positive ", ""),
            label.replace("progesterone-receptor negative ", ""),
            label.replace("Her2-receptor positive ", ""),
            label.replace("Her2-receptor negative ", ""),
            label.replace("triple-receptor negative ", ""),
            label.replace("luminal ", ""),
            label.replace(" A", ""),
            label.replace(" B", ""),
            label.replace(" C", ""),
            label.replace(" D", ""),
            label.replace(" E", ""),
            label.replace(" F", ""),
            label.replace(" G", ""),
            label.replace(" H", ""),
            label.replace(" I", ""),
            label.replace(" J", ""),
            label.replace(" K", ""),
            label.replace(" L", ""),
            label.replace(" M", ""),
            label.replace(" N", ""),
            label.replace(" O", ""),
            label.replace(" P", ""),
            label.replace(" Q", ""),
            label.replace(" R", ""),
            label.replace(" S", ""),
            label.replace(" T", ""),
            label.replace(" U", ""),
            label.replace(" V", ""),
            label.replace(" W", ""),
            label.replace(" X", ""),
            label.replace(" Y", ""),
            label.replace(" Z", ""),
        ]
        
        for var in variations:
            if var and var != label and len(var.strip()) > 0:
                synonyms.append(var.lower().strip())
        
        return list(set(synonyms))
    
    def expand_acronyms(self, disease_term: str) -> List[str]:
        """Expand acronyms in disease terms"""
        expanded = [disease_term]
        
        # Direct acronym matches
        for acronym, expansion in CANCER_ACRONYMS.items():
            if acronym.lower() in disease_term.lower():
                expanded.append(disease_term.replace(acronym, expansion))
                expanded.append(disease_term.replace(acronym, expansion).replace("(", "").replace(")", ""))
        
        # Handle complex terms with multiple acronyms
        complex_term = disease_term
        for acronym, expansion in CANCER_ACRONYMS.items():
            if acronym in complex_term:
                complex_term = complex_term.replace(acronym, expansion)
        
        if complex_term != disease_term:
            expanded.append(complex_term)
        
        return list(set(expanded))
    
    def fuzzy_match_disease(self, disease_term: str) -> List[Tuple[str, float, str]]:
        """Fuzzy match disease term to DOID labels"""
        matches = []
        
        # Expand acronyms first
        expanded_terms = self.expand_acronyms(disease_term)
        
        for term in expanded_terms:
            # Direct label matches
            for doid, label in self.doid_labels.items():
                if label.lower() == term.lower():
                    matches.append((doid, 1.0, label))
                
                # Synonym matches
                if doid in self.doid_synonyms:
                    for synonym in self.doid_synonyms[doid]:
                        if synonym == term.lower():
                            matches.append((doid, 0.95, label))
            
            # Fuzzy matches
            all_labels = list(self.doid_labels.values())
            fuzzy_matches = difflib.get_close_matches(
                term.lower(), 
                [label.lower() for label in all_labels], 
                n=10,  # Increased to get more candidates
                cutoff=CONFIG["FUZZY_CUTOFF"]
            )
            
            for match in fuzzy_matches:
                # Find the DOID for this label
                for doid, label in self.doid_labels.items():
                    if label.lower() == match:
                        score = difflib.SequenceMatcher(None, term.lower(), match).ratio()
                        
                        # Boost score for anatomical matches
                        boosted_score = score
                        # Check for lung cancer matches
                        if any(lung_term in term.lower() for lung_term in ["lung", "pulmonary", "non-small cell lung", "small cell lung"]) and "lung" in label.lower():
                            boosted_score = min(1.0, score + 0.15)
                        # Check for breast cancer matches
                        elif any(breast_term in term.lower() for breast_term in ["breast", "mammary"]) and "breast" in label.lower():
                            boosted_score = min(1.0, score + 0.15)
                        # Check for prostate cancer matches
                        elif "prostate" in term.lower() and "prostate" in label.lower():
                            boosted_score = min(1.0, score + 0.15)
                        # Check for colon cancer matches
                        elif any(colon_term in term.lower() for colon_term in ["colon", "colorectal", "large intestine"]) and "colon" in label.lower():
                            boosted_score = min(1.0, score + 0.15)
                        
                        matches.append((doid, boosted_score, label))
            
            # Also try matching against synonyms
            all_synonyms = []
            for doid, synonyms in self.doid_synonyms.items():
                for synonym in synonyms:
                    all_synonyms.append((doid, synonym))
            
            synonym_matches = difflib.get_close_matches(
                term.lower(),
                [synonym for _, synonym in all_synonyms],
                n=5,  # Increased to get more candidates
                cutoff=CONFIG["FUZZY_CUTOFF"]
            )
            
            for match in synonym_matches:
                for doid, synonym in all_synonyms:
                    if synonym == match:
                        score = difflib.SequenceMatcher(None, term.lower(), match).ratio()
                        label = self.doid_labels.get(doid, "")
                        
                        # Boost score for anatomical matches
                        boosted_score = score
                        # Check for lung cancer matches
                        if any(lung_term in term.lower() for lung_term in ["lung", "pulmonary", "non-small cell lung", "small cell lung"]) and "lung" in label.lower():
                            boosted_score = min(1.0, score + 0.15)
                        # Check for breast cancer matches
                        elif any(breast_term in term.lower() for breast_term in ["breast", "mammary"]) and "breast" in label.lower():
                            boosted_score = min(1.0, score + 0.15)
                        # Check for prostate cancer matches
                        elif "prostate" in term.lower() and "prostate" in label.lower():
                            boosted_score = min(1.0, score + 0.15)
                        # Check for colon cancer matches
                        elif any(colon_term in term.lower() for colon_term in ["colon", "colorectal", "large intestine"]) and "colon" in label.lower():
                            boosted_score = min(1.0, score + 0.15)
                        
                        matches.append((doid, boosted_score, label))
        
        # Remove duplicates and sort by score
        unique_matches = {}
        for doid, score, label in matches:
            if doid not in unique_matches or score > unique_matches[doid][1]:
                unique_matches[doid] = (doid, score, label)
        
        return sorted(unique_matches.values(), key=lambda x: x[1], reverse=True)
    
    def extract_disease_paths(self, doid: str) -> List[List[str]]:
        """Extract hierarchical paths for a DOID"""
        if doid in self.doid_hierarchy:
            return self.doid_hierarchy[doid].get("label_paths_to_root", [])
        return []
    
    def enrich_disease_data(self, data: List[Dict]) -> List[Dict]:
        """Enrich disease data in the JSON"""
        print("🔍 Enriching disease data...")
        
        # Extract all unique disease terms
        all_diseases = set()
        for entry in data:
            for drug in entry.get("extractedDrugs", []):
                diseases = drug.get("cancerIndication", [])
                if isinstance(diseases, list):
                    all_diseases.update(diseases)
                elif isinstance(diseases, str):
                    all_diseases.add(diseases)
        
        print(f"📊 Found {len(all_diseases)} unique disease terms")
        
        # Pre-compute disease matches
        disease_matches = {}
        for disease in tqdm(all_diseases, desc="Matching diseases"):
            if disease and disease.strip():
                matches = self.fuzzy_match_disease(disease.strip())
                disease_matches[disease] = matches
        
        # Enrich the data
        for entry in data:
            for drug in entry.get("extractedDrugs", []):
                diseases = drug.get("cancerIndication", [])
                if not diseases:
                    continue
                
                # Handle both string and list formats
                if isinstance(diseases, str):
                    diseases = [diseases]
                
                enriched_diseases = []
                for disease in diseases:
                    if not disease or not disease.strip():
                        continue
                    
                    matches = disease_matches.get(disease.strip(), [])
                    if matches:
                        best_match = matches[0]
                        doid, score, label = best_match
                        
                        # Extract paths
                        paths = self.extract_disease_paths(doid)
                        
                        enriched_disease = {
                            "input": disease,
                            "doid_id": doid,
                            "doid_label": label,
                            "match_score": score,
                            "match_status": self.get_match_status(score),
                            "hierarchy_paths": paths,
                            "expanded_terms": self.expand_acronyms(disease)
                        }
                        
                        enriched_diseases.append(enriched_disease)
                    else:
                        # No match found
                        enriched_disease = {
                            "input": disease,
                            "doid_id": None,
                            "doid_label": None,
                            "match_score": 0.0,
                            "match_status": "unknown",
                            "hierarchy_paths": [],
                            "expanded_terms": self.expand_acronyms(disease)
                        }
                        enriched_diseases.append(enriched_disease)
                
                # Add enriched disease data
                drug["diseaseOntology"] = enriched_diseases
        
        return data
    
    def get_match_status(self, score: float) -> str:
        """Get match status based on score"""
        if score >= 0.95:
            return "exact_match"
        elif score >= CONFIG["FUZZY_CUTOFF"]:
            return "fuzzy_match"
        elif score >= CONFIG["ACRONYM_CUTOFF"]:
            return "acronym_match"
        else:
            return "unknown"
    
    def generate_statistics(self, data: List[Dict]) -> Dict:
        """Generate enrichment statistics"""
        total_diseases = 0
        matched_diseases = 0
        match_status_counts = defaultdict(int)
        
        for entry in data:
            for drug in entry.get("extractedDrugs", []):
                disease_ontology = drug.get("diseaseOntology", [])
                for disease in disease_ontology:
                    total_diseases += 1
                    status = disease.get("match_status", "unknown")
                    match_status_counts[status] += 1
                    
                    if status != "unknown":
                        matched_diseases += 1
        
        return {
            "total_diseases": total_diseases,
            "matched_diseases": matched_diseases,
            "match_rate": matched_diseases / total_diseases if total_diseases > 0 else 0,
            "match_status_breakdown": dict(match_status_counts)
        }

def main():
    """Main execution function"""
    print("🚀 Enhanced Disease Ontology Enrichment")
    print("=" * 50)
    
    # Initialize enricher
    enricher = DiseaseEnricher()
    
    # Load DOID hierarchy
    if not enricher.load_doid_hierarchy():
        print("❌ Failed to load DOID hierarchy")
        return
    
    # Load input data
    print("📥 Loading input data...")
    with open(CONFIG["INPUT_JSON"], 'r') as f:
        data = json.load(f)
    
    print(f"✅ Loaded {len(data)} entries")
    
    # Enrich disease data
    enriched_data = enricher.enrich_disease_data(data)
    
    # Generate statistics
    stats = enricher.generate_statistics(enriched_data)
    
    # Save enriched data
    print("💾 Saving enriched data...")
    import os
    os.makedirs("dictionaries/disease", exist_ok=True)
    with open(CONFIG["OUTPUT_JSON"], 'w') as f:
        json.dump(enriched_data, f, indent=2)
    
    # Print statistics
    print("\n📊 Enrichment Statistics:")
    print(f"   • Total diseases: {stats['total_diseases']}")
    print(f"   • Matched diseases: {stats['matched_diseases']}")
    print(f"   • Match rate: {stats['match_rate']:.1%}")
    print("\n📋 Match Status Breakdown:")
    for status, count in stats['match_status_breakdown'].items():
        print(f"   • {status}: {count}")
    
    print(f"\n✅ Saved enriched data to {CONFIG['OUTPUT_JSON']}")

if __name__ == "__main__":
    main() 