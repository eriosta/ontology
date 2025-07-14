#!/usr/bin/env python3
"""
Company Name Cleaning and Standardization

This script cleans and standardizes company names in the ADC dataset.
It includes logic to extract company names from alphanumeric drug labels
when the company field is unknown.
"""

import json
import re
from pathlib import Path
from collections import defaultdict, Counter
import requests
from typing import Dict, List, Optional, Tuple

class CompanyCleaner:
    def __init__(self):
        self.company_dictionary = {}
        self.known_companies = set()
        self.drug_to_company_mapping = {}
        
    def load_data(self, file_path: str) -> List[Dict]:
        """Load the input JSON data"""
        with open(file_path, 'r') as f:
            return json.load(f)
    
    def save_data(self, data: List[Dict], file_path: str):
        """Save the enriched data"""
        with open(file_path, 'w') as f:
            json.dump(data, f, indent=2)
    
    def extract_company_from_drug_name(self, drug_name: str) -> Optional[str]:
        """
        Extract potential company name from alphanumeric drug labels.
        Common patterns: [Company][Number] or [Company][Letter][Number]
        """
        if not drug_name or drug_name == "unknown":
            return None
            
        # Common company prefixes that might be in drug names
        company_prefixes = [
            'ADC', 'ABT', 'AMG', 'BMS', 'GSK', 'JNJ', 'MRK', 'PFE', 'RHH', 'SNY',
            'AZ', 'RO', 'NVS', 'LLY', 'BMY', 'ABBV', 'TKY', 'DGN', 'IMD', 'GQ',
            'Affinity', 'MediLink', 'Heidelberg', 'ImmunoGen', 'Seattle', 'Genentech'
        ]
        
        # Look for company prefixes in drug name
        for prefix in company_prefixes:
            if drug_name.upper().startswith(prefix.upper()):
                return prefix
        
        # Pattern matching for alphanumeric codes
        # Look for patterns like: [Letters][Numbers] or [Letters][Letters][Numbers]
        patterns = [
            r'^([A-Z]{2,4})\d+',  # 2-4 letters followed by numbers
            r'^([A-Z]{2,4})[A-Z]\d+',  # 2-4 letters, 1 letter, numbers
        ]
        
        for pattern in patterns:
            match = re.match(pattern, drug_name.upper())
            if match:
                potential_company = match.group(1)
                # Check if this looks like a known company code
                if len(potential_company) >= 2:
                    return potential_company
        
        return None
    
    def clean_company_name(self, company_name: str) -> str:
        """Clean and standardize company names"""
        if not company_name or company_name == "unknown":
            return "unknown"
        
        # Remove common suffixes and clean up
        company = company_name.strip()
        
        # Common company name variations
        company_variations = {
            'affinity biopharma': 'Affinity Biopharma',
            'medilink': 'MediLink',
            'heidelberg pharma': 'Heidelberg Pharma',
            'immunogen': 'ImmunoGen',
            'seattle genetics': 'Seattle Genetics',
            'genentech': 'Genentech',
            'roche': 'Roche',
            'novartis': 'Novartis',
            'pfizer': 'Pfizer',
            'merck': 'Merck',
            'bristol-myers squibb': 'Bristol-Myers Squibb',
            'bms': 'Bristol-Myers Squibb',
            'gsk': 'GlaxoSmithKline',
            'glaxosmithkline': 'GlaxoSmithKline',
            'johnson & johnson': 'Johnson & Johnson',
            'jnj': 'Johnson & Johnson',
            'amgen': 'Amgen',
            'abbvie': 'AbbVie',
            'abbott': 'Abbott',
            'sanofi': 'Sanofi',
            'snyn': 'Sanofi',
            'astrazeneca': 'AstraZeneca',
            'az': 'AstraZeneca',
            'eli lilly': 'Eli Lilly',
            'lilly': 'Eli Lilly',
            'lly': 'Eli Lilly',
            'takeda': 'Takeda',
            'tky': 'Takeda',
            'daiichi sankyo': 'Daiichi Sankyo',
            'dgn': 'Daiichi Sankyo',
            'immunomedics': 'Immunomedics',
            'imd': 'Immunomedics'
        }
        
        # Check for exact matches in variations
        lower_company = company.lower()
        for variation, standard in company_variations.items():
            if lower_company == variation:
                return standard
        
        # Check for partial matches
        for variation, standard in company_variations.items():
            if variation in lower_company or lower_company in variation:
                return standard
        
        # If no match found, return the cleaned original
        return company.title()
    
    def build_company_dictionary(self, data: List[Dict]) -> Dict:
        """Build a dictionary of all unique company names and their frequencies"""
        company_counts = Counter()
        drug_company_mapping = {}
        
        for entry in data:
            for drug in entry.get("extractedDrugs", []):
                company = drug.get("company", "unknown")
                drug_name = drug.get("drugName", "")
                
                # Clean the company name
                cleaned_company = self.clean_company_name(company)
                
                # If company is unknown, try to extract from drug name
                if cleaned_company == "unknown" and drug_name:
                    extracted_company = self.extract_company_from_drug_name(drug_name)
                    if extracted_company:
                        cleaned_company = extracted_company
                        drug_company_mapping[drug_name] = extracted_company
                
                if cleaned_company != "unknown":
                    company_counts[cleaned_company] += 1
                    self.known_companies.add(cleaned_company)
                
                # Store the mapping for later use
                drug_company_mapping[drug_name] = cleaned_company
        
        self.drug_to_company_mapping = drug_company_mapping
        return dict(company_counts)
    
    def enrich_data(self, data: List[Dict]) -> List[Dict]:
        """Enrich the data with cleaned company names and additional metadata"""
        print("🏢 Enriching company data...")
        
        enriched_count = 0
        unknown_count = 0
        
        for entry in data:
            for drug in entry.get("extractedDrugs", []):
                original_company = drug.get("company", "unknown")
                drug_name = drug.get("drugName", "")
                
                # Clean the company name
                cleaned_company = self.clean_company_name(original_company)
                
                # If still unknown, try to extract from drug name
                if cleaned_company == "unknown" and drug_name:
                    extracted_company = self.extract_company_from_drug_name(drug_name)
                    if extracted_company:
                        cleaned_company = extracted_company
                        enriched_count += 1
                
                # Add enriched fields
                drug["companyCleaned"] = cleaned_company
                drug["companyOriginal"] = original_company
                drug["companyConfidence"] = 1 if cleaned_company != "unknown" else 0
                
                if cleaned_company == "unknown":
                    unknown_count += 1
                else:
                    enriched_count += 1
        
        print(f"   ✅ Enriched {enriched_count} company names")
        print(f"   ⚠️  {unknown_count} companies remain unknown")
        
        return data
    
    def generate_company_dictionary(self, data: List[Dict]) -> Dict:
        """Generate a comprehensive company dictionary"""
        company_dict = self.build_company_dictionary(data)
        
        # Add metadata
        dictionary = {
            "metadata": {
                "total_companies": len(company_dict),
                "known_companies": list(self.known_companies),
                "drug_to_company_mapping": self.drug_to_company_mapping
            },
            "companies": company_dict
        }
        
        return dictionary
    
    def run_pipeline(self, input_file: str, output_file: str):
        """Run the complete company cleaning pipeline"""
        print("🏢 Company Name Cleaning Pipeline")
        print("=" * 50)
        
        # Load data
        print(f"📂 Loading data from {input_file}...")
        data = self.load_data(input_file)
        print(f"   ✅ Loaded {len(data)} entries")
        
        # Build company dictionary
        print("📊 Building company dictionary...")
        company_dict = self.build_company_dictionary(data)
        print(f"   ✅ Found {len(company_dict)} unique companies")
        
        # Show top companies
        print("\n🏆 Top Companies:")
        for company, count in sorted(company_dict.items(), key=lambda x: x[1], reverse=True)[:10]:
            print(f"   • {company}: {count} drugs")
        
        # Enrich data
        enriched_data = self.enrich_data(data)
        
        # Save minimal individual output for unified pipeline debugging
        Path("individual_outputs").mkdir(parents=True, exist_ok=True)
        minimal_output = []
        for entry in data:
            minimal_entry = {"id": entry.get("id"), "extractedDrugs": []}
            for drug in entry.get("extractedDrugs", []):
                minimal_drug = {
                    "drugName": drug.get("drugName"),
                    "companyCleaned": drug.get("companyCleaned"),
                    "companyOriginal": drug.get("companyOriginal"),
                    "companyConfidence": drug.get("companyConfidence")
                }
                minimal_entry["extractedDrugs"].append(minimal_drug)
            minimal_output.append(minimal_entry)
        with open("individual_outputs/company_enriched.json", "w") as f:
            json.dump(minimal_output, f, indent=2)
        
        # Save enriched data (only if explicitly requested)
        if output_file != "aacrArticle_company_enriched.json":
            print(f"\n💾 Saving enriched data to {output_file}...")
            self.save_data(enriched_data, output_file)
        else:
            # Always save for the unified pipeline merge step
            print(f"\n💾 Saving enriched data to {output_file}...")
            self.save_data(enriched_data, output_file)
        
        # Save company dictionary
        dict_file = "dictionaries/company/company_dictionary.json"
        Path("dictionaries/company").mkdir(parents=True, exist_ok=True)
        
        full_dict = self.generate_company_dictionary(data)
        with open(dict_file, 'w') as f:
            json.dump(full_dict, f, indent=2)
        
        print(f"   ✅ Company dictionary saved to {dict_file}")
        
        print("\n🎉 Company cleaning pipeline completed!")

def main():
    """Main function"""
    cleaner = CompanyCleaner()
    
    input_file = "aacrArticle.json"
    output_file = "aacrArticle_company_enriched.json"
    
    if not Path(input_file).exists():
        print(f"❌ Input file not found: {input_file}")
        return
    
    cleaner.run_pipeline(input_file, output_file)

if __name__ == "__main__":
    main() 