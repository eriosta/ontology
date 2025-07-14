#!/usr/bin/env python3
"""
Trial Design Cleaning and Comprehensive Categorization

This script cleans and standardizes trial design information in the ADC dataset.
It categorizes trial designs into 7 comprehensive groups:
1. By Study Phase (preclinical, FIH, Phase 1-4)
2. By Randomization and Blinding (randomized, open-label, etc.)
3. By Structure and Setting (multicenter, parallel-group, etc.)
4. By Dose Finding and Expansion (dose-escalation, dose-expansion)
5. By Master Protocols and Precision Medicine (basket, umbrella, adaptive)
6. By Purpose or Population (superiority, pragmatic, observational)
7. By Analysis or Allocation (factorial, split-body, etc.)
"""

import json
import re
from pathlib import Path
from collections import defaultdict, Counter
import requests
from typing import Dict, List, Optional, Tuple
from difflib import SequenceMatcher

class TrialDesignCleaner:
    def __init__(self):
        self.trial_design_dictionary = {}
        self.known_designs = set()
        # Comprehensive trial design categories based on research
        self.design_categories = {
            # 1. By Study Phase
            "preclinical": {
                "keywords": ["preclinical", "pre-clinical", "preclinical study"],
                "description": "Preclinical studies in animals or in vitro models"
            },
            "first_in_human": {
                "keywords": ["first-in-human", "first in human", "fih", "first human"],
                "description": "First-in-human studies"
            },
            "phase_1": {
                "keywords": ["phase 1", "phase i", "phase1", "phase one"],
                "description": "Phase 1 clinical trials"
            },
            "phase_2": {
                "keywords": ["phase 2", "phase ii", "phase2", "phase two"],
                "description": "Phase 2 clinical trials"
            },
            "phase_3": {
                "keywords": ["phase 3", "phase iii", "phase3", "phase three"],
                "description": "Phase 3 clinical trials"
            },
            "phase_4": {
                "keywords": ["phase 4", "phase iv", "phase4", "post-marketing", "postmarketing"],
                "description": "Phase 4 post-marketing studies"
            },
            
            # 2. By Randomization and Blinding
            "randomized": {
                "keywords": ["randomized", "randomisation", "randomization", "randomised"],
                "description": "Randomized trials"
            },
            "non_randomized": {
                "keywords": ["non-randomized", "non-randomised", "nonrandomized", "nonrandomised"],
                "description": "Non-randomized trials"
            },
            "open_label": {
                "keywords": ["open-label", "open label", "open", "unblinded"],
                "description": "Open-label trials"
            },
            "single_blind": {
                "keywords": ["single-blind", "single blind", "single blinded"],
                "description": "Single-blind trials"
            },
            "double_blind": {
                "keywords": ["double-blind", "double blind", "double blinded", "blinded"],
                "description": "Double-blind trials"
            },
            "placebo_controlled": {
                "keywords": ["placebo-controlled", "placebo controlled", "placebo"],
                "description": "Placebo-controlled trials"
            },
            
            # 3. By Structure and Setting
            "single_center": {
                "keywords": ["single-center", "single center", "singlecentre", "single centre"],
                "description": "Single-center trials"
            },
            "multicenter": {
                "keywords": ["multicenter", "multi-center", "multi center", "multicentre", "multi-centre"],
                "description": "Multicenter trials"
            },
            "parallel_group": {
                "keywords": ["parallel-group", "parallel group", "parallel"],
                "description": "Parallel-group trials"
            },
            "crossover": {
                "keywords": ["crossover", "cross-over", "cross over"],
                "description": "Crossover trials"
            },
            "sequential": {
                "keywords": ["sequential", "sequential design"],
                "description": "Sequential trials"
            },
            "cluster_randomized": {
                "keywords": ["cluster-randomized", "cluster randomized", "cluster randomisation"],
                "description": "Cluster-randomized trials"
            },
            "stepped_wedge": {
                "keywords": ["stepped-wedge", "stepped wedge"],
                "description": "Stepped-wedge trials"
            },
            
            # 4. By Dose Finding and Expansion
            "dose_escalation": {
                "keywords": ["dose-escalation", "dose escalation", "escalation", "dose finding"],
                "description": "Dose-escalation studies"
            },
            "dose_expansion": {
                "keywords": ["dose-expansion", "dose expansion", "expansion"],
                "description": "Dose-expansion studies"
            },
            
            # 5. By Master Protocols and Precision Medicine
            "basket_trial": {
                "keywords": ["basket", "basket trial"],
                "description": "Basket trials"
            },
            "umbrella_trial": {
                "keywords": ["umbrella", "umbrella trial"],
                "description": "Umbrella trials"
            },
            "platform_trial": {
                "keywords": ["platform", "platform trial"],
                "description": "Platform trials"
            },
            "adaptive_trial": {
                "keywords": ["adaptive", "adaptive design", "adaptive trial"],
                "description": "Adaptive trials"
            },
            "enrichment_design": {
                "keywords": ["enrichment", "enrichment design"],
                "description": "Enrichment design trials"
            },
            "n_of_1_trial": {
                "keywords": ["n-of-1", "n of 1", "n=1"],
                "description": "N-of-1 trials"
            },
            
            # 6. By Purpose or Population
            "superiority": {
                "keywords": ["superiority", "superior"],
                "description": "Superiority trials"
            },
            "non_inferiority": {
                "keywords": ["non-inferiority", "non inferiority", "noninferiority"],
                "description": "Non-inferiority trials"
            },
            "equivalence": {
                "keywords": ["equivalence", "equivalent"],
                "description": "Equivalence trials"
            },
            "pragmatic": {
                "keywords": ["pragmatic", "real-world", "real world", "practical"],
                "description": "Pragmatic/real-world trials"
            },
            "registry_based": {
                "keywords": ["registry-based", "registry based", "registry"],
                "description": "Registry-based trials"
            },
            "observational": {
                "keywords": ["observational", "cohort", "case-control", "cross-sectional"],
                "description": "Observational studies"
            },
            "case_crossover": {
                "keywords": ["case-crossover", "case crossover"],
                "description": "Case-crossover studies"
            },
            
            # 7. By Analysis or Allocation
            "factorial": {
                "keywords": ["factorial", "factorial design"],
                "description": "Factorial trials"
            },
            "split_body": {
                "keywords": ["split-body", "split body", "split-body design"],
                "description": "Split-body trials"
            },
            "withdrawal_design": {
                "keywords": ["withdrawal", "withdrawal design"],
                "description": "Withdrawal design trials"
            },
            "delayed_start": {
                "keywords": ["delayed-start", "delayed start"],
                "description": "Delayed-start trials"
            }
        }
        
    def load_data(self, file_path: str) -> List[Dict]:
        """Load the input JSON data"""
        with open(file_path, 'r') as f:
            return json.load(f)
    
    def save_data(self, data: List[Dict], file_path: str):
        """Save the enriched data"""
        with open(file_path, 'w') as f:
            json.dump(data, f, indent=2)
    
    def clean_trial_design(self, trial_design: str) -> str:
        """Clean and standardize trial design descriptions"""
        if not trial_design or trial_design == "unknown":
            return "unknown"
        
        # Clean up the text
        design = trial_design.strip().lower()
        
        # Common trial design variations and standardizations
        design_variations = {
            'open-label, dose-escalation and expansion': 'Open-label, dose-escalation and expansion',
            'open-label dose-escalation': 'Open-label dose-escalation',
            'open label dose escalation': 'Open-label dose-escalation',
            'dose escalation study': 'Dose-escalation study',
            'phase 1 dose escalation': 'Phase 1 dose-escalation',
            'phase i dose escalation': 'Phase 1 dose-escalation',
            'first-in-human dose escalation': 'First-in-human dose-escalation',
            'fih dose escalation': 'First-in-human dose-escalation',
            'randomized phase 2': 'Randomized Phase 2',
            'randomized phase ii': 'Randomized Phase 2',
            'double-blind randomized': 'Double-blind randomized',
            'placebo-controlled randomized': 'Placebo-controlled randomized',
            'multicenter phase 1': 'Multicenter Phase 1',
            'basket trial': 'Basket trial',
            'umbrella trial': 'Umbrella trial',
            'adaptive design': 'Adaptive design',
            'crossover study': 'Crossover study',
            'parallel group': 'Parallel group study',
            'sequential design': 'Sequential design'
        }
        
        # Check for exact matches
        for variation, standard in design_variations.items():
            if design == variation:
                return standard
        
        # Check for partial matches
        for variation, standard in design_variations.items():
            if variation in design or design in variation:
                return standard
        
        # If no match found, return the cleaned original
        return trial_design.title()
    
    def categorize_trial_design(self, trial_design: str) -> Dict[str, bool]:
        """Categorize trial design into comprehensive research-based categories"""
        if not trial_design or trial_design == "unknown":
            return {category: False for category in self.design_categories.keys()}
        
        design_lower = trial_design.lower()
        categories = {}
        
        for category, category_info in self.design_categories.items():
            keywords = category_info["keywords"]
            categories[category] = any(keyword in design_lower for keyword in keywords)
        
        return categories
    
    def infer_trial_design_from_phase(self, phase: str) -> Optional[str]:
        """Infer trial design from phase information"""
        if not phase or phase == "unknown":
            return None
        
        phase_lower = phase.lower()
        
        # Common phase to design mappings
        phase_design_mapping = {
            "preclinical": "Preclinical study",
            "phase 1": "Dose-escalation study",
            "phase i": "Dose-escalation study",
            "phase1": "Dose-escalation study",
            "phase 2": "Expansion study",
            "phase ii": "Expansion study",
            "phase2": "Expansion study",
            "phase 3": "Randomized controlled trial",
            "phase iii": "Randomized controlled trial",
            "phase3": "Randomized controlled trial"
        }
        
        for phase_key, design in phase_design_mapping.items():
            if phase_key in phase_lower:
                return design
        
        return None
    
    def fetch_clinical_trials_info(self, drug_name: str, company: str) -> Optional[Dict]:
        """
        Fetch trial design information from ClinicalTrials.gov API
        This is a placeholder for future implementation
        """
        # TODO: Implement ClinicalTrials.gov API integration
        # For now, return None
        return None
    
    def build_trial_design_dictionary(self, data: List[Dict]) -> Dict:
        """Build a dictionary of all unique trial designs and their frequencies"""
        design_counts = Counter()
        design_categories = defaultdict(lambda: {"count": 0, "examples": []})
        
        for entry in data:
            for drug in entry.get("extractedDrugs", []):
                trial_design = drug.get("trialDesign", "unknown")
                phase = drug.get("phase", "unknown")
                
                # Clean the trial design
                cleaned_design = self.clean_trial_design(trial_design)
                
                # If trial design is unknown, try to infer from phase
                if cleaned_design == "unknown" and phase != "unknown":
                    inferred_design = self.infer_trial_design_from_phase(phase)
                    if inferred_design:
                        cleaned_design = inferred_design
                
                if cleaned_design != "unknown":
                    design_counts[cleaned_design] += 1
                    self.known_designs.add(cleaned_design)
                
                # Categorize the design
                categories = self.categorize_trial_design(cleaned_design)
                for category, is_present in categories.items():
                    if is_present:
                        design_categories[category]["count"] += 1
                        if cleaned_design not in design_categories[category]["examples"]:
                            design_categories[category]["examples"].append(cleaned_design)
        
        # Create organized category structure
        organized_categories = {
            "study_phase": {
                "preclinical": design_categories["preclinical"]["count"],
                "first_in_human": design_categories["first_in_human"]["count"],
                "phase_1": design_categories["phase_1"]["count"],
                "phase_2": design_categories["phase_2"]["count"],
                "phase_3": design_categories["phase_3"]["count"],
                "phase_4": design_categories["phase_4"]["count"]
            },
            "randomization_blinding": {
                "randomized": design_categories["randomized"]["count"],
                "non_randomized": design_categories["non_randomized"]["count"],
                "open_label": design_categories["open_label"]["count"],
                "single_blind": design_categories["single_blind"]["count"],
                "double_blind": design_categories["double_blind"]["count"],
                "placebo_controlled": design_categories["placebo_controlled"]["count"]
            },
            "structure_setting": {
                "single_center": design_categories["single_center"]["count"],
                "multicenter": design_categories["multicenter"]["count"],
                "parallel_group": design_categories["parallel_group"]["count"],
                "crossover": design_categories["crossover"]["count"],
                "sequential": design_categories["sequential"]["count"],
                "cluster_randomized": design_categories["cluster_randomized"]["count"],
                "stepped_wedge": design_categories["stepped_wedge"]["count"]
            },
            "dose_finding_expansion": {
                "dose_escalation": design_categories["dose_escalation"]["count"],
                "dose_expansion": design_categories["dose_expansion"]["count"]
            },
            "master_protocols_precision": {
                "basket_trial": design_categories["basket_trial"]["count"],
                "umbrella_trial": design_categories["umbrella_trial"]["count"],
                "platform_trial": design_categories["platform_trial"]["count"],
                "adaptive_trial": design_categories["adaptive_trial"]["count"],
                "enrichment_design": design_categories["enrichment_design"]["count"],
                "n_of_1_trial": design_categories["n_of_1_trial"]["count"]
            },
            "purpose_population": {
                "superiority": design_categories["superiority"]["count"],
                "non_inferiority": design_categories["non_inferiority"]["count"],
                "equivalence": design_categories["equivalence"]["count"],
                "pragmatic": design_categories["pragmatic"]["count"],
                "registry_based": design_categories["registry_based"]["count"],
                "observational": design_categories["observational"]["count"],
                "case_crossover": design_categories["case_crossover"]["count"]
            },
            "analysis_allocation": {
                "factorial": design_categories["factorial"]["count"],
                "split_body": design_categories["split_body"]["count"],
                "withdrawal_design": design_categories["withdrawal_design"]["count"],
                "delayed_start": design_categories["delayed_start"]["count"]
            }
        }
        
        return {
            "designs": dict(design_counts),
            "categories": dict(design_categories),
            "organized_categories": organized_categories
        }
    
    def enrich_data(self, data: List[Dict]) -> List[Dict]:
        """Enrich the data with cleaned trial designs and additional metadata"""
        print("🔬 Enriching trial design data...")
        
        enriched_count = 0
        unknown_count = 0
        
        for entry in data:
            for drug in entry.get("extractedDrugs", []):
                original_design = drug.get("trialDesign", "unknown")
                phase = drug.get("phase", "unknown")
                drug_name = drug.get("drugName", "")
                company = drug.get("company", "unknown")
                
                # Clean the trial design
                cleaned_design = self.clean_trial_design(original_design)
                
                # If still unknown, try to infer from phase
                if cleaned_design == "unknown" and phase != "unknown":
                    inferred_design = self.infer_trial_design_from_phase(phase)
                    if inferred_design:
                        cleaned_design = inferred_design
                        enriched_count += 1
                
                # Categorize the design
                categories = self.categorize_trial_design(cleaned_design)
                
                # Create organized categories structure for this drug
                organized_categories = {
                    "study_phase": {
                        "preclinical": categories.get("preclinical", False),
                        "first_in_human": categories.get("first_in_human", False),
                        "phase_1": categories.get("phase_1", False),
                        "phase_2": categories.get("phase_2", False),
                        "phase_3": categories.get("phase_3", False),
                        "phase_4": categories.get("phase_4", False)
                    },
                    "randomization_blinding": {
                        "randomized": categories.get("randomized", False),
                        "non_randomized": categories.get("non_randomized", False),
                        "open_label": categories.get("open_label", False),
                        "single_blind": categories.get("single_blind", False),
                        "double_blind": categories.get("double_blind", False),
                        "placebo_controlled": categories.get("placebo_controlled", False)
                    },
                    "structure_setting": {
                        "single_center": categories.get("single_center", False),
                        "multicenter": categories.get("multicenter", False),
                        "parallel_group": categories.get("parallel_group", False),
                        "crossover": categories.get("crossover", False),
                        "sequential": categories.get("sequential", False),
                        "cluster_randomized": categories.get("cluster_randomized", False),
                        "stepped_wedge": categories.get("stepped_wedge", False)
                    },
                    "dose_finding_expansion": {
                        "dose_escalation": categories.get("dose_escalation", False),
                        "dose_expansion": categories.get("dose_expansion", False)
                    },
                    "master_protocols_precision": {
                        "basket_trial": categories.get("basket_trial", False),
                        "umbrella_trial": categories.get("umbrella_trial", False),
                        "platform_trial": categories.get("platform_trial", False),
                        "adaptive_trial": categories.get("adaptive_trial", False),
                        "enrichment_design": categories.get("enrichment_design", False),
                        "n_of_1_trial": categories.get("n_of_1_trial", False)
                    },
                    "purpose_population": {
                        "superiority": categories.get("superiority", False),
                        "non_inferiority": categories.get("non_inferiority", False),
                        "equivalence": categories.get("equivalence", False),
                        "pragmatic": categories.get("pragmatic", False),
                        "registry_based": categories.get("registry_based", False),
                        "observational": categories.get("observational", False),
                        "case_crossover": categories.get("case_crossover", False)
                    },
                    "analysis_allocation": {
                        "factorial": categories.get("factorial", False),
                        "split_body": categories.get("split_body", False),
                        "withdrawal_design": categories.get("withdrawal_design", False),
                        "delayed_start": categories.get("delayed_start", False)
                    }
                }
                
                # Add enriched fields
                drug["trialDesignCleaned"] = cleaned_design
                drug["trialDesignOriginal"] = original_design
                drug["trialDesignOrganizedCategories"] = organized_categories
                drug["trialDesignConfidence"] = 1 if cleaned_design != "unknown" else 0
                
                if cleaned_design == "unknown":
                    unknown_count += 1
                else:
                    enriched_count += 1
        
        print(f"   ✅ Enriched {enriched_count} trial designs")
        print(f"   ⚠️  {unknown_count} trial designs remain unknown")
        
        return data
    
    def generate_trial_design_dictionary(self, data: List[Dict]) -> Dict:
        """Generate a comprehensive trial design dictionary"""
        design_dict = self.build_trial_design_dictionary(data)
        
        # Add metadata
        dictionary = {
            "metadata": {
                "total_designs": len(design_dict["designs"]),
                "known_designs": list(self.known_designs),
                "categories": list(self.design_categories.keys()),
                "category_descriptions": {category: info["description"] for category, info in self.design_categories.items()}
            },
            "designs": design_dict["designs"],
            "categories": design_dict["categories"],
            "organized_categories": design_dict["organized_categories"]
        }
        
        return dictionary
    
    def run_pipeline(self, input_file: str, output_file: str):
        """Run the complete trial design cleaning pipeline"""
        print("🔬 Trial Design Cleaning Pipeline")
        print("=" * 50)
        
        # Load data
        print(f"📂 Loading data from {input_file}...")
        data = self.load_data(input_file)
        print(f"   ✅ Loaded {len(data)} entries")
        
        # Build trial design dictionary
        print("📊 Building trial design dictionary...")
        design_dict = self.build_trial_design_dictionary(data)
        print(f"   ✅ Found {len(design_dict['designs'])} unique trial designs")
        
        # Show top designs
        print("\n🏆 Top Trial Designs:")
        for design, count in sorted(design_dict["designs"].items(), key=lambda x: x[1], reverse=True)[:10]:
            print(f"   • {design}: {count} drugs")
        
        # Show organized category breakdown
        print("\n📊 Trial Design Categories by Group:")
        
        group_names = {
            "study_phase": "1. By Study Phase",
            "randomization_blinding": "2. By Randomization and Blinding", 
            "structure_setting": "3. By Structure and Setting",
            "dose_finding_expansion": "4. By Dose Finding and Expansion",
            "master_protocols_precision": "5. By Master Protocols and Precision Medicine",
            "purpose_population": "6. By Purpose or Population",
            "analysis_allocation": "7. By Analysis or Allocation"
        }
        
        for group_key, group_name in group_names.items():
            if group_key in design_dict["organized_categories"]:
                print(f"\n   {group_name}:")
                for subcategory, count in design_dict["organized_categories"][group_key].items():
                    if count > 0:
                        description = self.design_categories[subcategory]["description"]
                        print(f"      • {subcategory}: {count} drugs ({description})")
        
        # Enrich data
        enriched_data = self.enrich_data(data)
        
        # Save enriched data (only if explicitly requested)
        if output_file != "aacrArticle_trial_design_enriched.json":
            print(f"\n💾 Saving enriched data to {output_file}...")
            self.save_data(enriched_data, output_file)
        else:
            # Always save for the unified pipeline merge step
            print(f"\n💾 Saving enriched data to {output_file}...")
            self.save_data(enriched_data, output_file)
        
        # Save trial design dictionary
        dict_file = "dictionaries/trial_design/trial_design_dictionary.json"
        Path("dictionaries/trial_design").mkdir(parents=True, exist_ok=True)
        
        full_dict = self.generate_trial_design_dictionary(data)
        with open(dict_file, 'w') as f:
            json.dump(full_dict, f, indent=2)
        
        print(f"   ✅ Trial design dictionary saved to {dict_file}")
        
        # Save minimal individual output for unified pipeline debugging
        Path("individual_outputs").mkdir(parents=True, exist_ok=True)
        minimal_output = []
        for entry in data:
            minimal_entry = {"id": entry.get("id"), "extractedDrugs": []}
            for drug in entry.get("extractedDrugs", []):
                minimal_drug = {
                    "drugName": drug.get("drugName"),
                    "trialDesignCleaned": drug.get("trialDesignCleaned"),
                    "trialDesignOriginal": drug.get("trialDesignOriginal"),
                    "trialDesignCategories": drug.get("trialDesignCategories"),
                    "trialDesignOrganizedCategories": drug.get("trialDesignOrganizedCategories"),
                    "trialDesignConfidence": drug.get("trialDesignConfidence")
                }
                minimal_entry["extractedDrugs"].append(minimal_drug)
            minimal_output.append(minimal_entry)
        with open("individual_outputs/trial_design_enriched.json", "w") as f:
            json.dump(minimal_output, f, indent=2)

        print("\n🎉 Trial design cleaning pipeline completed!")

def main():
    """Main function"""
    cleaner = TrialDesignCleaner()
    
    input_file = "aacrArticle.json"
    output_file = "aacrArticle_trial_design_enriched.json"
    
    if not Path(input_file).exists():
        print(f"❌ Input file not found: {input_file}")
        return
    
    cleaner.run_pipeline(input_file, output_file)

if __name__ == "__main__":
    main() 