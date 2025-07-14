#!/usr/bin/env python3
"""
Biomarker Strategy Cleaning and ADC-Specific Categorization

This script cleans and categorizes biomarker strategies in the ADC dataset based on
research-based categories including target antigen expression, pharmacokinetic profiling,
multigene expression signatures, resistance mechanisms, and molecular imaging.
"""

import json
import re
from pathlib import Path
from collections import defaultdict, Counter
from typing import Dict, List, Optional, Tuple
from difflib import SequenceMatcher

class BiomarkerStrategyCleaner:
    def __init__(self):
        self.biomarker_dictionary = {}
        self.known_strategies = set()
        
        # Define ADC biomarker strategy categories based on research
        self.biomarker_categories = {
            "target_antigen_expression": {
                "keywords": ["target antigen", "target expression", "antigen expression", "her2", "folate receptor", "trop2", "egfr", "cmet", "dll3", "ceacam5", "fgfr2b", "target verification", "target binding"],
                "description": "Target antigen expression assessment for patient selection and response prediction",
                "examples": ["HER2 expression", "folate receptor alpha expression", "TROP2 expression", "target antigen verification"]
            },
            "pharmacokinetic_profiling": {
                "keywords": ["pharmacokinetic", "pk", "total antibody", "conjugated adc", "unconjugated payload", "free payload", "plasma", "ligand binding assay", "lba", "lc-ms", "mass spectrometry", "exposure", "catabolism", "stability"],
                "description": "Multianalyte PK profiling including total antibody, conjugated ADC, and free payload",
                "examples": ["Total antibody quantification", "Conjugated ADC measurement", "Free payload analysis", "PK analyte profiling"]
            },
            "multigene_expression_signatures": {
                "keywords": ["multigene", "gene signature", "expression signature", "adc-trs", "adc treatment response score", "transcriptome", "gene expression", "proliferation", "adhesion", "tumor biology", "gene set", "gsea"],
                "description": "Multigene expression signatures integrating target expression with tumor biology features",
                "examples": ["ADC Treatment Response Score (ADC-TRS)", "multigene expression signature", "transcriptional profiling"]
            },
            "resistance_mechanisms": {
                "keywords": ["resistance", "resistant", "slc46a3", "noncleavable linker", "escape", "evasion", "refractory", "resistance mechanism"],
                "description": "Biomarkers for ADC resistance mechanisms and escape pathways",
                "examples": ["SLC46A3 expression", "resistance mechanism", "escape pathway"]
            },
            "molecular_imaging": {
                "keywords": ["molecular imaging", "pet", "spect", "imaging", "target quantification", "radiological", "scan", "tomography"],
                "description": "Molecular imaging for target quantification and visualization",
                "examples": ["PET imaging", "molecular imaging", "target quantification"]
            },
            "immunohistochemistry": {
                "keywords": ["ihc", "immunohistochemistry", "immunostaining", "tissue staining", "antibody staining", "tissue analysis"],
                "description": "IHC-based target antigen assessment",
                "examples": ["IHC analysis", "immunohistochemistry", "tissue staining"]
            },
            "in_vitro_assays": {
                "keywords": ["in vitro", "cell culture", "cell line", "assay", "binding assay", "affinity", "in vitro assay"],
                "description": "In vitro assays for target antigen assessment",
                "examples": ["In vitro binding assay", "cell culture assay", "affinity measurement"]
            },
            "exploratory_biomarkers": {
                "keywords": ["exploratory", "novel", "emerging", "investigational", "experimental", "pilot", "feasibility"],
                "description": "Exploratory and investigational biomarker approaches",
                "examples": ["Exploratory biomarker", "novel approach", "investigational biomarker"]
            }
        }
        
        # Common biomarker strategy variations
        self.strategy_variations = {
            "target engagement verified via elisa": "Target engagement verified via ELISA",
            "gene expression profiling using nanostring": "Gene expression profiling using Nanostring",
            "ihc analysis of tumor samples": "IHC analysis of tumor samples",
            "flow cytometry for cell surface markers": "Flow cytometry for cell surface markers",
            "circulating tumor dna analysis": "Circulating tumor DNA analysis",
            "pharmacodynamic biomarker assessment": "Pharmacodynamic biomarker assessment",
            "safety biomarker monitoring": "Safety biomarker monitoring",
            "predictive biomarker identification": "Predictive biomarker identification",
            "multiplex biomarker panel": "Multiplex biomarker panel"
        }
        
    def load_data(self, file_path: str) -> List[Dict]:
        """Load the input JSON data"""
        with open(file_path, 'r') as f:
            return json.load(f)
    
    def save_data(self, data: List[Dict], file_path: str):
        """Save the enriched data"""
        with open(file_path, 'w') as f:
            json.dump(data, f, indent=2)
    
    def clean_biomarker_strategy(self, strategy: str) -> str:
        """Clean and standardize biomarker strategy descriptions"""
        if not strategy or strategy == "unknown":
            return "unknown"
        
        # Clean up the text
        strategy_clean = strategy.strip()
        
        # Check for exact matches in variations
        strategy_lower = strategy_clean.lower()
        for variation, standard in self.strategy_variations.items():
            if strategy_lower == variation:
                return standard
        
        # Check for partial matches
        for variation, standard in self.strategy_variations.items():
            if variation in strategy_lower or strategy_lower in variation:
                return standard
        
        # If no match found, return the cleaned original
        return strategy_clean
    
    def categorize_biomarker_strategy(self, strategy: str) -> Dict[str, bool]:
        """Categorize biomarker strategy into ADC-specific categories based on research"""
        if not strategy or strategy == "unknown":
            return {category: False for category in self.biomarker_categories.keys()}
        
        strategy_lower = strategy.lower()
        categories = {}
        
        for category, category_info in self.biomarker_categories.items():
            keywords = category_info["keywords"]
            categories[category] = any(keyword in strategy_lower for keyword in keywords)
        
        return categories
    
    def extract_key_technologies(self, strategy: str) -> List[str]:
        """Extract key technologies mentioned in the strategy"""
        if not strategy or strategy == "unknown":
            return []
        
        # Common biomarker technologies
        technologies = [
            "ELISA", "IHC", "FACS", "Flow cytometry", "PCR", "qPCR", "RT-PCR",
            "Western blot", "Mass spectrometry", "Nanostring", "Microarray",
            "RNA-seq", "Next-generation sequencing", "NGS", "PET", "MRI", "CT",
            "Liquid biopsy", "CTC", "ctDNA", "Proteomics", "Transcriptomics"
        ]
        
        found_technologies = []
        strategy_lower = strategy.lower()
        
        for tech in technologies:
            if tech.lower() in strategy_lower:
                found_technologies.append(tech)
        
        return found_technologies
    
    def extract_key_biomarkers(self, strategy: str) -> List[str]:
        """Extract specific biomarker names mentioned in the strategy"""
        if not strategy or strategy == "unknown":
            return []
        
        # Common biomarker patterns
        biomarker_patterns = [
            r'\b[A-Z]{2,4}\d*\b',  # Abbreviations like EGFR, HER2, etc.
            r'\b[A-Z][a-z]+\s+[A-Z][a-z]+\b',  # Two-word names
            r'\b[A-Z][a-z]+\b'  # Single word names
        ]
        
        biomarkers = []
        for pattern in biomarker_patterns:
            matches = re.findall(pattern, strategy)
            biomarkers.extend(matches)
        
        return list(set(biomarkers))
    
    def calculate_strategy_complexity(self, strategy: str) -> int:
        """Calculate complexity score of biomarker strategy"""
        if not strategy or strategy == "unknown":
            return 0
        
        # Simple complexity scoring
        complexity = 0
        
        # Count different categories
        categories = self.categorize_biomarker_strategy(strategy)
        category_count = sum(categories.values())
        complexity += category_count * 2
        
        # Count technologies
        technologies = self.extract_key_technologies(strategy)
        complexity += len(technologies)
        
        # Count biomarkers
        biomarkers = self.extract_key_biomarkers(strategy)
        complexity += len(biomarkers)
        
        # Length factor
        complexity += len(strategy.split()) // 10
        
        return complexity
    
    def build_biomarker_dictionary(self, data: List[Dict]) -> Dict:
        """Build a dictionary of all unique biomarker strategies and their frequencies"""
        strategy_counts = Counter()
        category_counts = defaultdict(lambda: {"count": 0, "examples": []})
        technology_counts = Counter()
        
        for entry in data:
            for drug in entry.get("extractedDrugs", []):
                strategy = drug.get("biomarkerStrategy", "unknown")
                
                # Clean the strategy
                cleaned_strategy = self.clean_biomarker_strategy(strategy)
                
                if cleaned_strategy != "unknown":
                    strategy_counts[cleaned_strategy] += 1
                    self.known_strategies.add(cleaned_strategy)
                
                # Categorize by ADC-specific categories
                categories = self.categorize_biomarker_strategy(cleaned_strategy)
                for category, is_present in categories.items():
                    if is_present:
                        category_counts[category]["count"] += 1
                        if cleaned_strategy not in category_counts[category]["examples"]:
                            category_counts[category]["examples"].append(cleaned_strategy)
                
                # Count technologies
                technologies = self.extract_key_technologies(cleaned_strategy)
                for tech in technologies:
                    technology_counts[tech] += 1
        
        return {
            "strategies": dict(strategy_counts),
            "categories": dict(category_counts),
            "technologies": dict(technology_counts)
        }
    
    def enrich_data(self, data: List[Dict]) -> List[Dict]:
        """Enrich the data with cleaned biomarker strategies and additional metadata"""
        print("🧬 Enriching biomarker strategy data...")
        
        enriched_count = 0
        unknown_count = 0
        
        for entry in data:
            for drug in entry.get("extractedDrugs", []):
                original_strategy = drug.get("biomarkerStrategy", "unknown")
                
                # Clean the strategy
                cleaned_strategy = self.clean_biomarker_strategy(original_strategy)
                
                # Categorize the strategy
                categories = self.categorize_biomarker_strategy(cleaned_strategy)
                
                # Extract additional information
                technologies = self.extract_key_technologies(cleaned_strategy)
                biomarkers = self.extract_key_biomarkers(cleaned_strategy)
                complexity = self.calculate_strategy_complexity(cleaned_strategy)
                
                # Add enriched fields
                drug["biomarkerStrategyCleaned"] = cleaned_strategy
                drug["biomarkerStrategyOriginal"] = original_strategy
                drug["biomarkerStrategyCategories"] = categories
                drug["biomarkerTechnologies"] = technologies
                drug["biomarkerMolecules"] = biomarkers
                drug["biomarkerComplexity"] = complexity
                drug["biomarkerStrategyConfidence"] = 1 if cleaned_strategy != "unknown" else 0
                
                if cleaned_strategy == "unknown":
                    unknown_count += 1
                else:
                    enriched_count += 1
        
        print(f"   ✅ Enriched {enriched_count} biomarker strategies")
        print(f"   ⚠️  {unknown_count} biomarker strategies remain unknown")
        
        return data
    
    def generate_biomarker_dictionary(self, data: List[Dict]) -> Dict:
        """Generate a comprehensive biomarker strategy dictionary"""
        biomarker_dict = self.build_biomarker_dictionary(data)
        
        # Add metadata
        dictionary = {
            "metadata": {
                "total_strategies": len(biomarker_dict["strategies"]),
                "known_strategies": list(self.known_strategies),
                "categories": list(self.biomarker_categories.keys()),
                "category_descriptions": {category: info["description"] for category, info in self.biomarker_categories.items()},
                "category_examples": {category: info["examples"] for category, info in self.biomarker_categories.items()}
            },
            "strategies": biomarker_dict["strategies"],
            "categories": biomarker_dict["categories"],
            "technologies": biomarker_dict["technologies"]
        }
        
        return dictionary
    
    def run_pipeline(self, input_file: str, output_file: str):
        """Run the complete biomarker strategy cleaning pipeline"""
        print("🧬 Biomarker Strategy Cleaning Pipeline")
        print("=" * 50)
        
        # Load data
        print(f"📂 Loading data from {input_file}...")
        data = self.load_data(input_file)
        print(f"   ✅ Loaded {len(data)} entries")
        
        # Build biomarker dictionary
        print("📊 Building biomarker strategy dictionary...")
        biomarker_dict = self.build_biomarker_dictionary(data)
        print(f"   ✅ Found {len(biomarker_dict['strategies'])} unique strategies")
        
        # Show top strategies
        print("\n🏆 Top Biomarker Strategies:")
        for strategy, count in sorted(biomarker_dict["strategies"].items(), key=lambda x: x[1], reverse=True)[:10]:
            print(f"   • {strategy}: {count} drugs")
        
        # Show category breakdown
        print("\n📊 ADC Biomarker Strategy Categories:")
        for category, info in biomarker_dict["categories"].items():
            if info["count"] > 0:
                print(f"   • {category}: {info['count']} drugs")
        
        # Show top technologies
        print("\n🔬 Top Technologies:")
        for tech, count in sorted(biomarker_dict["technologies"].items(), key=lambda x: x[1], reverse=True)[:10]:
            print(f"   • {tech}: {count} mentions")
        
        # Enrich data
        enriched_data = self.enrich_data(data)
        
        # Save enriched data (only if explicitly requested)
        if output_file != "aacrArticle_biomarker_enriched.json":
            print(f"\n💾 Saving enriched data to {output_file}...")
            self.save_data(enriched_data, output_file)
        else:
            # Always save for the unified pipeline merge step
            print(f"\n💾 Saving enriched data to {output_file}...")
            self.save_data(enriched_data, output_file)
        
        # Save biomarker dictionary
        dict_file = "dictionaries/biomarker/biomarker_strategy_dictionary.json"
        Path("dictionaries/biomarker").mkdir(parents=True, exist_ok=True)
        
        full_dict = self.generate_biomarker_dictionary(data)
        with open(dict_file, 'w') as f:
            json.dump(full_dict, f, indent=2)
        
        print(f"   ✅ Biomarker strategy dictionary saved to {dict_file}")
        
        # Save minimal individual output for unified pipeline debugging
        Path("individual_outputs").mkdir(parents=True, exist_ok=True)
        minimal_output = []
        for entry in data:
            minimal_entry = {"id": entry.get("id"), "extractedDrugs": []}
            for drug in entry.get("extractedDrugs", []):
                minimal_drug = {
                    "drugName": drug.get("drugName"),
                    "biomarkerStrategyCleaned": drug.get("biomarkerStrategyCleaned"),
                    "biomarkerStrategyOriginal": drug.get("biomarkerStrategyOriginal"),
                    "biomarkerStrategyCategories": drug.get("biomarkerStrategyCategories"),
                    "biomarkerTechnologies": drug.get("biomarkerTechnologies"),
                    "biomarkerMolecules": drug.get("biomarkerMolecules"),
                    "biomarkerComplexity": drug.get("biomarkerComplexity"),
                    "biomarkerStrategyConfidence": drug.get("biomarkerStrategyConfidence")
                }
                minimal_entry["extractedDrugs"].append(minimal_drug)
            minimal_output.append(minimal_entry)
        with open("individual_outputs/biomarker_strategy_enriched.json", "w") as f:
            json.dump(minimal_output, f, indent=2)
        
        print("\n🎉 Biomarker strategy cleaning pipeline completed!")

def main():
    """Main function"""
    cleaner = BiomarkerStrategyCleaner()
    
    input_file = "aacrArticle.json"
    output_file = "aacrArticle_biomarker_enriched.json"
    
    if not Path(input_file).exists():
        print(f"❌ Input file not found: {input_file}")
        return
    
    cleaner.run_pipeline(input_file, output_file)

if __name__ == "__main__":
    main() 