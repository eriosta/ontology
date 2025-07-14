#!/usr/bin/env python3
"""
Simple runner script for the unified enrichment pipeline

Usage: python run_pipeline.py
"""

import sys
import json
from pathlib import Path
from unified_enrichment import UnifiedEnrichmentPipeline

def load_config():
    """Load configuration from JSON file"""
    config_path = Path("pipeline_config.json")
    if not config_path.exists():
        print("❌ Configuration file 'pipeline_config.json' not found!")
        print("Please ensure the configuration file exists in the current directory.")
        return None
    
    try:
        with open(config_path, 'r') as f:
            return json.load(f)
    except json.JSONDecodeError as e:
        print(f"❌ Error parsing configuration file: {e}")
        return None

def check_prerequisites():
    """Check if all required files and scripts exist"""
    print("🔍 Checking prerequisites...")
    
    required_files = [
        "aacrArticle.json",
        "hgnc_complete_set.tsv", 
        "taca.json"
    ]
    
    required_scripts = [
        "antigen.py",
        "disease.py",
        "drug.py", 
        "payload_linker.py",
        "company.py",
        "trial_design.py",
        "biomarker_strategy.py"
    ]
    
    missing_files = []
    missing_scripts = []
    
    for file_path in required_files:
        if not Path(file_path).exists():
            missing_files.append(file_path)
    
    for script_path in required_scripts:
        if not Path(script_path).exists():
            missing_scripts.append(script_path)
    
    if missing_files:
        print(f"❌ Missing required files: {missing_files}")
        return False
    
    if missing_scripts:
        print(f"❌ Missing required scripts: {missing_scripts}")
        return False
    
    print("✅ All prerequisites met!")
    return True

def main():
    """Main entry point"""
    print("🚀 ADC Data Enrichment Pipeline")
    print("=" * 50)
    
    # Load configuration
    config = load_config()
    if not config:
        sys.exit(1)
    
    # Check prerequisites
    if not check_prerequisites():
        print("\n❌ Pipeline cannot proceed due to missing prerequisites.")
        print("Please ensure all required files and scripts are present.")
        sys.exit(1)
    
    # Create and run pipeline
    try:
        pipeline = UnifiedEnrichmentPipeline()
        success = pipeline.run_pipeline()
        
        if success:
            print("\n" + "=" * 60)
            print("✅ PIPELINE COMPLETED SUCCESSFULLY!")
            print("=" * 60)
            print("📁 Output files:")
            print(f"   • Final enriched data: {config['output_files']['final_enriched']}")
            print(f"   • Dictionaries folder: {config['directories']['dictionaries']}/")
            print(f"   • Log file: {config['output_files']['log_file']}")
            print("\n📊 Data structure:")
            print("   • Each drug now has an 'ontology' field with:")
            print("     - drug: ChEMBL mappings")
            print("     - antigen: HGNC/TACA mappings") 
            print("     - disease: DOID/NCIT mappings")
            print("     - payload: ChEMBL small molecule mappings")
            print("     - linker: ChEMBL small molecule mappings")
            print("   • Additional enriched fields:")
            print("     - company: Cleaned company names")
            print("     - trial design: Categorized trial designs")
            print("     - biomarker strategy: ADC-specific categories")
            print("=" * 60)
        else:
            print("\n" + "=" * 60)
            print("❌ PIPELINE FAILED")
            print("=" * 60)
            print("Check the log file for detailed error information.")
            print("=" * 60)
            sys.exit(1)
            
    except KeyboardInterrupt:
        print("\n⚠️  Pipeline interrupted by user")
        sys.exit(1)
    except Exception as e:
        print(f"\n❌ Unexpected error: {e}")
        sys.exit(1)

if __name__ == "__main__":
    main() 