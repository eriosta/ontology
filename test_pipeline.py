#!/usr/bin/env python3
"""
Test script for the unified enrichment pipeline

This script validates the pipeline output and provides quality metrics.
"""

import json
import sys
from pathlib import Path
from collections import defaultdict, Counter

def load_json_safe(file_path):
    """Safely load JSON file with error handling"""
    try:
        with open(file_path, 'r') as f:
            return json.load(f)
    except FileNotFoundError:
        print(f"❌ File not found: {file_path}")
        return None
    except json.JSONDecodeError as e:
        print(f"❌ JSON decode error in {file_path}: {e}")
        return None

def analyze_enriched_data(data):
    """Analyze the enriched data and provide quality metrics"""
    print("📊 Analyzing enriched data...")
    
    total_entries = len(data)
    total_drugs = 0
    ontology_stats = defaultdict(lambda: {"matched": 0, "unknown": 0})
    
    # Collect statistics
    for entry in data:
        for drug in entry.get("extractedDrugs", []):
            total_drugs += 1
            ontology = drug.get("ontology", {})
            
            for ontology_type, ontology_data in ontology.items():
                match_status = ontology_data.get("match_status", "unknown")
                if match_status != "unknown":
                    ontology_stats[ontology_type]["matched"] += 1
                else:
                    ontology_stats[ontology_type]["unknown"] += 1
    
    # Print results
    print(f"📈 Summary Statistics:")
    print(f"   • Total entries: {total_entries}")
    print(f"   • Total drugs: {total_drugs}")
    print()
    
    print("🎯 Ontology Match Rates:")
    for ontology_type, stats in ontology_stats.items():
        total = stats["matched"] + stats["unknown"]
        if total > 0:
            match_rate = stats["matched"] / total * 100
            print(f"   • {ontology_type}: {match_rate:.1f}% ({stats['matched']}/{total})")
    
    return ontology_stats

def validate_dictionaries():
    """Validate that all dictionary files exist and are valid JSON"""
    print("🔍 Validating dictionary files...")
    
    dictionary_files = [
        "dictionaries/antigen/aacrArticle_hgnc.json",
        "dictionaries/disease/doid_cancer_leaf_paths.json", 
        "dictionaries/drug/chembl_drug_dictionary.json",
        "dictionaries/payload_linker/chembl_payload_dictionary.json",
        "dictionaries/payload_linker/chembl_linker_dictionary.json"
    ]
    
    valid_files = 0
    for file_path in dictionary_files:
        data = load_json_safe(file_path)
        if data is not None:
            print(f"   ✅ {file_path}")
            valid_files += 1
        else:
            print(f"   ❌ {file_path}")
    
    print(f"\n📁 Dictionary validation: {valid_files}/{len(dictionary_files)} files valid")
    return valid_files == len(dictionary_files)

def check_data_structure(data):
    """Check that the data structure is correct"""
    print("🔧 Validating data structure...")
    
    required_fields = ["id", "extractedDrugs"]
    ontology_fields = ["drug", "antigen", "disease", "payload", "linker"]
    
    structure_errors = []
    
    for i, entry in enumerate(data):
        # Check required fields
        for field in required_fields:
            if field not in entry:
                structure_errors.append(f"Entry {i}: Missing required field '{field}'")
        
        # Check drug ontology structure
        for j, drug in enumerate(entry.get("extractedDrugs", [])):
            ontology = drug.get("ontology", {})
            
            for ontology_type in ontology_fields:
                if ontology_type not in ontology:
                    structure_errors.append(f"Entry {i}, Drug {j}: Missing ontology field '{ontology_type}'")
                else:
                    ontology_data = ontology[ontology_type]
                    if "match_status" not in ontology_data:
                        structure_errors.append(f"Entry {i}, Drug {j}: Missing match_status in {ontology_type}")
    
    if structure_errors:
        print("   ❌ Structure validation failed:")
        for error in structure_errors[:5]:  # Show first 5 errors
            print(f"      • {error}")
        if len(structure_errors) > 5:
            print(f"      ... and {len(structure_errors) - 5} more errors")
        return False
    else:
        print("   ✅ Data structure is valid")
        return True

def sample_analysis(data):
    """Provide a sample analysis of the enriched data"""
    print("\n🔬 Sample Analysis:")
    
    # Find entries with good ontology coverage
    well_mapped_entries = []
    for entry in data:
        for drug in entry.get("extractedDrugs", []):
            ontology = drug.get("ontology", {})
            match_count = sum(1 for o in ontology.values() if o.get("match_status") != "unknown")
            if match_count >= 3:  # At least 3 ontology types matched
                well_mapped_entries.append((entry, drug, match_count))
    
    if well_mapped_entries:
        # Show best mapped example
        best_entry, best_drug, best_count = max(well_mapped_entries, key=lambda x: x[2])
        
        print(f"   📋 Best mapped drug: {best_drug.get('drugName', 'Unknown')}")
        print(f"      • Ontology matches: {best_count}/5")
        print(f"      • Antigen: {best_drug.get('targetAntigen', [])}")
        print(f"      • Disease: {best_drug.get('cancerIndication', [])}")
        
        # Show ontology details
        ontology = best_drug.get("ontology", {})
        for ontology_type, ontology_data in ontology.items():
            if ontology_data.get("match_status") != "unknown":
                print(f"      • {ontology_type}: {ontology_data.get('match_status')}")
    else:
        print("   ⚠️  No well-mapped entries found")

def main():
    """Main test function"""
    print("🧪 Testing Unified Enrichment Pipeline")
    print("=" * 50)
    
    # Check if enriched file exists
    enriched_file = "aacrArticle_fully_enriched.json"
    if not Path(enriched_file).exists():
        print(f"❌ Enriched file not found: {enriched_file}")
        print("Please run the pipeline first: python run_pipeline.py")
        sys.exit(1)
    
    # Load and analyze data
    data = load_json_safe(enriched_file)
    if data is None:
        sys.exit(1)
    
    print(f"✅ Loaded {len(data)} entries from {enriched_file}")
    
    # Run tests
    ontology_stats = analyze_enriched_data(data)
    dict_valid = validate_dictionaries()
    structure_valid = check_data_structure(data)
    sample_analysis(data)
    
    # Summary
    print("\n" + "=" * 50)
    print("📋 TEST SUMMARY")
    print("=" * 50)
    
    total_tests = 3
    passed_tests = 0
    
    if ontology_stats:
        print("✅ Data analysis completed")
        passed_tests += 1
    
    if dict_valid:
        print("✅ Dictionary validation passed")
        passed_tests += 1
    
    if structure_valid:
        print("✅ Data structure validation passed")
        passed_tests += 1
    
    print(f"\n🎯 Overall: {passed_tests}/{total_tests} tests passed")
    
    if passed_tests == total_tests:
        print("🎉 All tests passed! Pipeline is working correctly.")
    else:
        print("⚠️  Some tests failed. Check the output above for details.")
        sys.exit(1)

if __name__ == "__main__":
    main() 