#!/usr/bin/env python3
"""
Unified ADC Data Enrichment Pipeline

This script orchestrates all individual enrichment scripts to:
1. Generate all dictionaries in a 'dictionaries' folder
2. Enrich the input JSON with all ontology mappings
3. Create a comprehensive output with standardized data structure

Usage: python unified_enrichment.py
"""

import json
import os
import shutil
import subprocess
import sys
from pathlib import Path
from typing import Dict, List, Any, Optional
import logging

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler('enrichment_pipeline.log'),
        logging.StreamHandler()
    ]
)
logger = logging.getLogger(__name__)

# Configuration
CONFIG = {
    "INPUT_JSON": "aacrArticle.json",
    "OUTPUT_JSON": "aacrArticle_fully_enriched.json",
    "DICTIONARIES_FOLDER": "dictionaries",
    "SCRIPTS": {
        "antigen": "antigen.py",
        "disease": "disease_enhanced.py", 
        "drug": "drug.py",
        "payload_linker": "payload_linker.py"
    },
    "REQUIRED_FILES": {
        "HGNC_TSV": "hgnc_complete_set.tsv",
        "TACA_JSON": "taca.json"
    }
}

class UnifiedEnrichmentPipeline:
    """Main pipeline class that orchestrates all enrichment scripts"""
    
    def __init__(self):
        self.dictionaries_folder = Path(CONFIG["DICTIONARIES_FOLDER"])
        self.input_json = CONFIG["INPUT_JSON"]
        self.output_json = CONFIG["OUTPUT_JSON"]
        self.enriched_data = None
        
    def setup_directories(self):
        """Create necessary directories"""
        logger.info("Setting up directories...")
        
        # Create dictionaries folder
        self.dictionaries_folder.mkdir(exist_ok=True)
        
        # Create subdirectories for different dictionary types
        (self.dictionaries_folder / "antigen").mkdir(exist_ok=True)
        (self.dictionaries_folder / "disease").mkdir(exist_ok=True)
        (self.dictionaries_folder / "drug").mkdir(exist_ok=True)
        (self.dictionaries_folder / "payload_linker").mkdir(exist_ok=True)
        
        logger.info(f"✅ Created directory structure in {self.dictionaries_folder}")
    
    def check_required_files(self):
        """Verify all required input files exist"""
        logger.info("Checking required files...")
        
        missing_files = []
        for file_type, filename in CONFIG["REQUIRED_FILES"].items():
            if not Path(filename).exists():
                missing_files.append(filename)
        
        if missing_files:
            logger.error(f"❌ Missing required files: {missing_files}")
            return False
        
        logger.info("✅ All required files found")
        return True
    
    def run_script(self, script_name: str, script_path: str) -> bool:
        """Run an individual enrichment script"""
        logger.info(f"🔄 Running {script_name} script...")
        
        try:
            result = subprocess.run(
                [sys.executable, script_path],
                capture_output=True,
                text=True,
                check=True
            )
            logger.info(f"✅ {script_name} completed successfully")
            return True
        except subprocess.CalledProcessError as e:
            logger.error(f"❌ {script_name} failed: {e.stderr}")
            return False
    
    # Removed collect_dictionaries method - dictionaries are now saved directly to organized folders
    
    def load_enriched_data(self):
        """Load and merge all enriched data files"""
        logger.info("📥 Loading and merging enriched data...")
        
        # Start with the original input data
        with open(self.input_json, 'r') as f:
            self.enriched_data = json.load(f)
        logger.info(f"✅ Loaded base data from {self.input_json}")
        
        # Merge antigen enrichments
        antigen_file = "dictionaries/antigen/aacrArticle_hgnc.json"
        if Path(antigen_file).exists():
            with open(antigen_file, 'r') as f:
                antigen_data = json.load(f)
            self.merge_antigen_enrichments(antigen_data)
            logger.info(f"✅ Merged antigen enrichments from {antigen_file}")
        
        # Merge drug enrichments
        drug_file = self.dictionaries_folder / "drug" / "aacrArticle_chembl_enriched.json"
        if drug_file.exists():
            with open(drug_file, 'r') as f:
                drug_data = json.load(f)
            self.merge_drug_enrichments(drug_data)
            logger.info(f"✅ Merged drug enrichments from {drug_file}")
        
        # Merge payload/linker enrichments
        payload_linker_file = self.dictionaries_folder / "payload_linker" / "aacrArticle_chembl_payload_linker_enriched.json"
        if payload_linker_file.exists():
            with open(payload_linker_file, 'r') as f:
                payload_linker_data = json.load(f)
            self.merge_payload_linker_enrichments(payload_linker_data)
            logger.info(f"✅ Merged payload/linker enrichments from {payload_linker_file}")
        
        # Merge disease enrichments
        disease_file = self.dictionaries_folder / "disease" / "aacrArticle_disease_enriched.json"
        if disease_file.exists():
            with open(disease_file, 'r') as f:
                disease_data = json.load(f)
            self.merge_disease_enrichments(disease_data)
            logger.info(f"✅ Merged disease enrichments from {disease_file}")
        
        return True
    
    def merge_antigen_enrichments(self, antigen_data):
        """Merge antigen enrichments from antigen.py output"""
        # Create a mapping from entry ID to antigen enrichments
        antigen_map = {}
        for entry in antigen_data:
            entry_id = entry.get("id")
            if entry_id:
                antigen_map[entry_id] = {}
                for drug in entry.get("extractedDrugs", []):
                    drug_name = drug.get("drugName")
                    if drug_name:
                        antigen_map[entry_id][drug_name] = {
                            "targetOntology": drug.get("targetOntology", []),
                            "targetAntigenCanonicalized": drug.get("targetAntigenCanonicalized", [])
                        }
        
        # Merge into main data
        for entry in self.enriched_data:
            entry_id = entry.get("id")
            if entry_id in antigen_map:
                for drug in entry.get("extractedDrugs", []):
                    drug_name = drug.get("drugName")
                    if drug_name in antigen_map[entry_id]:
                        enrichments = antigen_map[entry_id][drug_name]
                        drug["targetOntology"] = enrichments.get("targetOntology", [])
                        drug["targetAntigenCanonicalized"] = enrichments.get("targetAntigenCanonicalized", [])
    
    def merge_drug_enrichments(self, drug_data):
        """Merge drug enrichments from drug.py output"""
        # Create a mapping from entry ID to drug enrichments
        drug_map = {}
        for entry in drug_data:
            entry_id = entry.get("id")
            if entry_id:
                drug_map[entry_id] = {}
                for drug in entry.get("extractedDrugs", []):
                    drug_name = drug.get("drugName")
                    if drug_name:
                        drug_map[entry_id][drug_name] = {
                            "drugNameChembl": drug.get("drugNameChembl"),
                            "mechanismOfActionChembl": drug.get("mechanismOfActionChembl", [])
                        }
        
        # Merge into main data
        for entry in self.enriched_data:
            entry_id = entry.get("id")
            if entry_id in drug_map:
                for drug in entry.get("extractedDrugs", []):
                    drug_name = drug.get("drugName")
                    if drug_name in drug_map[entry_id]:
                        enrichments = drug_map[entry_id][drug_name]
                        drug["drugNameChembl"] = enrichments.get("drugNameChembl")
                        drug["mechanismOfActionChembl"] = enrichments.get("mechanismOfActionChembl", [])
    
    def merge_payload_linker_enrichments(self, payload_linker_data):
        """Merge payload/linker enrichments from payload_linker.py output"""
        # Create a mapping from entry ID to payload/linker enrichments
        payload_linker_map = {}
        for entry in payload_linker_data:
            entry_id = entry.get("id")
            if entry_id:
                payload_linker_map[entry_id] = {}
                for drug in entry.get("extractedDrugs", []):
                    drug_name = drug.get("drugName")
                    if drug_name:
                        payload_linker_map[entry_id][drug_name] = {
                            "payloadNameChembl": drug.get("payloadNameChembl"),
                            "payloadOntology": drug.get("payloadOntology", []),
                            "payloadCanonicalized": drug.get("payloadCanonicalized", []),
                            "linkerNameChembl": drug.get("linkerNameChembl"),
                            "linkerOntology": drug.get("linkerOntology", []),
                            "linkerCanonicalized": drug.get("linkerCanonicalized", [])
                        }
        
        # Merge into main data
        for entry in self.enriched_data:
            entry_id = entry.get("id")
            if entry_id in payload_linker_map:
                for drug in entry.get("extractedDrugs", []):
                    drug_name = drug.get("drugName")
                    if drug_name in payload_linker_map[entry_id]:
                        enrichments = payload_linker_map[entry_id][drug_name]
                        drug["payloadNameChembl"] = enrichments.get("payloadNameChembl")
                        drug["payloadOntology"] = enrichments.get("payloadOntology", [])
                        drug["payloadCanonicalized"] = enrichments.get("payloadCanonicalized", [])
                        drug["linkerNameChembl"] = enrichments.get("linkerNameChembl")
                        drug["linkerOntology"] = enrichments.get("linkerOntology", [])
                        drug["linkerCanonicalized"] = enrichments.get("linkerCanonicalized", [])
    
    def merge_disease_enrichments(self, disease_data):
        """Merge disease enrichments from disease_enhanced.py output"""
        # Create a mapping from entry ID to disease enrichments
        disease_map = {}
        for entry in disease_data:
            entry_id = entry.get("id")
            if entry_id:
                disease_map[entry_id] = {}
                for drug in entry.get("extractedDrugs", []):
                    drug_name = drug.get("drugName")
                    if drug_name:
                        disease_map[entry_id][drug_name] = {
                            "diseaseOntology": drug.get("diseaseOntology", [])
                        }
        
        # Merge into main data
        for entry in self.enriched_data:
            entry_id = entry.get("id")
            if entry_id in disease_map:
                for drug in entry.get("extractedDrugs", []):
                    drug_name = drug.get("drugName")
                    if drug_name in disease_map[entry_id]:
                        enrichments = disease_map[entry_id][drug_name]
                        drug["diseaseOntology"] = enrichments.get("diseaseOntology", [])
    
    def create_comprehensive_enrichment(self):
        """Create the final comprehensive enriched JSON with standardized structure"""
        logger.info("🔧 Creating comprehensive enrichment...")
        
        if not self.enriched_data:
            logger.error("❌ No data loaded for enrichment")
            return False
        
        # Load all dictionaries for reference
        dictionaries = self.load_all_dictionaries()
        
        # Create comprehensive enrichment
        comprehensive_data = []
        
        for entry in self.enriched_data:
            enriched_entry = self.enrich_single_entry(entry, dictionaries)
            comprehensive_data.append(enriched_entry)
        
        # Save comprehensive output
        with open(self.output_json, 'w') as f:
            json.dump(comprehensive_data, f, indent=2)
        
        logger.info(f"✅ Saved comprehensive enrichment to {self.output_json}")
        return True
    
    def load_all_dictionaries(self) -> Dict[str, Any]:
        """Load all generated dictionaries"""
        dictionaries = {}
        
        # Load antigen dictionary (from HGNC results)
        antigen_file = self.dictionaries_folder / "antigen" / "aacrArticle_hgnc.json"
        if antigen_file.exists():
            with open(antigen_file, 'r') as f:
                dictionaries["antigen"] = json.load(f)
        
        # Load disease dictionary
        disease_file = self.dictionaries_folder / "disease" / "doid_cancer_leaf_paths.json"
        if disease_file.exists():
            with open(disease_file, 'r') as f:
                dictionaries["disease"] = json.load(f)
        
        # Load drug dictionary
        drug_file = self.dictionaries_folder / "drug" / "chembl_drug_dictionary.json"
        if drug_file.exists():
            with open(drug_file, 'r') as f:
                dictionaries["drug"] = json.load(f)
        
        # Load payload and linker dictionaries
        payload_file = self.dictionaries_folder / "payload_linker" / "chembl_payload_dictionary.json"
        linker_file = self.dictionaries_folder / "payload_linker" / "chembl_linker_dictionary.json"
        
        if payload_file.exists():
            with open(payload_file, 'r') as f:
                dictionaries["payload"] = json.load(f)
        
        if linker_file.exists():
            with open(linker_file, 'r') as f:
                dictionaries["linker"] = json.load(f)
        
        return dictionaries
    
    def enrich_single_entry(self, entry: Dict[str, Any], dictionaries: Dict[str, Any]) -> Dict[str, Any]:
        """Enrich a single entry with comprehensive ontology data"""
        
        enriched_entry = {
            "id": entry.get("id"),
            "createdAt": entry.get("createdAt"),
            "title": entry.get("title"),
            "abstract": entry.get("abstract"),
            "url": entry.get("url"),
            "extractedDrugs": []
        }
        
        for drug in entry.get("extractedDrugs", []):
            enriched_drug = self.enrich_single_drug(drug, dictionaries)
            enriched_entry["extractedDrugs"].append(enriched_drug)
        
        return enriched_entry
    
    def enrich_single_drug(self, drug: Dict[str, Any], dictionaries: Dict[str, Any]) -> Dict[str, Any]:
        """Enrich a single drug with all ontology mappings"""
        
        enriched_drug = {
            # Original fields
            "extractedAt": drug.get("extractedAt"),
            "drugName": drug.get("drugName"),
            "drugNameConfidence": drug.get("drugNameConfidence"),
            "company": drug.get("company"),
            "companyConfidence": drug.get("companyConfidence"),
            "cancerIndication": drug.get("cancerIndication"),
            "cancerIndicationConfidence": drug.get("cancerIndicationConfidence"),
            "targetAntigen": drug.get("targetAntigen"),
            "targetAntigenConfidence": drug.get("targetAntigenConfidence"),
            "mechanismOfAction": drug.get("mechanismOfAction"),
            "mechanismOfActionConfidence": drug.get("mechanismOfActionConfidence"),
            "payload": drug.get("payload"),
            "linker": drug.get("linker"),
            "phase": drug.get("phase"),
            
            # Enriched ontology fields
            "ontology": {
                "drug": self.get_drug_ontology(drug, dictionaries),
                "antigen": self.get_antigen_ontology(drug, dictionaries),
                "disease": self.get_disease_ontology(drug, dictionaries),
                "payload": self.get_payload_ontology(drug, dictionaries),
                "linker": self.get_linker_ontology(drug, dictionaries)
            }
        }
        
        return enriched_drug
    
    def get_drug_ontology(self, drug: Dict[str, Any], dictionaries: Dict[str, Any]) -> Dict[str, Any]:
        """Extract drug ontology information"""
        drug_ontology = {
            "chembl_id": None,
            "preferred_name": None,
            "max_phase": None,
            "mechanism_of_action": [],
            "targets": [],
            "indications": [],
            "match_status": "unknown"
        }
        
        # Check if we have ChEMBL enrichment
        if "drugNameChembl" in drug and drug.get("drugNameChembl"):
            drug_ontology["preferred_name"] = drug.get("drugNameChembl")
            drug_ontology["match_status"] = "chembl_match"
        
        if "mechanismOfActionChembl" in drug and drug.get("mechanismOfActionChembl"):
            drug_ontology["mechanism_of_action"] = drug.get("mechanismOfActionChembl", [])
        
        return drug_ontology
    
    def get_antigen_ontology(self, drug: Dict[str, Any], dictionaries: Dict[str, Any]) -> Dict[str, Any]:
        """Extract antigen ontology information"""
        antigen_ontology = {
            "hgnc_symbol": None,
            "hgnc_id": None,
            "ensembl_gene_id": None,
            "synonyms": [],
            "locus_type": None,
            "gene_group": [],
            "taca_subtype": None,
            "taca_family": None,
            "match_status": "unknown"
        }
        
        # Check if we have targetOntology enrichment
        if "targetOntology" in drug:
            target_ontology = drug.get("targetOntology", [])
            if target_ontology:
                # Take the first match
                first_match = target_ontology[0]
                antigen_ontology["match_status"] = first_match.get("match_type", "unknown")
                
                if first_match.get("HGNC"):
                    hgnc_data = first_match["HGNC"]
                    antigen_ontology.update({
                        "hgnc_symbol": hgnc_data.get("symbol"),
                        "hgnc_id": hgnc_data.get("hgnc_id"),
                        "ensembl_gene_id": hgnc_data.get("ensembl_gene_id"),
                        "synonyms": hgnc_data.get("synonyms", []),
                        "locus_type": hgnc_data.get("locus"),
                        "gene_group": hgnc_data.get("family", "").split(", ") if hgnc_data.get("family") else []
                    })
                
                if first_match.get("TACA"):
                    taca_data = first_match["TACA"]
                    antigen_ontology.update({
                        "taca_subtype": taca_data.get("subtype"),
                        "taca_family": taca_data.get("family")
                    })
        
        return antigen_ontology
    
    def get_disease_ontology(self, drug: Dict[str, Any], dictionaries: Dict[str, Any]) -> Dict[str, Any]:
        """Extract disease ontology information"""
        disease_ontology = {
            "doid_id": None,
            "doid_label": None,
            "ncit_id": None,
            "ncit_label": None,
            "synonyms": [],
            "hierarchy_path": [],
            "match_status": "unknown",
            "all_diseases": []  # Store all diseases for multiple indications
        }
        
        # Check if we have diseaseOntology enrichment
        if "diseaseOntology" in drug and drug.get("diseaseOntology"):
            disease_ontology_list = drug.get("diseaseOntology", [])
            if disease_ontology_list:
                # Store all diseases
                disease_ontology["all_diseases"] = disease_ontology_list
                
                # Take the first disease as primary (most relevant)
                first_disease = disease_ontology_list[0]
                disease_ontology.update({
                    "doid_id": first_disease.get("doid_id"),
                    "doid_label": first_disease.get("doid_label"),
                    "match_status": first_disease.get("match_status", "unknown"),
                    "hierarchy_path": first_disease.get("hierarchy_paths", [])[0] if first_disease.get("hierarchy_paths") else []
                })
                
                # Add expanded terms as synonyms
                expanded_terms = first_disease.get("expanded_terms", [])
                if expanded_terms:
                    disease_ontology["synonyms"] = expanded_terms
        
        return disease_ontology
    
    def get_payload_ontology(self, drug: Dict[str, Any], dictionaries: Dict[str, Any]) -> Dict[str, Any]:
        """Extract payload ontology information"""
        payload_ontology = {
            "chembl_id": None,
            "preferred_name": None,
            "max_phase": None,
            "molecule_type": None,
            "match_status": "unknown"
        }
        
        # Check if we have payload enrichment
        if "payloadNameChembl" in drug and drug.get("payloadNameChembl"):
            payload_ontology.update({
                "preferred_name": drug.get("payloadNameChembl"),
                "match_status": "chembl_match"
            })
        
        if "payloadOntology" in drug and drug.get("payloadOntology"):
            payload_ontology["match_status"] = "ontology_match"
        
        return payload_ontology
    
    def get_linker_ontology(self, drug: Dict[str, Any], dictionaries: Dict[str, Any]) -> Dict[str, Any]:
        """Extract linker ontology information"""
        linker_ontology = {
            "chembl_id": None,
            "preferred_name": None,
            "max_phase": None,
            "molecule_type": None,
            "match_status": "unknown"
        }
        
        # Check if we have linker enrichment
        if "linkerNameChembl" in drug and drug.get("linkerNameChembl"):
            linker_ontology.update({
                "preferred_name": drug.get("linkerNameChembl"),
                "match_status": "chembl_match"
            })
        
        if "linkerOntology" in drug and drug.get("linkerOntology"):
            linker_ontology["match_status"] = "ontology_match"
        
        return linker_ontology
    
    def run_pipeline(self):
        """Run the complete enrichment pipeline"""
        logger.info("🚀 Starting unified enrichment pipeline...")
        
        # Step 1: Setup
        self.setup_directories()
        
        if not self.check_required_files():
            logger.error("❌ Pipeline failed due to missing required files")
            return False
        
        # Step 2: Run individual scripts
        for script_name, script_path in CONFIG["SCRIPTS"].items():
            if not Path(script_path).exists():
                logger.warning(f"⚠️  Script {script_path} not found, skipping...")
                continue
            
            if not self.run_script(script_name, script_path):
                logger.error(f"❌ Pipeline failed at {script_name} script")
                return False
        
        # Step 3: Dictionaries are now saved directly to organized folders
        logger.info("📁 Dictionaries saved directly to organized folders")
        
        # Step 4: Load enriched data
        if not self.load_enriched_data():
            logger.error("❌ Failed to load enriched data")
            return False
        
        # Step 5: Create comprehensive enrichment
        if not self.create_comprehensive_enrichment():
            logger.error("❌ Failed to create comprehensive enrichment")
            return False
        
        logger.info("🎉 Unified enrichment pipeline completed successfully!")
        return True

def main():
    """Main entry point"""
    pipeline = UnifiedEnrichmentPipeline()
    success = pipeline.run_pipeline()
    
    if success:
        print("\n" + "="*60)
        print("✅ UNIFIED ENRICHMENT PIPELINE COMPLETED SUCCESSFULLY")
        print("="*60)
        print(f"📁 Dictionaries organized in: {CONFIG['DICTIONARIES_FOLDER']}/")
        print(f"📄 Final enriched output: {CONFIG['OUTPUT_JSON']}")
        print(f"📋 Log file: enrichment_pipeline.log")
        print("="*60)
    else:
        print("\n" + "="*60)
        print("❌ PIPELINE FAILED")
        print("="*60)
        print("Check the log file for details: enrichment_pipeline.log")
        print("="*60)
        sys.exit(1)

if __name__ == "__main__":
    main() 