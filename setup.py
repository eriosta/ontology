#!/usr/bin/env python3
"""
Setup script for the Unified ADC Data Enrichment Pipeline

This script installs dependencies and sets up the pipeline environment.
"""

import subprocess
import sys
import os
from pathlib import Path

def install_requirements():
    """Install required Python packages"""
    print("📦 Installing Python dependencies...")
    
    try:
        subprocess.check_call([
            sys.executable, "-m", "pip", "install", "-r", "requirements.txt"
        ])
        print("✅ Dependencies installed successfully")
        return True
    except subprocess.CalledProcessError as e:
        print(f"❌ Failed to install dependencies: {e}")
        return False

def check_dependencies():
    """Check if all required packages are installed"""
    print("🔍 Checking installed dependencies...")
    
    required_packages = [
        "pandas",
        "numpy", 
        "tqdm",
        "owlready2",
        "chembl_webresource_client"
    ]
    
    missing_packages = []
    
    for package in required_packages:
        try:
            __import__(package)
            print(f"   ✅ {package}")
        except ImportError:
            print(f"   ❌ {package}")
            missing_packages.append(package)
    
    if missing_packages:
        print(f"\n⚠️  Missing packages: {missing_packages}")
        return False
    else:
        print("\n✅ All dependencies are installed")
        return True

def create_directories():
    """Create necessary directories"""
    print("📁 Creating directories...")
    
    directories = [
        "dictionaries",
        "dictionaries/antigen",
        "dictionaries/disease", 
        "dictionaries/drug",
        "dictionaries/payload_linker"
    ]
    
    for directory in directories:
        Path(directory).mkdir(exist_ok=True)
        print(f"   ✅ Created {directory}/")

def check_input_files():
    """Check if required input files exist"""
    print("📄 Checking input files...")
    
    required_files = [
        "aacrArticle.json",
        "hgnc_complete_set.tsv",
        "taca.json"
    ]
    
    missing_files = []
    
    for file_path in required_files:
        if Path(file_path).exists():
            print(f"   ✅ {file_path}")
        else:
            print(f"   ❌ {file_path}")
            missing_files.append(file_path)
    
    if missing_files:
        print(f"\n⚠️  Missing input files: {missing_files}")
        print("Please ensure these files are present before running the pipeline.")
        return False
    else:
        print("\n✅ All input files found")
        return True

def main():
    """Main setup function"""
    print("🚀 Setting up Unified ADC Data Enrichment Pipeline")
    print("=" * 60)
    
    # Install dependencies
    if not install_requirements():
        print("❌ Setup failed at dependency installation")
        sys.exit(1)
    
    # Check dependencies
    if not check_dependencies():
        print("❌ Setup failed at dependency check")
        sys.exit(1)
    
    # Create directories
    create_directories()
    
    # Check input files
    input_files_ok = check_input_files()
    
    print("\n" + "=" * 60)
    print("📋 SETUP SUMMARY")
    print("=" * 60)
    
    if input_files_ok:
        print("✅ Setup completed successfully!")
        print("\n🎯 Next steps:")
        print("   1. Run the pipeline: python run_pipeline.py")
        print("   2. Test the output: python test_pipeline.py")
        print("   3. Check the documentation: PIPELINE_README.md")
    else:
        print("⚠️  Setup completed with warnings")
        print("\n📝 Please ensure all input files are present before running the pipeline")
    
    print("=" * 60)

if __name__ == "__main__":
    main() 