#!/usr/bin/env python3
"""
Script to visualize the percentage distribution of unique locus_type and gene_group values
from the HGNC-enriched AACR article data.
"""

import json
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from collections import Counter, defaultdict
import numpy as np

# Set style for better visualizations
plt.style.use('seaborn-v0_8')
sns.set_palette("husl")

def load_json_data(file_path):
    """Load the JSON data from the specified file."""
    with open(file_path, 'r') as f:
        return json.load(f)

def extract_unique_target_entries(data):
    """Extract unique target entries (drugs) from the JSON structure."""
    unique_targets = []
    
    for entry in data:
        for drug in entry.get("extractedDrugs", []):
            drug_name = drug.get("drugName", "Unknown")
            target_antigens = drug.get("targetAntigenCanonicalized", [])
            
            # Create a unique identifier for this drug-target combination
            if isinstance(target_antigens, list):
                for antigen in target_antigens:
                    unique_targets.append({
                        "drug_name": drug_name,
                        "target_antigen": antigen,
                        "hgnc_data": drug.get("HGNC", [])
                    })
            elif isinstance(target_antigens, str):
                unique_targets.append({
                    "drug_name": drug_name,
                    "target_antigen": target_antigens,
                    "hgnc_data": drug.get("HGNC", [])
                })
    
    return unique_targets

def analyze_locus_types(unique_targets):
    """Analyze the distribution of locus_type values for unique target entries."""
    locus_types = []
    
    for target_entry in unique_targets:
        hgnc_data = target_entry["hgnc_data"]
        if hgnc_data is not None:
            # Find the HGNC entry that matches this target antigen
            target_antigen = target_entry["target_antigen"]
            matching_hgnc = None
            
            for hgnc_entry in hgnc_data:
                if hgnc_entry and hgnc_entry.get("input") == target_antigen:
                    matching_hgnc = hgnc_entry
                    break
            
            if matching_hgnc:
                locus_types.append(matching_hgnc.get("locus_type"))
            else:
                locus_types.append(None)
        else:
            locus_types.append(None)
    
    # Count all values including None
    locus_type_counts = Counter(locus_types)
    
    # Calculate percentages
    total = len(locus_types)
    locus_type_percentages = {k: (v/total)*100 for k, v in locus_type_counts.items()}
    
    return locus_type_counts, locus_type_percentages, total

def analyze_gene_groups(unique_targets):
    """Analyze the distribution of gene_group values for unique target entries."""
    gene_groups = []
    empty_gene_groups = 0  # Count entries with empty gene_group lists
    
    for target_entry in unique_targets:
        hgnc_data = target_entry["hgnc_data"]
        if hgnc_data is not None:
            # Find the HGNC entry that matches this target antigen
            target_antigen = target_entry["target_antigen"]
            matching_hgnc = None
            
            for hgnc_entry in hgnc_data:
                if hgnc_entry and hgnc_entry.get("input") == target_antigen:
                    matching_hgnc = hgnc_entry
                    break
            
            if matching_hgnc:
                gene_group_list = matching_hgnc.get("gene_group", [])
                if isinstance(gene_group_list, list):
                    if len(gene_group_list) == 0:
                        empty_gene_groups += 1
                    else:
                        gene_groups.extend(gene_group_list)
                else:
                    empty_gene_groups += 1
            else:
                empty_gene_groups += 1
        else:
            empty_gene_groups += 1
    
    gene_group_counts = Counter(gene_groups)
    
    # Add empty gene groups to the counts
    if empty_gene_groups > 0:
        gene_group_counts["(Empty/None)"] = empty_gene_groups
    
    # Calculate percentages
    total = len(gene_groups) + empty_gene_groups
    gene_group_percentages = {k: (v/total)*100 for k, v in gene_group_counts.items()}
    
    return gene_group_counts, gene_group_percentages, total

def create_locus_type_visualization(locus_type_counts, locus_type_percentages, total):
    """Create visualization for locus_type distribution."""
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 8))
    
    # Sort by count for better visualization
    sorted_items = sorted(locus_type_counts.items(), key=lambda x: x[1], reverse=True)
    labels, counts = zip(*sorted_items)
    percentages = [locus_type_percentages[label] for label in labels]
    
    # Bar chart
    colors = ['red' if label is None else 'skyblue' for label in labels]
    bars = ax1.bar(range(len(labels)), counts, color=colors, alpha=0.7)
    ax1.set_title(f'Locus Type Distribution (Total: {total})', fontsize=14, fontweight='bold')
    ax1.set_xlabel('Locus Type')
    ax1.set_ylabel('Count')
    ax1.set_xticks(range(len(labels)))
    # Handle None labels
    display_labels = ['None' if label is None else label for label in labels]
    ax1.set_xticklabels(display_labels, rotation=45, ha='right')
    
    # Add value labels on bars
    for bar, count, pct in zip(bars, counts, percentages):
        height = bar.get_height()
        ax1.text(bar.get_x() + bar.get_width()/2., height + 0.5,
                f'{count}\n({pct:.1f}%)', ha='center', va='bottom', fontsize=10)
    
    # Pie chart
    if len(labels) <= 10:  # Only show pie chart if not too many categories
        pie_labels = ['None' if label is None else label for label in labels]
        ax2.pie(counts, labels=pie_labels, autopct='%1.1f%%', startangle=90)
        ax2.set_title('Locus Type Distribution (%)', fontsize=14, fontweight='bold')
    else:
        # Show top 10 in pie chart
        top_10_labels = labels[:10]
        top_10_counts = counts[:10]
        other_count = sum(counts[10:])
        
        if other_count > 0:
            pie_labels = ['None' if label is None else label for label in top_10_labels] + ['Others']
            pie_counts = list(top_10_counts) + [other_count]
        else:
            pie_labels = ['None' if label is None else label for label in top_10_labels]
            pie_counts = top_10_counts
            
        ax2.pie(pie_counts, labels=pie_labels, autopct='%1.1f%%', startangle=90)
        ax2.set_title('Top 10 Locus Types (%)', fontsize=14, fontweight='bold')
    
    plt.tight_layout()
    return fig

def create_gene_group_visualization(gene_group_counts, gene_group_percentages, total):
    """Create visualization for gene_group distribution."""
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(20, 8))
    
    # Sort by count for better visualization
    sorted_items = sorted(gene_group_counts.items(), key=lambda x: x[1], reverse=True)
    labels, counts = zip(*sorted_items)
    percentages = [gene_group_percentages[label] for label in labels]
    
    # Bar chart (top 20)
    top_20_labels = labels[:20]
    top_20_counts = counts[:20]
    top_20_percentages = percentages[:20]
    
    # Color empty/None values differently
    colors = ['red' if label == "(Empty/None)" else 'lightcoral' for label in top_20_labels]
    bars = ax1.bar(range(len(top_20_labels)), top_20_counts, color=colors, alpha=0.7)
    ax1.set_title(f'Top 20 Gene Groups (Total: {total})', fontsize=14, fontweight='bold')
    ax1.set_xlabel('Gene Group')
    ax1.set_ylabel('Count')
    ax1.set_xticks(range(len(top_20_labels)))
    ax1.set_xticklabels(top_20_labels, rotation=45, ha='right')
    
    # Add value labels on bars
    for bar, count, pct in zip(bars, top_20_counts, top_20_percentages):
        height = bar.get_height()
        ax1.text(bar.get_x() + bar.get_width()/2., height + 0.5,
                f'{count}\n({pct:.1f}%)', ha='center', va='bottom', fontsize=8)
    
    # Pie chart (top 10)
    top_10_labels = labels[:10]
    top_10_counts = counts[:10]
    other_count = sum(counts[10:])
    
    if other_count > 0:
        pie_labels = list(top_10_labels) + ['Others']
        pie_counts = list(top_10_counts) + [other_count]
    else:
        pie_labels = top_10_labels
        pie_counts = top_10_counts
        
    ax2.pie(pie_counts, labels=pie_labels, autopct='%1.1f%%', startangle=90)
    ax2.set_title('Top 10 Gene Groups (%)', fontsize=14, fontweight='bold')
    
    plt.tight_layout()
    return fig

def print_summary_statistics(locus_type_counts, locus_type_percentages, 
                           gene_group_counts, gene_group_percentages):
    """Print summary statistics."""
    print("=" * 80)
    print("HGNC DATA ANALYSIS SUMMARY (Unique Target Entries)")
    print("=" * 80)
    
    print(f"\n📊 LOCUS TYPE ANALYSIS:")
    print(f"Total unique locus types: {len(locus_type_counts)}")
    print(f"Total unique target entries: {sum(locus_type_counts.values())}")
    
    print(f"\nTop 10 Locus Types by Count:")
    for i, (locus_type, count) in enumerate(sorted(locus_type_counts.items(), 
                                                   key=lambda x: x[1], reverse=True)[:10], 1):
        percentage = locus_type_percentages[locus_type]
        display_name = "None" if locus_type is None else locus_type
        print(f"  {i:2d}. {display_name:<30} {count:>5} ({percentage:>5.1f}%)")
    
    print(f"\n📊 GENE GROUP ANALYSIS:")
    print(f"Total unique gene groups: {len(gene_group_counts)}")
    print(f"Total gene group entries: {sum(gene_group_counts.values())}")
    
    print(f"\nTop 10 Gene Groups by Count:")
    for i, (gene_group, count) in enumerate(sorted(gene_group_counts.items(), 
                                                   key=lambda x: x[1], reverse=True)[:10], 1):
        percentage = gene_group_percentages[gene_group]
        display_name = gene_group if gene_group != "(Empty/None)" else "(Empty/None)"
        print(f"  {i:2d}. {display_name:<50} {count:>5} ({percentage:>5.1f}%)")
    
    print("\n" + "=" * 80)

def main():
    """Main function to run the analysis."""
    # Configuration
    input_file = "aacrArticle_hgnc.json"
    
    print("🔄 Loading data...")
    data = load_json_data(input_file)
    
    print("📊 Extracting unique target entries...")
    unique_targets = extract_unique_target_entries(data)
    print(f"✅ Found {len(unique_targets)} unique target entries")
    
    print("🔍 Analyzing locus types...")
    locus_type_counts, locus_type_percentages, locus_total = analyze_locus_types(unique_targets)
    
    print("🔍 Analyzing gene groups...")
    gene_group_counts, gene_group_percentages, gene_group_total = analyze_gene_groups(unique_targets)
    
    # Print summary statistics
    print_summary_statistics(locus_type_counts, locus_type_percentages,
                           gene_group_counts, gene_group_percentages)
    
    # Create visualizations
    print("📈 Creating visualizations...")
    
    # Locus type visualization
    fig1 = create_locus_type_visualization(locus_type_counts, locus_type_percentages, locus_total)
    fig1.savefig('locus_type_distribution.png', dpi=300, bbox_inches='tight')
    print("✅ Saved locus_type_distribution.png")
    
    # Gene group visualization
    fig2 = create_gene_group_visualization(gene_group_counts, gene_group_percentages, gene_group_total)
    fig2.savefig('gene_group_distribution.png', dpi=300, bbox_inches='tight')
    print("✅ Saved gene_group_distribution.png")
    
    # Show plots
    plt.show()
    
    print("🎉 Analysis complete!")

if __name__ == "__main__":
    main() 