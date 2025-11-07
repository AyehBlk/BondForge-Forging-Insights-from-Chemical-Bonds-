"""
BondForge - Basic Demo Script
==============================

Forging Your First Insights

This script demonstrates how to use BondForge to analyze protein interactions.
Shows basic usage of the core analyzer (7 interaction types).

Usage:
    python example_demo.py path/to/protein.pdb

Author: BondForge Team
Date: 2025
"""

import sys
import os
from interaction_analyzer import InteractionAnalyzer
from interaction_visualizer import InteractionVisualizer


def run_demo(pdb_file):
    """
    Run complete demo analysis
    
    Parameters:
    -----------
    pdb_file : str
        Path to PDB file
    """
    
    print("="*70)
    print(" BONDFORGE - BASIC DEMO")
    print("Forging Insights from Chemical Bonds")
    print("="*70)
    print(f"\nAnalyzing: {pdb_file}\n")
    
    # Step 1: Initialize Analyzer
    print("Step 1: Loading structure...")
    try:
        analyzer = InteractionAnalyzer(pdb_file)
    except Exception as e:
        print(f"Error loading file: {e}")
        return
    
    # Step 2: Run Analysis
    print("\nStep 2: Running interaction analysis...")
    print("This may take a few minutes for large proteins...")
    analyzer.analyze_all()
    
    # Step 3: Display Summary
    print("\n" + "="*70)
    analyzer.print_summary()
    
    # Step 4: Show Top Results
    print("\n" + "="*70)
    print("TOP RESULTS")
    print("="*70)
    
    # Show top hydrogen bonds
    print("\n Top 5 Hydrogen Bonds:")
    hbonds = [i for i in analyzer.interactions if i.type == 'hydrogen_bond']
    hbonds_sorted = sorted(hbonds, key=lambda x: x.strength if x.strength else 0, reverse=True)
    
    for i, hb in enumerate(hbonds_sorted[:5], 1):
        res1 = f"{hb.residue1_id[2]}{hb.residue1_id[1]}"
        res2 = f"{hb.residue2_id[2]}{hb.residue2_id[1]}"
        print(f"  {i}. {res1} ↔ {res2}")
        print(f"     Distance: {hb.distance:.2f} Å, Strength: {hb.strength:.3f}")
    
    # Show hub residues
    if analyzer.hubs:
        print("\n Top 5 Hub Residues:")
        for i, hub in enumerate(analyzer.hubs[:5], 1):
            res = f"{hub.residue_id[2]}{hub.residue_id[1]}"
            print(f"  {i}. {res}")
            print(f"     Degree: {hub.degree}, Diversity: {hub.diversity}, Score: {hub.hub_score:.3f}")
            print(f"     Types: {hub.interaction_types}")
    
    # Show critical interactions
    if analyzer.critical_interactions:
        print("\n Top 5 Critical Interactions:")
        for i, critical in enumerate(analyzer.critical_interactions[:5], 1):
            interaction = critical['interaction']
            res1 = f"{interaction.residue1_id[2]}{interaction.residue1_id[1]}"
            res2 = f"{interaction.residue2_id[2]}{interaction.residue2_id[1]}"
            print(f"  {i}. {interaction.type}: {res1} - {res2}")
            print(f"     Score: {critical['criticality_score']}, Level: {critical['level']}")
            print(f"     Factors: {', '.join(critical['factors'])}")
    
    # Step 5: Export Results
    print("\n" + "="*70)
    print("Step 3: Exporting results...")
    
    output_dir = './demo_results'
    os.makedirs(output_dir, exist_ok=True)
    
    analyzer.export_to_csv(output_dir)
    print(f"✓ CSV files exported to {output_dir}/")
    
    # Step 6: Create Visualizations
    print("\nStep 4: Creating visualizations...")
    
    visualizer = InteractionVisualizer(analyzer)
    viz_dir = os.path.join(output_dir, 'visualizations')
    os.makedirs(viz_dir, exist_ok=True)
    
    try:
        visualizer.create_all_visualizations(output_dir=viz_dir)
        print(f"✓ Visualizations exported to {viz_dir}/")
    except Exception as e:
        print(f"Warning: Some visualizations could not be created: {e}")
    
    # Step 7: Summary
    print("\n" + "="*70)
    print("DEMO COMPLETE!")
    print("="*70)
    print(f"\nResults saved to: {output_dir}/")
    print("\nGenerated files:")
    print("  - interactions.csv          : All detected interactions")
    print("  - hubs.csv                  : Hub residue analysis")
    print("  - critical_interactions.csv : Critical interactions")
    print("  - visualizations/           : All plots and figures")
    print("\nNext steps:")
    print("  1. Review the CSV files in Excel or any spreadsheet software")
    print("  2. View the visualizations (PNG images)")
    print("  3. Use PyMOL/ChimeraX to visualize in 3D")
    print("  4. Integrate into your GUI application")
    print("\n" + "="*70 + "\n")


def print_usage():
    """Print usage instructions"""
    print("\nUsage: python example_demo.py <pdb_file>")
    print("\nExample:")
    print("  python example_demo.py protein.pdb")
    print("\nOr download a test structure from PDB:")
    print("  wget https://files.rcsb.org/download/1A2K.pdb")
    print("  python example_demo.py 1A2K.pdb")
    print()


if __name__ == "__main__":
    
    # Check command line arguments
    if len(sys.argv) != 2:
        print_usage()
        sys.exit(1)
    
    pdb_file = sys.argv[1]
    
    # Check if file exists
    if not os.path.exists(pdb_file):
        print(f"\nError: File '{pdb_file}' not found!")
        print_usage()
        sys.exit(1)
    
    # Run demo
    try:
        run_demo(pdb_file)
    except KeyboardInterrupt:
        print("\n\nDemo interrupted by user.")
        sys.exit(0)
    except Exception as e:
        print(f"\n\nError during analysis: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)
