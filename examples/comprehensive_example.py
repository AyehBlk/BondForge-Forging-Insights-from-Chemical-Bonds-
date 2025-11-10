#!/usr/bin/env python3
"""
BondForge - Comprehensive Example: Master Forge Demo
=====================================================

Forging Complete Insights from Chemical Bonds

This example demonstrates how to use both BondForge analyzers:
- Core Forge: 7 interaction types
- Master Forge: All 20 interaction types

Shows the full power of BondForge for protein interaction analysis.

Author: Ayeh Bolouki
Date: 2025-11-05
"""

import sys
import os

# Add src directory to path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))

from interaction_analyzer import ProteinInteractionAnalyzer
from extended_interaction_analyzer import ExtendedProteinInteractionAnalyzer
from interaction_visualizer import InteractionVisualizer


def example_core_analyzer(pdb_file):
    """
    Demonstrate the core forge with 7 interaction types.
    """
    print("\n" + "="*70)
    print(" BONDFORGE - CORE FORGE")
    print("7 Fundamental Interaction Types")
    print("="*70)
    
    # Initialize analyzer
    analyzer = ProteinInteractionAnalyzer(pdb_file)
    
    # Find hydrogen bonds
    print("\n1. Hydrogen Bonds:")
    h_bonds = analyzer.find_hydrogen_bonds()
    print(f"   Found {len(h_bonds)} hydrogen bonds")
    if h_bonds:
        print(f"   Example: {h_bonds[0]['donor']} → {h_bonds[0]['acceptor']}")
    
    # Find salt bridges
    print("\n2. Salt Bridges:")
    salt_bridges = analyzer.find_salt_bridges()
    print(f"   Found {len(salt_bridges)} salt bridges")
    if salt_bridges:
        print(f"   Example: {salt_bridges[0]['residue1']} ↔ {salt_bridges[0]['residue2']}")
    
    # Find disulfide bonds
    print("\n3. Disulfide Bonds:")
    disulfide = analyzer.find_disulfide_bonds()
    print(f"   Found {len(disulfide)} disulfide bonds")
    if disulfide:
        print(f"   Example: {disulfide[0]['residue1']} —S-S— {disulfide[0]['residue2']}")
    
    # Find hydrophobic interactions
    print("\n4. Hydrophobic Interactions:")
    hydrophobic = analyzer.find_hydrophobic_interactions()
    print(f"   Found {len(hydrophobic)} hydrophobic contacts")
    
    # Find pi-pi stacking
    print("\n5. Pi-Pi Stacking:")
    pi_pi = analyzer.find_pi_pi_stacking()
    print(f"   Found {len(pi_pi)} pi-pi stacking interactions")
    
    # Find cation-pi interactions
    print("\n6. Cation-Pi Interactions:")
    cation_pi = analyzer.find_cation_pi()
    print(f"   Found {len(cation_pi)} cation-pi interactions")
    
    # Find halogen bonds
    print("\n7. Halogen Bonds:")
    halogen = analyzer.find_halogen_bonds()
    print(f"   Found {len(halogen)} halogen bonds")
    
    # Identify hub residues
    print("\n8. Hub Residues (>5 interactions):")
    hubs = analyzer.identify_hub_residues(threshold=5)
    print(f"   Found {len(hubs)} hub residues")
    for hub in hubs[:5]:  # Show top 5
        print(f"   - {hub['residue']}: {hub['interaction_count']} interactions")
    
    # Export results
    output_file = 'core_interactions.csv'
    analyzer.export_to_csv(output_file)
    print(f"\n✓ Results exported to {output_file}")
    
    return analyzer


def example_extended_analyzer(pdb_file):
    """
    Demonstrate the master forge with all 20 interaction types.
    """
    print("\n" + "="*70)
    print(" BONDFORGE - MASTER FORGE")
    print("All 20 Interaction Types - Complete Analysis")
    print("="*70)
    
    # Initialize master forge
    analyzer = ExtendedProteinInteractionAnalyzer(pdb_file)
    
    # Forge complete analysis
    print("\n Forging all 20 interaction types...")
    results = analyzer.analyze_all_interactions()
    
    # Display summary
    print("\n" + "-"*70)
    print("INTERACTION SUMMARY")
    print("-"*70)
    
    interaction_types = [
        'hydrogen_bonds', 'salt_bridges', 'disulfide_bonds',
        'hydrophobic', 'pi_pi_stacking', 'cation_pi', 'halogen_bonds',
        'van_der_waals', 'anion_pi', 'sulfur_aromatic', 'ch_pi',
        'metal_coordination', 'carbonyl_pi', 'amide_aromatic',
        'sulfur_oxygen', 'backbone_carbonyl', 'aromatic_oxygen',
        'arginine_aromatic', 'backbone_amide_aromatic', 'backbone_carbonyl_charged'
    ]
    
    total = 0
    for i, itype in enumerate(interaction_types, 1):
        count = results['by_type'].get(itype, 0)
        total += count
        print(f"{i:2d}. {itype.replace('_', ' ').title():30s}: {count:5d}")
    
    print("-"*70)
    print(f"{'TOTAL INTERACTIONS':30s}: {total:5d}")
    print("-"*70)
    
    # Hub residues with extended analyzer
    print("\nHub Residues (highly connected):")
    hubs = analyzer.identify_hub_residues(threshold=10)
    print(f"Found {len(hubs)} hub residues with >10 interactions")
    for hub in hubs[:10]:  # Show top 10
        print(f"   {hub['residue']:8s}: {hub['interaction_count']:3d} interactions")
    
    # Critical interactions
    print("\nCritical Interactions:")
    critical = analyzer.identify_critical_interactions()
    print(f"Found {len(critical)} critical interactions")
    for interaction in critical[:5]:  # Show top 5
        print(f"   {interaction['residue1']} ↔ {interaction['residue2']}: "
              f"{interaction['type']}")
    
    # Export results
    output_file = 'extended_interactions.csv'
    analyzer.export_to_csv(output_file)
    print(f"\n✓ Complete results exported to {output_file}")
    
    return analyzer, results


def example_interface_analysis(pdb_file):
    """
    Demonstrate protein-protein interface analysis (for multi-chain structures).
    """
    print("\n" + "="*70)
    print("INTERFACE ANALYSIS - Protein-Protein Interactions")
    print("="*70)
    
    analyzer = ExtendedProteinInteractionAnalyzer(pdb_file)
    
    # Get available chains
    structure = analyzer.structure
    chains = [chain.id for chain in structure.get_chains()]
    
    print(f"\nAvailable chains: {', '.join(chains)}")
    
    if len(chains) >= 2:
        chain_A = chains[0]
        chain_B = chains[1]
        
        print(f"\nAnalyzing interface between chains {chain_A} and {chain_B}...")
        interface = analyzer.analyze_interface(chain_A, chain_B)
        
        print("\nInterface Statistics:")
        print(f"   Interface residues (chain {chain_A}): {len(interface['residues_A'])}")
        print(f"   Interface residues (chain {chain_B}): {len(interface['residues_B'])}")
        print(f"   Total interface interactions: {interface['interaction_count']}")
        print(f"   Interaction types at interface: {len(interface['interaction_types'])}")
        
        print("\n   Most common interactions at interface:")
        for itype, count in list(interface['interaction_types'].items())[:5]:
            print(f"      {itype}: {count}")
    else:
        print("\nNote: Interface analysis requires multi-chain structure")
        print("This structure has only one chain.")


def example_visualization(pdb_file):
    """
    Demonstrate visualization capabilities.
    """
    print("\n" + "="*70)
    print("VISUALIZATION - Network and Structure Plots")
    print("="*70)
    
    # Use extended analyzer for comprehensive visualization
    analyzer = ExtendedProteinInteractionAnalyzer(pdb_file)
    results = analyzer.analyze_all_interactions()
    
    # Create visualizer
    viz = InteractionVisualizer(analyzer)
    
    # Generate network plot
    print("\nGenerating interaction network plot...")
    try:
        viz.plot_interaction_network(results, output='interaction_network.png')
        print("✓ Network plot saved as 'interaction_network.png'")
    except Exception as e:
        print(f"   Note: Network plot generation requires matplotlib: {e}")
    
    # Generate PyMOL script
    print("\nGenerating PyMOL visualization script...")
    try:
        viz.generate_pymol_script('visualize_interactions.pml')
        print("✓ PyMOL script saved as 'visualize_interactions.pml'")
        print("   Load in PyMOL with: @visualize_interactions.pml")
    except Exception as e:
        print(f"   Note: PyMOL script generation failed: {e}")


def main():
    """
    Main function to run all BondForge examples.
    """
    print("\n" + "="*70)
    print("  BONDFORGE - COMPREHENSIVE DEMO")
    print(" Forging Insights from Chemical Bonds")
    print("="*70)
    
    # Check for PDB file argument
    if len(sys.argv) > 1:
        pdb_file = sys.argv[1]
    else:
        print("\nUsage: python comprehensive_example.py <pdb_file>")
        print("\nExample: python comprehensive_example.py 1abc.pdb")
        print("\nFor testing without a PDB file, you can download one from:")
        print("https://www.rcsb.org/structure/1ABC")
        return
    
    # Check if file exists
    if not os.path.exists(pdb_file):
        print(f"\nError: File '{pdb_file}' not found!")
        return
    
    print(f"\nAnalyzing structure: {pdb_file}")
    
    try:
        # Run all examples
        print("\n[1/4] Running core analyzer...")
        core_analyzer = example_core_analyzer(pdb_file)
        
        print("\n[2/4] Running extended analyzer...")
        extended_analyzer, results = example_extended_analyzer(pdb_file)
        
        print("\n[3/4] Running interface analysis...")
        example_interface_analysis(pdb_file)
        
        print("\n[4/4] Generating visualizations...")
        example_visualization(pdb_file)
        
        # Final summary
        print("\n" + "="*70)
        print("  FORGING COMPLETE!")
        print("="*70)
        print("\n Forged files:")
        print("   - core_interactions.csv        (Core forge: 7 types)")
        print("   - extended_interactions.csv    (Master forge: 20 types)")
        print("   - interaction_network.png      (Network visualization)")
        print("   - visualize_interactions.pml   (PyMOL script)")
        print("\n Next steps:")
        print("   1. Open CSV files in Excel for detailed analysis")
        print("   2. View network plot for interaction overview")
        print("   3. Use PyMOL script for 3D structure visualization")
        print("   4. Integrate BondForge into your own scripts")
        print("\n Keep forging insights from chemical bonds!")
        
    except Exception as e:
        print(f"\nError during analysis: {e}")
        import traceback
        traceback.print_exc()


if __name__ == '__main__':
    main()
