# 🔥 BondForge - Quick Start Guide

Forge your first insights in 5 minutes!

## Installation

### Option 1: Clone and Install (Recommended)

```bash
# Clone BondForge
git clone https://github.com/AyehBlk/BondForge.git
cd BondForge

# Create a virtual environment (optional but recommended)
python -m venv venv
source venv/bin/activate  # On Windows: venv\Scripts\activate

# Install dependencies
pip install -r requirements.txt

# Install BondForge
pip install -e .
```

### Option 2: Direct Installation

```bash
pip install git+https://github.com/AyehBlk/BondForge.git
```

## Your First Analysis

### Step 1: Get a PDB File

Download a sample protein structure:
```bash
# Using wget (Linux/Mac)
wget https://files.rcsb.org/download/1ABC.pdb

# Or using curl
curl -O https://files.rcsb.org/download/1ABC.pdb
```

Or use your own PDB file!

### Step 2: Basic Analysis

Create a file called `my_first_analysis.py`:

```python
from src.extended_interaction_analyzer import ExtendedProteinInteractionAnalyzer

# Load your protein structure
analyzer = ExtendedProteinInteractionAnalyzer('1ABC.pdb')

# Analyze all 20 interaction types
results = analyzer.analyze_all_interactions()

# Print summary
print(f"Total interactions found: {results['total_interactions']}")
print(f"\nInteraction breakdown:")
for interaction_type, count in results['by_type'].items():
    print(f"  {interaction_type}: {count}")

# Find important residues
hubs = analyzer.identify_hub_residues(threshold=5)
print(f"\nFound {len(hubs)} hub residues:")
for hub in hubs[:5]:
    print(f"  {hub['residue']}: {hub['interaction_count']} interactions")

# Export results
analyzer.export_to_csv('my_interactions.csv')
print("\n✓ Results saved to my_interactions.csv")
```

Run it:
```bash
python my_first_analysis.py
```

### Step 3: Visualize

```python
from src.extended_interaction_analyzer import ExtendedProteinInteractionAnalyzer
from src.interaction_visualizer import InteractionVisualizer

# Analyze
analyzer = ExtendedProteinInteractionAnalyzer('1ABC.pdb')
results = analyzer.analyze_all_interactions()

# Visualize
viz = InteractionVisualizer(analyzer)
viz.plot_interaction_network(results, output='network.png')
print("✓ Network visualization saved to network.png")
```

## Common Use Cases

### 1. Quick Protein Stability Check

```python
analyzer = ExtendedProteinInteractionAnalyzer('protein.pdb')

# Check stabilizing interactions
salt_bridges = analyzer.find_salt_bridges()
disulfide = analyzer.find_disulfide_bonds()
h_bonds = analyzer.find_hydrogen_bonds()

print(f"Stability factors:")
print(f"  Salt bridges: {len(salt_bridges)}")
print(f"  Disulfide bonds: {len(disulfide)}")
print(f"  H-bonds: {len(h_bonds)}")
```

### 2. Interface Analysis (Multi-chain)

```python
analyzer = ExtendedProteinInteractionAnalyzer('complex.pdb')

# Analyze interface between chains A and B
interface = analyzer.analyze_interface('A', 'B')

print(f"Interface interactions: {interface['interaction_count']}")
print(f"Key interaction types: {list(interface['interaction_types'].keys())}")
```

### 3. Binding Site Analysis

```python
# Analyze protein-ligand interactions
analyzer = ExtendedProteinInteractionAnalyzer('protein_ligand.pdb')
binding_site = analyzer.analyze_binding_site(ligand_chain='L')

print(f"Binding site residues: {len(binding_site['residues'])}")
print(f"Key interactions: {binding_site['interaction_types']}")
```

## Next Steps

1. **Read the docs**: Check `docs/` for detailed documentation
2. **Run examples**: See `examples/` for more complex usage
3. **Customize**: Adjust thresholds and parameters for your needs
4. **Contribute**: Found a bug or want to add features? See CONTRIBUTING.md

## Troubleshooting

### Common Issues

**Import Error**: Make sure you've installed all dependencies
```bash
pip install -r requirements.txt
```

**PDB Parsing Error**: Ensure your PDB file is properly formatted
```python
# Check if structure loaded correctly
print(f"Chains: {[c.id for c in analyzer.structure.get_chains()]}")
print(f"Residues: {len(list(analyzer.structure.get_residues()))}")
```

**No Interactions Found**: Check your distance thresholds
```python
# Adjust cutoffs if needed
h_bonds = analyzer.find_hydrogen_bonds(distance_cutoff=4.0)
```

## Getting Help

- **Documentation**: See `docs/` folder
- **Examples**: Check `examples/` folder  
- **Issues**: [Open an issue](https://github.com/AyehBlk/bondforge/issues)
- **Discussions**: [GitHub Discussions](https://github.com/AyehBlk/bondforge/discussions)

## Quick Reference

### Core Functions

```python
# Initialize
analyzer = ExtendedProteinInteractionAnalyzer('protein.pdb')

# Find specific interactions
h_bonds = analyzer.find_hydrogen_bonds()
salt_bridges = analyzer.find_salt_bridges()
pi_stacking = analyzer.find_pi_pi_stacking()

# Analyze all
results = analyzer.analyze_all_interactions()

# Find important residues
hubs = analyzer.identify_hub_residues()
critical = analyzer.identify_critical_interactions()

# Export
analyzer.export_to_csv('output.csv')
```

### Visualization

```python
from src.interaction_visualizer import InteractionVisualizer

viz = InteractionVisualizer(analyzer)
viz.plot_interaction_network(results, output='network.png')
viz.generate_pymol_script('structure.pml')
```

## What's Next?

Now that you've completed your first analysis:

1. Explore all 20 interaction types
2. Try different proteins and compare
3. Integrate into your research pipeline
4. Build custom analysis scripts
5. Share your findings!

Happy analyzing! 🧬✨
