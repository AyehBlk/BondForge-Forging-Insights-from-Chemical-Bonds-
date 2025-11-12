# Chemical Interaction Extraction Algorithms
## For Protein Analysis Software

---

## Table of Contents
1. [Core Interaction Detection Algorithms](#core-interaction-detection)
2. [Hub Interaction Detection](#hub-interaction-detection)
3. [Critical Interaction Identification](#critical-interaction-identification)
4. [Implementation Guidelines](#implementation-guidelines)

---

## 1. CORE INTERACTION DETECTION ALGORITHMS

### 1.1 Intra-Protein Interaction Detection

**Purpose**: Identify chemical interactions within a single protein structure

**Algorithm: INTRA_PROTEIN_INTERACTION_DETECTOR**

```
INPUT:
  - protein_structure: 3D coordinates of all atoms in protein
  - residue_list: List of all residues with their atoms
  - distance_cutoffs: Dictionary of interaction-specific distance thresholds
  
OUTPUT:
  - interaction_list: List of detected interactions with type and strength
  
PARAMETERS:
  distance_cutoffs = {
    'hydrogen_bond': 3.5 Å,
    'salt_bridge': 4.0 Å,
    'hydrophobic': 5.0 Å,
    'van_der_waals': 4.5 Å,
    'disulfide_bond': 2.5 Å,
    'pi_pi_stacking': 6.0 Å,
    'cation_pi': 6.0 Å
  }
  
  angle_cutoffs = {
    'hydrogen_bond_angle': 120-180°,  # Donor-H-Acceptor angle
    'salt_bridge_angle': 90-180°
  }

ALGORITHM:

1. INITIALIZE empty interaction_list

2. FOR each residue_i in residue_list:
   
   3. FOR each residue_j in residue_list where j > i:
      
      4. Calculate sequence_separation = |i - j|
      
      5. IF sequence_separation < 4:  # Skip nearby residues in sequence
         CONTINUE
      
      6. # Calculate minimum distance between residues
      7. min_distance = calculate_min_distance(residue_i, residue_j)
      
      8. IF min_distance > MAX(distance_cutoffs.values()):
         CONTINUE  # Too far for any interaction
      
      9. # Detect specific interaction types
      
      10. # A. HYDROGEN BONDS
      11. IF residue_i has H-bond donors OR acceptors:
          12. IF residue_j has H-bond acceptors OR donors:
              13. FOR each donor in [residue_i, residue_j]:
                  14. FOR each acceptor in [residue_i, residue_j]:
                      15. distance = calculate_distance(donor.H, acceptor)
                      16. IF distance <= distance_cutoffs['hydrogen_bond']:
                          17. angle = calculate_angle(donor-H-acceptor)
                          18. IF angle >= angle_cutoffs['hydrogen_bond_angle'][0]:
                              19. strength = calculate_hbond_strength(distance, angle)
                              20. ADD to interaction_list:
                                   {
                                     type: 'hydrogen_bond',
                                     residue1: residue_i,
                                     residue2: residue_j,
                                     atoms: [donor, acceptor],
                                     distance: distance,
                                     angle: angle,
                                     strength: strength
                                   }
      
      21. # B. SALT BRIDGES (Ionic interactions)
      22. IF residue_i is charged:
          23. IF residue_j is charged AND opposite_charge(residue_i, residue_j):
              24. distance = calculate_distance(charge_center_i, charge_center_j)
              25. IF distance <= distance_cutoffs['salt_bridge']:
                  26. strength = calculate_electrostatic_strength(distance, charges)
                  27. ADD to interaction_list:
                       {
                         type: 'salt_bridge',
                         residue1: residue_i,
                         residue2: residue_j,
                         distance: distance,
                         strength: strength
                       }
      
      28. # C. DISULFIDE BONDS
      29. IF residue_i.type == 'CYS' AND residue_j.type == 'CYS':
          30. SG_distance = calculate_distance(residue_i.SG, residue_j.SG)
          31. IF SG_distance <= distance_cutoffs['disulfide_bond']:
              32. ADD to interaction_list:
                   {
                     type: 'disulfide_bond',
                     residue1: residue_i,
                     residue2: residue_j,
                     distance: SG_distance,
                     strength: 'strong'  # Covalent bond
                   }
      
      33. # D. HYDROPHOBIC INTERACTIONS
      34. IF residue_i is hydrophobic AND residue_j is hydrophobic:
          35. distance = calculate_centroid_distance(residue_i, residue_j)
          36. IF distance <= distance_cutoffs['hydrophobic']:
              37. contact_area = calculate_contact_area(residue_i, residue_j)
              38. strength = contact_area / max_possible_area
              39. ADD to interaction_list:
                   {
                     type: 'hydrophobic',
                     residue1: residue_i,
                     residue2: residue_j,
                     distance: distance,
                     contact_area: contact_area,
                     strength: strength
                   }
      
      40. # E. PI-PI STACKING
      41. IF residue_i has aromatic_ring AND residue_j has aromatic_ring:
          42. distance = calculate_ring_centroid_distance(residue_i, residue_j)
          43. IF distance <= distance_cutoffs['pi_pi_stacking']:
              44. angle = calculate_ring_plane_angle(residue_i, residue_j)
              45. # Parallel stacking: angle ~0-30°, T-shaped: angle ~60-90°
              46. IF angle < 30 OR (angle > 60 AND angle < 90):
                  47. strength = calculate_pi_pi_strength(distance, angle)
                  48. ADD to interaction_list:
                       {
                         type: 'pi_pi_stacking',
                         residue1: residue_i,
                         residue2: residue_j,
                         distance: distance,
                         angle: angle,
                         geometry: 'parallel' if angle < 30 else 'T-shaped',
                         strength: strength
                       }
      
      49. # F. CATION-PI INTERACTIONS
      50. IF (residue_i is cation AND residue_j has aromatic) OR 
             (residue_j is cation AND residue_i has aromatic):
          51. cation = residue_i if residue_i is cation else residue_j
          52. aromatic = residue_j if residue_i is cation else residue_i
          53. distance = calculate_distance(cation.charge_center, aromatic.ring_center)
          54. IF distance <= distance_cutoffs['cation_pi']:
              55. strength = calculate_cation_pi_strength(distance)
              56. ADD to interaction_list:
                   {
                     type: 'cation_pi',
                     residue1: residue_i,
                     residue2: residue_j,
                     distance: distance,
                     strength: strength
                   }

57. RETURN interaction_list
```

---

### 1.2 Protein-Protein Interaction Detection

**Purpose**: Identify interactions at the interface between two protein chains

**Algorithm: PROTEIN_PROTEIN_INTERFACE_DETECTOR**

```
INPUT:
  - protein_chain1: Structure of first protein
  - protein_chain2: Structure of second protein
  - distance_cutoffs: Same as intra-protein
  
OUTPUT:
  - interface_interactions: List of interactions at interface
  - interface_residues: Set of residues at interface
  - interface_area: Buried surface area
  
ALGORITHM:

1. INITIALIZE empty interface_interactions list
2. INITIALIZE empty interface_residues_chain1 set
3. INITIALIZE empty interface_residues_chain2 set

4. # Step 1: Identify interface residues
5. FOR each residue_i in protein_chain1:
   6. FOR each residue_j in protein_chain2:
      7. min_distance = calculate_min_distance(residue_i, residue_j)
      8. IF min_distance <= 5.0 Å:  # Interface cutoff
         9. ADD residue_i to interface_residues_chain1
         10. ADD residue_j to interface_residues_chain2

11. # Step 2: Detect interactions between interface residues
12. FOR each residue_i in interface_residues_chain1:
    13. FOR each residue_j in interface_residues_chain2:
        
        14. # Use same interaction detection as intra-protein
        15. # But no sequence separation filter needed
        
        16. interactions = DETECT_ALL_INTERACTIONS(residue_i, residue_j, distance_cutoffs)
        
        17. FOR each interaction in interactions:
            18. ADD to interface_interactions

19. # Step 3: Calculate interface metrics
20. interface_area = calculate_buried_surface_area(
       protein_chain1, 
       protein_chain2, 
       interface_residues_chain1,
       interface_residues_chain2
    )

21. # Step 4: Classify interface hotspots
22. hotspot_residues = []
23. FOR each residue in interface_residues:
    24. interaction_count = COUNT interactions involving residue
    25. IF interaction_count >= 3:  # Hotspot threshold
        26. ADD residue to hotspot_residues

27. RETURN {
       interactions: interface_interactions,
       interface_residues: {chain1: interface_residues_chain1, 
                           chain2: interface_residues_chain2},
       interface_area: interface_area,
       hotspots: hotspot_residues
    }
```

---

### 1.3 Protein-Ligand Interaction Detection

**Purpose**: Identify interactions between protein and small molecule ligand

**Algorithm: PROTEIN_LIGAND_INTERACTION_DETECTOR**

```
INPUT:
  - protein_structure: 3D coordinates of protein
  - ligand_structure: 3D coordinates of ligand atoms
  - ligand_properties: Chemical properties of ligand (charges, groups, etc.)
  
OUTPUT:
  - binding_interactions: List of protein-ligand interactions
  - binding_residues: Residues in binding pocket
  - binding_affinity_estimate: Estimated binding strength
  
PARAMETERS:
  distance_cutoffs = {
    'hydrogen_bond': 3.5 Å,
    'hydrophobic': 4.5 Å,
    'pi_pi_stacking': 5.5 Å,
    'cation_pi': 5.5 Å,
    'halogen_bond': 4.0 Å,
    'metal_coordination': 3.0 Å,
    'van_der_waals': 4.0 Å
  }

ALGORITHM:

1. INITIALIZE empty binding_interactions list
2. INITIALIZE empty binding_residues set

3. # Step 1: Identify binding pocket residues
4. FOR each residue in protein_structure:
   5. min_distance = calculate_min_distance(residue, ligand_structure)
   6. IF min_distance <= 5.0 Å:  # Binding pocket cutoff
      7. ADD residue to binding_residues

8. # Step 2: Detect specific interactions
9. FOR each residue in binding_residues:
   
   10. # A. HYDROGEN BONDS
   11. FOR each protein_donor in residue.hbond_donors:
       12. FOR each ligand_acceptor in ligand_structure.acceptors:
           13. distance = calculate_distance(protein_donor.H, ligand_acceptor)
           14. IF distance <= distance_cutoffs['hydrogen_bond']:
               15. angle = calculate_angle(donor-H-acceptor)
               16. IF angle >= 120°:
                   17. strength = calculate_hbond_energy(distance, angle)
                   18. ADD to binding_interactions:
                        {
                          type: 'hydrogen_bond',
                          protein_residue: residue,
                          protein_atom: protein_donor,
                          ligand_atom: ligand_acceptor,
                          distance: distance,
                          angle: angle,
                          energy: strength
                        }
   
   19. # Similarly for ligand donors → protein acceptors
   
   20. # B. HYDROPHOBIC CONTACTS
   21. IF residue is hydrophobic:
       22. FOR each ligand_hydrophobic_atom in ligand_structure:
           23. distance = calculate_distance(residue.centroid, ligand_hydrophobic_atom)
           24. IF distance <= distance_cutoffs['hydrophobic']:
               25. contact_area = calculate_contact_surface(residue, ligand_hydrophobic_atom)
               26. ADD to binding_interactions:
                    {
                      type: 'hydrophobic',
                      protein_residue: residue,
                      ligand_atom: ligand_hydrophobic_atom,
                      distance: distance,
                      contact_area: contact_area
                    }
   
   27. # C. PI-PI STACKING
   28. IF residue has aromatic_ring:
       29. FOR each ligand_aromatic_ring in ligand_structure:
           30. distance = calculate_ring_distance(residue.ring, ligand_aromatic_ring)
           31. IF distance <= distance_cutoffs['pi_pi_stacking']:
               32. angle = calculate_ring_angle(residue.ring, ligand_aromatic_ring)
               33. IF angle < 30 OR (angle > 60 AND angle < 90):
                   34. ADD to binding_interactions
   
   35. # D. CATION-PI
   36. IF residue is cationic:
       37. FOR each ligand_aromatic_ring in ligand_structure:
           38. distance = calculate_distance(residue.charge_center, ligand_aromatic_ring.center)
           39. IF distance <= distance_cutoffs['cation_pi']:
               40. ADD to binding_interactions
   
   41. # E. HALOGEN BONDS (if ligand has halogens)
   42. IF ligand has halogen_atoms:
       43. FOR each halogen in ligand.halogen_atoms:
           44. FOR each acceptor in residue:
               45. distance = calculate_distance(halogen, acceptor)
               46. IF distance <= distance_cutoffs['halogen_bond']:
                   47. angle = calculate_halogen_bond_angle(halogen, acceptor)
                   48. IF angle >= 140°:  # Linear halogen bond
                       49. ADD to binding_interactions
   
   50. # F. METAL COORDINATION (if ligand or protein has metals)
   51. IF residue has metal_ion:
       52. FOR each ligand_donor in ligand_structure:
           53. distance = calculate_distance(metal_ion, ligand_donor)
           54. IF distance <= distance_cutoffs['metal_coordination']:
               55. ADD to binding_interactions

56. # Step 3: Calculate binding affinity estimate
57. binding_score = 0
58. FOR each interaction in binding_interactions:
    59. binding_score += interaction_energy_contribution(interaction)

60. # Step 4: Identify key binding residues
61. key_residues = []
62. FOR each residue in binding_residues:
    63. residue_contribution = SUM(energies of interactions involving residue)
    64. IF residue_contribution > threshold:
        65. ADD residue to key_residues

66. RETURN {
       interactions: binding_interactions,
       binding_residues: binding_residues,
       key_residues: key_residues,
       binding_score: binding_score,
       pocket_volume: calculate_pocket_volume(binding_residues)
    }
```

---

## 2. HUB INTERACTION DETECTION

**Purpose**: Identify residues that act as hubs with many connections

**Algorithm: HUB_RESIDUE_DETECTOR**

```
INPUT:
  - interaction_list: All detected interactions from above algorithms
  - protein_structure: Full protein structure
  - hub_threshold: Minimum connections to be considered a hub (default: 5)
  
OUTPUT:
  - hub_residues: List of hub residues with their connectivity
  - hub_network: Graph representation of hub interactions
  
ALGORITHM:

1. INITIALIZE residue_interaction_count dictionary
2. INITIALIZE residue_interaction_types dictionary

3. # Step 1: Count interactions per residue
4. FOR each interaction in interaction_list:
   5. residue1 = interaction.residue1
   6. residue2 = interaction.residue2
   7. interaction_type = interaction.type
   
   8. residue_interaction_count[residue1] += 1
   9. residue_interaction_count[residue2] += 1
   
   10. residue_interaction_types[residue1][interaction_type] += 1
   11. residue_interaction_types[residue2][interaction_type] += 1

12. # Step 2: Identify hub residues
13. hub_residues = []
14. FOR each residue in residue_interaction_count:
    15. IF residue_interaction_count[residue] >= hub_threshold:
        
        16. # Calculate hub metrics
        17. degree = residue_interaction_count[residue]
        18. diversity = COUNT(distinct interaction types for residue)
        19. avg_interaction_strength = AVERAGE(strengths of interactions)
        
        20. # Calculate hub centrality
        21. betweenness = calculate_betweenness_centrality(residue, interaction_network)
        22. closeness = calculate_closeness_centrality(residue, interaction_network)
        
        23. # Classify hub type
        24. IF diversity >= 3:
            25. hub_type = 'multi-functional_hub'
        26. ELSE IF residue is at protein surface:
            27. hub_type = 'interface_hub'
        28. ELSE:
            29. hub_type = 'structural_hub'
        
        30. ADD to hub_residues:
            {
              residue: residue,
              degree: degree,
              interaction_types: residue_interaction_types[residue],
              diversity: diversity,
              avg_strength: avg_interaction_strength,
              betweenness: betweenness,
              closeness: closeness,
              hub_type: hub_type,
              hub_score: calculate_hub_score(degree, diversity, centrality)
            }

31. # Step 3: Sort hubs by importance
32. SORT hub_residues by hub_score DESCENDING

33. # Step 4: Build hub network
34. hub_network = build_network_graph(hub_residues, interaction_list)

35. RETURN {
       hubs: hub_residues,
       network: hub_network,
       statistics: {
         total_hubs: LENGTH(hub_residues),
         avg_hub_degree: AVERAGE(degree of hubs),
         max_hub_degree: MAX(degree of hubs)
       }
    }

HELPER FUNCTION calculate_hub_score(degree, diversity, centrality):
    # Weighted score combining multiple metrics
    score = (0.4 * normalize(degree) + 
             0.3 * normalize(diversity) + 
             0.3 * centrality)
    RETURN score
```

---

## 3. CRITICAL INTERACTION IDENTIFICATION

**Purpose**: Identify interactions that are critical for protein stability or function

**Algorithm: CRITICAL_INTERACTION_DETECTOR**

```
INPUT:
  - interaction_list: All detected interactions
  - protein_structure: Full protein structure
  - secondary_structure: Alpha helices, beta sheets, loops
  - functional_sites: Active sites, binding sites (optional)
  
OUTPUT:
  - critical_interactions: Ranked list of critical interactions
  - stability_map: Map showing contribution to stability
  
ALGORITHM:

1. INITIALIZE critical_interactions list
2. INITIALIZE interaction_scores dictionary

3. # Step 1: Score each interaction by multiple criteria
4. FOR each interaction in interaction_list:
   
   5. score = 0
   6. criticality_factors = []
   
   7. # A. INTERACTION STRENGTH
   8. IF interaction.type == 'disulfide_bond':
       9. score += 10  # Covalent, very strong
       10. ADD 'covalent_bond' to criticality_factors
   
   11. ELSE IF interaction.type == 'salt_bridge':
       12. score += 7
       13. ADD 'strong_electrostatic' to criticality_factors
   
   14. ELSE IF interaction.type == 'hydrogen_bond':
       15. IF interaction.strength > 0.7:  # Strong H-bond
           16. score += 5
           17. ADD 'strong_hydrogen_bond' to criticality_factors
   
   18. # B. STRUCTURAL LOCATION
   19. IF interaction bridges two secondary structures:
       20. score += 8
       21. ADD 'structural_bridge' to criticality_factors
       22. # Example: Connects helix to beta sheet
   
   23. ELSE IF interaction is in protein core:
       24. score += 6
       25. ADD 'core_stabilization' to criticality_factors
   
   26. ELSE IF interaction is at domain interface:
       27. score += 7
       28. ADD 'domain_interface' to criticality_factors
   
   29. # C. SEQUENCE CONSERVATION
   30. IF both residues are highly conserved:  # From MSA if available
       31. conservation_score = calculate_conservation(interaction.residues)
       32. score += conservation_score * 5
       33. ADD 'evolutionary_conserved' to criticality_factors
   
   34. # D. BURIED SURFACE AREA
   35. IF interaction.contact_area > threshold:
       36. score += 4
       37. ADD 'large_interface' to criticality_factors
   
   38. # E. FUNCTIONAL IMPORTANCE
   39. IF interaction involves functional_site residues:
       40. score += 9
       41. ADD 'functional_site' to criticality_factors
   
   42. # F. NETWORK CENTRALITY
   43. IF either residue is a hub:
       44. score += 6
       45. ADD 'hub_interaction' to criticality_factors
   
   46. # G. STABILIZATION OF FOLD
   47. IF interaction reduces conformational entropy:
       48. entropy_reduction = calculate_entropy_reduction(interaction)
       49. score += entropy_reduction * 3
       50. ADD 'entropy_reduction' to criticality_factors
   
   51. # H. UNIQUE INTERACTIONS
   52. IF interaction.type is rare in protein:
       53. score += 3
       54. ADD 'unique_interaction' to criticality_factors
   
   55. interaction_scores[interaction] = score
   56. interaction.criticality_factors = criticality_factors

57. # Step 2: Rank interactions
58. ranked_interactions = SORT(interaction_scores) DESCENDING by score

59. # Step 3: Classify critical interactions
60. FOR each interaction in ranked_interactions:
    61. IF interaction_scores[interaction] >= 20:  # High criticality
        62. criticality_level = 'critical'
    63. ELSE IF interaction_scores[interaction] >= 12:  # Medium
        64. criticality_level = 'important'
    65. ELSE IF interaction_scores[interaction] >= 6:  # Low
        66. criticality_level = 'moderate'
    67. ELSE:
        68. criticality_level = 'minor'
    
    69. ADD to critical_interactions:
        {
          interaction: interaction,
          criticality_score: interaction_scores[interaction],
          criticality_level: criticality_level,
          factors: interaction.criticality_factors,
          predicted_effect_if_lost: predict_mutation_effect(interaction)
        }

70. # Step 4: Identify critical regions
71. critical_regions = []
72. clusters = cluster_interactions_by_location(critical_interactions)
73. FOR each cluster in clusters:
    74. IF cluster has >= 3 critical interactions:
        75. ADD cluster to critical_regions

76. # Step 5: Build stability map
77. stability_map = create_residue_stability_map(critical_interactions)

78. RETURN {
       critical_interactions: critical_interactions,
       critical_regions: critical_regions,
       stability_map: stability_map,
       statistics: {
         total_critical: COUNT(criticality_level == 'critical'),
         total_important: COUNT(criticality_level == 'important'),
         avg_criticality_score: AVERAGE(criticality_scores)
       }
    }

HELPER FUNCTION predict_mutation_effect(interaction):
    # Predict what happens if interaction is disrupted
    IF interaction.criticality_score >= 20:
        RETURN 'likely_destabilizing_or_loss_of_function'
    ELSE IF interaction.criticality_score >= 12:
        RETURN 'possibly_destabilizing'
    ELSE:
        RETURN 'minimal_effect'
```

---

## 4. IMPLEMENTATION GUIDELINES

### 4.1 Required Python Libraries

```python
# Core libraries
import numpy as np
import pandas as pd
from scipy.spatial import distance, cKDTree
from scipy.spatial.distance import cdist
import networkx as nx

# Bioinformatics libraries
from Bio.PDB import PDBParser, PDBIO, NeighborSearch
from Bio.PDB.DSSP import DSSP
from Bio.PDB.SASA import ShrakeRupley

# Molecular analysis
import MDAnalysis as mda  # Optional, for MD trajectories
from rdkit import Chem  # For ligand handling

# Visualization
import matplotlib.pyplot as plt
import seaborn as sns
```

### 4.2 Key Helper Functions to Implement

```python
def calculate_distance(atom1, atom2):
    """Calculate Euclidean distance between two atoms"""
    return np.linalg.norm(atom1.coord - atom2.coord)

def calculate_angle(atom1, atom2, atom3):
    """Calculate angle formed by three atoms"""
    v1 = atom1.coord - atom2.coord
    v2 = atom3.coord - atom2.coord
    cosine = np.dot(v1, v2) / (np.linalg.norm(v1) * np.linalg.norm(v2))
    return np.degrees(np.arccos(np.clip(cosine, -1, 1)))

def calculate_buried_surface_area(protein1, protein2):
    """Calculate interface buried surface area"""
    # Use Shrake-Rupley or FreeSASA
    sr = ShrakeRupley()
    sr.compute(protein1, level="R")
    sasa_alone = sum([r.sasa for r in protein1.get_residues()])
    
    # Now calculate together
    complex_structure = combine_structures(protein1, protein2)
    sr.compute(complex_structure, level="R")
    sasa_together = sum([r.sasa for r in complex_structure.get_residues()])
    
    buried_area = sasa_alone - sasa_together
    return buried_area

def is_charged_residue(residue):
    """Check if residue is charged"""
    charged_residues = {
        'ARG': '+', 'LYS': '+',  # Positive
        'ASP': '-', 'GLU': '-'   # Negative
    }
    return residue.resname in charged_residues

def is_hydrophobic_residue(residue):
    """Check if residue is hydrophobic"""
    hydrophobic_residues = {
        'ALA', 'VAL', 'ILE', 'LEU', 'MET', 
        'PHE', 'TRP', 'PRO', 'TYR'
    }
    return residue.resname in hydrophobic_residues

def has_aromatic_ring(residue):
    """Check if residue has aromatic ring"""
    aromatic_residues = {'PHE', 'TYR', 'TRP', 'HIS'}
    return residue.resname in aromatic_residues

def calculate_betweenness_centrality(residue, network_graph):
    """Calculate betweenness centrality for residue in network"""
    betweenness = nx.betweenness_centrality(network_graph)
    return betweenness.get(residue, 0)
```

### 4.3 Data Structures

```python
# Interaction object
class Interaction:
    def __init__(self):
        self.type = None  # 'hydrogen_bond', 'salt_bridge', etc.
        self.residue1 = None
        self.residue2 = None
        self.atoms = []
        self.distance = None
        self.angle = None
        self.strength = None
        self.energy = None
        self.criticality_score = 0
        
# Hub residue object
class HubResidue:
    def __init__(self):
        self.residue = None
        self.degree = 0
        self.interaction_types = {}
        self.diversity = 0
        self.hub_score = 0
        self.hub_type = None
        
# Critical interaction object  
class CriticalInteraction:
    def __init__(self):
        self.interaction = None
        self.criticality_score = 0
        self.criticality_level = None
        self.factors = []
        self.predicted_effect = None
```

### 4.4 Integration with GUI

```python
class InteractionAnalyzer:
    """Main class to be integrated with your GUI"""
    
    def __init__(self, pdb_file):
        self.structure = self.load_structure(pdb_file)
        self.interactions = []
        self.hubs = []
        self.critical_interactions = []
    
    def analyze_all(self):
        """Run all analyses"""
        self.detect_intra_protein_interactions()
        self.detect_protein_protein_interactions()
        self.detect_ligand_interactions()
        self.identify_hubs()
        self.identify_critical_interactions()
        
    def detect_intra_protein_interactions(self):
        """Implement INTRA_PROTEIN_INTERACTION_DETECTOR algorithm"""
        pass
    
    def detect_protein_protein_interactions(self):
        """Implement PROTEIN_PROTEIN_INTERFACE_DETECTOR algorithm"""
        pass
    
    def detect_ligand_interactions(self):
        """Implement PROTEIN_LIGAND_INTERACTION_DETECTOR algorithm"""
        pass
    
    def identify_hubs(self):
        """Implement HUB_RESIDUE_DETECTOR algorithm"""
        pass
    
    def identify_critical_interactions(self):
        """Implement CRITICAL_INTERACTION_DETECTOR algorithm"""
        pass
    
    def export_results(self, format='csv'):
        """Export results in various formats"""
        pass
    
    def visualize_interactions(self):
        """Create visualization of interactions"""
        pass
```

### 4.5 Performance Optimization Tips

1. **Use spatial indexing**: Use `scipy.spatial.cKDTree` for fast neighbor searches
2. **Distance cutoff prefiltering**: Skip distant residue pairs early
3. **Parallel processing**: Use `multiprocessing` for independent calculations
4. **Vectorization**: Use NumPy operations instead of loops where possible
5. **Caching**: Cache calculated distances and angles

```python
# Example optimization with cKDTree
from scipy.spatial import cKDTree

def find_neighbors_fast(atoms, cutoff=5.0):
    """Fast neighbor search using KD-tree"""
    coords = np.array([atom.coord for atom in atoms])
    tree = cKDTree(coords)
    pairs = tree.query_pairs(cutoff)
    return pairs
```

### 4.6 Output Format Recommendations

```python
# For CSV export
interaction_df = pd.DataFrame({
    'Type': [i.type for i in interactions],
    'Residue1': [i.residue1 for i in interactions],
    'Residue2': [i.residue2 for i in interactions],
    'Distance': [i.distance for i in interactions],
    'Strength': [i.strength for i in interactions],
    'Criticality': [i.criticality_score for i in interactions]
})

# For JSON export (for visualization)
output = {
    'interactions': [i.to_dict() for i in interactions],
    'hubs': [h.to_dict() for h in hubs],
    'critical_interactions': [c.to_dict() for c in critical_interactions],
    'statistics': calculate_statistics()
}
```

---

## 5. VALIDATION AND TESTING

### 5.1 Test Cases

```python
# Test with known structures
test_cases = {
    '1A2K': 'Barnase-Barstar complex (well-studied PPI)',
    '1PPE': 'Trypsin-BPTI complex (protease-inhibitor)',
    '3HFM': 'HIV protease (homodimer)',
    '1HVR': 'HIV protease with ligand'
}

# Expected results for validation
expected_1A2K = {
    'salt_bridges': 5-7,
    'hydrogen_bonds': 10-15,
    'interface_area': 1400-1600  # Å²
}
```

### 5.2 Quality Checks

1. Distance constraints satisfied for all interactions
2. Angle constraints satisfied for directional interactions
3. No duplicate interactions reported
4. Interaction counts within reasonable biological ranges
5. Hub residues match known literature (for test cases)

---

## NEXT STEPS FOR IMPLEMENTATION

1. Start with the simplest algorithm (hydrogen bond detection)
2. Test thoroughly with known structures
3. Add interaction types incrementally
4. Implement hub detection once basic interactions work
5. Add critical interaction identification last
6. Integrate visualization
7. Build user-friendly GUI wrapper

This design gives you a complete blueprint for implementation. Each algorithm is modular, so you can build and test them independently!
