"""
BondForge - Core Interaction Analysis Module
=============================================

Forging Insights from Chemical Bonds

This module implements the core algorithms for detecting and analyzing
7 types of chemical interactions in protein structures.

Part of BondForge: Comprehensive protein interaction analysis toolkit

Author: Your Name
Date: 2025
"""

import numpy as np
import pandas as pd
from scipy.spatial import distance, cKDTree
from scipy.spatial.distance import cdist
import networkx as nx
from typing import List, Dict, Tuple, Set
from dataclasses import dataclass
from Bio.PDB import PDBParser, PDBIO, NeighborSearch, Selection
from Bio.PDB.SASA import ShrakeRupley
import warnings
warnings.filterwarnings('ignore')


# ============================================================================
# DATA STRUCTURES
# ============================================================================

@dataclass
class Interaction:
    """Represents a chemical interaction between residues"""
    type: str  # 'hydrogen_bond', 'salt_bridge', etc.
    residue1_id: tuple  # (chain, resseq, resname)
    residue2_id: tuple
    atoms: List[str]  # Atom names involved
    distance: float
    angle: float = None
    strength: float = None
    energy: float = None
    criticality_score: float = 0.0
    criticality_factors: List[str] = None
    
    def to_dict(self):
        """Convert to dictionary for export"""
        return {
            'type': self.type,
            'residue1': f"{self.residue1_id[0]}:{self.residue1_id[2]}{self.residue1_id[1]}",
            'residue2': f"{self.residue2_id[0]}:{self.residue2_id[2]}{self.residue2_id[1]}",
            'atoms': ','.join(self.atoms),
            'distance': round(self.distance, 2),
            'angle': round(self.angle, 2) if self.angle else None,
            'strength': round(self.strength, 3) if self.strength else None,
            'criticality_score': round(self.criticality_score, 2)
        }


@dataclass
class HubResidue:
    """Represents a hub residue with many interactions"""
    residue_id: tuple
    degree: int
    interaction_types: Dict[str, int]
    diversity: int
    hub_score: float
    hub_type: str
    betweenness: float = 0.0
    closeness: float = 0.0
    
    def to_dict(self):
        return {
            'residue': f"{self.residue_id[0]}:{self.residue_id[2]}{self.residue_id[1]}",
            'degree': self.degree,
            'diversity': self.diversity,
            'hub_score': round(self.hub_score, 3),
            'hub_type': self.hub_type,
            'betweenness': round(self.betweenness, 3),
            'interaction_types': self.interaction_types
        }


# ============================================================================
# RESIDUE CLASSIFICATION HELPERS
# ============================================================================

CHARGED_RESIDUES = {
    'ARG': '+', 'LYS': '+',  # Positive
    'ASP': '-', 'GLU': '-'   # Negative
}

HYDROPHOBIC_RESIDUES = {
    'ALA', 'VAL', 'ILE', 'LEU', 'MET', 
    'PHE', 'TRP', 'PRO', 'CYS'
}

AROMATIC_RESIDUES = {'PHE', 'TYR', 'TRP', 'HIS'}

POLAR_RESIDUES = {'SER', 'THR', 'ASN', 'GLN', 'TYR', 'CYS'}

HBOND_DONORS = {
    'ARG': ['NH1', 'NH2', 'NE'],
    'ASN': ['ND2'],
    'GLN': ['NE2'],
    'HIS': ['ND1', 'NE2'],
    'LYS': ['NZ'],
    'SER': ['OG'],
    'THR': ['OG1'],
    'TRP': ['NE1'],
    'TYR': ['OH'],
    'CYS': ['SG'],
    'BACKBONE': ['N']  # All residues have backbone N
}

HBOND_ACCEPTORS = {
    'ASP': ['OD1', 'OD2'],
    'GLU': ['OE1', 'OE2'],
    'ASN': ['OD1'],
    'GLN': ['OE1'],
    'HIS': ['ND1', 'NE2'],
    'SER': ['OG'],
    'THR': ['OG1'],
    'TYR': ['OH'],
    'CYS': ['SG'],
    'BACKBONE': ['O']  # All residues have backbone O
}


def is_charged_residue(residue):
    """Check if residue is charged"""
    return residue.get_resname() in CHARGED_RESIDUES


def get_charge(residue):
    """Get charge of residue"""
    return CHARGED_RESIDUES.get(residue.get_resname(), None)


def is_hydrophobic_residue(residue):
    """Check if residue is hydrophobic"""
    return residue.get_resname() in HYDROPHOBIC_RESIDUES


def is_aromatic_residue(residue):
    """Check if residue has aromatic ring"""
    return residue.get_resname() in AROMATIC_RESIDUES


def is_polar_residue(residue):
    """Check if residue is polar"""
    return residue.get_resname() in POLAR_RESIDUES


# ============================================================================
# GEOMETRY CALCULATION HELPERS
# ============================================================================

def calculate_distance(atom1, atom2):
    """Calculate Euclidean distance between two atoms"""
    return np.linalg.norm(atom1.coord - atom2.coord)


def calculate_angle(atom1, atom2, atom3):
    """
    Calculate angle formed by three atoms (in degrees)
    atom2 is the vertex of the angle
    """
    v1 = atom1.coord - atom2.coord
    v2 = atom3.coord - atom2.coord
    
    # Normalize vectors
    v1_norm = v1 / np.linalg.norm(v1)
    v2_norm = v2 / np.linalg.norm(v2)
    
    # Calculate angle
    cosine = np.dot(v1_norm, v2_norm)
    cosine = np.clip(cosine, -1, 1)  # Handle numerical errors
    angle = np.degrees(np.arccos(cosine))
    
    return angle


def calculate_dihedral(atom1, atom2, atom3, atom4):
    """Calculate dihedral angle formed by four atoms"""
    v1 = atom2.coord - atom1.coord
    v2 = atom3.coord - atom2.coord
    v3 = atom4.coord - atom3.coord
    
    n1 = np.cross(v1, v2)
    n2 = np.cross(v2, v3)
    
    n1_norm = n1 / np.linalg.norm(n1)
    n2_norm = n2 / np.linalg.norm(n2)
    
    cosine = np.dot(n1_norm, n2_norm)
    cosine = np.clip(cosine, -1, 1)
    
    return np.degrees(np.arccos(cosine))


def get_residue_centroid(residue):
    """Calculate geometric center of residue"""
    coords = [atom.coord for atom in residue.get_atoms()]
    return np.mean(coords, axis=0)


def get_residue_id(residue):
    """Get unique identifier for residue"""
    return (residue.get_parent().id, residue.id[1], residue.get_resname())


# ============================================================================
# MAIN INTERACTION ANALYZER CLASS
# ============================================================================

class InteractionAnalyzer:
    """
    Main class for analyzing protein interactions
    """
    
    def __init__(self, pdb_file: str):
        """
        Initialize analyzer with PDB file
        
        Parameters:
        -----------
        pdb_file : str
            Path to PDB file
        """
        self.pdb_file = pdb_file
        self.structure = None
        self.interactions = []
        self.hubs = []
        self.critical_interactions = []
        
        # Distance cutoffs (in Angstroms)
        self.distance_cutoffs = {
            'hydrogen_bond': 3.5,
            'salt_bridge': 4.0,
            'hydrophobic': 5.0,
            'van_der_waals': 4.5,
            'disulfide_bond': 2.5,
            'pi_pi_stacking': 6.0,
            'cation_pi': 6.0,
            'halogen_bond': 4.0
        }
        
        # Angle cutoffs (in degrees)
        self.angle_cutoffs = {
            'hydrogen_bond_min': 120,
            'hydrogen_bond_max': 180,
            'halogen_bond_min': 140
        }
        
        self.load_structure()
    
    
    def load_structure(self):
        """Load PDB structure"""
        parser = PDBParser(QUIET=True)
        self.structure = parser.get_structure('protein', self.pdb_file)
        print(f"Loaded structure with {len(list(self.structure.get_residues()))} residues")
    
    
    def analyze_all(self):
        """Run all analyses"""
        print("\n=== Starting Complete Interaction Analysis ===\n")
        
        print("1. Detecting intra-protein interactions...")
        self.detect_intra_protein_interactions()
        print(f"   Found {len(self.interactions)} interactions")
        
        print("\n2. Identifying hub residues...")
        self.identify_hubs()
        print(f"   Found {len(self.hubs)} hub residues")
        
        print("\n3. Identifying critical interactions...")
        self.identify_critical_interactions()
        print(f"   Found {len([i for i in self.critical_interactions if i['level'] == 'critical'])} critical interactions")
        
        print("\n=== Analysis Complete ===\n")
    
    
    def detect_intra_protein_interactions(self):
        """
        Detect all interactions within protein structure
        Implements INTRA_PROTEIN_INTERACTION_DETECTOR algorithm
        """
        self.interactions = []
        
        # Get all residues
        residues = list(self.structure.get_residues())
        
        # Build spatial index for fast neighbor search
        atoms = list(self.structure.get_atoms())
        ns = NeighborSearch(atoms)
        
        # Analyze each residue pair
        for i, residue_i in enumerate(residues):
            
            # Find nearby residues using neighbor search
            for residue_j in residues[i+1:]:
                
                # Skip nearby residues in sequence
                if residue_i.get_parent() == residue_j.get_parent():  # Same chain
                    seq_separation = abs(residue_i.id[1] - residue_j.id[1])
                    if seq_separation < 4:
                        continue
                
                # Check if residues are close enough
                min_dist = self._get_min_distance(residue_i, residue_j)
                if min_dist > max(self.distance_cutoffs.values()):
                    continue
                
                # Detect specific interaction types
                self._detect_hydrogen_bonds(residue_i, residue_j)
                self._detect_salt_bridges(residue_i, residue_j)
                self._detect_disulfide_bonds(residue_i, residue_j)
                self._detect_hydrophobic_interactions(residue_i, residue_j)
                self._detect_pi_pi_stacking(residue_i, residue_j)
                self._detect_cation_pi(residue_i, residue_j)
    
    
    def _get_min_distance(self, residue_i, residue_j):
        """Calculate minimum distance between two residues"""
        min_dist = float('inf')
        for atom_i in residue_i.get_atoms():
            for atom_j in residue_j.get_atoms():
                dist = calculate_distance(atom_i, atom_j)
                if dist < min_dist:
                    min_dist = dist
        return min_dist
    
    
    def _detect_hydrogen_bonds(self, residue_i, residue_j):
        """Detect hydrogen bonds between two residues"""
        
        # Get potential donors and acceptors
        donors_i = self._get_hbond_donors(residue_i)
        acceptors_i = self._get_hbond_acceptors(residue_i)
        donors_j = self._get_hbond_donors(residue_j)
        acceptors_j = self._get_hbond_acceptors(residue_j)
        
        # Check all donor-acceptor pairs
        for donor_residue, acceptor_residue, donors, acceptors in [
            (residue_i, residue_j, donors_i, acceptors_j),
            (residue_j, residue_i, donors_j, acceptors_i)
        ]:
            for donor_atom in donors:
                for acceptor_atom in acceptors:
                    
                    # Check distance
                    dist = calculate_distance(donor_atom, acceptor_atom)
                    if dist <= self.distance_cutoffs['hydrogen_bond']:
                        
                        # For proper H-bond, check angle if hydrogen atom available
                        # Simplified: just use donor heavy atom
                        angle = self._estimate_hbond_angle(donor_atom, acceptor_atom, donor_residue)
                        
                        if angle >= self.angle_cutoffs['hydrogen_bond_min']:
                            strength = self._calculate_hbond_strength(dist, angle)
                            
                            interaction = Interaction(
                                type='hydrogen_bond',
                                residue1_id=get_residue_id(donor_residue),
                                residue2_id=get_residue_id(acceptor_residue),
                                atoms=[donor_atom.name, acceptor_atom.name],
                                distance=dist,
                                angle=angle,
                                strength=strength
                            )
                            self.interactions.append(interaction)
    
    
    def _get_hbond_donors(self, residue):
        """Get hydrogen bond donor atoms from residue"""
        donors = []
        resname = residue.get_resname()
        
        # Get sidechain donors
        if resname in HBOND_DONORS:
            for atom_name in HBOND_DONORS[resname]:
                if atom_name in residue:
                    donors.append(residue[atom_name])
        
        # Add backbone N
        if 'N' in residue:
            donors.append(residue['N'])
        
        return donors
    
    
    def _get_hbond_acceptors(self, residue):
        """Get hydrogen bond acceptor atoms from residue"""
        acceptors = []
        resname = residue.get_resname()
        
        # Get sidechain acceptors
        if resname in HBOND_ACCEPTORS:
            for atom_name in HBOND_ACCEPTORS[resname]:
                if atom_name in residue:
                    acceptors.append(residue[atom_name])
        
        # Add backbone O
        if 'O' in residue:
            acceptors.append(residue['O'])
        
        return acceptors
    
    
    def _estimate_hbond_angle(self, donor, acceptor, donor_residue):
        """
        Estimate hydrogen bond angle
        This is simplified - ideally would use actual H atom position
        """
        # Use CA as reference point for angle calculation
        if 'CA' in donor_residue:
            ca = donor_residue['CA']
            angle = calculate_angle(ca, donor, acceptor)
            return angle
        return 180  # Default to linear
    
    
    def _calculate_hbond_strength(self, distance, angle):
        """
        Calculate hydrogen bond strength based on distance and angle
        Returns value between 0 and 1
        """
        # Distance component (optimal at ~2.8 Å)
        optimal_distance = 2.8
        distance_score = np.exp(-((distance - optimal_distance) ** 2) / 0.5)
        
        # Angle component (optimal at 180°)
        angle_score = (angle - 120) / 60  # Normalize 120-180 to 0-1
        
        strength = 0.6 * distance_score + 0.4 * angle_score
        return max(0, min(1, strength))
    
    
    def _detect_salt_bridges(self, residue_i, residue_j):
        """Detect salt bridges (ionic interactions)"""
        
        # Check if residues are charged with opposite charges
        if not (is_charged_residue(residue_i) and is_charged_residue(residue_j)):
            return
        
        charge_i = get_charge(residue_i)
        charge_j = get_charge(residue_j)
        
        if charge_i == charge_j:  # Same charge, no salt bridge
            return
        
        # Get charged atoms
        charged_atoms_i = self._get_charged_atoms(residue_i)
        charged_atoms_j = self._get_charged_atoms(residue_j)
        
        # Find closest pair
        min_dist = float('inf')
        best_pair = None
        
        for atom_i in charged_atoms_i:
            for atom_j in charged_atoms_j:
                dist = calculate_distance(atom_i, atom_j)
                if dist < min_dist:
                    min_dist = dist
                    best_pair = (atom_i, atom_j)
        
        # If within cutoff, add salt bridge
        if min_dist <= self.distance_cutoffs['salt_bridge']:
            strength = self._calculate_salt_bridge_strength(min_dist)
            
            interaction = Interaction(
                type='salt_bridge',
                residue1_id=get_residue_id(residue_i),
                residue2_id=get_residue_id(residue_j),
                atoms=[best_pair[0].name, best_pair[1].name],
                distance=min_dist,
                strength=strength
            )
            self.interactions.append(interaction)
    
    
    def _get_charged_atoms(self, residue):
        """Get charged atoms from residue"""
        resname = residue.get_resname()
        charged_atoms = []
        
        charge_atom_map = {
            'ARG': ['NH1', 'NH2', 'NE'],
            'LYS': ['NZ'],
            'ASP': ['OD1', 'OD2'],
            'GLU': ['OE1', 'OE2']
        }
        
        if resname in charge_atom_map:
            for atom_name in charge_atom_map[resname]:
                if atom_name in residue:
                    charged_atoms.append(residue[atom_name])
        
        return charged_atoms
    
    
    def _calculate_salt_bridge_strength(self, distance):
        """Calculate salt bridge strength based on distance"""
        # Coulombic interaction: strength ~ 1/r
        optimal_distance = 3.0
        strength = optimal_distance / distance
        return min(1, strength)
    
    
    def _detect_disulfide_bonds(self, residue_i, residue_j):
        """Detect disulfide bonds between cysteines"""
        
        # Check if both are cysteines
        if residue_i.get_resname() != 'CYS' or residue_j.get_resname() != 'CYS':
            return
        
        # Check SG-SG distance
        if 'SG' not in residue_i or 'SG' not in residue_j:
            return
        
        sg_i = residue_i['SG']
        sg_j = residue_j['SG']
        dist = calculate_distance(sg_i, sg_j)
        
        if dist <= self.distance_cutoffs['disulfide_bond']:
            interaction = Interaction(
                type='disulfide_bond',
                residue1_id=get_residue_id(residue_i),
                residue2_id=get_residue_id(residue_j),
                atoms=['SG', 'SG'],
                distance=dist,
                strength=1.0  # Covalent bond - always strong
            )
            self.interactions.append(interaction)
    
    
    def _detect_hydrophobic_interactions(self, residue_i, residue_j):
        """Detect hydrophobic interactions"""
        
        # Check if both residues are hydrophobic
        if not (is_hydrophobic_residue(residue_i) and is_hydrophobic_residue(residue_j)):
            return
        
        # Calculate centroid distance
        centroid_i = get_residue_centroid(residue_i)
        centroid_j = get_residue_centroid(residue_j)
        dist = np.linalg.norm(centroid_i - centroid_j)
        
        if dist <= self.distance_cutoffs['hydrophobic']:
            # Estimate contact area (simplified)
            contact_area = self._estimate_contact_area(residue_i, residue_j, dist)
            strength = min(1.0, contact_area / 20.0)  # Normalize by typical area
            
            interaction = Interaction(
                type='hydrophobic',
                residue1_id=get_residue_id(residue_i),
                residue2_id=get_residue_id(residue_j),
                atoms=['centroid', 'centroid'],
                distance=dist,
                strength=strength
            )
            self.interactions.append(interaction)
    
    
    def _estimate_contact_area(self, residue_i, residue_j, distance):
        """Estimate contact area between residues (simplified)"""
        # This is a simplified estimate
        # Real implementation would use SASA calculations
        max_area = 30.0  # Typical maximum contact area
        contact = max_area * np.exp(-(distance - 3.5) ** 2)
        return contact
    
    
    def _detect_pi_pi_stacking(self, residue_i, residue_j):
        """Detect pi-pi stacking interactions"""
        
        # Check if both have aromatic rings
        if not (is_aromatic_residue(residue_i) and is_aromatic_residue(residue_j)):
            return
        
        # Get ring centers (simplified - using specific atoms)
        ring_center_i = self._get_aromatic_center(residue_i)
        ring_center_j = self._get_aromatic_center(residue_j)
        
        if ring_center_i is None or ring_center_j is None:
            return
        
        dist = np.linalg.norm(ring_center_i - ring_center_j)
        
        if dist <= self.distance_cutoffs['pi_pi_stacking']:
            # Simplified angle calculation
            # Real implementation would calculate ring plane angles
            angle = 30  # Placeholder
            
            geometry = 'parallel' if angle < 30 else 'T-shaped'
            strength = 1.0 - (dist / self.distance_cutoffs['pi_pi_stacking'])
            
            interaction = Interaction(
                type='pi_pi_stacking',
                residue1_id=get_residue_id(residue_i),
                residue2_id=get_residue_id(residue_j),
                atoms=['ring', 'ring'],
                distance=dist,
                angle=angle,
                strength=strength
            )
            self.interactions.append(interaction)
    
    
    def _get_aromatic_center(self, residue):
        """Get center of aromatic ring"""
        resname = residue.get_resname()
        
        # Define aromatic atoms for each residue type
        aromatic_atoms = {
            'PHE': ['CG', 'CD1', 'CD2', 'CE1', 'CE2', 'CZ'],
            'TYR': ['CG', 'CD1', 'CD2', 'CE1', 'CE2', 'CZ'],
            'TRP': ['CD2', 'CE2', 'CE3', 'CZ2', 'CZ3', 'CH2'],
            'HIS': ['CG', 'ND1', 'CD2', 'CE1', 'NE2']
        }
        
        if resname not in aromatic_atoms:
            return None
        
        coords = []
        for atom_name in aromatic_atoms[resname]:
            if atom_name in residue:
                coords.append(residue[atom_name].coord)
        
        if len(coords) < 3:
            return None
        
        return np.mean(coords, axis=0)
    
    
    def _detect_cation_pi(self, residue_i, residue_j):
        """Detect cation-pi interactions"""
        
        # Check if one is cation and other is aromatic
        cation = None
        aromatic = None
        
        if is_charged_residue(residue_i) and get_charge(residue_i) == '+':
            if is_aromatic_residue(residue_j):
                cation = residue_i
                aromatic = residue_j
        elif is_charged_residue(residue_j) and get_charge(residue_j) == '+':
            if is_aromatic_residue(residue_i):
                cation = residue_j
                aromatic = residue_i
        
        if cation is None or aromatic is None:
            return
        
        # Get cation center and aromatic center
        cation_center = get_residue_centroid(cation)
        aromatic_center = self._get_aromatic_center(aromatic)
        
        if aromatic_center is None:
            return
        
        dist = np.linalg.norm(cation_center - aromatic_center)
        
        if dist <= self.distance_cutoffs['cation_pi']:
            strength = 1.0 - (dist / self.distance_cutoffs['cation_pi'])
            
            interaction = Interaction(
                type='cation_pi',
                residue1_id=get_residue_id(cation),
                residue2_id=get_residue_id(aromatic),
                atoms=['cation', 'aromatic'],
                distance=dist,
                strength=strength
            )
            self.interactions.append(interaction)
    
    
    def identify_hubs(self, hub_threshold=5):
        """
        Identify hub residues with many interactions
        Implements HUB_RESIDUE_DETECTOR algorithm
        """
        self.hubs = []
        
        # Count interactions per residue
        residue_interactions = {}
        residue_types = {}
        
        for interaction in self.interactions:
            for res_id in [interaction.residue1_id, interaction.residue2_id]:
                if res_id not in residue_interactions:
                    residue_interactions[res_id] = []
                    residue_types[res_id] = {}
                
                residue_interactions[res_id].append(interaction)
                
                int_type = interaction.type
                residue_types[res_id][int_type] = residue_types[res_id].get(int_type, 0) + 1
        
        # Identify hubs
        for res_id, interactions in residue_interactions.items():
            degree = len(interactions)
            
            if degree >= hub_threshold:
                diversity = len(residue_types[res_id])
                avg_strength = np.mean([i.strength for i in interactions if i.strength])
                
                # Classify hub type
                if diversity >= 3:
                    hub_type = 'multi-functional'
                else:
                    hub_type = 'structural'
                
                # Calculate hub score
                hub_score = self._calculate_hub_score(degree, diversity, avg_strength)
                
                hub = HubResidue(
                    residue_id=res_id,
                    degree=degree,
                    interaction_types=residue_types[res_id],
                    diversity=diversity,
                    hub_score=hub_score,
                    hub_type=hub_type
                )
                
                self.hubs.append(hub)
        
        # Sort by hub score
        self.hubs.sort(key=lambda x: x.hub_score, reverse=True)
    
    
    def _calculate_hub_score(self, degree, diversity, avg_strength):
        """Calculate hub score from multiple metrics"""
        # Normalize values
        degree_norm = min(1.0, degree / 20.0)
        diversity_norm = min(1.0, diversity / 6.0)
        
        score = 0.4 * degree_norm + 0.3 * diversity_norm + 0.3 * avg_strength
        return score
    
    
    def identify_critical_interactions(self):
        """
        Identify critical interactions for stability/function
        Implements CRITICAL_INTERACTION_DETECTOR algorithm
        """
        self.critical_interactions = []
        
        for interaction in self.interactions:
            score = 0
            factors = []
            
            # A. Interaction strength
            if interaction.type == 'disulfide_bond':
                score += 10
                factors.append('covalent_bond')
            elif interaction.type == 'salt_bridge':
                score += 7
                factors.append('strong_electrostatic')
            elif interaction.type == 'hydrogen_bond' and interaction.strength > 0.7:
                score += 5
                factors.append('strong_hydrogen_bond')
            
            # B. Structural location (simplified - would need secondary structure info)
            # For now, use distance from surface as proxy
            
            # C. Unique interactions
            type_count = sum(1 for i in self.interactions if i.type == interaction.type)
            if type_count < 5:  # Rare interaction type
                score += 3
                factors.append('unique_interaction')
            
            # Classify criticality
            if score >= 15:
                level = 'critical'
            elif score >= 10:
                level = 'important'
            elif score >= 5:
                level = 'moderate'
            else:
                level = 'minor'
            
            self.critical_interactions.append({
                'interaction': interaction,
                'criticality_score': score,
                'level': level,
                'factors': factors
            })
        
        # Sort by score
        self.critical_interactions.sort(key=lambda x: x['criticality_score'], reverse=True)
    
    
    def export_to_csv(self, output_dir='.'):
        """Export all results to CSV files"""
        import os
        
        # Export interactions
        interactions_df = pd.DataFrame([i.to_dict() for i in self.interactions])
        interactions_df.to_csv(os.path.join(output_dir, 'interactions.csv'), index=False)
        print(f"Exported {len(self.interactions)} interactions to interactions.csv")
        
        # Export hubs
        if self.hubs:
            hubs_df = pd.DataFrame([h.to_dict() for h in self.hubs])
            hubs_df.to_csv(os.path.join(output_dir, 'hubs.csv'), index=False)
            print(f"Exported {len(self.hubs)} hub residues to hubs.csv")
        
        # Export critical interactions
        if self.critical_interactions:
            critical_df = pd.DataFrame([
                {
                    **c['interaction'].to_dict(),
                    'criticality_score': c['criticality_score'],
                    'criticality_level': c['level'],
                    'factors': ','.join(c['factors'])
                }
                for c in self.critical_interactions
            ])
            critical_df.to_csv(os.path.join(output_dir, 'critical_interactions.csv'), index=False)
            print(f"Exported {len(self.critical_interactions)} critical interactions to critical_interactions.csv")
    
    
    def get_statistics(self):
        """Get summary statistics"""
        stats = {
            'total_interactions': len(self.interactions),
            'total_hubs': len(self.hubs),
            'interaction_types': {},
            'critical_count': len([i for i in self.critical_interactions if i['level'] == 'critical']),
            'important_count': len([i for i in self.critical_interactions if i['level'] == 'important'])
        }
        
        # Count by type
        for interaction in self.interactions:
            itype = interaction.type
            stats['interaction_types'][itype] = stats['interaction_types'].get(itype, 0) + 1
        
        return stats
    
    
    def print_summary(self):
        """Print analysis summary"""
        stats = self.get_statistics()
        
        print("\n" + "="*60)
        print("INTERACTION ANALYSIS SUMMARY")
        print("="*60)
        print(f"Total Interactions: {stats['total_interactions']}")
        print(f"Hub Residues: {stats['total_hubs']}")
        print(f"Critical Interactions: {stats['critical_count']}")
        print(f"Important Interactions: {stats['important_count']}")
        print("\nInteraction Types:")
        for itype, count in sorted(stats['interaction_types'].items()):
            print(f"  {itype}: {count}")
        print("="*60 + "\n")


# ============================================================================
# EXAMPLE USAGE
# ============================================================================

if __name__ == "__main__":
    
    # Example: Analyze a protein structure
    pdb_file = "example.pdb"  # Replace with your PDB file
    
    # Create analyzer
    analyzer = InteractionAnalyzer(pdb_file)
    
    # Run complete analysis
    analyzer.analyze_all()
    
    # Print summary
    analyzer.print_summary()
    
    # Export results
    analyzer.export_to_csv(output_dir='./results')
    
    # Access specific results
    print("\nTop 5 Hub Residues:")
    for i, hub in enumerate(analyzer.hubs[:5], 1):
        print(f"{i}. {hub.residue_id[2]}{hub.residue_id[1]} - Degree: {hub.degree}, Score: {hub.hub_score:.3f}")
    
    print("\nTop 5 Critical Interactions:")
    for i, critical in enumerate(analyzer.critical_interactions[:5], 1):
        interaction = critical['interaction']
        print(f"{i}. {interaction.type}: {interaction.residue1_id[2]}{interaction.residue1_id[1]} - "
              f"{interaction.residue2_id[2]}{interaction.residue2_id[1]} "
              f"(Score: {critical['criticality_score']})")
