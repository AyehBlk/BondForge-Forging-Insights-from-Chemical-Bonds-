"""
BondForge - Advanced Energy Calculator
======================================

Quantum mechanical energy corrections and advanced force field calculations.

Author: Ayeh Bolouki
Email: ayehbolouki1988@gmail.com | ayehgeek@gmail.com
"""

import numpy as np
from dataclasses import dataclass
from typing import Dict, Optional
import logging

logger = logging.getLogger(__name__)

# Physical constants
COULOMB_CONSTANT = 332.0636  # kcal·Å/(mol·e²)
KB = 0.001987204  # kcal/(mol·K)
TEMPERATURE = 298.15  # K
EPSILON_VACUUM = 1.0
EPSILON_WATER = 78.5
EPSILON_PROTEIN = 4.0
IONIC_STRENGTH_DEFAULT = 0.15  # M


@dataclass
class EnergyComponents:
    """Container for energy decomposition"""
    electrostatic: float
    polarization: float
    charge_transfer: float
    dispersion: float
    exchange_repulsion: float
    total: float
    
    def __str__(self):
        return (f"Total: {self.total:.3f} kcal/mol\n"
                f"  Electrostatic: {self.electrostatic:.3f}\n"
                f"  Polarization: {self.polarization:.3f}\n"
                f"  Charge Transfer: {self.charge_transfer:.3f}\n"
                f"  Dispersion: {self.dispersion:.3f}\n"
                f"  Exchange-Repulsion: {self.exchange_repulsion:.3f}")


class AdvancedEnergyCalculator:
    """
    Advanced energy calculations with QM corrections
    
    Based on:
    - Morokuma energy decomposition
    - Symmetry-Adapted Perturbation Theory (SAPT)
    - DFT-SAPT methods
    
    References:
    - Morozov & Kortemme (2006) PNAS
    - Hohenstein & Sherrill (2012) WIREs Comput Mol Sci
    - Stone (2013) "The Theory of Intermolecular Forces"
    """
    
    def __init__(self, force_field: str = 'amber'):
        """Initialize calculator with force field parameters"""
        self.force_field = force_field
        self._load_parameters()
    
    def _load_parameters(self):
        """Load force field parameters"""
        self.vdw_params = {
            'C': {'epsilon': 0.070, 'sigma': 1.908},
            'N': {'epsilon': 0.200, 'sigma': 1.824},
            'O': {'epsilon': 0.200, 'sigma': 1.661},
            'S': {'epsilon': 0.250, 'sigma': 2.000},
            'H': {'epsilon': 0.015, 'sigma': 1.459}
        }
        
        self.polarizabilities = {
            'C': 1.76, 'N': 1.10, 'O': 0.80,
            'S': 2.90, 'H': 0.67
        }
        
        self.c6_coefficients = {
            ('C', 'C'): 46.6, ('C', 'N'): 35.4,
            ('C', 'O'): 26.1, ('N', 'N'): 26.8,
            ('N', 'O'): 19.8, ('O', 'O'): 14.6
        }
    
    def calculate_hydrogen_bond_energy(
        self,
        donor_atom: str,
        acceptor_atom: str,
        distance: float,
        angle: float,
        donor_charge: float = -0.4,
        acceptor_charge: float = -0.5
    ) -> EnergyComponents:
        """
        Calculate H-bond energy with QM corrections
        
        E = E_elec + E_pol + E_ct + E_disp + E_rep
        """
        # 1. Electrostatic energy
        epsilon_eff = self._distance_dependent_dielectric(distance)
        E_elec = COULOMB_CONSTANT * donor_charge * acceptor_charge / (epsilon_eff * distance)
        
        # Angle correction
        angle_factor = np.cos(np.radians(angle - 180)) ** 2
        E_elec *= angle_factor
        
        # 2. Polarization energy
        alpha_acceptor = self.polarizabilities.get(acceptor_atom, 1.0)
        electric_field = abs(donor_charge) / (distance ** 2)
        E_pol = -0.5 * alpha_acceptor * (electric_field ** 2)
        E_pol *= angle_factor
        
        # 3. Charge transfer energy
        overlap = np.exp(-0.5 * (distance - 2.8))
        E_ct = -3.0 * overlap * angle_factor
        
        # 4. Dispersion energy
        pair = tuple(sorted([donor_atom, acceptor_atom]))
        C6 = self.c6_coefficients.get(pair, 20.0)
        C8 = C6 * 50
        
        beta = 1.5
        E_disp = -(C6 / distance**6) * self._damping(6, beta, distance)
        E_disp += -(C8 / distance**8) * self._damping(8, beta, distance)
        E_disp *= 0.0006275  # Convert to kcal/mol
        
        # 5. Exchange-repulsion
        A, B = 1000.0, 3.5
        E_rep = A * np.exp(-B * distance)
        
        E_total = E_elec + E_pol + E_ct + E_disp + E_rep
        
        return EnergyComponents(
            electrostatic=E_elec,
            polarization=E_pol,
            charge_transfer=E_ct,
            dispersion=E_disp,
            exchange_repulsion=E_rep,
            total=E_total
        )
    
    def calculate_salt_bridge_energy(
        self,
        distance: float,
        geometry: str = 'monodentate',
        ionic_strength: float = IONIC_STRENGTH_DEFAULT,
        buried_fraction: float = 0.5
    ) -> float:
        """Calculate salt bridge energy with screening"""
        q1, q2 = +1.0, -1.0
        
        # Distance-dependent dielectric
        lambda_param, d_offset = 0.003627, 2.0
        epsilon_eff = EPSILON_PROTEIN + (EPSILON_WATER - EPSILON_PROTEIN) * (
            1 - np.exp(-lambda_param * (distance - d_offset))
        )
        epsilon_eff = EPSILON_PROTEIN * buried_fraction + epsilon_eff * (1 - buried_fraction)
        
        # Coulomb energy
        E_coulomb = COULOMB_CONSTANT * q1 * q2 / (epsilon_eff * distance)
        
        # Debye-Hückel screening
        kappa = 0.316 * np.sqrt(ionic_strength)
        screening = np.exp(-kappa * distance)
        
        E_screened = E_coulomb * screening
        
        # Geometry correction
        geometry_factor = 1.5 if geometry == 'bidentate' else 1.0
        
        return E_screened * geometry_factor
    
    def calculate_vdw_energy_quantum(
        self,
        atom1_type: str,
        atom2_type: str,
        distance: float,
        is_aromatic: bool = False
    ) -> float:
        """Lennard-Jones potential with quantum corrections"""
        # Get LJ parameters
        epsilon1 = self.vdw_params[atom1_type]['epsilon']
        epsilon2 = self.vdw_params[atom2_type]['epsilon']
        sigma1 = self.vdw_params[atom1_type]['sigma']
        sigma2 = self.vdw_params[atom2_type]['sigma']
        
        epsilon = np.sqrt(epsilon1 * epsilon2)
        sigma = (sigma1 + sigma2) / 2.0
        
        # Standard LJ
        r_ratio = sigma / distance
        E_lj = 4 * epsilon * (r_ratio**12 - r_ratio**6)
        
        # D3 dispersion correction
        pair = tuple(sorted([atom1_type, atom2_type]))
        C6 = self.c6_coefficients.get(pair, 20.0) * 0.0006275
        C8 = C6 * 50
        
        s6, s8, a1, a2 = 1.0, 1.0, 0.40, 5.0
        R0 = sigma
        
        def fermi_damping(r, R0, n):
            return 1.0 / (1.0 + 6.0 * (r / (s6 * R0))**(-a2 * n))
        
        damp6 = fermi_damping(distance, R0, 1)
        damp8 = fermi_damping(distance, R0, 2)
        
        E_d3 = -(s6 * C6 / distance**6) * damp6 - (s8 * C8 / distance**8) * damp8
        
        if is_aromatic:
            E_d3 *= 1.2
        
        return E_lj + E_d3
    
    def calculate_pi_pi_energy(
        self,
        distance: float,
        angle: float,
        geometry: str = 'parallel'
    ) -> float:
        """Pi-pi stacking energy with orientation dependence"""
        if geometry == 'parallel':
            r_opt, E_opt = 3.8, -2.4
            angle_penalty = (angle / 30.0) ** 2
        elif geometry == 't-shaped':
            r_opt, E_opt = 4.5, -2.6
            angle_penalty = ((angle - 90) / 30.0) ** 2
        else:
            r_opt, E_opt = 5.0, -1.5
            angle_penalty = 0
        
        alpha = 2.0
        distance_term = np.exp(-2 * alpha * (distance - r_opt)) - \
                       2 * np.exp(-alpha * (distance - r_opt))
        
        return E_opt * distance_term * np.exp(-angle_penalty)
    
    @staticmethod
    def _distance_dependent_dielectric(distance: float) -> float:
        """Mehler-Solmajer distance-dependent dielectric"""
        lambda_param, d_offset = 0.003627, 2.0
        return EPSILON_PROTEIN + (EPSILON_WATER - EPSILON_PROTEIN) * (
            1 - np.exp(-lambda_param * (distance - d_offset))
        )
    
    @staticmethod
    def _damping(n, beta, r):
        """Tang-Toennies damping function"""
        b_r = beta * r
        damp = 1.0 - np.exp(-b_r) * sum(
            b_r**k / np.math.factorial(k) for k in range(n + 1)
        )
        return damp
