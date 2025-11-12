"""
BondForge - Statistical Validator
=================================

Statistical validation of interactions against PDB database.

Author: Ayeh Bolouki
Email: ayehbolouki1988@gmail.com | ayehgeek@gmail.com
"""

import numpy as np
from scipy import stats
from typing import Dict, Optional, Tuple
import logging

logger = logging.getLogger(__name__)


class StatisticalValidator:
    """
    Statistical validation of interactions against PDB database
    
    Methods:
    - Z-score calculation
    - Percentile ranking
    - Confidence intervals
    - Outlier detection
    """
    
    def __init__(self, pdb_statistics: Optional[Dict] = None):
        """
        Initialize validator
        
        Parameters:
        -----------
        pdb_statistics : dict
            Pre-computed statistics from PDB analysis
            If None, uses built-in statistics
        """
        if pdb_statistics is None:
            self.statistics = self._load_default_statistics()
        else:
            self.statistics = pdb_statistics
    
    def _load_default_statistics(self) -> Dict:
        """
        Load default statistics from PDB analysis
        
        Based on:
        - McDonald & Thornton (1994) J Mol Biol
        - Torshin et al. (2002) Protein Engineering
        """
        return {
            'hydrogen_bond': {
                'mean_distance': 2.9,
                'std_distance': 0.3,
                'mean_angle': 155,
                'std_angle': 15,
                'frequency': 0.85
            },
            'salt_bridge': {
                'mean_distance': 3.5,
                'std_distance': 0.5,
                'frequency': 0.15
            },
            'hydrophobic': {
                'mean_distance': 4.5,
                'std_distance': 0.6,
                'frequency': 2.5
            },
            'pi_pi_stacking': {
                'mean_distance': 4.5,
                'std_distance': 0.5,
                'mean_angle_parallel': 10,
                'std_angle_parallel': 8,
                'mean_angle_tshaped': 85,
                'std_angle_tshaped': 10,
                'frequency': 0.05
            },
            'cation_pi': {
                'mean_distance': 5.0,
                'std_distance': 0.7,
                'frequency': 0.03
            },
            'disulfide_bond': {
                'mean_distance': 2.05,
                'std_distance': 0.05,
                'frequency': 0.02
            }
        }
    
    def calculate_zscore(
        self,
        interaction_type: str,
        observed_distance: float,
        observed_angle: Optional[float] = None
    ) -> Dict:
        """
        Calculate Z-score for interaction
        
        Z = (X - μ) / σ
        
        Parameters:
        -----------
        interaction_type : str
            Type of interaction
        observed_distance : float
            Observed distance (Å)
        observed_angle : float, optional
            Observed angle (degrees)
        
        Returns:
        --------
        dict
            Z-scores, percentiles, and classification
        """
        if interaction_type not in self.statistics:
            return {
                'z_score_distance': None,
                'z_score_angle': None,
                'is_typical': True,
                'confidence': 'unknown'
            }
        
        stats_data = self.statistics[interaction_type]
        
        # Distance Z-score
        mu_dist = stats_data['mean_distance']
        sigma_dist = stats_data['std_distance']
        z_distance = (observed_distance - mu_dist) / sigma_dist
        percentile_dist = stats.norm.cdf(z_distance)
        
        result = {
            'z_score_distance': z_distance,
            'percentile_distance': percentile_dist,
            'is_typical_distance': abs(z_distance) < 2.0
        }
        
        # Angle Z-score (if applicable)
        if observed_angle is not None and 'mean_angle' in stats_data:
            mu_angle = stats_data['mean_angle']
            sigma_angle = stats_data['std_angle']
            z_angle = (observed_angle - mu_angle) / sigma_angle
            percentile_angle = stats.norm.cdf(z_angle)
            
            result['z_score_angle'] = z_angle
            result['percentile_angle'] = percentile_angle
            result['is_typical_angle'] = abs(z_angle) < 2.0
            result['is_typical'] = result['is_typical_distance'] and result['is_typical_angle']
        else:
            result['is_typical'] = result['is_typical_distance']
        
        # Confidence classification
        max_z = max(abs(z_distance), abs(result.get('z_score_angle', 0)))
        if max_z < 1.0:
            confidence = 'high'
        elif max_z < 2.0:
            confidence = 'medium'
        elif max_z < 3.0:
            confidence = 'low'
        else:
            confidence = 'very_low'
        
        result['confidence'] = confidence
        
        return result
    
    def calculate_interaction_confidence(
        self,
        interaction_type: str,
        distance: float,
        angle: Optional[float] = None,
        energy: Optional[float] = None,
        structural_context: Optional[Dict] = None
    ) -> Tuple[float, Dict]:
        """
        Calculate overall confidence score for interaction
        
        Combines:
        - Geometric criteria (distance, angle)
        - Energetic favorability
        - Statistical typicality
        - Structural context
        
        Parameters:
        -----------
        interaction_type : str
            Type of interaction
        distance : float
            Distance (Å)
        angle : float, optional
            Angle (degrees)
        energy : float, optional
            Interaction energy (kcal/mol)
        structural_context : dict, optional
            Additional context
        
        Returns:
        --------
        tuple
            (overall_confidence, component_scores)
            confidence ∈ [0, 1]
        """
        scores = {}
        
        # 1. Geometric score (Z-score based)
        z_result = self.calculate_zscore(interaction_type, distance, angle)
        if z_result['z_score_distance'] is not None:
            z_dist = abs(z_result['z_score_distance'])
            scores['geometric'] = np.exp(-0.5 * z_dist)
        else:
            scores['geometric'] = 0.5
        
        # 2. Energetic score
        if energy is not None:
            if energy < -1.0:
                scores['energetic'] = 1.0
            elif energy < 0:
                scores['energetic'] = 0.5 + 0.5 * abs(energy)
            else:
                scores['energetic'] = 0.2
        else:
            scores['energetic'] = 0.5
        
        # 3. Statistical score
        if interaction_type in self.statistics:
            freq = self.statistics[interaction_type]['frequency']
            scores['statistical'] = min(1.0, freq / 1.0)
        else:
            scores['statistical'] = 0.5
        
        # 4. Structural context score
        if structural_context:
            context_score = 0.5
            if 'buried_fraction' in structural_context:
                context_score += 0.2 * structural_context['buried_fraction']
            if 'in_secondary_structure' in structural_context:
                if structural_context['in_secondary_structure']:
                    context_score += 0.2
            scores['structural'] = min(1.0, context_score)
        else:
            scores['structural'] = 0.5
        
        # Weighted average
        weights = [0.35, 0.30, 0.20, 0.15]
        confidence = np.average(list(scores.values()), weights=weights)
        
        return confidence, scores
    
    def is_outlier(
        self,
        interaction_type: str,
        distance: float,
        angle: Optional[float] = None,
        threshold: float = 3.0
    ) -> bool:
        """
        Check if interaction is a statistical outlier
        
        Parameters:
        -----------
        interaction_type : str
            Type of interaction
        distance : float
            Observed distance
        angle : float, optional
            Observed angle
        threshold : float
            Z-score threshold (default: 3.0)
        
        Returns:
        --------
        bool
            True if outlier
        """
        z_result = self.calculate_zscore(interaction_type, distance, angle)
        
        if z_result['z_score_distance'] is None:
            return False
        
        is_distance_outlier = abs(z_result['z_score_distance']) > threshold
        
        if angle is not None and z_result.get('z_score_angle') is not None:
            is_angle_outlier = abs(z_result['z_score_angle']) > threshold
            return is_distance_outlier or is_angle_outlier
        
        return is_distance_outlier
