"""
BondForge - Configuration Management
====================================

Hierarchical configuration system for BondForge.

Author: Ayeh Bolouki
Email: ayehbolouki1988@gmail.com | ayehgeek@gmail.com
"""

import os
import yaml
from pathlib import Path
from typing import Any, Dict, Optional
import logging

logger = logging.getLogger(__name__)


class Config:
    """
    Configuration manager with hierarchical loading
    
    Priority order (highest to lowest):
    1. CLI arguments
    2. Project config (./bondforge.yaml)
    3. User config (~/.bondforge/config.yaml)
    4. System config (/etc/bondforge/config.yaml)
    5. Default config
    """
    
    DEFAULT_CONFIG = {
        'analysis': {
            'interaction_types': [
                'hydrogen_bonds', 'salt_bridges', 'disulfide_bonds',
                'hydrophobic', 'pi_pi_stacking', 'cation_pi', 'halogen_bonds',
                'van_der_waals', 'anion_pi', 'sulfur_aromatic', 'ch_pi',
                'metal_coordination', 'carbonyl_pi', 'amide_aromatic',
                'sulfur_oxygen', 'backbone_carbonyl', 'aromatic_oxygen',
                'arginine_aromatic', 'backbone_amide_aromatic', 
                'backbone_carbonyl_charged'
            ],
            'distance_cutoffs': {
                'hydrogen_bond': 3.5,
                'salt_bridge': 4.0,
                'disulfide_bond': 2.5,
                'hydrophobic': 5.0,
                'pi_pi_stacking': 6.0,
                'cation_pi': 6.0,
                'halogen_bond': 4.0,
                'van_der_waals': 4.5
            },
            'angle_cutoffs': {
                'hydrogen_bond_min': 120,
                'hydrogen_bond_max': 180,
                'halogen_bond_min': 140
            },
            'energy_calculation': {
                'enabled': True,
                'method': 'advanced',  # 'simple' or 'advanced'
                'force_field': 'amber'
            }
        },
        'performance': {
            'n_jobs': -1,  # Use all cores
            'chunk_size': 1000,
            'cache_enabled': True,
            'cache_dir': '.bondforge_cache'
        },
        'output': {
            'format': ['csv'],
            'precision': 3,
            'include_metadata': True,
            'compress': False
        },
        'visualization': {
            'dpi': 300,
            'format': 'png',
            'style': 'publication',
            'color_scheme': 'viridis'
        }
    }
    
    def __init__(self, config_file: Optional[str] = None):
        """
        Initialize configuration
        
        Parameters:
        -----------
        config_file : str, optional
            Path to specific config file
        """
        self.config = self.DEFAULT_CONFIG.copy()
        self._load_config_files()
        
        if config_file:
            self._load_file(config_file)
    
    def _load_config_files(self):
        """Load configuration from files (system → user → project)"""
        config_locations = [
            Path('/etc/bondforge/config.yaml'),
            Path.home() / '.bondforge' / 'config.yaml',
            Path.cwd() / 'bondforge.yaml'
        ]
        
        for config_file in config_locations:
            if config_file.exists():
                self._load_file(str(config_file))
    
    def _load_file(self, config_file: str):
        """Load a specific config file"""
        try:
            with open(config_file) as f:
                user_config = yaml.safe_load(f)
                if user_config:
                    self._merge_config(user_config)
                    logger.debug(f"Loaded config from {config_file}")
        except Exception as e:
            logger.warning(f"Failed to load {config_file}: {e}")
    
    def _merge_config(self, user_config: Dict):
        """Deep merge user config into default config"""
        def deep_merge(base: Dict, update: Dict):
            for key, value in update.items():
                if key in base and isinstance(base[key], dict) and isinstance(value, dict):
                    deep_merge(base[key], value)
                else:
                    base[key] = value
        
        deep_merge(self.config, user_config)
    
    def get(self, key_path: str, default: Any = None) -> Any:
        """
        Get config value using dot notation
        
        Example: config.get('analysis.distance_cutoffs.hydrogen_bond')
        """
        keys = key_path.split('.')
        value = self.config
        for key in keys:
            if isinstance(value, dict):
                value = value.get(key)
            else:
                return default
        return value if value is not None else default
    
    def set(self, key_path: str, value: Any):
        """
        Set config value using dot notation
        
        Example: config.set('analysis.distance_cutoffs.hydrogen_bond', 3.2)
        """
        keys = key_path.split('.')
        config = self.config
        for key in keys[:-1]:
            if key not in config:
                config[key] = {}
            config = config[key]
        config[keys[-1]] = value
    
    def save(self, filepath: Optional[str] = None):
        """
        Save current configuration to file
        
        Parameters:
        -----------
        filepath : str, optional
            Path to save config. If None, saves to user config location
        """
        if filepath is None:
            filepath = Path.home() / '.bondforge' / 'config.yaml'
        else:
            filepath = Path(filepath)
        
        # Create directory if needed
        filepath.parent.mkdir(parents=True, exist_ok=True)
        
        with open(filepath, 'w') as f:
            yaml.dump(self.config, f, default_flow_style=False)
        
        logger.info(f"Configuration saved to {filepath}")
    
    def to_dict(self) -> Dict:
        """Return configuration as dictionary"""
        return self.config.copy()


# Global config instance
_config = None

def get_config() -> Config:
    """Get global configuration instance"""
    global _config
    if _config is None:
        _config = Config()
    return _config

def reset_config():
    """Reset configuration to defaults"""
    global _config
    _config = Config()
