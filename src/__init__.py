"""
BondForge - Forging Insights from Chemical Bonds
=================================================

A comprehensive, production-ready toolkit for analyzing 20 types of chemical 
interactions in protein structures.

Author: Ayeh Bolouki
Email: ayehbolouki1988@gmail.com | ayehgeek@gmail.com
GitHub: https://github.com/AyehBlk/BondForge

Modules:
    interaction_analyzer: Core analyzer with 7 interaction types
    extended_interaction_analyzer: Extended analyzer with 20 interaction types
    interaction_visualizer: Visualization and plotting tools
    cli: Command-line interface
    config: Configuration management
    energy_calculator: Advanced energy calculations
    statistical_validator: Statistical validation tools
"""

__version__ = '2.0.0'
__author__ = 'Ayeh Bolouki'
__email__ = 'ayehbolouki1988@gmail.com'
__package_name__ = 'BondForge'
__tagline__ = 'Forging Insights from Chemical Bonds'

# Try to import main classes
try:
    from .interaction_analyzer import ProteinInteractionAnalyzer
    from .extended_interaction_analyzer import ExtendedProteinInteractionAnalyzer
    from .interaction_visualizer import InteractionVisualizer
    
    __all__ = [
        'ProteinInteractionAnalyzer',
        'ExtendedProteinInteractionAnalyzer',
        'InteractionVisualizer',
    ]
except ImportError:
    # If imports fail, just expose version info
    __all__ = []

def get_version():
    """Get BondForge version"""
    return __version__

def get_info():
    """Get package information"""
    return {
        'name': __package_name__,
        'version': __version__,
        'author': __author__,
        'email': __email__,
        'tagline': __tagline__,
        'url': 'https://github.com/AyehBlk/BondForge',
    }
