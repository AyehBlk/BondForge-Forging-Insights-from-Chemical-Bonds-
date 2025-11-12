#!/usr/bin/env python3
"""
BondForge - Professional Protein Interaction Analysis Toolkit
==============================================================

Setup and installation script for BondForge 2.0

Author: Ayeh Bolouki
Email: ayehbolouki1988@gmail.com | ayehgeek@gmail.com
GitHub: https://github.com/AyehBlk/BondForge
"""

from setuptools import setup, find_packages
from pathlib import Path
import sys

# Ensure Python 3.7+
if sys.version_info < (3, 7):
    sys.exit('BondForge requires Python 3.7 or higher')

# Read the contents of README file
this_directory = Path(__file__).parent
readme_file = this_directory / 'README.md'
if readme_file.exists():
    long_description = readme_file.read_text(encoding='utf-8')
else:
    long_description = 'BondForge: Professional Protein Interaction Analysis Toolkit'

# Read requirements
requirements_file = this_directory / 'requirements.txt'
if requirements_file.exists():
    with open(requirements_file, encoding='utf-8') as f:
        requirements = [
            line.strip() 
            for line in f 
            if line.strip() and not line.startswith('#') and not line.startswith('--')
        ]
else:
    # Minimal requirements if file not found
    requirements = [
        'numpy>=1.21.0',
        'pandas>=1.3.0',
        'scipy>=1.7.0',
        'biopython>=1.79',
        'matplotlib>=3.4.0',
        'seaborn>=0.11.0',
        'networkx>=2.6.0',
        'click>=8.0.0',
        'rich>=10.0.0',
        'PyYAML>=5.4.0',
        'joblib>=1.1.0',
    ]

setup(
    # Package metadata
    name='bondforge',
    version='2.0.0',
    author='Ayeh Bolouki',
    author_email='ayehbolouki1988@gmail.com',
    maintainer='Ayeh Bolouki',
    maintainer_email='ayehbolouki1988@gmail.com',
    
    # Description
    description='🔥 BondForge: Professional toolkit for analyzing 20 types of chemical interactions in protein structures',
    long_description=long_description,
    long_description_content_type='text/markdown',
    
    # URLs
    url='https://github.com/AyehBlk/BondForge',
    project_urls={
        'Documentation': 'https://github.com/AyehBlk/BondForge/tree/main/docs',
        'Source': 'https://github.com/AyehBlk/BondForge',
        'Bug Reports': 'https://github.com/AyehBlk/BondForge/issues',
        'Discussions': 'https://github.com/AyehBlk/BondForge/discussions',
    },
    
    # Package configuration
    packages=find_packages(where='src', exclude=['tests', 'tests.*']),
    package_dir={'': 'src'},
    include_package_data=True,
    
    # Package data
    package_data={
        'bondforge': [
            'data/*.yaml',
            'data/*.json',
            'templates/*',
        ],
    },
    
    # Dependencies
    python_requires='>=3.7',
    install_requires=requirements,
    
    # Optional dependencies for different use cases
    extras_require={
        # Development tools
        'dev': [
            'pytest>=6.2.0',
            'pytest-cov>=3.0.0',
            'pytest-xdist>=2.4.0',
            'black>=21.0',
            'flake8>=3.9.0',
            'mypy>=0.910',
            'pylint>=2.11.0',
            'pre-commit>=2.15.0',
            'hypothesis>=6.23.0',
        ],
        
        # Documentation
        'docs': [
            'sphinx>=4.0.0',
            'sphinx-rtd-theme>=0.5.0',
            'sphinx-click>=3.0.0',
            'myst-parser>=0.15.0',
            'sphinx-autodoc-typehints>=1.12.0',
        ],
        
        # Machine learning features
        'ml': [
            'scikit-learn>=1.0.0',
            'xgboost>=1.5.0',
            'tensorflow>=2.7.0',
        ],
        
        # Advanced visualization
        'viz': [
            'plotly>=5.3.0',
            'py3Dmol>=1.8.0',
            'ipywidgets>=7.6.0',
        ],
        
        # Molecular dynamics analysis
        'md': [
            'MDAnalysis>=2.0.0',
        ],
        
        # Ligand analysis
        'ligand': [
            'rdkit>=2021.03.1',
        ],
        
        # Performance optimization
        'performance': [
            'ray>=1.9.0',
            'joblib>=1.1.0',
        ],
        
        # Complete installation (all features)
        'all': [
            'pytest>=6.2.0',
            'pytest-cov>=3.0.0',
            'black>=21.0',
            'sphinx>=4.0.0',
            'sphinx-rtd-theme>=0.5.0',
            'scikit-learn>=1.0.0',
            'xgboost>=1.5.0',
            'plotly>=5.3.0',
            'py3Dmol>=1.8.0',
            'MDAnalysis>=2.0.0',
            'ray>=1.9.0',
        ],
    },
    
    # Console scripts (CLI commands)
    entry_points={
        'console_scripts': [
            # Main CLI command
            'bondforge=bondforge.cli:main',
            
            # Alternative command names
            'forge=bondforge.cli:main',
            'bf=bondforge.cli:main',
            
            # Specific analysis commands
            'bondforge-analyze=bondforge.cli:analyze',
            'bondforge-batch=bondforge.cli:batch',
            'bondforge-compare=bondforge.cli:compare',
            'bondforge-visualize=bondforge.cli:visualize',
        ],
    },
    
    # PyPI classifiers
    classifiers=[
        # Development status
        'Development Status :: 4 - Beta',
        
        # Intended audience
        'Intended Audience :: Science/Research',
        'Intended Audience :: Developers',
        'Intended Audience :: Healthcare Industry',
        
        # Topics
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'Topic :: Scientific/Engineering :: Chemistry',
        'Topic :: Scientific/Engineering :: Medical Science Apps.',
        'Topic :: Software Development :: Libraries :: Python Modules',
        
        # License
        'License :: OSI Approved :: MIT License',
        
        # Python versions
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
        'Programming Language :: Python :: 3.9',
        'Programming Language :: Python :: 3.10',
        'Programming Language :: Python :: 3.11',
        'Programming Language :: Python :: 3 :: Only',
        
        # Other classifiers
        'Natural Language :: English',
        'Operating System :: OS Independent',
        'Environment :: Console',
        'Typing :: Typed',
    ],
    
    # Keywords for PyPI search
    keywords=[
        'bondforge',
        'protein',
        'bioinformatics',
        'structural-biology',
        'chemical-interactions',
        'bonds',
        'molecular-interactions',
        'protein-structure',
        'pdb',
        'computational-biology',
        'drug-discovery',
        'protein-engineering',
        'network-analysis',
        'hydrogen-bonds',
        'salt-bridges',
        'hydrophobic-interactions',
        'pi-stacking',
        'force-field',
        'energy-calculations',
    ],
    
    # Zip safety
    zip_safe=False,
    
    # Platform
    platforms=['any'],
    
    # License
    license='MIT',
)
