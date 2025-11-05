#!/usr/bin/env python3
"""
BondForge - Forging Insights from Chemical Bonds
Setup and installation script
"""

from setuptools import setup, find_packages
import os

# Read the contents of README file
this_directory = os.path.abspath(os.path.dirname(__file__))
with open(os.path.join(this_directory, 'README.md'), encoding='utf-8') as f:
    long_description = f.read()

# Read requirements
with open(os.path.join(this_directory, 'requirements.txt'), encoding='utf-8') as f:
    requirements = [line.strip() for line in f if line.strip() and not line.startswith('#')]

setup(
    name='bondforge',
    version='1.0.0',
    author='Your Name',
    author_email='your.email@example.com',
    description='BondForge: Forge insights from 20 types of chemical bonds. Comprehensive protein interaction analysis toolkit.',
    long_description=long_description,
    long_description_content_type='text/markdown',
    url='https://github.com/yourusername/BondForge',
    packages=find_packages(where='src'),
    package_dir={'': 'src'},
    classifiers=[
        'Development Status :: 4 - Beta',
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'Topic :: Scientific/Engineering :: Chemistry',
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
        'Programming Language :: Python :: 3.9',
        'Programming Language :: Python :: 3.10',
        'Programming Language :: Python :: 3.11',
    ],
    keywords='bondforge protein bioinformatics structural-biology chemical-interactions bonds forge',
    python_requires='>=3.7',
    install_requires=requirements,
    extras_require={
        'dev': [
            'pytest>=6.0',
            'pytest-cov>=2.0',
            'black>=21.0',
            'flake8>=3.9',
            'mypy>=0.910',
        ],
        'docs': [
            'sphinx>=4.0',
            'sphinx-rtd-theme>=0.5',
        ],
    },
    entry_points={
        'console_scripts': [
            'bondforge=src.extended_interaction_analyzer:main',
            'forge-analyze=src.extended_interaction_analyzer:main',
        ],
    },
    include_package_data=True,
    zip_safe=False,
)
