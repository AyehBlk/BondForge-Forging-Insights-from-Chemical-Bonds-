# Contributing to BondForge

First off, thank you for considering contributing to BondForge! It's people like you that make this tool better for the entire structural biology community.

## Code of Conduct

This project and everyone participating in it is governed by our commitment to fostering an open and welcoming environment. We expect all contributors to be respectful and professional.

## How Can I Contribute?

### Reporting Bugs

Before creating bug reports, please check the existing issues to avoid duplicates. When you create a bug report, include as many details as possible:

- **Use a clear and descriptive title**
- **Describe the exact steps to reproduce the problem**
- **Provide specific examples** (code snippets, PDB files if possible)
- **Describe the behavior you observed and expected**
- **Include your environment details** (Python version, OS, BioPython version)

### Suggesting Enhancements

Enhancement suggestions are tracked as GitHub issues. When creating an enhancement suggestion:

- **Use a clear and descriptive title**
- **Provide a detailed description** of the suggested enhancement
- **Explain why this enhancement would be useful** to most users
- **List any scientific references** if applicable

### Pull Requests

1. Fork the repo and create your branch from `main`
2. If you've added code, add tests
3. Ensure the test suite passes
4. Make sure your code follows the existing style
5. Write a clear commit message
6. Include relevant issue numbers in your PR description

## Development Setup

```bash
# Clone your fork
git clone https://github.com/AyehBlk/BondForge.git
cd BondForge

# Create a virtual environment
python -m venv venv
source venv/bin/activate  # On Windows: venv\Scripts\activate

# Install in development mode with dev dependencies
pip install -e .[dev]

# Run tests
pytest tests/
```

## Coding Standards

### Python Style Guide

- Follow [PEP 8](https://www.python.org/dev/peps/pep-0008/)
- Use type hints where appropriate
- Maximum line length: 100 characters
- Use meaningful variable names

### Documentation

- All public functions must have docstrings
- Use Google-style docstrings
- Include parameter types and return values
- Add usage examples for complex functions

Example:
```python
def find_hydrogen_bonds(self, distance_cutoff=3.5, angle_cutoff=120):
    """
    Identify hydrogen bonds in the protein structure.
    
    Args:
        distance_cutoff (float): Maximum D-A distance in Angstroms. Default: 3.5
        angle_cutoff (float): Minimum D-H-A angle in degrees. Default: 120
    
    Returns:
        list: List of dictionaries containing hydrogen bond information
        
    Example:
        >>> analyzer = ProteinInteractionAnalyzer('protein.pdb')
        >>> h_bonds = analyzer.find_hydrogen_bonds()
        >>> print(f"Found {len(h_bonds)} hydrogen bonds")
    """
```

### Testing

- Write unit tests for new features
- Aim for >80% code coverage
- Test edge cases
- Use pytest framework

```python
def test_hydrogen_bond_detection():
    analyzer = ProteinInteractionAnalyzer('test_structure.pdb')
    h_bonds = analyzer.find_hydrogen_bonds()
    assert len(h_bonds) > 0
    assert 'donor' in h_bonds[0]
    assert 'acceptor' in h_bonds[0]
```

## Scientific Accuracy

### Adding New Interaction Types

When proposing a new interaction type:

1. **Provide scientific references** - Include at least 2-3 peer-reviewed papers
2. **Define clear criteria** - Distance cutoffs, angle requirements, etc.
3. **Explain biological relevance** - Why is this interaction important?
4. **Include energy estimates** - If available from literature
5. **Test on real structures** - Validate against known examples

### Algorithm Modifications

- Document the scientific basis for changes
- Compare with existing implementations if available
- Validate against experimental data when possible
- Include references in comments

## Commit Messages

- Use present tense ("Add feature" not "Added feature")
- Use imperative mood ("Move cursor to..." not "Moves cursor to...")
- Limit first line to 72 characters
- Reference issues and pull requests

Examples:
```
Add anion-pi interaction detection

Implements detection of anion-pi interactions based on
McGaughey et al. (1998). Includes distance and angle criteria.

Fixes #123
```

## Branch Naming

- `feature/` - New features
- `bugfix/` - Bug fixes
- `docs/` - Documentation updates
- `test/` - Test additions or modifications

Examples:
- `feature/metal-coordination`
- `bugfix/hydrogen-bond-angles`
- `docs/api-reference`

## Review Process

1. All submissions require review
2. Reviewers will check:
   - Code quality and style
   - Scientific accuracy
   - Test coverage
   - Documentation completeness
3. Address review comments promptly
4. Once approved, maintainers will merge

## License

By contributing, you agree that your contributions will be licensed under the MIT License.

## Questions?

Feel free to open an issue with the `question` label, or contact the maintainers directly.

## Recognition

Contributors will be acknowledged in:
- README.md
- CHANGELOG.md
- Release notes

Thank you for making protein interaction analysis better for everyone! 🧬✨
