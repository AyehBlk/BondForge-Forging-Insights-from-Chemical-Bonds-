# Tests

Unit tests and integration tests for the BondForge.

## Status

⚠️ **Test suite is under development**

We're working on comprehensive test coverage. Contributions welcome!

## Planned Test Coverage

### Unit Tests

- [ ] `test_interaction_analyzer.py` - Core analyzer functions
- [ ] `test_extended_analyzer.py` - Extended analyzer functions
- [ ] `test_visualizer.py` - Visualization functions
- [ ] `test_distance_calculations.py` - Geometric calculations
- [ ] `test_angle_calculations.py` - Angular measurements
- [ ] `test_energy_estimates.py` - Energy calculations

### Integration Tests

- [ ] `test_full_analysis.py` - End-to-end analysis pipeline
- [ ] `test_export.py` - Export functionality
- [ ] `test_large_structures.py` - Performance on large proteins
- [ ] `test_multi_chain.py` - Multi-chain structure handling

### Test Data

- [ ] Sample PDB files for testing
- [ ] Known interaction examples
- [ ] Edge case structures
- [ ] Invalid input handling

## Running Tests

Once tests are implemented, run with:

```bash
# Run all tests
pytest tests/

# Run with coverage
pytest tests/ --cov=src --cov-report=html

# Run specific test file
pytest tests/test_interaction_analyzer.py

# Run with verbose output
pytest tests/ -v
```

## Contributing Tests

We especially need help with:

1. **Test Data**: Finding good test structures with known interactions
2. **Edge Cases**: Unusual structures, missing atoms, etc.
3. **Performance Tests**: Benchmarking on large structures
4. **Validation Tests**: Comparing with other tools

See [CONTRIBUTING.md](../CONTRIBUTING.md) for guidelines.

## Test Structure Template

```python
import pytest
from src.interaction_analyzer import ProteinInteractionAnalyzer

class TestHydrogenBonds:
    @pytest.fixture
    def analyzer(self):
        return ProteinInteractionAnalyzer('tests/data/test_structure.pdb')
    
    def test_find_hydrogen_bonds(self, analyzer):
        h_bonds = analyzer.find_hydrogen_bonds()
        assert len(h_bonds) > 0
        assert 'donor' in h_bonds[0]
        assert 'acceptor' in h_bonds[0]
    
    def test_hydrogen_bond_distance(self, analyzer):
        h_bonds = analyzer.find_hydrogen_bonds(distance_cutoff=3.0)
        for bond in h_bonds:
            assert bond['distance'] <= 3.0
```

## CI/CD

Tests will run automatically on:
- Pull requests
- Pushes to main/develop branches
- Release tags

See `.github/workflows/ci.yml` for configuration.

---

**Help Wanted!** 🙏

If you'd like to contribute to test development, please:
1. Open an issue to discuss
2. Submit a PR with new tests
3. Help validate existing functionality

Thank you for helping make this software more reliable!
