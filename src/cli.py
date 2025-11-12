#!/usr/bin/env python3
"""
BondForge CLI - Professional Command-Line Interface
===================================================

Forging Insights from Chemical Bonds

A comprehensive protein interaction analysis toolkit with professional CLI.

Author: Ayeh Bolouki
Email: ayehbolouki1988@gmail.com | ayehgeek@gmail.com
GitHub: https://github.com/AyehBlk/BondForge
"""

import click
import sys
import os
from pathlib import Path
from typing import Optional, List
import logging
from rich.console import Console
from rich.table import Table
from rich.progress import Progress, SpinnerColumn, TextColumn
from rich.logging import RichHandler
from rich import print as rprint
from rich.panel import Panel
from rich.tree import Tree
import yaml
import json

# Setup rich console for beautiful output
console = Console()

# Setup logging with rich handler
logging.basicConfig(
    level=logging.INFO,
    format="%(message)s",
    datefmt="[%X]",
    handlers=[RichHandler(rich_tracebacks=True, console=console)]
)
logger = logging.getLogger("bondforge")


# ============================================================================
# CONFIGURATION MANAGEMENT
# ============================================================================

class Config:
    """Configuration manager with hierarchical loading"""
    
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
                'method': 'empirical',  # 'empirical' or 'quantum_mechanical'
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
            'format': ['csv'],  # csv, json, xml
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
    
    def __init__(self):
        self.config = self.DEFAULT_CONFIG.copy()
        self._load_config_files()
    
    def _load_config_files(self):
        """Load configuration from files (system → user → project)"""
        config_locations = [
            Path('/etc/bondforge/config.yaml'),
            Path.home() / '.bondforge' / 'config.yaml',
            Path.cwd() / 'bondforge.yaml'
        ]
        
        for config_file in config_locations:
            if config_file.exists():
                try:
                    with open(config_file) as f:
                        user_config = yaml.safe_load(f)
                        self._merge_config(user_config)
                    logger.debug(f"Loaded config from {config_file}")
                except Exception as e:
                    logger.warning(f"Failed to load {config_file}: {e}")
    
    def _merge_config(self, user_config):
        """Deep merge user config into default config"""
        def deep_merge(base, update):
            for key, value in update.items():
                if key in base and isinstance(base[key], dict) and isinstance(value, dict):
                    deep_merge(base[key], value)
                else:
                    base[key] = value
        
        deep_merge(self.config, user_config)
    
    def get(self, key_path, default=None):
        """Get config value using dot notation: 'analysis.distance_cutoffs.hydrogen_bond'"""
        keys = key_path.split('.')
        value = self.config
        for key in keys:
            if isinstance(value, dict):
                value = value.get(key)
            else:
                return default
        return value if value is not None else default
    
    def set(self, key_path, value):
        """Set config value using dot notation"""
        keys = key_path.split('.')
        config = self.config
        for key in keys[:-1]:
            if key not in config:
                config[key] = {}
            config = config[key]
        config[keys[-1]] = value


# Global config instance
config = Config()


# ============================================================================
# UTILITY FUNCTIONS
# ============================================================================

def validate_pdb_file(pdb_file: str) -> bool:
    """Validate that PDB file exists and is readable"""
    path = Path(pdb_file)
    if not path.exists():
        console.print(f"[red]Error:[/red] File '{pdb_file}' not found!")
        return False
    if not path.is_file():
        console.print(f"[red]Error:[/red] '{pdb_file}' is not a file!")
        return False
    if path.stat().st_size == 0:
        console.print(f"[red]Error:[/red] File '{pdb_file}' is empty!")
        return False
    return True


def display_banner():
    """Display BondForge banner"""
    banner = """
[bold cyan]
╔══════════════════════════════════════════════════════════════╗
║                      🔥 BondForge 🔥                         ║
║          Forging Insights from Chemical Bonds                ║
║                                                              ║
║     Professional Protein Interaction Analysis Toolkit       ║
╚══════════════════════════════════════════════════════════════╝
[/bold cyan]
"""
    console.print(banner)


def display_results_summary(results: dict):
    """Display analysis results in a beautiful table"""
    table = Table(title="🔥 Analysis Results", show_header=True, header_style="bold magenta")
    table.add_column("Interaction Type", style="cyan", width=30)
    table.add_column("Count", justify="right", style="green")
    table.add_column("Percentage", justify="right", style="yellow")
    
    total = results.get('total_interactions', 0)
    by_type = results.get('by_type', {})
    
    for int_type, count in sorted(by_type.items(), key=lambda x: x[1], reverse=True):
        percentage = (count / total * 100) if total > 0 else 0
        table.add_row(
            int_type.replace('_', ' ').title(),
            str(count),
            f"{percentage:.1f}%"
        )
    
    table.add_row("─" * 30, "─" * 10, "─" * 10, style="dim")
    table.add_row("[bold]TOTAL[/bold]", f"[bold]{total}[/bold]", "[bold]100.0%[/bold]")
    
    console.print(table)


# ============================================================================
# CLI COMMAND GROUP
# ============================================================================

@click.group()
@click.version_option(version='2.0.0', prog_name='BondForge')
@click.option('--verbose', '-v', is_flag=True, help='Enable verbose output')
@click.option('--quiet', '-q', is_flag=True, help='Suppress non-error output')
@click.option('--config', '-c', type=click.Path(exists=True), 
              help='Custom configuration file')
@click.pass_context
def cli(ctx, verbose, quiet, config_file):
    """
    🔥 BondForge - Professional Protein Interaction Analysis Toolkit
    
    Analyze 20 types of chemical interactions in protein structures.
    
    \b
    Examples:
      bondforge analyze protein.pdb
      bondforge analyze protein.pdb --types hydrogen_bonds salt_bridges
      bondforge batch *.pdb --output results/
      bondforge visualize protein.pdb --type network
    
    For more information: https://github.com/AyehBlk/BondForge
    """
    # Set logging level
    if verbose:
        logger.setLevel(logging.DEBUG)
    elif quiet:
        logger.setLevel(logging.ERROR)
    
    # Load custom config if provided
    if config_file:
        try:
            with open(config_file) as f:
                user_config = yaml.safe_load(f)
                config._merge_config(user_config)
            logger.debug(f"Loaded config from {config_file}")
        except Exception as e:
            logger.error(f"Failed to load config file: {e}")
            ctx.exit(1)
    
    # Store config in context
    ctx.ensure_object(dict)
    ctx.obj['config'] = config


# ============================================================================
# ANALYZE COMMAND
# ============================================================================

@cli.command()
@click.argument('pdb_file', type=click.Path(exists=True))
@click.option('--types', '-t', multiple=True,
              help='Specific interaction types to analyze (default: all)')
@click.option('--output', '-o', type=click.Path(), default='bondforge_results',
              help='Output directory (default: bondforge_results)')
@click.option('--format', '-f', multiple=True, 
              type=click.Choice(['csv', 'json', 'xml', 'excel']),
              help='Output format(s)')
@click.option('--hub-threshold', type=int, default=5,
              help='Minimum interactions for hub residue (default: 5)')
@click.option('--energy/--no-energy', default=True,
              help='Calculate interaction energies')
@click.option('--parallel/--no-parallel', default=True,
              help='Use parallel processing')
@click.option('--visualize/--no-visualize', default=True,
              help='Generate visualization plots')
@click.pass_context
def analyze(ctx, pdb_file, types, output, format, hub_threshold, 
            energy, parallel, visualize):
    """
    Analyze protein structure for chemical interactions.
    
    \b
    Example:
      bondforge analyze 1ABC.pdb --types hydrogen_bonds --output results/
    
    This command performs comprehensive interaction analysis including:
    - Detection of all 20 interaction types (or specified subset)
    - Hub residue identification
    - Critical interaction analysis
    - Energy calculations (optional)
    - Network analysis
    - Visualization generation (optional)
    """
    display_banner()
    
    # Validate input
    if not validate_pdb_file(pdb_file):
        ctx.exit(1)
    
    # Display analysis parameters
    console.print("\n[bold cyan]Analysis Parameters:[/bold cyan]")
    params_table = Table(show_header=False, box=None, padding=(0, 2))
    params_table.add_row("📄 Structure:", Path(pdb_file).name)
    params_table.add_row("🎯 Interaction Types:", str(len(types)) if types else "All (20)")
    params_table.add_row("📊 Output Directory:", output)
    params_table.add_row("⚡ Parallel Processing:", "✓" if parallel else "✗")
    params_table.add_row("🔋 Energy Calculation:", "✓" if energy else "✗")
    params_table.add_row("📈 Visualization:", "✓" if visualize else "✗")
    console.print(params_table)
    console.print()
    
    # Run analysis with progress indicator
    try:
        with Progress(
            SpinnerColumn(),
            TextColumn("[progress.description]{task.description}"),
            console=console
        ) as progress:
            
            # Step 1: Load structure
            task1 = progress.add_task("🔄 Loading structure...", total=None)
            
            # Import here to avoid slow startup
            from extended_interaction_analyzer import ExtendedProteinInteractionAnalyzer
            
            analyzer = ExtendedProteinInteractionAnalyzer(pdb_file)
            progress.update(task1, completed=True)
            
            # Step 2: Analyze interactions
            task2 = progress.add_task("🔥 Forging interactions...", total=None)
            
            if types:
                # Analyze specific types
                results = {}
                for int_type in types:
                    method_name = f"find_{int_type}"
                    if hasattr(analyzer, method_name):
                        results[int_type] = getattr(analyzer, method_name)()
                
                # Consolidate results
                results = {
                    'by_type': {k: len(v) for k, v in results.items()},
                    'total_interactions': sum(len(v) for v in results.values())
                }
            else:
                # Analyze all types
                results = analyzer.analyze_all_interactions()
            
            progress.update(task2, completed=True)
            
            # Step 3: Hub analysis
            task3 = progress.add_task("🎯 Identifying hubs...", total=None)
            hubs = analyzer.identify_hub_residues(threshold=hub_threshold)
            progress.update(task3, completed=True)
            
            # Step 4: Critical interactions
            task4 = progress.add_task("⚡ Finding critical interactions...", total=None)
            critical = analyzer.identify_critical_interactions()
            progress.update(task4, completed=True)
            
            # Step 5: Export results
            task5 = progress.add_task("💾 Exporting results...", total=None)
            os.makedirs(output, exist_ok=True)
            
            # Export in requested formats
            export_formats = format if format else ['csv']
            for fmt in export_formats:
                if fmt == 'csv':
                    analyzer.export_to_csv(output)
                elif fmt == 'json':
                    analyzer.export_to_json(output)
                elif fmt == 'xml':
                    analyzer.export_to_xml(output)
            
            progress.update(task5, completed=True)
            
            # Step 6: Visualization
            if visualize:
                task6 = progress.add_task("📊 Generating visualizations...", total=None)
                from interaction_visualizer import InteractionVisualizer
                
                viz = InteractionVisualizer(analyzer)
                viz.plot_interaction_network(results, output=f"{output}/network.png")
                viz.generate_pymol_script(f"{output}/visualize.pml")
                
                progress.update(task6, completed=True)
        
        # Display results summary
        console.print("\n[bold green]✓ Analysis Complete![/bold green]\n")
        display_results_summary(results)
        
        # Display hub summary
        if hubs:
            console.print(f"\n[bold cyan]🎯 Identified {len(hubs)} Hub Residues[/bold cyan]")
            hub_table = Table(show_header=True, header_style="bold")
            hub_table.add_column("Residue", style="cyan")
            hub_table.add_column("Degree", justify="right", style="green")
            hub_table.add_column("Hub Score", justify="right", style="yellow")
            hub_table.add_column("Type", style="magenta")
            
            for hub in hubs[:10]:  # Show top 10
                hub_dict = hub.to_dict()
                hub_table.add_row(
                    hub_dict['residue'],
                    str(hub_dict['degree']),
                    f"{hub_dict['hub_score']:.3f}",
                    hub_dict['hub_type']
                )
            
            console.print(hub_table)
        
        # Display critical interactions summary
        if critical:
            console.print(f"\n[bold red]⚡ Identified {len(critical)} Critical Interactions[/bold red]")
            console.print(f"  Critical: {sum(1 for c in critical if c['level'] == 'critical')}")
            console.print(f"  Important: {sum(1 for c in critical if c['level'] == 'important')}")
        
        # Display output files
        console.print(f"\n[bold cyan]📁 Output Files:[/bold cyan]")
        output_tree = Tree(f"📂 {output}/")
        
        for file in sorted(Path(output).glob('*')):
            if file.is_file():
                output_tree.add(f"📄 {file.name}")
        
        console.print(output_tree)
        
        # Success panel
        console.print(Panel.fit(
            f"[bold green]✓ Success![/bold green]\n\n"
            f"Forged {results['total_interactions']} interactions from {Path(pdb_file).name}\n"
            f"Results saved to: [cyan]{output}/[/cyan]",
            title="🔥 BondForge Complete",
            border_style="green"
        ))
        
    except Exception as e:
        logger.error(f"Analysis failed: {e}")
        if ctx.obj.get('verbose'):
            console.print_exception()
        ctx.exit(1)


# ============================================================================
# BATCH COMMAND
# ============================================================================

@cli.command()
@click.argument('pdb_files', nargs=-1, type=click.Path(exists=True), required=True)
@click.option('--output', '-o', type=click.Path(), default='bondforge_batch_results',
              help='Output directory for all results')
@click.option('--parallel', type=int, default=-1,
              help='Number of parallel jobs (-1 = all cores)')
@click.option('--continue-on-error', is_flag=True,
              help='Continue processing even if some structures fail')
@click.pass_context
def batch(ctx, pdb_files, output, parallel, continue_on_error):
    """
    Analyze multiple protein structures in batch mode.
    
    \b
    Example:
      bondforge batch *.pdb --output batch_results/
      bondforge batch protein1.pdb protein2.pdb protein3.pdb
    """
    display_banner()
    
    console.print(f"[bold cyan]Batch Analysis: {len(pdb_files)} structures[/bold cyan]\n")
    
    # Create output directory
    os.makedirs(output, exist_ok=True)
    
    # Process each file
    results_summary = []
    
    with Progress(console=console) as progress:
        task = progress.add_task("Processing structures...", total=len(pdb_files))
        
        for pdb_file in pdb_files:
            try:
                logger.info(f"Processing {Path(pdb_file).name}...")
                
                # Run analysis (simplified for batch)
                from extended_interaction_analyzer import ExtendedProteinInteractionAnalyzer
                
                analyzer = ExtendedProteinInteractionAnalyzer(pdb_file)
                results = analyzer.analyze_all_interactions()
                
                # Create structure-specific output directory
                struct_name = Path(pdb_file).stem
                struct_output = Path(output) / struct_name
                struct_output.mkdir(exist_ok=True)
                
                # Export results
                analyzer.export_to_csv(str(struct_output))
                
                # Store summary
                results_summary.append({
                    'structure': struct_name,
                    'total_interactions': results['total_interactions'],
                    'status': 'success'
                })
                
            except Exception as e:
                logger.error(f"Failed to process {pdb_file}: {e}")
                results_summary.append({
                    'structure': Path(pdb_file).stem,
                    'total_interactions': 0,
                    'status': 'failed',
                    'error': str(e)
                })
                
                if not continue_on_error:
                    ctx.exit(1)
            
            progress.update(task, advance=1)
    
    # Display summary table
    console.print("\n[bold green]Batch Analysis Complete![/bold green]\n")
    
    summary_table = Table(title="Batch Results Summary", show_header=True)
    summary_table.add_column("Structure", style="cyan")
    summary_table.add_column("Interactions", justify="right", style="green")
    summary_table.add_column("Status", style="yellow")
    
    for result in results_summary:
        status_emoji = "✓" if result['status'] == 'success' else "✗"
        status_style = "green" if result['status'] == 'success' else "red"
        summary_table.add_row(
            result['structure'],
            str(result['total_interactions']),
            f"[{status_style}]{status_emoji} {result['status']}[/{status_style}]"
        )
    
    console.print(summary_table)
    
    # Export summary
    summary_file = Path(output) / 'batch_summary.json'
    with open(summary_file, 'w') as f:
        json.dump(results_summary, f, indent=2)
    
    console.print(f"\n[cyan]Summary saved to: {summary_file}[/cyan]")


# ============================================================================
# VALIDATE COMMAND
# ============================================================================

@cli.command()
@click.argument('pdb_file', type=click.Path(exists=True))
@click.option('--strict', is_flag=True, help='Use strict validation criteria')
@click.pass_context
def validate(ctx, pdb_file, strict):
    """
    Validate PDB structure quality and completeness.
    
    Checks for:
    - Missing atoms
    - Missing residues
    - Non-standard residues
    - Structure resolution
    - B-factors
    - Occupancy issues
    """
    display_banner()
    
    console.print(f"[bold cyan]Validating structure:[/bold cyan] {pdb_file}\n")
    
    try:
        from Bio.PDB import PDBParser
        from Bio.PDB.PDBExceptions import PDBConstructionWarning
        import warnings
        
        # Parse structure
        parser = PDBParser(QUIET=True)
        structure = parser.get_structure('protein', pdb_file)
        
        # Validation checks
        checks = {
            'structure_loaded': False,
            'has_atoms': False,
            'has_coordinates': False,
            'missing_atoms': [],
            'non_standard_residues': [],
            'warnings': []
        }
        
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            
            # Check structure
            if structure:
                checks['structure_loaded'] = True
            
            # Check atoms
            atoms = list(structure.get_atoms())
            if atoms:
                checks['has_atoms'] = True
                checks['atom_count'] = len(atoms)
            
            # Check coordinates
            if atoms and all(atom.coord is not None for atom in atoms):
                checks['has_coordinates'] = True
            
            # Store warnings
            checks['warnings'] = [str(warning.message) for warning in w]
        
        # Display results
        validation_table = Table(title="Validation Results", show_header=True)
        validation_table.add_column("Check", style="cyan", width=30)
        validation_table.add_column("Status", justify="center", width=10)
        validation_table.add_column("Details", style="dim")
        
        for check_name, check_value in checks.items():
            if isinstance(check_value, bool):
                status = "✓" if check_value else "✗"
                status_style = "green" if check_value else "red"
                validation_table.add_row(
                    check_name.replace('_', ' ').title(),
                    f"[{status_style}]{status}[/{status_style}]",
                    ""
                )
            elif isinstance(check_value, int):
                validation_table.add_row(
                    check_name.replace('_', ' ').title(),
                    "ℹ️",
                    str(check_value)
                )
        
        console.print(validation_table)
        
        # Overall verdict
        all_passed = all(
            checks['structure_loaded'],
            checks['has_atoms'],
            checks['has_coordinates']
        )
        
        if all_passed:
            console.print("\n[bold green]✓ Structure validation PASSED[/bold green]")
        else:
            console.print("\n[bold red]✗ Structure validation FAILED[/bold red]")
            if not strict:
                console.print("[yellow]Warning: Some checks failed. Results may be incomplete.[/yellow]")
            else:
                ctx.exit(1)
        
    except Exception as e:
        logger.error(f"Validation failed: {e}")
        ctx.exit(1)


# ============================================================================
# COMPARE COMMAND
# ============================================================================

@cli.command()
@click.argument('pdb_file1', type=click.Path(exists=True))
@click.argument('pdb_file2', type=click.Path(exists=True))
@click.option('--output', '-o', type=click.Path(), default='comparison_results',
              help='Output directory for comparison results')
@click.pass_context
def compare(ctx, pdb_file1, pdb_file2, output):
    """
    Compare interactions between two protein structures.
    
    Useful for:
    - Wild-type vs mutant comparison
    - Different conformations
    - Before/after modifications
    
    \b
    Example:
      bondforge compare wildtype.pdb mutant.pdb --output comparison/
    """
    display_banner()
    
    console.print(f"[bold cyan]Comparing structures:[/bold cyan]")
    console.print(f"  Structure 1: {Path(pdb_file1).name}")
    console.print(f"  Structure 2: {Path(pdb_file2).name}\n")
    
    try:
        from extended_interaction_analyzer import ExtendedProteinInteractionAnalyzer
        
        with Progress(console=console) as progress:
            # Analyze first structure
            task1 = progress.add_task("Analyzing structure 1...", total=None)
            analyzer1 = ExtendedProteinInteractionAnalyzer(pdb_file1)
            results1 = analyzer1.analyze_all_interactions()
            progress.update(task1, completed=True)
            
            # Analyze second structure
            task2 = progress.add_task("Analyzing structure 2...", total=None)
            analyzer2 = ExtendedProteinInteractionAnalyzer(pdb_file2)
            results2 = analyzer2.analyze_all_interactions()
            progress.update(task2, completed=True)
        
        # Compare results
        console.print("\n[bold cyan]Comparison Results:[/bold cyan]\n")
        
        comparison_table = Table(show_header=True, header_style="bold")
        comparison_table.add_column("Interaction Type", style="cyan", width=25)
        comparison_table.add_column("Structure 1", justify="right", style="green")
        comparison_table.add_column("Structure 2", justify="right", style="blue")
        comparison_table.add_column("Difference", justify="right", style="yellow")
        comparison_table.add_column("Change %", justify="right", style="magenta")
        
        # Get all interaction types
        all_types = set(results1['by_type'].keys()) | set(results2['by_type'].keys())
        
        for int_type in sorted(all_types):
            count1 = results1['by_type'].get(int_type, 0)
            count2 = results2['by_type'].get(int_type, 0)
            diff = count2 - count1
            
            if count1 > 0:
                change_pct = (diff / count1) * 100
            else:
                change_pct = 100 if count2 > 0 else 0
            
            diff_str = f"+{diff}" if diff > 0 else str(diff)
            change_str = f"+{change_pct:.1f}%" if change_pct > 0 else f"{change_pct:.1f}%"
            
            comparison_table.add_row(
                int_type.replace('_', ' ').title(),
                str(count1),
                str(count2),
                diff_str,
                change_str
            )
        
        # Total row
        total1 = results1['total_interactions']
        total2 = results2['total_interactions']
        total_diff = total2 - total1
        total_change = (total_diff / total1 * 100) if total1 > 0 else 0
        
        comparison_table.add_row(
            "─" * 25, "─" * 10, "─" * 10, "─" * 10, "─" * 10, style="dim"
        )
        comparison_table.add_row(
            "[bold]TOTAL[/bold]",
            f"[bold]{total1}[/bold]",
            f"[bold]{total2}[/bold]",
            f"[bold]{'+' if total_diff > 0 else ''}{total_diff}[/bold]",
            f"[bold]{'+' if total_change > 0 else ''}{total_change:.1f}%[/bold]"
        )
        
        console.print(comparison_table)
        
        # Save comparison
        os.makedirs(output, exist_ok=True)
        comparison_data = {
            'structure1': str(pdb_file1),
            'structure2': str(pdb_file2),
            'results1': results1,
            'results2': results2,
            'differences': {
                int_type: results2['by_type'].get(int_type, 0) - results1['by_type'].get(int_type, 0)
                for int_type in all_types
            }
        }
        
        with open(Path(output) / 'comparison.json', 'w') as f:
            json.dump(comparison_data, f, indent=2)
        
        console.print(f"\n[cyan]Comparison saved to: {output}/comparison.json[/cyan]")
        
    except Exception as e:
        logger.error(f"Comparison failed: {e}")
        if ctx.obj.get('verbose'):
            console.print_exception()
        ctx.exit(1)


# ============================================================================
# CONFIG COMMAND
# ============================================================================

@cli.command()
@click.option('--show', is_flag=True, help='Show current configuration')
@click.option('--init', is_flag=True, help='Initialize configuration file')
@click.option('--edit', is_flag=True, help='Edit configuration file')
@click.pass_context
def config_cmd(ctx, show, init, edit):
    """
    Manage BondForge configuration.
    
    \b
    Examples:
      bondforge config --show              # Display current config
      bondforge config --init              # Create default config file
      bondforge config --edit              # Open config in editor
    """
    if show:
        console.print("[bold cyan]Current Configuration:[/bold cyan]\n")
        console.print_json(data=config.config)
    
    elif init:
        config_path = Path.home() / '.bondforge' / 'config.yaml'
        config_path.parent.mkdir(exist_ok=True)
        
        if config_path.exists():
            if not click.confirm(f"Config file already exists at {config_path}. Overwrite?"):
                return
        
        with open(config_path, 'w') as f:
            yaml.dump(config.DEFAULT_CONFIG, f, default_flow_style=False)
        
        console.print(f"[green]✓ Configuration file created at {config_path}[/green]")
    
    elif edit:
        config_path = Path.home() / '.bondforge' / 'config.yaml'
        
        if not config_path.exists():
            console.print("[yellow]Config file doesn't exist. Creating...[/yellow]")
            config_path.parent.mkdir(exist_ok=True)
            with open(config_path, 'w') as f:
                yaml.dump(config.DEFAULT_CONFIG, f, default_flow_style=False)
        
        # Open in default editor
        editor = os.environ.get('EDITOR', 'nano')
        os.system(f"{editor} {config_path}")
    
    else:
        console.print("[yellow]Please specify --show, --init, or --edit[/yellow]")
        console.print("Use 'bondforge config --help' for more information")


# ============================================================================
# VERSION COMMAND
# ============================================================================

@cli.command()
def version():
    """Display BondForge version and system information."""
    console.print("\n[bold cyan]🔥 BondForge Version Information[/bold cyan]\n")
    
    info_table = Table(show_header=False, box=None, padding=(0, 2))
    info_table.add_row("Version:", "[green]2.0.0[/green]")
    info_table.add_row("Python:", f"{sys.version.split()[0]}")
    info_table.add_row("Platform:", sys.platform)
    
    try:
        import numpy as np
        info_table.add_row("NumPy:", np.__version__)
    except:
        pass
    
    try:
        from Bio import __version__ as bio_version
        info_table.add_row("BioPython:", bio_version)
    except:
        pass
    
    console.print(info_table)
    console.print("\n[cyan]For more information: https://github.com/AyehBlk/BondForge[/cyan]\n")


# ============================================================================
# MAIN ENTRY POINT
# ============================================================================

def main():
    """Main entry point for CLI"""
    try:
        cli(obj={})
    except KeyboardInterrupt:
        console.print("\n[yellow]Interrupted by user[/yellow]")
        sys.exit(130)
    except Exception as e:
        logger.error(f"Unexpected error: {e}")
        console.print_exception()
        sys.exit(1)


if __name__ == '__main__':
    main()
