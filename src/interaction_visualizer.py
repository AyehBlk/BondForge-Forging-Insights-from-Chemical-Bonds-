"""
BondForge - Visualization Module
=================================

Forging Insights from Chemical Bonds - Visual Representation

This module provides visualization functions for BondForge interaction analysis.
Creates publication-ready figures, network graphs, and 3D visualizations.

Part of BondForge: Comprehensive protein interaction analysis toolkit

Author: Your Name
Date: 2025
"""

import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import pandas as pd
import networkx as nx
from matplotlib.patches import Circle, FancyBboxPatch
from matplotlib.collections import LineCollection
import warnings
warnings.filterwarnings('ignore')


class InteractionVisualizer:
    """
    Class for visualizing protein interaction analysis results
    """
    
    def __init__(self, analyzer):
        """
        Initialize with InteractionAnalyzer object
        
        Parameters:
        -----------
        analyzer : InteractionAnalyzer
            Analyzer object containing interaction data
        """
        self.analyzer = analyzer
        self.interactions = analyzer.interactions
        self.hubs = analyzer.hubs
        
        # Set style
        sns.set_style("whitegrid")
        self.colors = {
            'hydrogen_bond': '#3498db',  # Blue
            'salt_bridge': '#e74c3c',     # Red
            'disulfide_bond': '#f39c12',  # Orange
            'hydrophobic': '#95a5a6',     # Gray
            'pi_pi_stacking': '#9b59b6',  # Purple
            'cation_pi': '#16a085',       # Teal
            'halogen_bond': '#27ae60'     # Green
        }
    
    
    def plot_interaction_summary(self, save_path=None):
        """
        Create summary plots of all interactions
        """
        fig, axes = plt.subplots(2, 2, figsize=(15, 12))
        
        # 1. Interaction type distribution
        self._plot_interaction_types(axes[0, 0])
        
        # 2. Interaction distance distribution
        self._plot_distance_distribution(axes[0, 1])
        
        # 3. Hub degree distribution
        self._plot_hub_distribution(axes[1, 0])
        
        # 4. Criticality distribution
        self._plot_criticality_distribution(axes[1, 1])
        
        plt.tight_layout()
        
        if save_path:
            plt.savefig(save_path, dpi=300, bbox_inches='tight')
            print(f"Summary plot saved to {save_path}")
        
        return fig
    
    
    def _plot_interaction_types(self, ax):
        """Plot distribution of interaction types"""
        type_counts = {}
        for interaction in self.interactions:
            itype = interaction.type
            type_counts[itype] = type_counts.get(itype, 0) + 1
        
        types = list(type_counts.keys())
        counts = list(type_counts.values())
        colors = [self.colors.get(t, '#34495e') for t in types]
        
        ax.bar(range(len(types)), counts, color=colors, alpha=0.7, edgecolor='black')
        ax.set_xticks(range(len(types)))
        ax.set_xticklabels(types, rotation=45, ha='right')
        ax.set_ylabel('Number of Interactions')
        ax.set_title('Interaction Type Distribution', fontsize=14, fontweight='bold')
        ax.grid(axis='y', alpha=0.3)
    
    
    def _plot_distance_distribution(self, ax):
        """Plot distribution of interaction distances"""
        distances_by_type = {}
        
        for interaction in self.interactions:
            itype = interaction.type
            if itype not in distances_by_type:
                distances_by_type[itype] = []
            distances_by_type[itype].append(interaction.distance)
        
        # Create box plot
        data = []
        labels = []
        colors_list = []
        
        for itype, distances in distances_by_type.items():
            data.append(distances)
            labels.append(itype)
            colors_list.append(self.colors.get(itype, '#34495e'))
        
        bp = ax.boxplot(data, labels=labels, patch_artist=True)
        
        for patch, color in zip(bp['boxes'], colors_list):
            patch.set_facecolor(color)
            patch.set_alpha(0.7)
        
        ax.set_xticklabels(labels, rotation=45, ha='right')
        ax.set_ylabel('Distance (Å)')
        ax.set_title('Interaction Distance Distribution', fontsize=14, fontweight='bold')
        ax.grid(axis='y', alpha=0.3)
    
    
    def _plot_hub_distribution(self, ax):
        """Plot hub degree distribution"""
        if not self.hubs:
            ax.text(0.5, 0.5, 'No hub residues found', 
                   ha='center', va='center', fontsize=12)
            ax.set_title('Hub Degree Distribution', fontsize=14, fontweight='bold')
            return
        
        degrees = [hub.degree for hub in self.hubs]
        
        ax.hist(degrees, bins=20, color='#3498db', alpha=0.7, edgecolor='black')
        ax.set_xlabel('Interaction Degree')
        ax.set_ylabel('Number of Hubs')
        ax.set_title('Hub Degree Distribution', fontsize=14, fontweight='bold')
        ax.grid(axis='y', alpha=0.3)
        
        # Add mean line
        mean_degree = np.mean(degrees)
        ax.axvline(mean_degree, color='red', linestyle='--', linewidth=2, 
                  label=f'Mean: {mean_degree:.1f}')
        ax.legend()
    
    
    def _plot_criticality_distribution(self, ax):
        """Plot criticality level distribution"""
        if not self.analyzer.critical_interactions:
            ax.text(0.5, 0.5, 'No critical interactions analyzed', 
                   ha='center', va='center', fontsize=12)
            ax.set_title('Criticality Distribution', fontsize=14, fontweight='bold')
            return
        
        levels = [c['level'] for c in self.analyzer.critical_interactions]
        level_counts = pd.Series(levels).value_counts()
        
        colors_map = {
            'critical': '#e74c3c',
            'important': '#f39c12',
            'moderate': '#3498db',
            'minor': '#95a5a6'
        }
        
        colors = [colors_map.get(level, '#34495e') for level in level_counts.index]
        
        ax.bar(range(len(level_counts)), level_counts.values, 
              color=colors, alpha=0.7, edgecolor='black')
        ax.set_xticks(range(len(level_counts)))
        ax.set_xticklabels(level_counts.index, rotation=45, ha='right')
        ax.set_ylabel('Number of Interactions')
        ax.set_title('Criticality Level Distribution', fontsize=14, fontweight='bold')
        ax.grid(axis='y', alpha=0.3)
    
    
    def plot_interaction_network(self, top_n=50, save_path=None):
        """
        Create network graph of interactions
        
        Parameters:
        -----------
        top_n : int
            Number of top interactions to show (by criticality)
        """
        # Create network graph
        G = nx.Graph()
        
        # Select top interactions
        if self.analyzer.critical_interactions:
            interactions_to_plot = [
                c['interaction'] 
                for c in self.analyzer.critical_interactions[:top_n]
            ]
        else:
            interactions_to_plot = self.interactions[:top_n]
        
        # Add nodes and edges
        for interaction in interactions_to_plot:
            res1 = f"{interaction.residue1_id[2]}{interaction.residue1_id[1]}"
            res2 = f"{interaction.residue2_id[2]}{interaction.residue2_id[1]}"
            
            G.add_node(res1)
            G.add_node(res2)
            G.add_edge(res1, res2, 
                      type=interaction.type,
                      distance=interaction.distance,
                      strength=interaction.strength or 0.5)
        
        # Create figure
        fig, ax = plt.subplots(figsize=(16, 12))
        
        # Layout
        pos = nx.spring_layout(G, k=2, iterations=50)
        
        # Calculate node sizes based on degree
        degrees = dict(G.degree())
        node_sizes = [degrees[node] * 200 + 300 for node in G.nodes()]
        
        # Identify hub nodes
        hub_residues = set()
        for hub in self.hubs:
            res_label = f"{hub.residue_id[2]}{hub.residue_id[1]}"
            if res_label in G.nodes():
                hub_residues.add(res_label)
        
        # Color nodes
        node_colors = ['#e74c3c' if node in hub_residues else '#3498db' 
                      for node in G.nodes()]
        
        # Draw network
        nx.draw_networkx_nodes(G, pos, 
                              node_size=node_sizes,
                              node_color=node_colors,
                              alpha=0.7,
                              ax=ax)
        
        # Draw edges by type
        for itype, color in self.colors.items():
            edges = [(u, v) for u, v, d in G.edges(data=True) if d['type'] == itype]
            if edges:
                nx.draw_networkx_edges(G, pos, edges,
                                      edge_color=color,
                                      width=2,
                                      alpha=0.6,
                                      ax=ax)
        
        # Draw labels
        nx.draw_networkx_labels(G, pos, 
                               font_size=8,
                               font_weight='bold',
                               ax=ax)
        
        # Add legend
        legend_elements = [
            plt.Line2D([0], [0], color=color, linewidth=3, label=itype)
            for itype, color in self.colors.items()
        ]
        legend_elements.extend([
            plt.Line2D([0], [0], marker='o', color='w', 
                      markerfacecolor='#e74c3c', markersize=10, label='Hub'),
            plt.Line2D([0], [0], marker='o', color='w', 
                      markerfacecolor='#3498db', markersize=10, label='Regular')
        ])
        
        ax.legend(handles=legend_elements, loc='upper left', fontsize=10)
        ax.set_title(f'Protein Interaction Network (Top {top_n} Interactions)', 
                    fontsize=16, fontweight='bold')
        ax.axis('off')
        
        plt.tight_layout()
        
        if save_path:
            plt.savefig(save_path, dpi=300, bbox_inches='tight')
            print(f"Network plot saved to {save_path}")
        
        return fig
    
    
    def plot_residue_heatmap(self, save_path=None):
        """
        Create heatmap showing interaction frequency between residue types
        """
        # Create matrix of residue-residue interactions
        residue_types = ['ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', 'GLY', 
                        'HIS', 'ILE', 'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER', 
                        'THR', 'TRP', 'TYR', 'VAL']
        
        matrix = np.zeros((len(residue_types), len(residue_types)))
        res_to_idx = {res: i for i, res in enumerate(residue_types)}
        
        for interaction in self.interactions:
            res1 = interaction.residue1_id[2]
            res2 = interaction.residue2_id[2]
            
            if res1 in res_to_idx and res2 in res_to_idx:
                i, j = res_to_idx[res1], res_to_idx[res2]
                matrix[i, j] += 1
                matrix[j, i] += 1  # Symmetric
        
        # Create heatmap
        fig, ax = plt.subplots(figsize=(14, 12))
        
        im = ax.imshow(matrix, cmap='YlOrRd', aspect='auto')
        
        # Set ticks
        ax.set_xticks(range(len(residue_types)))
        ax.set_yticks(range(len(residue_types)))
        ax.set_xticklabels(residue_types, rotation=90)
        ax.set_yticklabels(residue_types)
        
        # Add colorbar
        cbar = plt.colorbar(im, ax=ax)
        cbar.set_label('Number of Interactions', rotation=270, labelpad=20)
        
        # Add title
        ax.set_title('Residue-Residue Interaction Frequency', 
                    fontsize=16, fontweight='bold', pad=20)
        
        plt.tight_layout()
        
        if save_path:
            plt.savefig(save_path, dpi=300, bbox_inches='tight')
            print(f"Heatmap saved to {save_path}")
        
        return fig
    
    
    def plot_hub_analysis(self, top_n=10, save_path=None):
        """
        Create detailed visualization of top hub residues
        """
        if not self.hubs:
            print("No hub residues to visualize")
            return None
        
        top_hubs = self.hubs[:top_n]
        
        fig, axes = plt.subplots(1, 2, figsize=(16, 6))
        
        # 1. Hub scores
        hub_labels = [f"{h.residue_id[2]}{h.residue_id[1]}" for h in top_hubs]
        hub_scores = [h.hub_score for h in top_hubs]
        
        axes[0].barh(range(len(hub_labels)), hub_scores, color='#3498db', alpha=0.7)
        axes[0].set_yticks(range(len(hub_labels)))
        axes[0].set_yticklabels(hub_labels)
        axes[0].set_xlabel('Hub Score')
        axes[0].set_title('Top Hub Residues by Score', fontsize=14, fontweight='bold')
        axes[0].invert_yaxis()
        axes[0].grid(axis='x', alpha=0.3)
        
        # 2. Interaction diversity
        diversities = [h.diversity for h in top_hubs]
        degrees = [h.degree for h in top_hubs]
        
        scatter = axes[1].scatter(diversities, degrees, 
                                 s=[h.hub_score * 1000 for h in top_hubs],
                                 c=range(len(top_hubs)),
                                 cmap='viridis',
                                 alpha=0.6,
                                 edgecolor='black',
                                 linewidth=1.5)
        
        # Add labels
        for i, hub in enumerate(top_hubs):
            label = f"{hub.residue_id[2]}{hub.residue_id[1]}"
            axes[1].annotate(label, (diversities[i], degrees[i]),
                           xytext=(5, 5), textcoords='offset points',
                           fontsize=8, fontweight='bold')
        
        axes[1].set_xlabel('Interaction Diversity (# of types)')
        axes[1].set_ylabel('Degree (# of interactions)')
        axes[1].set_title('Hub Diversity vs. Degree', fontsize=14, fontweight='bold')
        axes[1].grid(alpha=0.3)
        
        # Add colorbar
        cbar = plt.colorbar(scatter, ax=axes[1])
        cbar.set_label('Hub Rank', rotation=270, labelpad=20)
        
        plt.tight_layout()
        
        if save_path:
            plt.savefig(save_path, dpi=300, bbox_inches='tight')
            print(f"Hub analysis plot saved to {save_path}")
        
        return fig
    
    
    def plot_criticality_analysis(self, top_n=20, save_path=None):
        """
        Visualize critical interactions
        """
        if not self.analyzer.critical_interactions:
            print("No critical interactions to visualize")
            return None
        
        top_critical = self.analyzer.critical_interactions[:top_n]
        
        fig, ax = plt.subplots(figsize=(14, 8))
        
        # Prepare data
        labels = []
        scores = []
        colors = []
        
        color_map = {
            'critical': '#e74c3c',
            'important': '#f39c12',
            'moderate': '#3498db',
            'minor': '#95a5a6'
        }
        
        for c in top_critical:
            interaction = c['interaction']
            res1 = f"{interaction.residue1_id[2]}{interaction.residue1_id[1]}"
            res2 = f"{interaction.residue2_id[2]}{interaction.residue2_id[1]}"
            label = f"{res1}-{res2}\n({interaction.type})"
            
            labels.append(label)
            scores.append(c['criticality_score'])
            colors.append(color_map[c['level']])
        
        # Create horizontal bar chart
        y_pos = np.arange(len(labels))
        ax.barh(y_pos, scores, color=colors, alpha=0.7, edgecolor='black')
        
        ax.set_yticks(y_pos)
        ax.set_yticklabels(labels, fontsize=8)
        ax.set_xlabel('Criticality Score')
        ax.set_title(f'Top {top_n} Critical Interactions', 
                    fontsize=16, fontweight='bold')
        ax.invert_yaxis()
        ax.grid(axis='x', alpha=0.3)
        
        # Add legend
        legend_elements = [
            plt.Rectangle((0, 0), 1, 1, facecolor=color, edgecolor='black', label=level)
            for level, color in color_map.items()
        ]
        ax.legend(handles=legend_elements, loc='lower right')
        
        plt.tight_layout()
        
        if save_path:
            plt.savefig(save_path, dpi=300, bbox_inches='tight')
            print(f"Criticality analysis plot saved to {save_path}")
        
        return fig
    
    
    def create_all_visualizations(self, output_dir='./visualizations'):
        """
        Create all visualizations and save to directory
        """
        import os
        os.makedirs(output_dir, exist_ok=True)
        
        print(f"\nGenerating all visualizations in {output_dir}/...\n")
        
        # 1. Summary
        self.plot_interaction_summary(
            save_path=os.path.join(output_dir, 'summary.png')
        )
        
        # 2. Network
        self.plot_interaction_network(
            save_path=os.path.join(output_dir, 'network.png')
        )
        
        # 3. Heatmap
        self.plot_residue_heatmap(
            save_path=os.path.join(output_dir, 'heatmap.png')
        )
        
        # 4. Hub analysis
        if self.hubs:
            self.plot_hub_analysis(
                save_path=os.path.join(output_dir, 'hubs.png')
            )
        
        # 5. Criticality analysis
        if self.analyzer.critical_interactions:
            self.plot_criticality_analysis(
                save_path=os.path.join(output_dir, 'criticality.png')
            )
        
        print("\nAll visualizations created successfully!\n")
    
    
    def plot_3d_structure_with_interactions(self, view_interactions=True):
        """
        Create 3D visualization using py3Dmol (for Jupyter/web display)
        This requires py3Dmol library
        """
        try:
            import py3Dmol
        except ImportError:
            print("py3Dmol not installed. Install with: pip install py3Dmol")
            return None
        
        view = py3Dmol.view(width=800, height=600)
        
        # Add protein structure
        with open(self.analyzer.pdb_file, 'r') as f:
            pdb_data = f.read()
        
        view.addModel(pdb_data, 'pdb')
        view.setStyle({'cartoon': {'color': 'spectrum'}})
        
        if view_interactions:
            # Highlight interaction residues
            for interaction in self.interactions[:50]:  # Show top 50
                res1_id = interaction.residue1_id[1]
                res2_id = interaction.residue2_id[1]
                
                # Color by interaction type
                color = self.colors.get(interaction.type, 'gray')
                
                view.addStyle({'resi': [res1_id, res2_id]}, 
                            {'stick': {'color': color}})
        
        view.zoomTo()
        return view


# ============================================================================
# INTEGRATION WITH GUI
# ============================================================================

class VisualizationPanel:
    """
    Class to integrate visualizations into PyQt5 GUI
    """
    
    def __init__(self, parent_widget):
        """
        Initialize visualization panel
        
        Parameters:
        -----------
        parent_widget : QWidget
            Parent Qt widget to embed matplotlib figures
        """
        self.parent = parent_widget
        
    
    def embed_matplotlib_figure(self, fig):
        """
        Embed matplotlib figure into Qt widget
        """
        from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg
        
        canvas = FigureCanvasQTAgg(fig)
        return canvas


# ============================================================================
# EXAMPLE USAGE
# ============================================================================

if __name__ == "__main__":
    # Assuming you have already run the analyzer
    from interaction_analyzer import InteractionAnalyzer
    
    # Load and analyze
    analyzer = InteractionAnalyzer("example.pdb")
    analyzer.analyze_all()
    
    # Create visualizer
    visualizer = InteractionVisualizer(analyzer)
    
    # Generate all visualizations
    visualizer.create_all_visualizations(output_dir='./results/visualizations')
    
    # Or create individual plots
    # visualizer.plot_interaction_summary()
    # visualizer.plot_interaction_network()
    # plt.show()
