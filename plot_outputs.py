#!/usr/bin/env python3
"""
GAPpy Model Output Visualization Script

This script creates comprehensive plots of the GAPpy forest model outputs
including forest dynamics, soil biogeochemistry, and environmental conditions.

Features:
- Forest composition and biomass analysis
- Soil carbon and nitrogen dynamics
- Environmental stress indicators
- Comprehensive dashboard view
- Customizable output formats and styling

Usage:
    python plot_outputs.py [--output-dir OUTPUT_DIR] [--plots-dir PLOTS_DIR]
                          [--format FORMAT] [--show/--no-show]
"""

import argparse
import os
import sys
from pathlib import Path
from typing import Tuple, List, Optional, Dict, Any
import warnings

import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import numpy as np

# Suppress matplotlib warnings
warnings.filterwarnings('ignore', category=UserWarning, module='matplotlib')

# Configuration
DEFAULT_OUTPUT_DIR = "output_data"
DEFAULT_PLOTS_DIR = "output_data/plots"
DEFAULT_FORMAT = "png"
DEFAULT_DPI = 150
DEFAULT_STYLE = 'default'

# Enhanced color palettes for better differentiation
# Use a combination of highly distinguishable colors
SPECIES_COLORS = [
    '#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd',  # Classic 5
    '#8c564b', '#e377c2', '#7f7f7f', '#bcbd22', '#17becf',  # Next 5
    '#1a9850', '#d73027', '#4575b4', '#fee08b', '#313695',  # Additional 5
    '#a50026', '#006837', '#542788', '#f46d43', '#74add1'   # More contrast
]

# Line styles for additional differentiation when colors repeat
LINE_STYLES = ['-', '--', '-.', ':']
MARKERS = ['o', 's', '^', 'v', 'D', 'p', '*', 'h']

SOIL_COLORS = {"surface": "#8B4513", "mineral": "#D2691E", "available": "#228B22"}
CLIMATE_COLORS = {"temperature": "#DC143C", "water": "#1E90FF", "stress": "#FF8C00"}


class GAPpyPlotter:
    """
    Comprehensive plotting class for GAPpy model outputs.

    Handles data loading, processing, and visualization of forest ecosystem
    model outputs including biomass, species composition, soil biogeochemistry,
    and environmental conditions.
    """

    def __init__(self, output_dir: str = DEFAULT_OUTPUT_DIR,
                 plots_dir: str = DEFAULT_PLOTS_DIR,
                 style: str = DEFAULT_STYLE):
        """
        Initialize the plotter with configuration.

        Args:
            output_dir: Directory containing GAPpy output CSV files
            plots_dir: Directory to save generated plots
            style: Matplotlib style to use for plots
        """
        self.output_dir = Path(output_dir)
        self.plots_dir = Path(plots_dir)
        self.style = style

        # Data containers
        self.site_data: Optional[pd.DataFrame] = None
        self.species_data: Optional[pd.DataFrame] = None
        self.soil_data: Optional[pd.DataFrame] = None
        self.genus_data: Optional[pd.DataFrame] = None

        # Create plots directory
        self.plots_dir.mkdir(exist_ok=True, parents=True)

        # Set plotting style
        plt.style.use(self.style)

        # Configure matplotlib for cleaner plots
        plt.rcParams['font.size'] = 10
        plt.rcParams['axes.labelsize'] = 11
        plt.rcParams['axes.titlesize'] = 12
        plt.rcParams['legend.fontsize'] = 9
        plt.rcParams['figure.titlesize'] = 14

    def load_data(self) -> bool:
        """
        Load all GAPpy output data files.

        Returns:
            True if data loaded successfully, False otherwise
        """
        required_files = {
            "site_data": "site_data.csv",
            "species_data": "species_data.csv",
            "soil_data": "soil_data.csv"
        }

        optional_files = {
            "genus_data": "genus_data.csv"
        }

        print(f"Loading data from: {self.output_dir}")

        try:
            # Load required files
            for attr_name, filename in required_files.items():
                filepath = self.output_dir / filename
                if not filepath.exists():
                    print(f"Error: Required file {filepath} not found")
                    return False

                df = pd.read_csv(filepath)

                # Fix malformed scientific notation (e.g., "0.000000e+" -> 0.0)
                if attr_name in ['species_data', 'genus_data']:
                    numeric_cols = ['max_diam', 'max_hgt', 'leaf_area_ind', 'basal_area',
                                  'total_biomC', 'pl_biomC_std', 'total_biomN', 'pl_biomN_std']
                    for col in numeric_cols:
                        if col in df.columns:
                            df[col] = pd.to_numeric(df[col], errors='coerce').fillna(0.0)

                # Rename columns to match expected names
                if attr_name == 'species_data':
                    column_mapping = {
                        'total_biomC': 'biomass_c',
                        'total_biomN': 'biomass_n',
                        'pl_biomC_std': 'biomass_c_std',
                        'pl_biomN_std': 'biomass_n_std'
                    }
                    df = df.rename(columns=column_mapping)

                setattr(self, attr_name, df)
                print(f"  ✓ {filename}: {len(df)} records")

            # Load optional files
            for attr_name, filename in optional_files.items():
                filepath = self.output_dir / filename
                if filepath.exists():
                    df = pd.read_csv(filepath)

                    # Fix malformed scientific notation (e.g., "0.000000e+" -> 0.0)
                    if attr_name in ['species_data', 'genus_data']:
                        numeric_cols = ['max_diam', 'max_hgt', 'leaf_area_ind', 'basal_area',
                                      'total_biomC', 'pl_biomC_std', 'total_biomN', 'pl_biomN_std']
                        for col in numeric_cols:
                            if col in df.columns:
                                df[col] = pd.to_numeric(df[col], errors='coerce').fillna(0.0)

                    # Rename columns to match expected names
                    if attr_name == 'genus_data':
                        column_mapping = {
                            'total_biomC': 'biomass_c',
                            'total_biomN': 'biomass_n',
                            'pl_biomC_std': 'biomass_c_std',
                            'pl_biomN_std': 'biomass_n_std'
                        }
                        df = df.rename(columns=column_mapping)

                    setattr(self, attr_name, df)
                    print(f"  ✓ {filename}: {len(df)} records")
                else:
                    print(f"  - {filename}: not found (optional)")

            return True

        except Exception as e:
            print(f"Error loading data: {e}")
            return False

    def validate_data(self) -> bool:
        """
        Validate loaded data for required columns and basic consistency.

        Returns:
            True if data is valid, False otherwise
        """
        if self.species_data is None or self.site_data is None or self.soil_data is None:
            print("Error: Required data not loaded")
            return False

        # Check for required columns
        required_columns = {
            'species_data': ['year', 'genus', 'species', 'biomass_c', 'biomass_n', 'basal_area'],
            'site_data': ['year', 'deg_days', 'grow_days', 'rain'],
            'soil_data': ['year', 'a0c0', 'a0n0', 'ac0', 'an0']
        }

        for data_name, columns in required_columns.items():
            df = getattr(self, data_name)
            missing_cols = [col for col in columns if col not in df.columns]
            if missing_cols:
                print(f"Error: Missing columns in {data_name}: {missing_cols}")
                return False

        # Check data consistency
        # Note: Species data is typically written at print_interval (e.g., every 10 years)
        # while site/soil data is written every year, so different counts are expected
        years_species = set(self.species_data['year'].unique())
        years_site = set(self.site_data['year'].unique())
        years_soil = set(self.soil_data['year'].unique())

        # Check if species years are a subset of site/soil years (expected behavior)
        if not years_species.issubset(years_site) or not years_species.issubset(years_soil):
            print("Warning: Species data years not consistent with site/soil data")
            print(f"  Species: {len(years_species)} years")
            print(f"  Site: {len(years_site)} years")
            print(f"  Soil: {len(years_soil)} years")
        elif len(years_species) < len(years_site):
            print(f"Note: Species data written at intervals ({len(years_species)} years), "
                  f"site/soil data written annually ({len(years_site)} years)")

        return True

    def prepare_species_summary(self) -> pd.DataFrame:
        """
        Prepare aggregated species data for plotting.

        Returns:
            DataFrame with species aggregated by year and genus
        """
        if self.species_data is None:
            raise ValueError("Species data not loaded")

        # Use genus if available, otherwise use species
        genus_col = 'genus' if 'genus' in self.species_data.columns else 'species'

        # Calculate tree counts if not present
        if 'n_trees' not in self.species_data.columns:
            # Estimate from basal area or use constant
            self.species_data['n_trees'] = np.maximum(1,
                self.species_data['basal_area'] / 100)  # Rough estimate

        summary = self.species_data.groupby(['year', genus_col]).agg({
            'biomass_c': 'sum',
            'biomass_n': 'sum',
            'basal_area': 'sum',
            'n_trees': 'sum'
        }).reset_index()

        # Add derived metrics
        summary['cn_ratio'] = summary['biomass_c'] / summary['biomass_n'].replace(0, np.nan)
        summary['biomass_per_tree'] = summary['biomass_c'] / summary['n_trees'].replace(0, np.nan)

        return summary

    def plot_forest_dynamics(self, figsize: Tuple[int, int] = (14, 10)) -> plt.Figure:
        """Plot forest composition and biomass over time."""
        fig, axes = plt.subplots(2, 2, figsize=figsize)
        fig.suptitle('Forest Dynamics Over Time', fontsize=14, fontweight='bold')

        species_summary = self.prepare_species_summary()
        genus_col = 'genus' if 'genus' in species_summary.columns else 'species'

        # Get unique genera/species for consistent colors
        unique_taxa = species_summary[genus_col].unique()
        color_map = dict(zip(unique_taxa, SPECIES_COLORS[:len(unique_taxa)]))

        # Plot 1: Stacked Area Plot for Biomass Carbon by Species
        ax1 = axes[0, 0]

        # Pivot data for stacked area plot
        biomass_pivot = species_summary.pivot_table(
            index='year',
            columns=genus_col,
            values='biomass_c',
            fill_value=0
        )

        # Create stacked area plot
        ax1.stackplot(biomass_pivot.index,
                     *[biomass_pivot[taxon] for taxon in biomass_pivot.columns],
                     labels=biomass_pivot.columns,
                     colors=[color_map[taxon] for taxon in biomass_pivot.columns],
                     alpha=0.7)

        ax1.set_xlabel('Year')
        ax1.set_ylabel('Biomass Carbon (kg/m²)')
        ax1.set_title('Biomass Carbon by Species')
        ax1.legend(bbox_to_anchor=(1.05, 1), loc='upper left', fontsize=8, frameon=False)
        ax1.grid(True, alpha=0.2)

        # Plot 2: Tree Population by Species
        ax2 = axes[0, 1]
        for idx, taxon in enumerate(unique_taxa):
            subset = species_summary[species_summary[genus_col] == taxon]
            linestyle = LINE_STYLES[idx % len(LINE_STYLES)]
            marker = MARKERS[idx % len(MARKERS)] if len(unique_taxa) <= 10 else None
            ax2.plot(subset['year'], subset['n_trees'],
                    linewidth=2.5, label=taxon, color=color_map[taxon],
                    linestyle=linestyle, marker=marker, markersize=4,
                    markevery=max(1, len(subset) // 10), alpha=0.9)

        ax2.set_xlabel('Year')
        ax2.set_ylabel('Number of Trees')
        ax2.set_title('Tree Population by Species')
        ax2.legend(bbox_to_anchor=(1.05, 1), loc='upper left', frameon=False, fontsize=8)
        ax2.grid(True, alpha=0.2)

        # Plot 3: Basal Area by Species
        ax3 = axes[1, 0]
        for idx, taxon in enumerate(unique_taxa):
            subset = species_summary[species_summary[genus_col] == taxon]
            linestyle = LINE_STYLES[idx % len(LINE_STYLES)]
            marker = MARKERS[idx % len(MARKERS)] if len(unique_taxa) <= 10 else None
            ax3.plot(subset['year'], subset['basal_area'],
                    linewidth=2.5, label=taxon, color=color_map[taxon],
                    linestyle=linestyle, marker=marker, markersize=4,
                    markevery=max(1, len(subset) // 10), alpha=0.9)

        ax3.set_xlabel('Year')
        ax3.set_ylabel('Basal Area (cm²)')
        ax3.set_title('Basal Area by Species')
        ax3.legend(bbox_to_anchor=(1.05, 1), loc='upper left', frameon=False, fontsize=8)
        ax3.grid(True, alpha=0.2)

        # Plot 4: Species Diversity (Shannon Index)
        ax4 = axes[1, 1]
        diversity_data = []

        for year in species_summary['year'].unique():
            year_data = species_summary[species_summary['year'] == year]
            total_biomass = year_data['biomass_c'].sum()

            if total_biomass > 0:
                proportions = year_data['biomass_c'] / total_biomass
                # Shannon diversity index
                shannon = -np.sum(proportions * np.log(proportions + 1e-10))
                diversity_data.append({'year': year, 'shannon': shannon})

        if diversity_data:
            div_df = pd.DataFrame(diversity_data)
            ax4.plot(div_df['year'], div_df['shannon'],
                    linewidth=2, color='#7570B3')
            ax4.set_xlabel('Year')
            ax4.set_ylabel('Shannon Diversity Index')
            ax4.set_title('Species Diversity')
            ax4.grid(True, alpha=0.2)

        plt.tight_layout()
        return fig

    def plot_soil_biogeochemistry(self, figsize: Tuple[int, int] = (14, 10)) -> plt.Figure:
        """Plot soil carbon and nitrogen dynamics."""
        fig, axes = plt.subplots(2, 2, figsize=figsize)
        fig.suptitle('Soil Biogeochemistry Over Time', fontsize=14, fontweight='bold')

        # Plot 1: Soil Carbon Pools
        ax1 = axes[0, 0]
        ax1.plot(self.soil_data['year'], self.soil_data['a0c0'],
                linewidth=3, label='Surface C (A0)',
                color=SOIL_COLORS['surface'], linestyle='-', marker='o',
                markevery=max(1, len(self.soil_data) // 20), markersize=5, alpha=0.9)
        ax1.plot(self.soil_data['year'], self.soil_data['ac0'],
                linewidth=3, label='Mineral C (A)',
                color=SOIL_COLORS['mineral'], linestyle='--', marker='s',
                markevery=max(1, len(self.soil_data) // 20), markersize=5, alpha=0.9)

        # Add total soil carbon
        total_c = self.soil_data['a0c0'] + self.soil_data['ac0']
        ax1.plot(self.soil_data['year'], total_c,
                linewidth=3.5, label='Total C',
                color='#000000', linestyle='-.', marker='D',
                markevery=max(1, len(self.soil_data) // 20), markersize=5, alpha=0.9)

        ax1.set_xlabel('Year')
        ax1.set_ylabel('Carbon Pool (kg C/m²)')
        ax1.set_title('Soil Carbon Pools')
        ax1.legend(frameon=False)
        ax1.grid(True, alpha=0.2)

        # Plot 2: Soil Nitrogen Pools
        ax2 = axes[0, 1]
        ax2.plot(self.soil_data['year'], self.soil_data['a0n0'],
                linewidth=3, label='Surface N (A0)',
                color='#1E90FF', linestyle='-', marker='o',
                markevery=max(1, len(self.soil_data) // 20), markersize=5, alpha=0.9)
        ax2.plot(self.soil_data['year'], self.soil_data['an0'],
                linewidth=3, label='Mineral N (A)',
                color='#4169E1', linestyle='--', marker='s',
                markevery=max(1, len(self.soil_data) // 20), markersize=5, alpha=0.9)

        # Add total soil nitrogen
        total_n = self.soil_data['a0n0'] + self.soil_data['an0']
        ax2.plot(self.soil_data['year'], total_n,
                linewidth=3.5, label='Total N',
                color='#00008B', linestyle='-.', marker='D',
                markevery=max(1, len(self.soil_data) // 20), markersize=5, alpha=0.9)

        ax2.set_xlabel('Year')
        ax2.set_ylabel('Nitrogen Pool (kg N/m²)')
        ax2.set_title('Soil Nitrogen Pools')
        ax2.legend(frameon=False)
        ax2.grid(True, alpha=0.2)

        # Plot 3: Available Nitrogen (if present)
        ax3 = axes[1, 0]
        if 'avail_n' in self.soil_data.columns:
            ax3.plot(self.soil_data['year'], self.soil_data['avail_n'] * 1000,
                    linewidth=2, color=SOIL_COLORS['available'])
            ax3.set_ylabel('Available N (g N/m²)')
        else:
            # Estimate available N from total N
            estimated_avail = total_n * 0.05  # Assume 5% is available
            ax3.plot(self.soil_data['year'], estimated_avail * 1000,
                    linewidth=2, color=SOIL_COLORS['available'],
                    linestyle=':', alpha=0.6)
            ax3.set_ylabel('Estimated Avail. N (g N/m²)')

        ax3.set_xlabel('Year')
        ax3.set_title('Available Nitrogen')
        ax3.grid(True, alpha=0.2)

        # Plot 4: Soil C:N Ratios
        ax4 = axes[1, 1]
        # Handle division by zero
        cn_ratio_A0 = self.soil_data['a0c0'] / self.soil_data['a0n0'].replace(0, np.nan)
        cn_ratio_A = self.soil_data['ac0'] / self.soil_data['an0'].replace(0, np.nan)
        cn_ratio_total = total_c / total_n.replace(0, np.nan)

        ax4.plot(self.soil_data['year'], cn_ratio_A0,
                linewidth=3, label='Surface C:N', color='#FF8C00',
                linestyle='-', marker='o', markevery=max(1, len(self.soil_data) // 20),
                markersize=5, alpha=0.9)
        ax4.plot(self.soil_data['year'], cn_ratio_A,
                linewidth=3, label='Mineral C:N', color='#DC143C',
                linestyle='--', marker='s', markevery=max(1, len(self.soil_data) // 20),
                markersize=5, alpha=0.9)
        ax4.plot(self.soil_data['year'], cn_ratio_total,
                linewidth=3.5, label='Total C:N', color='#8B0000',
                linestyle='-.', marker='D', markevery=max(1, len(self.soil_data) // 20),
                markersize=5, alpha=0.9)

        ax4.set_xlabel('Year')
        ax4.set_ylabel('C:N Ratio')
        ax4.set_title('Soil C:N Ratios')
        ax4.legend(frameon=False)
        ax4.grid(True, alpha=0.2)

        plt.tight_layout()
        return fig

    def plot_environmental_conditions(self, figsize: Tuple[int, int] = (14, 10)) -> plt.Figure:
        """Plot environmental and climate conditions."""
        fig, axes = plt.subplots(2, 2, figsize=figsize)
        fig.suptitle('Environmental Conditions Over Time', fontsize=14, fontweight='bold')

        # Plot 1: Temperature-related variables
        ax1 = axes[0, 0]
        ax1.plot(self.site_data['year'], self.site_data['deg_days'],
                linewidth=2, color=CLIMATE_COLORS['temperature'],
                label='Degree Days')

        if 'grow_days' in self.site_data.columns:
            ax1_twin = ax1.twinx()
            ax1_twin.plot(self.site_data['year'], self.site_data['grow_days'],
                         linewidth=2, color='#2CA02C', label='Growing Days')
            ax1_twin.set_ylabel('Growing Days', color='#2CA02C')
            ax1_twin.tick_params(axis='y', labelcolor='#2CA02C')

        ax1.set_xlabel('Year')
        ax1.set_ylabel('Degree Days', color=CLIMATE_COLORS['temperature'])
        ax1.set_title('Temperature Conditions')
        ax1.tick_params(axis='y', labelcolor=CLIMATE_COLORS['temperature'])
        ax1.grid(True, alpha=0.2)

        # Plot 2: Water balance
        ax2 = axes[0, 1]
        ax2.plot(self.site_data['year'], self.site_data['rain'],
                linewidth=2, label='Rainfall',
                color=CLIMATE_COLORS['water'])

        # Add other water variables if present
        water_vars = ['pot_evap', 'act_evap', 'runoff']
        colors = ['#FFA500', '#E74C3C', '#9B59B6']

        for var, color in zip(water_vars, colors):
            if var in self.site_data.columns:
                ax2.plot(self.site_data['year'], self.site_data[var],
                        linewidth=2, label=var.replace('_', ' ').title(),
                        color=color)

        ax2.set_xlabel('Year')
        ax2.set_ylabel('Water Flux (m/year)')
        ax2.set_title('Water Balance')
        ax2.legend(frameon=False)
        ax2.grid(True, alpha=0.2)

        # Plot 3: Drought stress indicators
        ax3 = axes[1, 0]
        drought_vars = ['dry_days_upper', 'dry_days_base', 'drought_days']
        drought_colors = ['#DC143C', '#FF8C00', '#FFD700']
        drought_styles = ['-', '--', '-.']
        drought_markers = ['o', 's', '^']
        drought_found = False

        for idx, var in enumerate(drought_vars):
            if var in self.site_data.columns:
                ax3.plot(self.site_data['year'], self.site_data[var],
                        linewidth=3, label=var.replace('_', ' ').title(),
                        color=drought_colors[idx], linestyle=drought_styles[idx],
                        marker=drought_markers[idx], markevery=max(1, len(self.site_data) // 20),
                        markersize=5, alpha=0.9)
                drought_found = True

        if not drought_found:
            # Create estimated drought stress from temperature and precipitation
            drought_estimate = np.maximum(0, self.site_data['deg_days'] / 50 - self.site_data['rain'] * 10)
            ax3.plot(self.site_data['year'], drought_estimate,
                    linewidth=2, label='Estimated Drought Stress',
                    color=CLIMATE_COLORS['stress'], alpha=0.6)

        ax3.set_xlabel('Year')
        ax3.set_ylabel('Stress Indicator')
        ax3.set_title('Drought Stress')
        ax3.legend(frameon=False)
        ax3.grid(True, alpha=0.2)

        # Plot 4: Additional environmental variables
        ax4 = axes[1, 1]
        env_vars = ['flood_days', 'wind_days', 'freeze_days']
        env_colors = ['#1E90FF', '#9370DB', '#708090']
        env_styles = ['-', '--', '-.']
        env_markers = ['o', 's', 'D']
        env_found = False

        for idx, var in enumerate(env_vars):
            if var in self.site_data.columns:
                ax4.plot(self.site_data['year'], self.site_data[var],
                        linewidth=3, label=var.replace('_', ' ').title(),
                        color=env_colors[idx], linestyle=env_styles[idx],
                        marker=env_markers[idx], markevery=max(1, len(self.site_data) // 20),
                        markersize=5, alpha=0.9)
                env_found = True

        if not env_found:
            ax4.text(0.5, 0.5, 'No additional\nenvironmental\nvariables found',
                    ha='center', va='center', transform=ax4.transAxes,
                    fontsize=10, color='#666666')

        ax4.set_xlabel('Year')
        ax4.set_ylabel('Days')
        ax4.set_title('Other Environmental Stress')
        if env_found:
            ax4.legend(frameon=False)
        ax4.grid(True, alpha=0.2)

        plt.tight_layout()
        return fig

    def create_summary_dashboard(self, figsize: Tuple[int, int] = (16, 12)) -> plt.Figure:
        """Create a comprehensive dashboard view."""
        fig = plt.figure(figsize=figsize)
        fig.suptitle('GAPpy Model Output Dashboard', fontsize=14, fontweight='bold')

        # Create complex grid layout
        gs = gridspec.GridSpec(3, 4, hspace=0.3, wspace=0.3, figure=fig)

        species_summary = self.prepare_species_summary()
        genus_col = 'genus' if 'genus' in species_summary.columns else 'species'

        # Total ecosystem biomass
        ax1 = fig.add_subplot(gs[0, :2])
        total_biomass = species_summary.groupby('year')['biomass_c'].sum()
        ax1.plot(total_biomass.index, total_biomass.values,
                linewidth=2, color='#2E7D32')
        ax1.fill_between(total_biomass.index, total_biomass.values, alpha=0.2, color='#4CAF50')
        ax1.set_xlabel('Year')
        ax1.set_ylabel('Total Biomass Carbon (kg/m²)')
        ax1.set_title('Total Ecosystem Biomass', fontweight='bold')
        ax1.grid(True, alpha=0.2)

        # Species composition pie chart (final year)
        ax2 = fig.add_subplot(gs[0, 2])
        final_year = species_summary['year'].max()
        final_biomass = species_summary[species_summary['year'] == final_year].groupby(genus_col)['biomass_c'].sum()
        final_biomass = final_biomass[final_biomass > 0]

        if len(final_biomass) > 0:
            # Limit to top 8 species for readability
            if len(final_biomass) > 8:
                top_species = final_biomass.nlargest(7)
                other_biomass = final_biomass.iloc[7:].sum()
                if other_biomass > 0:
                    top_species['Others'] = other_biomass
                final_biomass = top_species

            wedges, texts, autotexts = ax2.pie(final_biomass.values, labels=final_biomass.index,
                                              autopct='%1.1f%%', startangle=90)
            # Improve text readability
            for autotext in autotexts:
                autotext.set_color('white')
                autotext.set_fontweight('bold')
            ax2.set_title(f'Species Composition\n(Year {final_year})', fontweight='bold')
        else:
            ax2.text(0.5, 0.5, 'No Living\nBiomass', ha='center', va='center',
                    transform=ax2.transAxes, fontsize=12, fontweight='bold')
            ax2.set_title(f'Species Composition\n(Year {final_year})', fontweight='bold')

        # Total tree count
        ax3 = fig.add_subplot(gs[0, 3])
        total_trees = species_summary.groupby('year')['n_trees'].sum()
        ax3.plot(total_trees.index, total_trees.values,
                linewidth=2, color='#7B1FA2')
        ax3.fill_between(total_trees.index, total_trees.values, alpha=0.2, color='#9C27B0')
        ax3.set_xlabel('Year')
        ax3.set_ylabel('Total Trees')
        ax3.set_title('Tree Population', fontweight='bold')
        ax3.grid(True, alpha=0.2)

        # Soil carbon accumulation
        ax4 = fig.add_subplot(gs[1, :2])
        total_soil_c = self.soil_data['a0c0'] + self.soil_data['ac0']
        ax4.plot(self.soil_data['year'], total_soil_c,
                linewidth=2, color='#6D4C41')
        ax4.fill_between(self.soil_data['year'], total_soil_c, alpha=0.2, color='#8D6E63')
        ax4.set_xlabel('Year')
        ax4.set_ylabel('Total Soil Carbon (kg C/m²)')
        ax4.set_title('Soil Carbon Accumulation', fontweight='bold')
        ax4.grid(True, alpha=0.2)

        # Environmental stress indicators
        ax5 = fig.add_subplot(gs[1, 2:])
        if 'dry_days_upper' in self.site_data.columns:
            ax5.plot(self.site_data['year'], self.site_data['dry_days_upper'],
                    linewidth=3, label='Drought Days', color='#DC143C',
                    linestyle='-', marker='o', markevery=max(1, len(self.site_data) // 20),
                    markersize=5, alpha=0.9)
        if 'flood_days' in self.site_data.columns:
            ax5.plot(self.site_data['year'], self.site_data['flood_days'],
                    linewidth=3, label='Flood Days', color='#1E90FF',
                    linestyle='--', marker='s', markevery=max(1, len(self.site_data) // 20),
                    markersize=5, alpha=0.9)

        # If no stress data, create estimate
        if 'dry_days_upper' not in self.site_data.columns and 'flood_days' not in self.site_data.columns:
            stress_estimate = self.site_data['deg_days'] / 100
            ax5.plot(self.site_data['year'], stress_estimate,
                    linewidth=2, label='Temperature Stress',
                    color=CLIMATE_COLORS['stress'])

        ax5.set_xlabel('Year')
        ax5.set_ylabel('Stress Indicator')
        ax5.set_title('Environmental Stress', fontweight='bold')
        ax5.legend(frameon=False)
        ax5.grid(True, alpha=0.2)

        # Forest productivity (biomass per tree)
        ax6 = fig.add_subplot(gs[2, :2])
        productivity = total_biomass / total_trees.replace(0, np.nan)
        productivity = productivity.dropna()
        if len(productivity) > 0:
            ax6.plot(productivity.index, productivity.values,
                    linewidth=2, color='#F57C00')
            ax6.fill_between(productivity.index, productivity.values, alpha=0.2, color='#FF9800')
        ax6.set_xlabel('Year')
        ax6.set_ylabel('Biomass per Tree (kg C/tree)')
        ax6.set_title('Forest Productivity', fontweight='bold')
        ax6.grid(True, alpha=0.2)

        # Ecosystem C:N ratio
        ax7 = fig.add_subplot(gs[2, 2:])
        total_biomass_n = species_summary.groupby('year')['biomass_n'].sum()
        ecosystem_cn = total_biomass / total_biomass_n.replace(0, np.nan)
        ecosystem_cn = ecosystem_cn.dropna()
        if len(ecosystem_cn) > 0:
            ax7.plot(ecosystem_cn.index, ecosystem_cn.values,
                    linewidth=2, color='#C62828')
            ax7.fill_between(ecosystem_cn.index, ecosystem_cn.values, alpha=0.2, color='#E57373')
        ax7.set_xlabel('Year')
        ax7.set_ylabel('C:N Ratio')
        ax7.set_title('Ecosystem C:N Ratio', fontweight='bold')
        ax7.grid(True, alpha=0.2)

        return fig

    def save_plots(self, figures: List[plt.Figure],
                  format: str = DEFAULT_FORMAT,
                  dpi: int = DEFAULT_DPI) -> None:
        """Save all plots to files."""
        plot_names = [
            "forest_dynamics",
            "soil_biogeochemistry",
            "environmental_conditions",
            "summary_dashboard"
        ]

        print(f"\nSaving plots to: {self.plots_dir}")

        for fig, name in zip(figures, plot_names):
            filepath = self.plots_dir / f"{name}.{format}"
            fig.savefig(filepath, dpi=dpi, bbox_inches='tight',
                       facecolor='white', edgecolor='none')
            print(f"  ✓ {filepath.name}")

    def print_summary_statistics(self) -> None:
        """Display comprehensive summary statistics."""
        print("\n" + "=" * 60)
        print("GAPpy SIMULATION SUMMARY")
        print("=" * 60)

        # Basic simulation info
        years = self.site_data['year'].max() - self.site_data['year'].min() + 1
        print(f"Simulation period: {years} years ({self.site_data['year'].min()}-{self.site_data['year'].max()})")

        # Forest metrics
        species_summary = self.prepare_species_summary()
        final_year_data = species_summary[species_summary['year'] == species_summary['year'].max()]

        final_biomass = final_year_data['biomass_c'].sum()
        final_biomass_n = final_year_data['biomass_n'].sum()
        final_trees = final_year_data['n_trees'].sum()
        n_species = len(final_year_data)

        print(f"\nFOREST METRICS (Final Year):")
        print(f"  Total biomass C: {final_biomass:.2f} kg C/m²")
        print(f"  Total biomass N: {final_biomass_n:.3f} kg N/m²")
        print(f"  Total trees: {final_trees:.0f}")
        print(f"  Number of species: {n_species}")
        if final_trees > 0:
            print(f"  Average biomass per tree: {final_biomass/final_trees:.3f} kg C/tree")
        if final_biomass_n > 0:
            print(f"  Forest C:N ratio: {final_biomass/final_biomass_n:.1f}")

        # Soil metrics
        final_soil = self.soil_data[self.soil_data['year'] == self.soil_data['year'].max()].iloc[0]
        total_soil_c = final_soil['a0c0'] + final_soil['ac0']
        total_soil_n = final_soil['a0n0'] + final_soil['an0']

        print(f"\nSOIL METRICS (Final Year):")
        print(f"  Total soil C: {total_soil_c:.2f} kg C/m²")
        print(f"  Total soil N: {total_soil_n:.3f} kg N/m²")
        if total_soil_n > 0:
            print(f"  Soil C:N ratio: {total_soil_c/total_soil_n:.1f}")

        # Climate metrics
        final_site = self.site_data[self.site_data['year'] == self.site_data['year'].max()].iloc[0]
        print(f"\nCLIMATE METRICS (Final Year):")
        print(f"  Degree days: {final_site['deg_days']:.0f}")
        print(f"  Annual rainfall: {final_site['rain']:.2f} m")

        # Change over time
        initial_biomass = species_summary[species_summary['year'] == species_summary['year'].min()]['biomass_c'].sum()
        initial_soil_c = self.soil_data[self.soil_data['year'] == self.soil_data['year'].min()].iloc[0]['a0c0'] + \
                        self.soil_data[self.soil_data['year'] == self.soil_data['year'].min()].iloc[0]['ac0']

        print(f"\nCHANGE OVER TIME:")
        print(f"  Forest biomass change: {final_biomass - initial_biomass:+.2f} kg C/m²")
        print(f"  Soil carbon change: {total_soil_c - initial_soil_c:+.2f} kg C/m²")
        print(f"  Total ecosystem C change: {(final_biomass + total_soil_c) - (initial_biomass + initial_soil_c):+.2f} kg C/m²")

        print("=" * 60)

    def create_all_plots(self, save_plots: bool = True,
                        show_plots: bool = True,
                        format: str = DEFAULT_FORMAT) -> List[plt.Figure]:
        """Create all visualization plots."""
        if not self.load_data():
            return []

        if not self.validate_data():
            return []

        print("\nCreating visualizations...")
        figures = []

        # Create plots
        try:
            fig1 = self.plot_forest_dynamics()
            figures.append(fig1)
            print("  ✓ Forest dynamics")

            fig2 = self.plot_soil_biogeochemistry()
            figures.append(fig2)
            print("  ✓ Soil biogeochemistry")

            fig3 = self.plot_environmental_conditions()
            figures.append(fig3)
            print("  ✓ Environmental conditions")

            fig4 = self.create_summary_dashboard()
            figures.append(fig4)
            print("  ✓ Summary dashboard")

        except Exception as e:
            import traceback
            print(f"Error creating plots: {e}")
            print("\nFull traceback:")
            traceback.print_exc()
            return figures

        # Save plots
        if save_plots and figures:
            self.save_plots(figures, format=format)

        # Print summary statistics
        self.print_summary_statistics()

        # Show plots
        if show_plots:
            plt.show()

        return figures


def main():
    """Main function with command line interface."""
    parser = argparse.ArgumentParser(
        description="Visualize GAPpy forest model outputs",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python plot_outputs.py                           # Use default settings
  python plot_outputs.py --output-dir data/       # Custom output directory
  python plot_outputs.py --format pdf --no-show   # Save as PDF, don't display
  python plot_outputs.py --plots-dir figures/     # Custom plots directory
        """
    )

    parser.add_argument('--output-dir', '-o',
                       default=DEFAULT_OUTPUT_DIR,
                       help=f'Directory containing GAPpy output files (default: {DEFAULT_OUTPUT_DIR})')

    parser.add_argument('--plots-dir', '-p',
                       default=DEFAULT_PLOTS_DIR,
                       help=f'Directory to save plot files (default: {DEFAULT_PLOTS_DIR})')

    parser.add_argument('--format', '-f',
                       choices=['png', 'pdf', 'svg', 'jpg'],
                       default=DEFAULT_FORMAT,
                       help=f'Output format for plots (default: {DEFAULT_FORMAT})')

    parser.add_argument('--dpi', '-d',
                       type=int, default=DEFAULT_DPI,
                       help=f'DPI for saved plots (default: {DEFAULT_DPI})')

    parser.add_argument('--style', '-s',
                       default=DEFAULT_STYLE,
                       help=f'Matplotlib style (default: {DEFAULT_STYLE})')

    show_group = parser.add_mutually_exclusive_group()
    show_group.add_argument('--show', action='store_true', default=True,
                           help='Display plots interactively (default)')
    show_group.add_argument('--no-show', dest='show', action='store_false',
                           help='Do not display plots interactively')

    args = parser.parse_args()

    print("GAPpy Output Visualization")
    print("=" * 50)
    print(f"Output directory: {args.output_dir}")
    print(f"Plots directory: {args.plots_dir}")
    print(f"Format: {args.format}")

    try:
        # Create plotter and generate visualizations
        plotter = GAPpyPlotter(
            output_dir=args.output_dir,
            plots_dir=args.plots_dir,
            style=args.style
        )

        figures = plotter.create_all_plots(
            save_plots=True,
            show_plots=args.show,
            format=args.format
        )

        if figures:
            print(f"\n✓ Successfully created {len(figures)} visualizations")
            print(f"✓ Plots saved to: {plotter.plots_dir}")
        else:
            print("\n✗ No plots were created")
            sys.exit(1)

    except KeyboardInterrupt:
        print("\n\nVisualization cancelled by user")
        sys.exit(1)
    except Exception as e:
        print(f"\nError: {e}")
        sys.exit(1)

    print("\nVisualization complete!")


if __name__ == "__main__":
    main()