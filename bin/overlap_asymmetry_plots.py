"""
Step 3: Visualization Functions for Overlapping Frame Asymmetry

This module provides plotting functions to visualize asymmetric evolution
in overlapping reading frames, similar to Szklarczyk et al. 2007 Figure 5.

Author: Analysis for mitochondrial microproteins
"""

import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import pandas as pd
from typing import Dict, List
import warnings
from adjustText import adjust_text

# Set plotting style
sns.set_style("whitegrid")
sns.set_context("notebook", font_scale=1.1)


def plot_asymmetry_scatter(results_df: pd.DataFrame,
                           output_file: str = None,
                           figsize: tuple = (10, 10),
                           title = 'Asymmetric Evolution in Overlapping Reading Frames\n' +
                'Points above diagonal: Alternative frame evolving faster',
                           include_labels: bool = True, 
                           new_labels: Dict[str, str] = None,
                ):
    """
    Create scatter plot of dN values for canonical vs alternative frames.
    
    This recreates the style of Szklarczyk et al. 2007 Figure 5A, showing
    asymmetric evolution when points deviate from the diagonal.
    
    Parameters:
    -----------
    results_df : pd.DataFrame
        Results from compare_multiple_overlaps()
    output_file : str, optional
        Path to save figure
    figsize : tuple
        Figure size (width, height)
    title : str
        Title of the plot
    include_labels : bool
        Whether to include labels for each point
    new_labels : Dict[str, str]
        Dictionary of new labels for points
    """
    fig, ax = plt.subplots(figsize=figsize)
    
    # Extract data
    dn_canonical = results_df['dN_canonical'].values
    dn_alternative = results_df['dN_alternative'].values
    significant = results_df['Significant_Asymmetry'].values
    
    # Color code by significance
    colors = ['red' if sig else 'gray' for sig in significant]
    
    # Create scatter plot
    scatter = ax.scatter(dn_canonical, dn_alternative,
                        c=colors, s=100, alpha=0.6, edgecolors='black', linewidth=1.5)

    # Add diagonal line (y=x) representing symmetric evolution
    max_val = max(dn_canonical.max(), dn_alternative.max()) * 1.1
    diagonal_line = ax.plot([0, max_val], [0, max_val], 'k--', alpha=0.5, linewidth=2,
                            label='Symmetric evolution (dN_alt = dN_can)')

    if include_labels:
        labels = results_df['Region'].values
        if new_labels:
            labels = [new_labels.get(label, label) for label in labels]

        # Add labels for each point - collect text objects for auto-adjustment
        # Calculate base offset (percentage of axis range)
        x_range = dn_canonical.max() - dn_canonical.min()
        y_range = dn_alternative.max() - dn_alternative.min()
        base_offset_x = x_range * 0.09  # 2% of x range
        base_offset_y = y_range * -0.09  # 2% of y range

        texts = []
        for idx, row in results_df.iterrows():
            # Apply base offset so labels don't start on the points
            text = ax.text(row['dN_canonical'] + base_offset_x,
                          row['dN_alternative'] + base_offset_y,
                          labels[idx],
                          fontsize=10, alpha=1, fontweight='bold')
            texts.append(text)

        # Automatically adjust text positions to avoid overlaps with labels, points, and the diagonal line
        adjust_text(texts,
                   x=dn_canonical, y=dn_alternative,  # Avoid scatter points
                   ax=ax,
                   expand_points=(2.0, 2.0),    # Expand point boundaries
                   expand_text=(0.2, 0.2),      # Expand text boundaries
                   force_points=(0.2, 0.2),     # Force pushing away from points
                   force_text=(0.2, 0.2),       # Force pushing labels apart
                   add_objects=diagonal_line,   # Avoid the diagonal line
                   lim=500,                     # Increase iteration limit
                   autoalign='xy',              # Auto-align labels
                   arrowprops=dict(arrowstyle='-', color='gray', lw=1, alpha=0))
    
    # Formatting
    ax.set_xlabel('dN Canonical Frame', fontsize=12, fontweight='bold')
    ax.set_ylabel('dN Alternative Frame', fontsize=12, fontweight='bold')
    ax.set_title(title,
                fontsize=14, fontweight='bold', pad=20)
    
    # Add legend
    from matplotlib.patches import Patch
    legend_elements = [
        Patch(facecolor='red', alpha=0.6, label='Significant asymmetry (p<0.05)'),
        Patch(facecolor='gray', alpha=0.6, label='Not significant'),
        plt.Line2D([0], [0], color='k', linestyle='--', linewidth=2, 
                   label='Symmetric evolution')
    ]
    # Legend outside
    ax.legend(handles=legend_elements, loc='upper left', fontsize=10, bbox_to_anchor=(1, 1))    
    # Set equal aspect ratio
    ax.set_aspect('equal', adjustable='box')
    
    # Grid
    ax.grid(True, alpha=0.3)
    sns.despine()
    plt.tight_layout()
    
    if output_file:
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        print(f"Scatter plot saved to {output_file}")
    
    return fig, ax


def plot_omega_comparison(results_df: pd.DataFrame,
                         output_file: str = None,
                         figsize: tuple = (12, 6)):
    """
    Create bar plot comparing omega (dN/dS) values between frames.
    
    Parameters:
    -----------
    results_df : pd.DataFrame
        Results from compare_multiple_overlaps()
    output_file : str, optional
        Path to save figure
    figsize : tuple
        Figure size (width, height)
    """
    fig, ax = plt.subplots(figsize=figsize)
    
    # Prepare data
    regions = results_df['Region'].values
    omega_can = results_df['omega_canonical'].values
    omega_alt = results_df['omega_alternative'].values
    significant = results_df['Significant_Asymmetry'].values
    
    # Set up bar positions
    x = np.arange(len(regions))
    width = 0.35
    
    # Create bars
    bars1 = ax.bar(x - width/2, omega_can, width, label='Canonical Frame',
                   color='steelblue', alpha=0.8, edgecolor='black')
    bars2 = ax.bar(x + width/2, omega_alt, width, label='Alternative Frame',
                   color='coral', alpha=0.8, edgecolor='black')
    
    # Mark significant asymmetries with asterisks
    for i, sig in enumerate(significant):
        if sig:
            max_height = max(omega_can[i], omega_alt[i])
            ax.text(i, max_height * 1.1, '*', ha='center', va='bottom',
                   fontsize=12, fontweight='bold', color='red')
    
    # Add horizontal line at omega=1 (neutral evolution)
    ax.axhline(y=1, color='black', linestyle='--', linewidth=2, alpha=0.5,
              label='Neutral evolution (ω=1)')
    
    # Formatting
    ax.set_xlabel('Overlapping Region', fontsize=12, fontweight='bold')
    ax.set_ylabel('ω (dN/dS)', fontsize=12, fontweight='bold')
    ax.set_title('Selection Pressure Comparison in Overlapping Frames\n' +
                'ω < 1: Purifying selection | ω > 1: Positive/relaxed selection',
                fontsize=14, fontweight='bold', pad=20)
    ax.set_xticks(x)
    ax.set_xticklabels(regions, rotation=45, ha='center')
    ax.legend(fontsize=10, loc='upper left')
    
    # Add note about significance
    ax.text(0.02, 0.98, '* = Significant asymmetry (p<0.05)',
           transform=ax.transAxes, fontsize=10, verticalalignment='top',
           bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
    
    # Grid
    ax.grid(True, alpha=0.3, axis='y')
    
    plt.tight_layout()
    
    if output_file:
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        print(f"Omega comparison plot saved to {output_file}")
    
    return fig, ax


def plot_asymmetry_ratios(results_df: pd.DataFrame,
                          output_file: str = None,
                          figsize: tuple = (10, 6)):
    """
    Create bar plot of asymmetry ratios (dN_alt / dN_can).
    
    Ratios > 1 indicate alternative frame evolving faster.
    Ratios < 1 indicate canonical frame evolving faster.
    
    Parameters:
    -----------
    results_df : pd.DataFrame
        Results from compare_multiple_overlaps()
    output_file : str, optional
        Path to save figure
    figsize : tuple
        Figure size (width, height)
    """
    fig, ax = plt.subplots(figsize=figsize)
    
    # Prepare data
    regions = results_df['Region'].values
    ratios = results_df['Asymmetry_Ratio_dN'].values
    significant = results_df['Significant_Asymmetry'].values
    
    # Color code by significance
    colors = ['red' if sig else 'gray' for sig in significant]
    
    # Create bar plot
    bars = ax.bar(range(len(regions)), ratios, color=colors, 
                  alpha=0.7, edgecolor='black', linewidth=1.5)
    
    # Add horizontal line at ratio=1 (symmetric evolution)
    ax.axhline(y=1, color='black', linestyle='--', linewidth=2, alpha=0.7,
              label='Symmetric evolution (ratio=1)')
    
    # Formatting
    ax.set_xlabel('Overlapping Region', fontsize=12, fontweight='bold')
    ax.set_ylabel('Asymmetry Ratio (dN_alt / dN_can)', fontsize=12, fontweight='bold')
    ax.set_title('Asymmetric Evolution Ratios\n' +
                'Ratio > 1: Alternative frame evolving faster',
                fontsize=14, fontweight='bold', pad=20)
    ax.set_xticks(range(len(regions)))
    ax.set_xticklabels(regions, rotation=45, ha='center')
    
    # Add legend
    from matplotlib.patches import Patch
    legend_elements = [
        Patch(facecolor='red', alpha=0.7, label='Significant (p<0.05)'),
        Patch(facecolor='gray', alpha=0.7, label='Not significant'),
        plt.Line2D([0], [0], color='k', linestyle='--', linewidth=2,
                   label='Symmetric (ratio=1)')
    ]
    ax.legend(handles=legend_elements, loc='upper left', fontsize=10)
    
    # Grid
    ax.grid(True, alpha=0.3, axis='y')
    
    # Add values on bars
    for i, (bar, ratio) in enumerate(zip(bars, ratios)):
        height = bar.get_height()
        ax.text(bar.get_x() + bar.get_width()/2., height,
               f'{ratio:.2f}',
               ha='center', va='bottom' if height > 1 else 'top',
               fontsize=9, fontweight='bold')
    
    plt.tight_layout()
    
    if output_file:
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        print(f"Asymmetry ratio plot saved to {output_file}")
    
    return fig, ax


def plot_substitution_counts(results_df: pd.DataFrame,
                             output_file: str = None,
                             figsize: tuple = (12, 6),
                             title : str = 'Nonsynonymous vs Synonymous Substitutions\n' +
                'in Overlapping Reading Frames',
                            new_labels : Dict[str, str] = None):
    """
    Create stacked bar plot showing nonsynonymous vs synonymous substitutions as percentages.

    Parameters:
    -----------
    results_df : pd.DataFrame
        Results from compare_multiple_overlaps()
    output_file : str, optional
        Path to save figure
    figsize : tuple
        Figure size (width, height)
    """
    fig, ax = plt.subplots(figsize=figsize)

    # Prepare data
    regions = results_df['Region'].values
    n_subs_can = results_df['N_subs_canonical'].values
    s_subs_can = results_df['S_subs_canonical'].values
    n_subs_alt = results_df['N_subs_alternative'].values
    s_subs_alt = results_df['S_subs_alternative'].values

    # Calculate totals and percentages for each frame
    total_can = n_subs_can + s_subs_can
    total_alt = n_subs_alt + s_subs_alt

    # Avoid division by zero
    n_pct_can = np.where(total_can > 0, (n_subs_can / total_can) * 100, 0)
    s_pct_can = np.where(total_can > 0, (s_subs_can / total_can) * 100, 0)
    n_pct_alt = np.where(total_alt > 0, (n_subs_alt / total_alt) * 100, 0)
    s_pct_alt = np.where(total_alt > 0, (s_subs_alt / total_alt) * 100, 0)

    # Set up bar positions
    x = np.arange(len(regions))
    width = 0.35

    # Create stacked bars for canonical frame (percentages)
    p1 = ax.bar(x - width/2, n_pct_can, width, label='N subs (canonical)',
               color='#d62728', alpha=0.8)
    p2 = ax.bar(x - width/2, s_pct_can, width, bottom=n_pct_can,
               label='S subs (canonical)', color='#1f77b4', alpha=0.8)

    # Create stacked bars for alternative frame (percentages)
    p3 = ax.bar(x + width/2, n_pct_alt, width, label='N subs (alternative)',
               color='#ff7f0e', alpha=0.8)
    p4 = ax.bar(x + width/2, s_pct_alt, width, bottom=n_pct_alt,
               label='S subs (alternative)', color='#2ca02c', alpha=0.8)

    # Add percentage labels with absolute counts on bars
    for i in range(len(regions)):
        # Canonical frame - N subs (bottom segment)
        if n_pct_can[i] > 5:  # Only show label if segment is large enough
            ax.text(x[i] - width/2, n_pct_can[i]/2,
                   f'{n_pct_can[i]:.1f}%\n(n={int(n_subs_can[i])})',
                   ha='center', va='center', fontsize=8, fontweight='bold', color='white')

        # Canonical frame - S subs (top segment)
        if s_pct_can[i] > 5:
            ax.text(x[i] - width/2, n_pct_can[i] + s_pct_can[i]/2,
                   f'{s_pct_can[i]:.1f}%\n(n={int(s_subs_can[i])})',
                   ha='center', va='center', fontsize=8, fontweight='bold', color='white')

        # Alternative frame - N subs (bottom segment)
        if n_pct_alt[i] > 5:
            ax.text(x[i] + width/2, n_pct_alt[i]/2,
                   f'{n_pct_alt[i]:.1f}%\n(n={int(n_subs_alt[i])})',
                   ha='center', va='center', fontsize=8, fontweight='bold', color='white')

        # Alternative frame - S subs (top segment)
        if s_pct_alt[i] > 5:
            ax.text(x[i] + width/2, n_pct_alt[i] + s_pct_alt[i]/2,
                   f'{s_pct_alt[i]:.1f}%\n(n={int(s_subs_alt[i])})',
                   ha='center', va='center', fontsize=8, fontweight='bold', color='white')

    # Formatting
    ax.set_xlabel('Overlapping Region', fontsize=12, fontweight='bold')
    ax.set_ylabel('Percentage of Substitutions', fontsize=12, fontweight='bold')
    ax.set_title(title,
                fontsize=14, fontweight='bold', pad=20)
    ax.set_xticks(x)
    if new_labels:
        labels = [new_labels.get(label, label) for label in regions]
    else:
        labels = regions
    ax.set_xticklabels(labels, rotation=45, ha='center')

    # Set y-axis to 0-100%
    ax.set_ylim(0, 100)

    # Make legend smaller, transparent, and auto-position to avoid bars
    ax.legend(fontsize=8, loc='best', ncol=2, fancybox=True,
             framealpha=0.7, edgecolor='gray')

    # Grid
    ax.grid(True, alpha=0.3, axis='y')
    sns.despine()
    plt.tight_layout()

    if output_file:
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        print(f"Substitution counts plot saved to {output_file}")

    return fig, ax


def plot_beta_stop(results_df: pd.DataFrame,
                   output_file: str = None,
                   figsize: tuple = (12, 6),
                   title: str = 'Selection Against Stop Codons (βSTOP)\n' +
                   'Lower values indicate stronger functional constraint',
                   new_labels : Dict[str, str] = None,
                   alpha: float = 0.05,
                   multiple_testing_correction: str = None):
    """
    Create bar plot of Beta_STOP values with significance asterisks.

    Beta_STOP < 1 indicates selection against stop codons in alternative frame.
    Asterisks are added for significant values based on Beta_STOP_P_value.

    Parameters:
    -----------
    results_df : pd.DataFrame
        Results from compare_multiple_overlaps()
    output_file : str, optional
        Path to save figure
    figsize : tuple
        Figure size (width, height)
    title : str
        Title for the plot
    new_labels : Dict[str, str], optional
        Dictionary to replace region labels
    alpha : float
        Significance threshold (default: 0.05)
    multiple_testing_correction : str, optional
        Method for multiple testing correction. Options:
        - None: No correction
        - 'bonferroni': Bonferroni correction
        - 'fdr_bh': Benjamini-Hochberg FDR correction
        - 'sidak': Sidak correction
    """
    from statsmodels.stats.multitest import multipletests

    fig, ax = plt.subplots(figsize=figsize)

    # Check if Beta_STOP is in columns
    if 'Beta_STOP' not in results_df.columns:
        warnings.warn("Beta_STOP column not found in results")
        return fig, ax

    # Check if p-values are available
    if 'Beta_STOP_P_value' not in results_df.columns:
        warnings.warn("Beta_STOP_P_value column not found in results. Asterisks will not be displayed.")
        p_values = None
        significant = np.array([False] * len(results_df))
    else:
        p_values = results_df['Beta_STOP_P_value'].values

        # Apply multiple testing correction if requested
        if multiple_testing_correction is not None:
            if multiple_testing_correction == 'bonferroni':
                corrected_alpha = alpha / len(p_values)
                significant = p_values < corrected_alpha
                correction_label = f'Bonferroni-corrected α={corrected_alpha:.4f}'
            elif multiple_testing_correction == 'fdr_bh':
                reject, pvals_corrected, _, _ = multipletests(p_values, alpha=alpha, method='fdr_bh')
                significant = reject
                correction_label = f'FDR (Benjamini-Hochberg) α={alpha}'
            elif multiple_testing_correction == 'sidak':
                corrected_alpha = 1 - (1 - alpha) ** (1 / len(p_values))
                significant = p_values < corrected_alpha
                correction_label = f'Šidák-corrected α={corrected_alpha:.4f}'
            else:
                warnings.warn(f"Unknown correction method: {multiple_testing_correction}. Using uncorrected p-values.")
                significant = p_values < alpha
                correction_label = f'Uncorrected α={alpha}'
        else:
            # No correction
            significant = p_values < alpha
            correction_label = f'Uncorrected α={alpha}'

    # Prepare data
    regions = results_df['Region'].values
    beta_stop = results_df['Beta_STOP'].values

    # Set up bar positions
    x = np.arange(len(regions))

    # Create bars
    bars = ax.bar(x, beta_stop, color='purple', alpha=0.7, edgecolor='black')

    # Add horizontal line at Beta=1 (Neutral expectation if normalized correctly)
    ax.axhline(y=1, color='black', linestyle='--', linewidth=2, alpha=0.5,
              label='Neutral expectation')

    # Check if observed stop counts are available
    if 'Observed_STOP_F1' in results_df.columns:
        observed_stops = results_df['Observed_STOP_F1'].values
    else:
        observed_stops = None

    # Add significance asterisks above bars
    if p_values is not None:
        for i, (bar, val, sig, pval) in enumerate(zip(bars, beta_stop, significant, p_values)):
            height = bar.get_height()
            # Add value on bar
            ax.text(bar.get_x() + bar.get_width()/2., height,
                   f'{val:.2f}',
                   ha='center', va='bottom',
                   fontsize=10, fontweight='bold')
            # Add asterisk if significant
            if sig:
                ax.text(bar.get_x() + bar.get_width()/2., height * 1.15,
                       '*', ha='center', va='bottom',
                       fontsize=15, fontweight='bold', color='red')
            # Add observed stop count below x-axis
            if observed_stops is not None:
                ax.text(bar.get_x() + bar.get_width()/2., -0.02,
                       f'n={int(observed_stops[i])}',
                       ha='center', va='top',
                       fontsize=9, style='italic', color='darkblue',
                       transform=ax.get_xaxis_transform())
    else:
        # Just add values without significance
        for i, (bar, val) in enumerate(zip(bars, beta_stop)):
            height = bar.get_height()
            ax.text(bar.get_x() + bar.get_width()/2., height,
                   f'{val:.2f}',
                   ha='center', va='bottom',
                   fontsize=10, fontweight='bold')
            # Add observed stop count below x-axis
            if observed_stops is not None:
                ax.text(bar.get_x() + bar.get_width()/2., -0.02,
                       f'n={int(observed_stops[i])}',
                       ha='center', va='top',
                       fontsize=9, style='italic', color='darkblue',
                       transform=ax.get_xaxis_transform())

    # Formatting
    ax.set_xlabel('Overlapping Region', fontsize=10, fontweight='bold')
    ax.set_ylabel('βSTOP (Observed / Baseline)', fontsize=10, fontweight='bold')
    ax.set_title(title,
                fontsize=14, fontweight='bold', pad=20)
    ax.set_xticks(x)
    if new_labels:
        labels = [new_labels.get(label, label) for label in regions]
    else:
        labels = regions
    ax.set_xticklabels(labels, rotation=45, ha='center')

    # Update legend to include significance info
    legend_elements = [plt.Line2D([0], [0], color='black', linestyle='--', linewidth=2,
                                  label='Neutral expectation')]

    ax.legend(handles=legend_elements, fontsize=8, loc='upper right', bbox_to_anchor=(1.1, 1))

    # Grid
    ax.grid(True, alpha=0.3, axis='y')
    sns.despine()
    plt.tight_layout()

    if output_file:
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        print(f"Beta STOP plot saved to {output_file}")

    return fig, ax


def plot_raw_substitution_rates(results_df: pd.DataFrame,
                                 output_file: str = None,
                                 figsize: tuple = (14, 8),
                                 title: str = 'Raw Substitution Rates Across Overlapping Regions\n' +
                                 'SS/SN/NS/NN: Syn/Nonsyn in Canonical and Alternative frames',
                                 new_labels: Dict[str, str] = None,
                                 include_stop: bool = True):
    """
    Create grouped bar plot comparing raw substitution rates across all overlapping regions.

    This function plots the raw rates returned from ARFomeNormalized model:
    - SS: Synonymous in both frames
    - SN: Synonymous in canonical, Nonsynonymous in alternative
    - NS: Nonsynonymous in canonical, Synonymous in alternative
    - NN: Nonsynonymous in both frames
    - STOP_F1: Stop codon introduced in alternative frame

    Parameters:
    -----------
    results_df : pd.DataFrame
        Results from compare_multiple_overlaps()
    output_file : str, optional
        Path to save figure
    figsize : tuple
        Figure size (width, height)
    title : str
        Title of the plot
    new_labels : Dict[str, str], optional
        Dictionary to replace region labels
    include_stop : bool
        Whether to include STOP_F1 rate in the plot (default True)
    """
    # Check if required columns exist
    required_cols = ['Rate_SS', 'Rate_SN', 'Rate_NS', 'Rate_NN']
    if not all(col in results_df.columns for col in required_cols):
        warnings.warn("Raw rate columns not found in results_df. Make sure to use updated compare_multiple_overlaps.")
        return None, None

    fig, ax = plt.subplots(figsize=figsize)

    # Prepare data
    regions = results_df['Region'].values
    rate_ss = results_df['Rate_SS'].values
    rate_sn = results_df['Rate_SN'].values
    rate_ns = results_df['Rate_NS'].values
    rate_nn = results_df['Rate_NN'].values

    # Set up bar positions
    x = np.arange(len(regions))

    # Determine number of bars per group
    if include_stop and 'Rate_STOP_F1' in results_df.columns:
        rate_stop = results_df['Rate_STOP_F1'].values
        n_bars = 5
        width = 0.15  # Width of each bar
    else:
        rate_stop = None
        n_bars = 4
        width = 0.18

    # Create grouped bars
    offset = width * (n_bars - 1) / 2
    bars1 = ax.bar(x - offset + 0*width, rate_ss, width,
                   label='SS (Syn-Syn)', color='#2ca02c', alpha=0.8, edgecolor='black')
    bars2 = ax.bar(x - offset + 1*width, rate_sn, width,
                   label='SN (Syn-Nonsyn)', color='#ff7f0e', alpha=0.8, edgecolor='black')
    bars3 = ax.bar(x - offset + 2*width, rate_ns, width,
                   label='NS (Nonsyn-Syn)', color='#1f77b4', alpha=0.8, edgecolor='black')
    bars4 = ax.bar(x - offset + 3*width, rate_nn, width,
                   label='NN (Nonsyn-Nonsyn)', color='#d62728', alpha=0.8, edgecolor='black')

    if rate_stop is not None and include_stop:
        bars5 = ax.bar(x - offset + 4*width, rate_stop, width,
                       label='STOP_F1 (Stop in Alt)', color='#9467bd', alpha=0.8, edgecolor='black')

    # Formatting
    ax.set_xlabel('Overlapping Region', fontsize=12, fontweight='bold')
    ax.set_ylabel('Substitution Rate (Observed / Potential)', fontsize=12, fontweight='bold')
    ax.set_title(title, fontsize=14, fontweight='bold', pad=20)
    ax.set_xticks(x)

    if new_labels:
        labels = [new_labels.get(label, label) for label in regions]
    else:
        labels = regions
    ax.set_xticklabels(labels, rotation=45, ha='center')

    # Add legend with better positioning
    ax.legend(fontsize=10, loc='upper left', ncol=2, fancybox=True,
              framealpha=0.9, edgecolor='gray')

    # Grid
    ax.grid(True, alpha=0.3, axis='y')
    sns.despine()
    plt.tight_layout()

    if output_file:
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        print(f"Raw substitution rates plot saved to {output_file}")

    return fig, ax


def create_comprehensive_asymmetry_figure(results_df: pd.DataFrame,
                                         output_file: str = None,
                                         figsize: tuple = (18, 12)):
    """
    Create a comprehensive multi-panel figure showing all asymmetry analyses.
    
    This creates a figure similar to publication-quality multi-panel plots,
    combining scatter, bar, ratio plots, and Beta_STOP.
    
    Parameters:
    -----------
    results_df : pd.DataFrame
        Results from compare_multiple_overlaps()
    output_file : str, optional
        Path to save figure
    figsize : tuple
        Figure size (width, height)
    """
    fig = plt.figure(figsize=figsize)
    
    # Create grid for subplots (2x3 grid to fit everything, or 2x2 and combine?)
    # Let's do 2 rows, 3 columns layout
    gs = fig.add_gridspec(2, 3, hspace=0.3, wspace=0.3)
    
    # Panel A: Scatter plot (dN canonical vs alternative)
    ax1 = fig.add_subplot(gs[0, 0])
    dn_canonical = results_df['dN_canonical'].values
    dn_alternative = results_df['dN_alternative'].values
    significant = results_df['Significant_Asymmetry'].values
    colors = ['red' if sig else 'gray' for sig in significant]
    
    ax1.scatter(dn_canonical, dn_alternative, c=colors, s=100, 
               alpha=0.6, edgecolors='black', linewidth=1.5)
    max_val = max(dn_canonical.max(), dn_alternative.max()) * 1.1 if len(dn_canonical) > 0 else 1.0
    ax1.plot([0, max_val], [0, max_val], 'k--', alpha=0.5, linewidth=2)
    
    for idx, row in results_df.iterrows():
        ax1.annotate(row['Region'], 
                    (row['dN_canonical'], row['dN_alternative']),
                    xytext=(3, 3), textcoords='offset points',
                    fontsize=8, alpha=0.7)
    
    ax1.set_xlabel('dN Canonical', fontweight='bold')
    ax1.set_ylabel('dN Alternative', fontweight='bold')
    ax1.set_title('A. Asymmetric dN Evolution', fontweight='bold', loc='left')
    ax1.grid(True, alpha=0.3)
    ax1.set_aspect('equal', adjustable='box')
    
    # Panel B: Omega comparison
    ax2 = fig.add_subplot(gs[0, 1])
    regions = results_df['Region'].values
    omega_can = results_df['omega_canonical'].values
    omega_alt = results_df['omega_alternative'].values
    
    x = np.arange(len(regions))
    width = 0.35
    ax2.bar(x - width/2, omega_can, width, label='Canonical',
           color='steelblue', alpha=0.8, edgecolor='black')
    ax2.bar(x + width/2, omega_alt, width, label='Alternative',
           color='coral', alpha=0.8, edgecolor='black')
    ax2.axhline(y=1, color='black', linestyle='--', linewidth=2, alpha=0.5)
    
    for i, sig in enumerate(significant):
        if sig:
            max_height = max(omega_can[i], omega_alt[i])
            ax2.text(i, max_height * 1.1, '*', ha='center', va='bottom',
                    fontsize=16, fontweight='bold', color='red')
    
    ax2.set_xlabel('Region', fontweight='bold')
    ax2.set_ylabel('ω (dN/dS)', fontweight='bold')
    ax2.set_title('B. Selection Pressure (ω)', fontweight='bold', loc='left')
    ax2.set_xticks(x)
    ax2.set_xticklabels(regions, rotation=45, ha='center', fontsize=9)
    ax2.legend(fontsize=9)
    ax2.grid(True, alpha=0.3, axis='y')
    
    # Panel C: Beta STOP
    ax3 = fig.add_subplot(gs[0, 2])
    if 'Beta_STOP' in results_df.columns:
        beta_stop = results_df['Beta_STOP'].values
        ax3.bar(x, beta_stop, color='purple', alpha=0.7, edgecolor='black')
        ax3.axhline(y=1, color='black', linestyle='--', linewidth=2, alpha=0.5)
        
        for i, val in enumerate(beta_stop):
            ax3.text(i, val, f'{val:.2f}', ha='center', va='bottom', fontsize=8)
            
    ax3.set_xlabel('Region', fontweight='bold')
    ax3.set_ylabel('βSTOP', fontweight='bold')
    ax3.set_title('C. Stop Codon Selection', fontweight='bold', loc='left')
    ax3.set_xticks(x)
    ax3.set_xticklabels(regions, rotation=45, ha='center', fontsize=9)
    ax3.grid(True, alpha=0.3, axis='y')

    # Panel D: Asymmetry ratios
    ax4 = fig.add_subplot(gs[1, 0])
    ratios = results_df['Asymmetry_Ratio_dN'].values
    
    bars = ax4.bar(range(len(regions)), ratios, color=colors,
                  alpha=0.7, edgecolor='black', linewidth=1.5)
    ax4.axhline(y=1, color='black', linestyle='--', linewidth=2, alpha=0.7)
    
    for i, (bar, ratio) in enumerate(zip(bars, ratios)):
        height = bar.get_height()
        ax4.text(bar.get_x() + bar.get_width()/2., height,
                f'{ratio:.2f}',
                ha='center', va='bottom' if height > 1 else 'top',
                fontsize=8, fontweight='bold')
    
    ax4.set_xlabel('Region', fontweight='bold')
    ax4.set_ylabel('Asymmetry Ratio', fontweight='bold')
    ax4.set_title('D. Asymmetry Ratio (dN_alt/dN_can)', fontweight='bold', loc='left')
    ax4.set_xticks(range(len(regions)))
    ax4.set_xticklabels(regions, rotation=45, ha='center', fontsize=9)
    ax4.grid(True, alpha=0.3, axis='y')
    
    # Panel E: Substitution counts
    ax5 = fig.add_subplot(gs[1, 1:]) # Span 2 columns
    n_subs_can = results_df['N_subs_canonical'].values
    s_subs_can = results_df['S_subs_canonical'].values
    n_subs_alt = results_df['N_subs_alternative'].values
    s_subs_alt = results_df['S_subs_alternative'].values
    
    ax5.bar(x - width/2, n_subs_can, width, label='N (can)',
           color='#d62728', alpha=0.8)
    ax5.bar(x - width/2, s_subs_can, width, bottom=n_subs_can,
           label='S (can)', color='#1f77b4', alpha=0.8)
    ax5.bar(x + width/2, n_subs_alt, width, label='N (alt)',
           color='#ff7f0e', alpha=0.8)
    ax5.bar(x + width/2, s_subs_alt, width, bottom=n_subs_alt,
           label='S (alt)', color='#2ca02c', alpha=0.8)
    
    ax5.set_xlabel('Region', fontweight='bold')
    ax5.set_ylabel('Substitution Count', fontweight='bold')
    ax5.set_title('E. Substitution Types', fontweight='bold', loc='left')
    ax5.set_xticks(x)
    ax5.set_xticklabels(regions, rotation=45, ha='center', fontsize=9)
    ax5.legend(fontsize=8, ncol=2, loc='upper left')
    ax5.grid(True, alpha=0.3, axis='y')
    
    # Add main title
    fig.suptitle('Asymmetric Evolution in Overlapping Reading Frames',
                fontsize=16, fontweight='bold', y=0.98)
    
    if output_file:
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        print(f"Comprehensive figure saved to {output_file}")
    
    return fig


def plot_rrna_mdp_titv(results: Dict,
                       output_file: str = None,
                       figsize: tuple = (14, 6),
                       title: str = 'Ti/Tv Ratios: MDPs vs Random rRNA Regions',
                       new_labels: Dict[str, str] = None,
                       show_raw: bool = False,
                       alpha_threshold: float = 0.05,
                       show_global: bool = False,
                       test_results: Dict = None):
    """
    Plot Ti/Tv ratios for MDPs/genes vs random rRNA regions at each codon position.

    Shows violin plots of random region distributions with MDP values as markers
    for rRNA MDP analysis. For canonical gene analysis, shows only bar plots
    without random region background.

    Supports both single results dictionaries and combined dictionaries:
    - Single: output from analyze_rrna_mdp_asymmetry() or analyze_canonical_gene_titv()
    - Combined: {'rnr1_results': results1, 'rnr2_results': results2, ...}

    Parameters:
    -----------
    results : Dict
        Output from analyze_rrna_mdp_asymmetry(), analyze_canonical_gene_titv(),
        or a combined dictionary with multiple result sets
    output_file : str, optional
        Path to save figure
    figsize : tuple
        Figure size (width, height)
    title : str
        Plot title
    new_labels : Dict[str, str], optional
        Dictionary to rename MDP/gene labels
    show_raw : bool
        If True, show raw Ti/Tv ratios (can be inf); if False, show pseudocount ratios
    alpha_threshold : float
        Significance threshold for marking with asterisks
    show_global : bool
        Whether to include global statistics as a separate entry
    test_results : Dict, optional
        Statistical test results to display on the plot:
        - For MDPs: output from test_mdp_diversity_permutation()
        Displays p-values and markers for significant results.
        Also supports combined dictionaries matching the structure of results.

    Returns:
    --------
    fig, axes
    """
    # Detect if this is a combined dictionary
    is_combined = not any(key in results for key in ['mdp_results', 'gene_results', 'comparison', 'random_results'])

    # Flatten combined results into a single structure
    all_gene_data = {}
    all_random_data = {}
    all_comparison_data = {}
    all_test_results = {}
    global_mdp_stats = None
    global_random_stats = None
    global_gene_stats = None
    is_canonical_mode = False

    # Flatten test_results if it's also a combined dictionary
    if test_results is not None:
        test_is_combined = not any(isinstance(v, dict) and ('p_value_permutation' in v or 'p_value' in v)
                                   for v in test_results.values())
        if test_is_combined:
            # Flatten nested test results
            for source_key, source_test_results in test_results.items():
                if isinstance(source_test_results, dict):
                    for mdp_name, test_data in source_test_results.items():
                        all_test_results[mdp_name] = test_data
        else:
            # Already flat
            all_test_results = test_results
    else:
        all_test_results = None

    # Override test_results with flattened version
    if all_test_results is not None:
        test_results = all_test_results

    if is_combined:
        # Combined dictionary structure
        for source_key, source_results in results.items():
            if 'mdp_results' in source_results:
                # rRNA MDP results
                for gene_name, gene_data in source_results['mdp_results'].items():
                    all_gene_data[gene_name] = gene_data
                if 'random_results' in source_results:
                    for gene_name, random_data in source_results['random_results'].items():
                        all_random_data[gene_name] = random_data
                if 'comparison' in source_results:
                    for gene_name, comp_data in source_results['comparison'].items():
                        all_comparison_data[gene_name] = comp_data
            elif 'gene_results' in source_results:
                # Canonical gene results
                is_canonical_mode = True
                for gene_name, gene_data in source_results['gene_results'].items():
                    all_gene_data[gene_name] = gene_data
    else:
        # Single dictionary structure
        if 'mdp_results' in results:
            all_gene_data = results['mdp_results']
            all_random_data = results.get('random_results', {})
            all_comparison_data = results.get('comparison', {})
        elif 'gene_results' in results:
            is_canonical_mode = True
            all_gene_data = results['gene_results']
            global_gene_stats = results.get('global_gene_stats')
        else:
            raise ValueError("Invalid results dictionary. Expected output from analyze_rrna_mdp_asymmetry() or analyze_canonical_gene_titv()")

    # Get gene/MDP names
    gene_names = list(all_gene_data.keys())
    if len(gene_names) == 0:
        raise ValueError("No genes/MDPs found in results")

    # Add global statistics as a separate entry if available
    if not is_canonical_mode and global_mdp_stats is not None and show_global:
        gene_names.append('GLOBAL')
    elif is_canonical_mode and global_gene_stats is not None and show_global:
        gene_names.append('GLOBAL')


    # Apply label mapping if provided
    if new_labels:
        display_names = [new_labels.get(name, name) for name in gene_names if name != 'GLOBAL']
        if show_global:
            if is_canonical_mode:
                display_names.append('ALL GENES')
            else:
                display_names.append('ALL MDPs')
    else:
        display_names = [name if name != 'GLOBAL' else ('ALL GENES' if is_canonical_mode else 'ALL MDPs')
                        for name in gene_names]

    # Create figure with 3 subplots (one per codon position)
    fig, axes = plt.subplots(1, 3, figsize=figsize, sharey=False)

    positions = ['first', 'second', 'third']
    position_labels = ['1st Position', '2nd Position', '3rd Position']
    colors = ['#3498db', '#e74c3c', '#2ecc71']  # Blue, Red, Green

    for pos_idx, (pos, pos_label, color) in enumerate(zip(positions, position_labels, colors)):
        ax = axes[pos_idx]

        # Prepare data
        gene_values = []
        gene_positions = []
        significant_genes = []
        violin_data = []

        for gene_idx, gene_name in enumerate(gene_names):
            # Handle GLOBAL entry specially
            if gene_name == 'GLOBAL':
                if not is_canonical_mode and global_mdp_stats is not None:
                    # Use global MDP stats
                    gene_val = global_mdp_stats[f'{pos}_pos_titv']
                    # Add global random stats for violin plot
                    if global_random_stats is not None:
                        # Create a mock distribution centered on global random value
                        global_random_val = global_random_stats[f'{pos}_pos_titv']
                        violin_data.append([global_random_val])
                    else:
                        violin_data.append([np.nan])
                elif is_canonical_mode and global_gene_stats is not None:
                    # Use global gene stats
                    gene_val = global_gene_stats[f'{pos}_pos_titv']
                    violin_data.append([np.nan])
                else:
                    continue

                gene_values.append(gene_val)
                gene_positions.append(gene_idx + 1)
                significant_genes.append(False)  # Global is not individually tested
            else:
                gene_data = all_gene_data[gene_name]

                # Get gene Ti/Tv value
                if show_raw:
                    gene_val = gene_data[f'{pos}_pos']['ti_tv_ratio_raw']
                else:
                    gene_val = gene_data[f'{pos}_pos']['ti_tv_ratio']

                gene_values.append(gene_val)
                gene_positions.append(gene_idx + 1)

                # Check if significant (only for MDP mode)
                if not is_canonical_mode:
                    # First check if permutation test results are provided
                    if test_results is not None and gene_name in test_results:
                        is_sig = test_results[gene_name].get('significant_permutation', False)
                        significant_genes.append(is_sig)
                    # Fall back to comparison data from analyze_rrna_mdp_asymmetry
                    elif gene_name in all_comparison_data:
                        is_sig = all_comparison_data[gene_name].get('significant', False)
                        is_sig = is_sig == 'True' or is_sig is True
                        significant_genes.append(is_sig)
                    else:
                        significant_genes.append(False)
                else:
                    significant_genes.append(False)

                # Get random region data (only for MDP mode)
                if not is_canonical_mode and gene_name in all_random_data:
                    random_vals = all_random_data[gene_name].get(f'{pos}_pos_titv', [])
                    if len(random_vals) > 0:
                        violin_data.append(random_vals)
                    else:
                        violin_data.append([np.nan])
                else:
                    violin_data.append([np.nan])

        # Plot based on mode
        if not is_canonical_mode and len(violin_data) > 0 and any(len(v) > 0 for v in violin_data):
            # MDP mode: violin plot + scatter
            parts = ax.violinplot(violin_data, positions=range(1, len(gene_names) + 1),
                                 showmeans=False, showmedians=True, showextrema=False,
                                 widths=0.7)

            # Style violins
            for pc in parts['bodies']:
                pc.set_facecolor(color)
                pc.set_alpha(0.3)
                pc.set_edgecolor('black')
                pc.set_linewidth(1)

            # Style median lines
            parts['cmedians'].set_color('black')
            parts['cmedians'].set_linewidth(2)
            parts['cmedians'].set_alpha(0.7)

            # Overlay gene values as scatter points
            # Use gold star for GLOBAL entries, red diamond for significant, gray diamond for others
            for idx, (pos, val, sig, name) in enumerate(zip(gene_positions, gene_values, significant_genes, gene_names)):
                if name == 'GLOBAL':
                    ax.scatter([pos], [val], c='gold', s=150, marker='*',
                             edgecolors='black', linewidth=2, zorder=11)
                elif sig:
                    ax.scatter([pos], [val], c='red', s=80, marker='D',
                             edgecolors='black', linewidth=2, zorder=10)
                else:
                    ax.scatter([pos], [val], c='darkgray', s=80, marker='D',
                             edgecolors='black', linewidth=2, zorder=10)

            # Add significance asterisks
            for gene_idx, (gene_val, is_sig) in enumerate(zip(gene_values, significant_genes)):
                if is_sig and np.isfinite(gene_val):
                    if len(violin_data[gene_idx]) > 0:
                        max_violin = np.max([v for v in violin_data[gene_idx] if np.isfinite(v)])
                        y_pos = max(max_violin, gene_val) * 1.05
                    else:
                        y_pos = gene_val * 1.05

                    ax.text(gene_idx + 1, y_pos, '*', ha='center', va='bottom',
                           fontsize=12, fontweight='bold', color='red')
        else:
            # Canonical gene mode: bar plot only
            # Use different color for GLOBAL entry
            bar_colors = ['gold' if name == 'GLOBAL' else color for name in gene_names]
            ax.bar(gene_positions, gene_values, color=bar_colors, alpha=0.7, edgecolor='black', linewidth=1.5)

            # Add values on bars
            for gene_idx, gene_val in enumerate(gene_values):
                if np.isfinite(gene_val):
                    ax.text(gene_idx + 1, gene_val * 1.02, f'{gene_val:.2f}',
                           ha='center', va='bottom', fontsize=9, fontweight='bold')

        # Formatting
        label_text = 'Gene' if is_canonical_mode else 'MDP'
        ax.set_xlabel(label_text, fontsize=16, fontweight='bold')
        ax.set_ylabel('Ti/Tv Ratio' + (' (with +1 pseudocount)' if not show_raw else ' (raw)'),
                     fontsize=16, fontweight='bold')
        ax.set_title(pos_label, fontsize=16, fontweight='bold')
        ax.set_xticks(range(1, len(gene_names) + 1))
        ax.set_xticklabels(display_names, rotation=45, ha='center', fontsize=16)
        ax.grid(True, alpha=0.3, axis='y')

        # Add horizontal line at Ti/Tv = 1 (equal transitions and transversions)
        ax.axhline(y=1, color='gray', linestyle='--', linewidth=1.5, alpha=0.5)

    # Add legend to first subplot
    from matplotlib.patches import Patch
    if not is_canonical_mode:
        legend_elements = [
            plt.Line2D([0], [0], marker='D', color='w', markerfacecolor='red',
                      markersize=8, markeredgecolor='black', markeredgewidth=1.5,
                      label='MDP (Significant)', linestyle='None'),
            plt.Line2D([0], [0], marker='D', color='w', markerfacecolor='darkgray',
                      markersize=8, markeredgecolor='black', markeredgewidth=1.5,
                      label='MDP (N.S)', linestyle='None'),
            Patch(facecolor=colors[0], alpha=0.3, edgecolor='black',
                  label='Random rRNA regions'),
            plt.Line2D([0], [0], color='black', linewidth=2, alpha=0.7,
                      label='Median (random)')
        ]
        if show_global:
            legend_elements.append(
                plt.Line2D([0], [0], marker='*', color='w', markerfacecolor='gold',
                          markersize=12, markeredgecolor='black', markeredgewidth=1.5,
                          label='All MDPs combined', linestyle='None')
            )
        axes[0].legend(handles=legend_elements, loc='upper left', fontsize=12,
                       framealpha=0.9)
    else:
        pass

    # Add main title
    fig.suptitle(title, fontsize=14, fontweight='bold', y=0.98)

    # Add statistical test results annotation for MDP mode (permutation test)
    if not is_canonical_mode and test_results is not None:
        # Count significant MDPs
        n_significant = sum(1 for mdp_name, res in test_results.items()
                          if res.get('significant_permutation', False))
        n_total = len(test_results)

        # Get average p-value and n_permutations from first result
        avg_p_perm = np.mean([res.get('p_value_permutation', 1.0)
                             for res in test_results.values()])
        n_perm = list(test_results.values())[0].get('n_permutations', 'N/A') if test_results else 'N/A'

        # Create annotation text
        annotation_text = f"Permutation test (n={n_perm}): {n_significant}/{n_total} MDPs significant, avg p = {avg_p_perm:.4f}"

        # Add annotation at bottom of figure
        fig.text(0.5, 0.01, annotation_text,
                ha='center', va='bottom', fontsize=10, style='italic',
                bbox=dict(boxstyle='round,pad=0.5', facecolor='lightblue', alpha=0.3))

    plt.tight_layout(rect=[0, 0.03, 1, 0.96])

    # Despine
    for ax in axes:
        sns.despine(ax=ax)

    if output_file:
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        print(f"Ti/Tv plot saved to {output_file}")

    return fig, axes


def plot_rrna_mdp_diversity(results: Dict,
                            output_file: str = None,
                            figsize: tuple = (10, 6),
                            title: str = 'Nucleotide Diversity: MDPs vs Random rRNA Regions',
                            new_labels: Dict[str, str] = None,
                            alpha_threshold: float = 0.05,
                            show_global: bool = False,
                            test_results: Dict = None):
    """
    Plot nucleotide diversity for MDPs/genes vs random rRNA regions.

    Shows violin plots of random region distributions with MDP values as markers
    for rRNA MDP analysis. For canonical gene analysis, shows only bar plots
    without random region background.

    Supports both single results dictionaries and combined dictionaries:
    - Single: output from analyze_rrna_mdp_asymmetry() or analyze_canonical_gene_titv()
    - Combined: {'rnr1_results': results1, 'rnr2_results': results2, ...}

    Parameters:
    -----------
    results : Dict
        Output from analyze_rrna_mdp_asymmetry(), analyze_canonical_gene_titv(),
        or a combined dictionary with multiple result sets
    output_file : str, optional
        Path to save figure
    figsize : tuple
        Figure size (width, height)
    title : str
        Plot title
    new_labels : Dict[str, str], optional
        Dictionary to rename MDP/gene labels
    alpha_threshold : float
        Significance threshold for marking with asterisks
    show_global : bool
        Whether to include global statistics as a separate entry
    test_results : Dict, optional
        Statistical test results to display on the plot:
        - For canonical genes: output from test_codon_position_diversity_difference()
        - For MDPs: output from test_mdp_diversity_permutation()
        Displays p-values and markers for significant results.
        Also supports combined dictionaries matching the structure of results.

    Returns:
    --------
    fig, axes
    """
    # Detect if this is a combined dictionary
    is_combined = not any(key in results for key in ['mdp_results', 'gene_results', 'comparison', 'random_results'])

    # Flatten combined results into a single structure
    all_gene_data = {}
    all_random_data = {}
    all_comparison_data = {}
    all_test_results = {}
    global_mdp_stats = None
    global_random_stats = None
    global_gene_stats = None
    is_canonical_mode = False

    # Flatten test_results if it's also a combined dictionary
    if test_results is not None:
        test_is_combined = not any(isinstance(v, dict) and ('p_value_permutation' in v or 'p_value' in v)
                                   for v in test_results.values())
        if test_is_combined:
            # Flatten nested test results
            for source_key, source_test_results in test_results.items():
                if isinstance(source_test_results, dict):
                    for mdp_name, test_data in source_test_results.items():
                        all_test_results[mdp_name] = test_data
        else:
            # Already flat
            all_test_results = test_results
    else:
        all_test_results = None

    # Override test_results with flattened version
    if all_test_results is not None:
        test_results = all_test_results

    if is_combined:
        # Combined dictionary structure
        for source_key, source_results in results.items():
            if 'mdp_results' in source_results:
                # rRNA MDP results
                for gene_name, gene_data in source_results['mdp_results'].items():
                    all_gene_data[gene_name] = gene_data
                if 'random_results' in source_results:
                    for gene_name, random_data in source_results['random_results'].items():
                        all_random_data[gene_name] = random_data
                if 'comparison' in source_results:
                    for gene_name, comp_data in source_results['comparison'].items():
                        all_comparison_data[gene_name] = comp_data
                # Get global stats (from last source)
                global_mdp_stats = source_results.get('global_mdp_stats')
                global_random_stats = source_results.get('global_random_stats')
            elif 'gene_results' in source_results:
                # Canonical gene results
                is_canonical_mode = True
                for gene_name, gene_data in source_results['gene_results'].items():
                    all_gene_data[gene_name] = gene_data
                # Get global stats
                global_gene_stats = source_results.get('global_gene_stats')
    else:
        # Single dictionary structure
        if 'mdp_results' in results:
            all_gene_data = results['mdp_results']
            all_random_data = results.get('random_results', {})
            all_comparison_data = results.get('comparison', {})
            global_mdp_stats = results.get('global_mdp_stats')
            global_random_stats = results.get('global_random_stats')
        elif 'gene_results' in results:
            is_canonical_mode = True
            all_gene_data = results['gene_results']
            global_gene_stats = results.get('global_gene_stats')
        else:
            raise ValueError("Invalid results dictionary. Expected output from analyze_rrna_mdp_asymmetry() or analyze_canonical_gene_titv()")

    # Get gene/MDP names
    gene_names = list(all_gene_data.keys())
    if len(gene_names) == 0:
        raise ValueError("No genes/MDPs found in results")

    # Add global statistics as a separate entry if available

    if not is_canonical_mode and global_mdp_stats is not None and show_global:
        gene_names.append('GLOBAL')
    elif is_canonical_mode and global_gene_stats is not None and show_global:
        gene_names.append('GLOBAL')
    # Apply label mapping if provided
    if new_labels:
        display_names = [new_labels.get(name, name) for name in gene_names if name != 'GLOBAL']
        if show_global:
            if is_canonical_mode:
                display_names.append('ALL GENES')
            else:
                display_names.append('ALL MDPs')
    else:
        display_names = [name if name != 'GLOBAL' else ('ALL GENES' if is_canonical_mode else 'ALL MDPs')
                        for name in gene_names]

    # Create figure with 2 subplots
    fig, axes = plt.subplots(1, 2, figsize=figsize, sharey=True)

    colors = ['#9b59b6', '#1abc9c']  # Purple, Teal

    for plot_idx, (gene_metric, random_key, label) in enumerate([
        ('combined_12_diversity', 'combined_12_diversity', '1st+2nd Position'),
        ('third_pos', 'third_diversity', '3rd Position')
    ]):
        ax = axes[plot_idx]
        color = colors[plot_idx]

        # Prepare data
        gene_values = []
        gene_positions = []
        significant_genes = []
        violin_data = []

        for gene_idx, gene_name in enumerate(gene_names):
            # Handle GLOBAL entry specially
            if gene_name == 'GLOBAL':
                if not is_canonical_mode and global_mdp_stats is not None:
                    # Use global MDP stats
                    if plot_idx == 0:
                        gene_val = global_mdp_stats['roa_diversity_12']
                    else:
                        gene_val = global_mdp_stats['roa_diversity_third']

                    # Add global random stats for violin plot
                    if global_random_stats is not None:
                        if plot_idx == 0:
                            global_random_val = global_random_stats['roa_diversity_12']
                        else:
                            global_random_val = global_random_stats['roa_diversity_third']
                        violin_data.append([global_random_val])
                    else:
                        violin_data.append([np.nan])
                elif is_canonical_mode and global_gene_stats is not None:
                    # Use global gene stats
                    if plot_idx == 0:
                        gene_val = global_gene_stats['roa_diversity_12']
                    else:
                        gene_val = global_gene_stats['roa_diversity_third']
                    violin_data.append([np.nan])
                else:
                    continue

                gene_values.append(gene_val)
                gene_positions.append(gene_idx + 1)
                significant_genes.append(False)  # Global is not individually tested
            else:
                gene_data = all_gene_data[gene_name]

                # Get gene diversity value
                if plot_idx == 0:
                    gene_val = gene_data['combined_12_diversity']
                else:
                    gene_val = gene_data['third_pos']['nucleotide_diversity']

                gene_values.append(gene_val)
                gene_positions.append(gene_idx + 1)

                # Check if significant (only for MDP mode)
                if not is_canonical_mode:
                    # First check if permutation test results are provided
                    if test_results is not None and gene_name in test_results:
                        is_sig = test_results[gene_name].get('significant_empirical', False)
                        significant_genes.append(is_sig)
                    # Fall back to comparison data from analyze_rrna_mdp_asymmetry
                    elif gene_name in all_comparison_data:
                        is_sig = all_comparison_data[gene_name].get('significant', False)
                        is_sig = is_sig == 'True' or is_sig is True
                        significant_genes.append(is_sig)
                    else:
                        significant_genes.append(False)
                else:
                    significant_genes.append(False)

                # Get random region data (only for MDP mode)
                if not is_canonical_mode and gene_name in all_random_data:
                    random_vals = all_random_data[gene_name].get(random_key, [])
                    if len(random_vals) > 0:
                        violin_data.append(random_vals)
                    else:
                        violin_data.append([np.nan])
                else:
                    violin_data.append([np.nan])

        # Plot based on mode
        if not is_canonical_mode and len(violin_data) > 0 and any(len(v) > 0 for v in violin_data):
            # MDP mode: violin plot + scatter
            parts = ax.violinplot(violin_data, positions=range(1, len(gene_names) + 1),
                                 showmeans=False, showmedians=True, showextrema=False,
                                 widths=0.7)

            # Style violins
            for pc in parts['bodies']:
                pc.set_facecolor(color)
                pc.set_alpha(0.3)
                pc.set_edgecolor('black')
                pc.set_linewidth(1)

            # Style median lines
            parts['cmedians'].set_color('black')
            parts['cmedians'].set_linewidth(2)
            parts['cmedians'].set_alpha(0.7)

            # Overlay gene values as scatter points
            # Use gold star for GLOBAL entries, red diamond for significant, gray diamond for others
            for idx, (pos, val, sig, name) in enumerate(zip(gene_positions, gene_values, significant_genes, gene_names)):
                if name == 'GLOBAL':
                    ax.scatter([pos], [val], c='gold', s=150, marker='*',
                             edgecolors='black', linewidth=2, zorder=11)
                elif sig:
                    ax.scatter([pos], [val], c='red', s=80, marker='D',
                             edgecolors='black', linewidth=2, zorder=10)
                else:
                    ax.scatter([pos], [val], c='darkgray', s=80, marker='D',
                             edgecolors='black', linewidth=2, zorder=10)

            # Add significance asterisks
            for gene_idx, (gene_val, is_sig) in enumerate(zip(gene_values, significant_genes)):
                if is_sig and np.isfinite(gene_val):
                    if len(violin_data[gene_idx]) > 0:
                        max_violin = np.max([v for v in violin_data[gene_idx] if np.isfinite(v)])
                        y_pos = max(max_violin, gene_val) * 1.05
                    else:
                        y_pos = gene_val * 1.05

                    ax.text(gene_idx + 1, y_pos, '*', ha='center', va='bottom',
                           fontsize=12, fontweight='bold', color='red')
        else:
            # Canonical gene mode: bar plot only
            ax.bar(gene_positions, gene_values, color=color, alpha=0.7, edgecolor='black', linewidth=1.5)

            # Add values on bars
            for gene_idx, gene_val in enumerate(gene_values):
                if np.isfinite(gene_val):
                    ax.text(gene_idx + 1, gene_val * 1.02, f'{gene_val:.4f}',
                           ha='center', va='bottom', fontsize=9, fontweight='bold')

        # Formatting
        label_text = 'Gene' if is_canonical_mode else 'MDP'
        ax.set_xlabel(label_text, fontsize=16, fontweight='bold')
        if plot_idx == 0:
            ax.set_ylabel('Average nucleotide diversity (π)', fontsize=16, fontweight='bold')
        ax.set_title(label, fontsize=16, fontweight='bold')
        ax.set_xticks(range(1, len(gene_names) + 1))
        ax.set_xticklabels(display_names, rotation=45, ha='center', fontsize=16)
        ax.grid(True, alpha=0.3, axis='y')

    # Add legend to first subplot
    from matplotlib.patches import Patch
    if not is_canonical_mode:
        legend_elements = [
            plt.Line2D([0], [0], marker='D', color='w', markerfacecolor='red',
                      markersize=7, markeredgecolor='black', markeredgewidth=1.5,
                      label='MDP (Significant)', linestyle='None'),
            plt.Line2D([0], [0], marker='D', color='w', markerfacecolor='darkgray',
                      markersize=7, markeredgecolor='black', markeredgewidth=1.5,
                      label='MDP (N.S)', linestyle='None'),
            Patch(facecolor='gray', alpha=0.3, edgecolor='black',
                  label='Random rRNA regions'),
            plt.Line2D([0], [0], color='black', linewidth=2, alpha=0.7,
                      label='Median (random)')
        ]
        axes[0].legend(handles=legend_elements, loc='upper left', fontsize=12,
                       framealpha=0.9)
    else:
        # Canonical mode legend
        pass

    # Add main title
    fig.suptitle(title, fontsize=14, fontweight='bold', y=0.98)

    # Add statistical test results annotation for canonical mode
    if is_canonical_mode and test_results is not None and 'p_value_empirical' in test_results:
        p_val = test_results['p_value_empirical']
        cohens_d = test_results.get('cohens_d', np.nan)

        # Format p-value
        if p_val < 0.001:
            p_str = "p < 0.001"
        elif p_val < 0.01:
            p_str = f"p = {p_val:.3f}"
        else:
            p_str = f"p = {p_val:.4f}"

        # Format Cohen's d
        if np.isfinite(cohens_d):
            d_str = f"Cohen's d = {cohens_d:.2f}"
            annotation_text = f"{test_results['test_name']}: {p_str}, {d_str}"
        else:
            annotation_text = f"{test_results['test_name']}: {p_str}"

        # Add significance marker if applicable
        if test_results.get('significant_empirical', False):
            annotation_text += " *"

        # Add annotation at bottom of figure
        fig.text(0.5, 0.01, annotation_text,
                ha='center', va='bottom', fontsize=10, style='italic',
                bbox=dict(boxstyle='round,pad=0.5', facecolor='wheat', alpha=0.3))

    # Add statistical test results annotation for MDP mode (permutation test)
    if not is_canonical_mode and test_results is not None:
        # Count significant MDPs
        n_significant = sum(1 for mdp_name, res in test_results.items()
                          if res.get('significant_empirical', False))
        n_total = len(test_results)

        # Get average p-value and n_permutations from first result
        avg_p_perm = np.mean([res.get('p_value_empirical', 1.0)
                             for res in test_results.values()])
        n_perm = list(test_results.values())[0].get('n_permutations', 'N/A') if test_results else 'N/A'

        # Create annotation text
        annotation_text = f"Permutation test (n={n_perm}): {n_significant}/{n_total} MDPs significant, avg p = {avg_p_perm:.4f}"

        # Add annotation at bottom of figure
        fig.text(0.5, 0.01, annotation_text,
                ha='center', va='bottom', fontsize=10, style='italic',
                bbox=dict(boxstyle='round,pad=0.5', facecolor='lightblue', alpha=0.3))

    # Despine
    for ax in axes:
        sns.despine(ax=ax)

    plt.tight_layout(rect=[0, 0.03, 1, 0.96])

    if output_file:
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        print(f"Diversity plot saved to {output_file}")

    return fig, axes


def plot_rrna_mdp_titv_difference(results: Dict,
                                   output_file: str = None,
                                   figsize: tuple = (10, 6),
                                   title: str = 'Ti/Tv Difference (3rd - 1st+2nd): MDPs vs Random rRNA',
                                   new_labels: Dict[str, str] = None,
                                   alpha_threshold: float = 0.05,
                                   show_global: bool = False):
    """
    Plot Ti/Tv DIFFERENCE (3rd position - 1st+2nd positions) for MDPs/genes vs random regions.

    For rRNA MDP analysis: shows violin plots of random region distributions with MDP values.
    For canonical gene analysis: shows only bar plots without random region background.

    This plot directly visualizes what the statistical test is measuring: whether the
    codon-position-dependent asymmetry (protein-like pattern) is stronger in MDPs
    than in random rRNA regions.

    A positive difference indicates Ti/Tv(3rd) > Ti/Tv(1st+2nd), which is the
    protein-coding pattern (synonymous 3rd codon positions tolerate more transitions).

    Supports both single results dictionaries and combined dictionaries:
    - Single: output from analyze_rrna_mdp_asymmetry() or analyze_canonical_gene_titv()
    - Combined: {'rnr1_results': results1, 'rnr2_results': results2, ...}

    Parameters:
    -----------
    results : Dict
        Output from analyze_rrna_mdp_asymmetry(), analyze_canonical_gene_titv(),
        or a combined dictionary with multiple result sets
    output_file : str, optional
        Path to save figure
    figsize : tuple
        Figure size (width, height)
    title : str
        Plot title
    new_labels : Dict[str, str], optional
        Dictionary to rename MDP/gene labels
    alpha_threshold : float
        Significance threshold for marking with asterisks
    show_global : bool
        Whether to include global statistics as a separate entry

    Returns:
    --------
    fig, ax
    """
    # Detect if this is a combined dictionary
    is_combined = not any(key in results for key in ['mdp_results', 'gene_results', 'comparison', 'random_results'])

    # Flatten combined results into a single structure
    all_gene_data = {}
    all_random_data = {}
    all_comparison_data = {}
    global_mdp_stats = None
    global_random_stats = None
    global_gene_stats = None
    is_canonical_mode = False

    if is_combined:
        # Combined dictionary structure
        for source_key, source_results in results.items():
            if 'mdp_results' in source_results:
                # rRNA MDP results
                for gene_name, gene_data in source_results['mdp_results'].items():
                    all_gene_data[gene_name] = gene_data
                if 'random_results' in source_results:
                    for gene_name, random_data in source_results['random_results'].items():
                        all_random_data[gene_name] = random_data
                if 'comparison' in source_results:
                    for gene_name, comp_data in source_results['comparison'].items():
                        all_comparison_data[gene_name] = comp_data
                # Get global stats (from last source)
                global_mdp_stats = source_results.get('global_mdp_stats')
                global_random_stats = source_results.get('global_random_stats')
            elif 'gene_results' in source_results:
                # Canonical gene results
                is_canonical_mode = True
                for gene_name, gene_data in source_results['gene_results'].items():
                    all_gene_data[gene_name] = gene_data
                # Get global stats
                global_gene_stats = source_results.get('global_gene_stats')
    else:
        # Single dictionary structure
        if 'mdp_results' in results:
            all_gene_data = results['mdp_results']
            all_random_data = results.get('random_results', {})
            all_comparison_data = results.get('comparison', {})
            global_mdp_stats = results.get('global_mdp_stats')
            global_random_stats = results.get('global_random_stats')
        elif 'gene_results' in results:
            is_canonical_mode = True
            all_gene_data = results['gene_results']
            global_gene_stats = results.get('global_gene_stats')
        else:
            raise ValueError("Invalid results dictionary. Expected output from analyze_rrna_mdp_asymmetry() or analyze_canonical_gene_titv()")

    # Get gene/MDP names
    gene_names = list(all_gene_data.keys())
    if len(gene_names) == 0:
        raise ValueError("No genes/MDPs found in results")

    # Add global statistics as a separate entry if available
    if not is_canonical_mode and global_mdp_stats is not None and show_global:
        gene_names.append('GLOBAL')
    elif is_canonical_mode and global_gene_stats is not None and show_global:
        gene_names.append('GLOBAL')


    # Apply label mapping if provided
    if new_labels:
        display_names = [new_labels.get(name, name) for name in gene_names if name != 'GLOBAL']
        if show_global:
            if is_canonical_mode:
                display_names.append('ALL GENES')
            else:
                display_names.append('ALL MDPs')
    else:
        display_names = [name if name != 'GLOBAL' else ('ALL GENES' if is_canonical_mode else 'ALL MDPs')
                        for name in gene_names]

    # Create figure
    fig, ax = plt.subplots(figsize=figsize)

    # Prepare data
    gene_diffs = []
    significant_genes = []
    p_values = []
    violin_data = []

    for gene_idx, gene_name in enumerate(gene_names):
        # Handle GLOBAL entry specially
        if gene_name == 'GLOBAL':
            if not is_canonical_mode and global_mdp_stats is not None:
                # Calculate global MDP Ti/Tv difference
                global_titv_12 = global_mdp_stats['combined_12_titv']
                global_titv_3 = global_mdp_stats['third_pos_titv']
                gene_diff = global_titv_3 - global_titv_12 if (np.isfinite(global_titv_12) and np.isfinite(global_titv_3)) else np.nan

                # Add global random stats for violin plot
                if global_random_stats is not None:
                    global_random_titv_12 = global_random_stats['combined_12_titv']
                    global_random_titv_3 = global_random_stats['third_pos_titv']
                    global_random_diff = global_random_titv_3 - global_random_titv_12 if (
                        np.isfinite(global_random_titv_12) and np.isfinite(global_random_titv_3)
                    ) else np.nan
                    violin_data.append([global_random_diff])
                else:
                    violin_data.append([np.nan])
            elif is_canonical_mode and global_gene_stats is not None:
                # Calculate global gene Ti/Tv difference
                global_titv_12 = global_gene_stats['combined_12_titv']
                global_titv_3 = global_gene_stats['third_pos_titv']
                gene_diff = global_titv_3 - global_titv_12 if (np.isfinite(global_titv_12) and np.isfinite(global_titv_3)) else np.nan
                violin_data.append([np.nan])
            else:
                continue

            gene_diffs.append(gene_diff)
            significant_genes.append(False)  # Global is not individually tested
            p_values.append(np.nan)
        else:
            gene_data = all_gene_data[gene_name]

            # Calculate gene Ti/Tv difference
            if is_canonical_mode:
                # For canonical genes, use titv_diff if available
                gene_diff = gene_data.get('titv_diff', np.nan)
                if np.isnan(gene_diff):
                    # Calculate it if not pre-calculated
                    titv_12 = gene_data.get('combined_12_titv', np.nan)
                    titv_3 = gene_data.get('third_pos', {}).get('ti_tv_ratio', np.nan)
                    gene_diff = titv_3 - titv_12 if (np.isfinite(titv_12) and np.isfinite(titv_3)) else np.nan
            else:
                # For MDP results, get from comparison
                if gene_name in all_comparison_data:
                    gene_diff = all_comparison_data[gene_name].get('mdp_titv_diff', np.nan)
                else:
                    # Calculate it if not in comparison
                    titv_12 = gene_data.get('combined_12_titv', np.nan)
                    titv_3 = gene_data.get('third_pos', {}).get('ti_tv_ratio', np.nan)
                    gene_diff = titv_3 - titv_12 if (np.isfinite(titv_12) and np.isfinite(titv_3)) else np.nan

            gene_diffs.append(gene_diff)

            # Get significance (only for MDP mode)
            if not is_canonical_mode and gene_name in all_comparison_data:
                comp = all_comparison_data[gene_name]
                is_sig = comp.get('significant', False)
                is_sig = is_sig == 'True' or is_sig is True
                significant_genes.append(is_sig)
                p_values.append(comp.get('p_value', np.nan))
            else:
                significant_genes.append(False)
                p_values.append(np.nan)

            # Get random region data (only for MDP mode)
            if not is_canonical_mode and gene_name in all_random_data:
                random = all_random_data[gene_name]
                random_titv_12 = random.get('combined_12_titv', [])
                random_titv_3 = random.get('third_pos_titv', [])

                # Calculate differences for each random region
                random_diffs = []
                for i in range(min(len(random_titv_12), len(random_titv_3))):
                    if np.isfinite(random_titv_12[i]) and np.isfinite(random_titv_3[i]):
                        diff = random_titv_3[i] - random_titv_12[i]
                        random_diffs.append(diff)

                violin_data.append(random_diffs if len(random_diffs) > 0 else [np.nan])
            else:
                violin_data.append([np.nan])

    # Plot based on mode
    positions = range(1, len(gene_names) + 1)

    if not is_canonical_mode and len(violin_data) > 0 and any(len(v) > 0 for v in violin_data):
        # MDP mode: violin plot + scatter
        parts = ax.violinplot(violin_data, positions=positions,
                             showmeans=False, showmedians=True, showextrema=False,
                             widths=0.7)

        # Style violins
        for pc in parts['bodies']:
            pc.set_facecolor('#3498db')
            pc.set_alpha(0.3)
            pc.set_edgecolor('black')
            pc.set_linewidth(1)

        # Style median lines
        parts['cmedians'].set_color('black')
        parts['cmedians'].set_linewidth(2)
        parts['cmedians'].set_alpha(0.7)

        # Overlay gene differences as scatter points
        # Use gold star for GLOBAL entries, red diamond for significant, gray diamond for others
        for idx, (pos, diff, sig, name) in enumerate(zip(positions, gene_diffs, significant_genes, gene_names)):
            if name == 'GLOBAL':
                ax.scatter([pos], [diff], c='gold', s=150, marker='*',
                         edgecolors='black', linewidth=2, zorder=11)
            elif sig:
                ax.scatter([pos], [diff], c='red', s=80, marker='D',
                         edgecolors='black', linewidth=2, zorder=10)
            else:
                ax.scatter([pos], [diff], c='darkgray', s=80, marker='D',
                         edgecolors='black', linewidth=2, zorder=10)

        # Add significance asterisks
        for gene_idx, (gene_diff, is_sig) in enumerate(zip(gene_diffs, significant_genes)):
            if is_sig and np.isfinite(gene_diff):
                if len(violin_data[gene_idx]) > 0:
                    valid_random = [v for v in violin_data[gene_idx] if np.isfinite(v)]
                    if len(valid_random) > 0:
                        max_violin = np.max(valid_random)
                        y_pos = max(max_violin, gene_diff) * 1.05
                    else:
                        y_pos = gene_diff * 1.05
                else:
                    y_pos = gene_diff * 1.05

                ax.text(gene_idx + 1, y_pos, '*', ha='center', va='bottom',
                       fontsize=12, fontweight='bold', color='red')
    else:
        # Canonical gene mode: bar plot only
        # Color bars: gold for GLOBAL, green for positive diff, red for negative diff
        bar_colors = []
        for name, diff in zip(gene_names, gene_diffs):
            if name == 'GLOBAL':
                bar_colors.append('gold')
            elif diff > 0:
                bar_colors.append('#2ecc71')
            else:
                bar_colors.append('#e74c3c')
        ax.bar(positions, gene_diffs, color=bar_colors, alpha=0.7, edgecolor='black', linewidth=1.5)

        # Add values on bars
        for gene_idx, gene_diff in enumerate(gene_diffs):
            if np.isfinite(gene_diff) and gene_diff > 0:
                y_offset = gene_diff * 1.02 if gene_diff > 0 else gene_diff * 0.98
                va = 'bottom' if gene_diff > 0 else 'top'
                ax.text(gene_idx + 1, y_offset, f'{gene_diff:.3f}',
                       ha='center', va=va, fontsize=9, fontweight='bold')

        # Add horizontal line at 0
        ax.axhline(y=0, color='black', linestyle='-', linewidth=1.5, alpha=0.7)

    # Formatting
    label_text = 'Gene' if is_canonical_mode else 'MDP'
    ax.set_xlabel(label_text, fontsize=16, fontweight='bold')
    ax.set_ylabel('Ti/Tv Difference\n[Ti/Tv(3rd) - Ti/Tv(1st+2nd)]', fontsize=16, fontweight='bold')
    ax.set_title(title, fontsize=14, fontweight='bold', pad=20)
    ax.set_xticks(positions)
    ax.set_xticklabels(display_names, rotation=45, ha='center', fontsize=16)
    ax.grid(True, alpha=0.3, axis='y')

    # Add legend
    from matplotlib.patches import Patch
    if not is_canonical_mode:
        legend_elements = [
            plt.Line2D([0], [0], marker='D', color='w', markerfacecolor='red',
                      markersize=7, markeredgecolor='black', markeredgewidth=1.5,
                      label='MDP (Significant)', linestyle='None'),
            plt.Line2D([0], [0], marker='D', color='w', markerfacecolor='darkgray',
                      markersize=7, markeredgecolor='black', markeredgewidth=1.5,
                      label='MDP (N.S)', linestyle='None'),
            Patch(facecolor='#3498db', alpha=0.3, edgecolor='black',
                  label='Random rRNA regions'),
            plt.Line2D([0], [0], color='black', linewidth=2, alpha=0.7,
                      label='Median (random)'),
        ]
        if show_global:
            legend_elements.append(
                plt.Line2D([0], [0], marker='*', color='w', markerfacecolor='gold',
                          markersize=12, markeredgecolor='black', markeredgewidth=1.5,
                          label='All MDPs combined', linestyle='None')
            )
        ax.legend(handles=legend_elements, loc='upper left', fontsize=8,
                 framealpha=0.9)
    else:
        # Canonical mode legend
        legend_elements = [
            Patch(facecolor='#2ecc71', alpha=0.7, edgecolor='black',
                  label='Protein pattern (positive)'),
            Patch(facecolor='#e74c3c', alpha=0.7, edgecolor='black',
                  label='No protein pattern (negative)')
        ]
        ax.legend(handles=legend_elements, loc='upper left', fontsize=12,
                 framealpha=0.9)

    sns.despine(ax=ax)
    plt.tight_layout()

    if output_file:
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        print(f"Ti/Tv difference plot saved to {output_file}")

    return fig, ax


def plot_mcnemar_asymmetry(results_df: pd.DataFrame,
                           output_file: str = None,
                           figsize: tuple = (16, 5),
                           title: str = 'McNemar Test: Frame Asymmetry Analysis',
                           new_labels: Dict[str, str] = None):
    """
    Create visualization of McNemar test results for overlapping frame asymmetry.

    Shows:
    1. Heatmaps of contingency tables for each region
    2. Bar chart comparing off-diagonal counts (SN vs NS) - the key test
    3. P-values and significance markers

    The McNemar test evaluates whether mutations have asymmetric effects on the
    two overlapping frames. If SN ≠ NS, there is significant asymmetry.

    Contingency table structure:
                    Alternative Frame
                    Syn         Nonsyn
    Canonical  Syn  SS          SN
               Nonsyn NS        NN

    Parameters:
    -----------
    results_df : pd.DataFrame
        Results from compare_multiple_overlaps() with McNemar test results
    output_file : str, optional
        Path to save figure
    figsize : tuple
        Figure size (width, height)
    title : str
        Main title for the plot
    new_labels : Dict[str, str], optional
        Dictionary to replace region labels

    Returns:
    --------
    fig, axes
    """
    # Check if McNemar results are present
    required_cols = ['Observed_SS', 'Observed_SN', 'Observed_NS', 'Observed_NN', 'Fisher_P_value']
    if not all(col in results_df.columns for col in required_cols):
        warnings.warn("McNemar test columns not found in results_df. Required: Obs_SS, Obs_SN, Obs_NS, Obs_NN, McNemar_P_value")
        return None, None

    # Prepare data
    regions = results_df['Region'].values
    ss = results_df['Observed_SS'].values
    sn = results_df['Observed_SN'].values
    ns = results_df['Observed_NS'].values
    nn = results_df['Observed_NN'].values
    p_values = results_df['Fisher_P_value'].values
    significant = results_df['Significant_Asymmetry'].values

    # Apply label mapping
    if new_labels:
        display_names = [new_labels.get(name, name) for name in regions]
    else:
        display_names = list(regions)

    # Create figure with 2 subplots
    fig, axes = plt.subplots(1, 2, figsize=figsize)

    # ===== Panel A: Heatmaps of contingency tables =====
    ax1 = axes[0]

    # Create a grid for multiple small heatmaps
    n_regions = len(regions)
    n_cols = min(4, n_regions)  # Max 4 heatmaps per row
    n_rows = int(np.ceil(n_regions / n_cols))

    # Create subgrid within ax1
    from matplotlib.gridspec import GridSpecFromSubplotSpec
    gs = GridSpecFromSubplotSpec(n_rows, n_cols, subplot_spec=ax1.get_subplotspec(),
                                  hspace=0.4, wspace=0.3)
    ax1.axis('off')  # Turn off the parent axis

    # Plot each contingency table as a small heatmap
    for idx, (region, ss_val, sn_val, ns_val, nn_val, pval, sig) in enumerate(
        zip(display_names, ss, sn, ns, nn, p_values, significant)):

        row = idx // n_cols
        col = idx % n_cols
        ax_sub = fig.add_subplot(gs[row, col])

        # Create contingency table
        table = np.array([[ss_val, sn_val],
                         [ns_val, nn_val]])

        # Plot heatmap
        im = ax_sub.imshow(table, cmap='YlOrRd', aspect='auto')

        # Add text annotations
        for i in range(2):
            for j in range(2):
                text = ax_sub.text(j, i, f'{int(table[i, j])}',
                                 ha="center", va="center", color="black",
                                 fontsize=12, fontweight='bold')

        # Set labels
        ax_sub.set_xticks([0, 1])
        ax_sub.set_yticks([0, 1])
        ax_sub.set_xticklabels(['Syn', 'Nonsyn'], fontsize=9)
        ax_sub.set_yticklabels(['Syn', 'Nonsyn'], fontsize=9)
        ax_sub.set_xlabel('Alternative', fontsize=9, fontweight='bold')
        ax_sub.set_ylabel('Canonical', fontsize=9, fontweight='bold')

        # Title with region name and significance
        sig_marker = '*' if sig else ''
        title_color = 'red' if sig else 'black'
        ax_sub.set_title(f'{region}{sig_marker}\np={pval:.3f}',
                        fontsize=9, fontweight='bold', color=title_color)

        # Highlight off-diagonal cells (key for McNemar test)
        from matplotlib.patches import Rectangle
        # SN cell (top-right)
        rect1 = Rectangle((0.5, -0.5), 1, 1, linewidth=2, edgecolor='blue',
                         facecolor='none', linestyle='--')
        ax_sub.add_patch(rect1)
        # NS cell (bottom-left)
        rect2 = Rectangle((-0.5, 0.5), 1, 1, linewidth=2, edgecolor='blue',
                         facecolor='none', linestyle='--')
        ax_sub.add_patch(rect2)

    # ===== Panel B: Bar chart of off-diagonal counts =====
    ax2 = axes[1]

    x = np.arange(len(regions))
    width = 0.35

    # Create bars for SN and NS
    bars1 = ax2.bar(x - width/2, sn, width, label='SN (Syn→Nonsyn)',
                   color='#ff7f0e', alpha=0.8, edgecolor='black', linewidth=1.5)
    bars2 = ax2.bar(x + width/2, ns, width, label='NS (Nonsyn→Syn)',
                   color='#1f77b4', alpha=0.8, edgecolor='black', linewidth=1.5)

    # Add value labels on bars
    for bar in bars1:
        height = bar.get_height()
        if height > 0:
            ax2.text(bar.get_x() + bar.get_width()/2., height,
                   f'{int(height)}',
                   ha='center', va='bottom', fontsize=9, fontweight='bold')

    for bar in bars2:
        height = bar.get_height()
        if height > 0:
            ax2.text(bar.get_x() + bar.get_width()/2., height,
                   f'{int(height)}',
                   ha='center', va='bottom', fontsize=9, fontweight='bold')

    # Add significance asterisks
    for i, (sig, pval) in enumerate(zip(significant, p_values)):
        if sig:
            max_height = max(sn[i], ns[i])
            ax2.text(i, max_height * 1.15, '***' if pval < 0.001 else '**' if pval < 0.01 else '*',
                   ha='center', va='bottom', fontsize=14, fontweight='bold', color='red')

    # Formatting
    ax2.set_xlabel('Overlapping Region', fontsize=12, fontweight='bold')
    ax2.set_ylabel('Mutation Count', fontsize=12, fontweight='bold')
    ax2.set_title('Off-Diagonal Counts (McNemar Test Focus)\n' +
                 'If SN ≠ NS: Asymmetric selection pressure',
                 fontsize=12, fontweight='bold', pad=10)
    ax2.set_xticks(x)
    ax2.set_xticklabels(display_names, rotation=45, ha='right', fontsize=10)
    ax2.legend(fontsize=10, loc='upper left', framealpha=0.9)
    ax2.grid(True, alpha=0.3, axis='y')

    # Add annotation explaining the test
    annotation_text = (
        "McNemar's Test compares SN vs NS counts.\n"
        "SN = Syn in canonical, Nonsyn in alternative\n"
        "NS = Nonsyn in canonical, Syn in alternative\n"
        "Significant difference indicates asymmetric evolution"
    )
    ax2.text(0.02, 0.98, annotation_text,
           transform=ax2.transAxes, fontsize=8, verticalalignment='top',
           bbox=dict(boxstyle='round', facecolor='lightblue', alpha=0.3))

    # Main title
    fig.suptitle(title, fontsize=14, fontweight='bold', y=0.98)

    plt.tight_layout(rect=[0, 0, 1, 0.96])
    sns.despine(ax=ax2)

    if output_file:
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        print(f"McNemar asymmetry plot saved to {output_file}")

    return fig, axes


def plot_mutation_frequency_distributions(mutations_df: pd.DataFrame,
                                          frequency_threshold: float = 1.0,
                                          output_file: str = None,
                                          figsize: tuple = (16, 12),
                                          title: str = 'Mutation Frequency Distributions by Category'):
    """
    Create comprehensive visualization of mutation frequency distributions across categories.

    This function displays:
    1. Mutation count by category
    2. Frequency distribution (violin + swarm plot) for each category
    3. Histogram of frequencies highlighting mutations above threshold
    4. Summary statistics table

    Parameters:
    -----------
    mutations_df : pd.DataFrame
        DataFrame with columns: Region, Mutation_Category, Mutation, Count, Frequency_Pct
    frequency_threshold : float
        Threshold (in %) to highlight high-frequency mutations (default: 1.0%)
    output_file : str, optional
        Path to save figure
    figsize : tuple
        Figure size (width, height)
    title : str
        Main title for the figure

    Returns:
    --------
    fig, axes : matplotlib figure and axes objects
    """
    if mutations_df.empty:
        print("Warning: Empty mutations DataFrame provided")
        return None, None

    # Define category order and colors
    categories = ['SS', 'SN', 'NS', 'NN', 'STOP_F1']
    colors = {
        'SS': '#2ecc71',      # Green - neutral in both
        'SN': '#e74c3c',      # Red - deleterious in alt
        'NS': '#3498db',      # Blue - deleterious in can
        'NN': '#9b59b6',      # Purple - deleterious in both
        'STOP_F1': '#e67e22'  # Orange - stop in alt
    }

    # Filter to only categories present in data
    present_categories = [cat for cat in categories if cat in mutations_df['Mutation_Category'].unique()]

    # Create figure with subplots
    fig = plt.figure(figsize=figsize)
    gs = fig.add_gridspec(3, 2, hspace=0.3, wspace=0.3, height_ratios=[1, 1.5, 1])

    # ==================== Panel 1: Mutation Counts by Category ====================
    ax1 = fig.add_subplot(gs[0, :])

    # Count mutations per category
    category_counts = mutations_df.groupby('Mutation_Category').size().reindex(present_categories, fill_value=0)
    high_freq_counts = mutations_df[mutations_df['Frequency_Pct'] >= frequency_threshold].groupby('Mutation_Category').size().reindex(present_categories, fill_value=0)

    x = np.arange(len(present_categories))
    width = 0.35

    bars1 = ax1.bar(x - width/2, category_counts, width, label='All mutations',
                   color=[colors[cat] for cat in present_categories], alpha=0.7)
    bars2 = ax1.bar(x + width/2, high_freq_counts, width, label=f'≥{frequency_threshold}% frequency',
                   color=[colors[cat] for cat in present_categories], alpha=1.0, edgecolor='black', linewidth=2)

    # Add count labels on bars
    for bars in [bars1, bars2]:
        for bar in bars:
            height = bar.get_height()
            if height > 0:
                ax1.text(bar.get_x() + bar.get_width()/2., height,
                        f'{int(height)}',
                        ha='center', va='bottom', fontsize=9, fontweight='bold')

    ax1.set_xlabel('Mutation Category', fontsize=12, fontweight='bold')
    ax1.set_ylabel('Number of Mutations', fontsize=12, fontweight='bold')
    ax1.set_title('Mutation Count by Category', fontsize=13, fontweight='bold', pad=10)
    ax1.set_xticks(x)
    ax1.set_xticklabels(present_categories)
    ax1.legend(fontsize=10, loc='upper right')
    ax1.grid(True, alpha=0.3, axis='y')

    # ==================== Panel 2: Frequency Distribution (Violin + Swarm) ====================
    ax2 = fig.add_subplot(gs[1, :])

    # Prepare data for plotting
    plot_data = mutations_df[mutations_df['Mutation_Category'].isin(present_categories)].copy()

    # Violin plot
    parts = ax2.violinplot([plot_data[plot_data['Mutation_Category'] == cat]['Frequency_Pct'].values
                           for cat in present_categories],
                          positions=range(len(present_categories)),
                          widths=0.7,
                          showmeans=False,
                          showmedians=True)

    # Color the violin plots
    for i, pc in enumerate(parts['bodies']):
        pc.set_facecolor(colors[present_categories[i]])
        pc.set_alpha(0.3)

    # Add swarm plot for individual points
    for i, cat in enumerate(present_categories):
        cat_data = plot_data[plot_data['Mutation_Category'] == cat]
        y_vals = cat_data['Frequency_Pct'].values

        if len(y_vals) > 0:
            # Add jitter to x positions
            x_vals = np.random.normal(i, 0.04, size=len(y_vals))

            # Separate high and low frequency mutations
            high_freq_mask = y_vals >= frequency_threshold

            # Plot low frequency mutations
            ax2.scatter(x_vals[~high_freq_mask], y_vals[~high_freq_mask],
                       alpha=0.6, s=30, color=colors[cat], edgecolors='gray', linewidths=0.5)

            # Highlight high frequency mutations
            if high_freq_mask.sum() > 0:
                ax2.scatter(x_vals[high_freq_mask], y_vals[high_freq_mask],
                           alpha=1.0, s=80, color=colors[cat], edgecolors='black', linewidths=2,
                           marker='D', label=f'{cat} ≥{frequency_threshold}%' if i == 0 else '')

    # Add horizontal line at threshold
    ax2.axhline(y=frequency_threshold, color='red', linestyle='--', linewidth=2, alpha=0.7,
               label=f'{frequency_threshold}% threshold')

    ax2.set_xlabel('Mutation Category', fontsize=12, fontweight='bold')
    ax2.set_ylabel('Frequency in Population (%)', fontsize=12, fontweight='bold')
    ax2.set_title('Frequency Distribution of Individual Mutations', fontsize=13, fontweight='bold', pad=10)
    ax2.set_xticks(range(len(present_categories)))
    ax2.set_xticklabels(present_categories)
    ax2.set_yscale('log')
    ax2.grid(True, alpha=0.3, axis='y')
    ax2.legend(fontsize=9, loc='upper right')

    # ==================== Panel 3: Frequency Histogram ====================
    ax3 = fig.add_subplot(gs[2, 0])

    # Create histogram
    all_freqs = mutations_df['Frequency_Pct'].values
    bins = np.logspace(np.log10(all_freqs.min()), np.log10(all_freqs.max()), 30)

    ax3.hist(all_freqs[all_freqs < frequency_threshold], bins=bins, alpha=0.6,
            color='gray', label=f'<{frequency_threshold}%', edgecolor='black')
    ax3.hist(all_freqs[all_freqs >= frequency_threshold], bins=bins, alpha=0.9,
            color='red', label=f'≥{frequency_threshold}%', edgecolor='black')

    ax3.axvline(x=frequency_threshold, color='red', linestyle='--', linewidth=2, alpha=0.7)
    ax3.set_xlabel('Frequency (%)', fontsize=11, fontweight='bold')
    ax3.set_ylabel('Count', fontsize=11, fontweight='bold')
    ax3.set_title('Frequency Distribution (All Categories)', fontsize=12, fontweight='bold')
    ax3.set_xscale('log')
    ax3.legend(fontsize=9)
    ax3.grid(True, alpha=0.3)

    # ==================== Panel 4: Summary Statistics Table ====================
    ax4 = fig.add_subplot(gs[2, 1])
    ax4.axis('off')

    # Calculate summary statistics
    stats_data = []
    for cat in present_categories:
        cat_data = mutations_df[mutations_df['Mutation_Category'] == cat]
        high_freq_data = cat_data[cat_data['Frequency_Pct'] >= frequency_threshold]

        stats_data.append([
            cat,
            len(cat_data),
            len(high_freq_data),
            f"{cat_data['Frequency_Pct'].mean():.2f}",
            f"{cat_data['Frequency_Pct'].median():.2f}",
            f"{cat_data['Frequency_Pct'].max():.2f}"
        ])

    # Create table
    table_headers = ['Category', 'Total\nMutations', f'Mutations\n≥{frequency_threshold}%',
                    'Mean\nFreq (%)', 'Median\nFreq (%)', 'Max\nFreq (%)']

    table = ax4.table(cellText=stats_data,
                     colLabels=table_headers,
                     cellLoc='center',
                     loc='center',
                     colWidths=[0.15, 0.15, 0.18, 0.17, 0.17, 0.17])

    table.auto_set_font_size(False)
    table.set_fontsize(9)
    table.scale(1, 2)

    # Color code the category column
    for i, cat in enumerate(present_categories):
        table[(i+1, 0)].set_facecolor(colors[cat])
        table[(i+1, 0)].set_alpha(0.6)

    # Style header
    for j in range(len(table_headers)):
        table[(0, j)].set_facecolor('#34495e')
        table[(0, j)].set_text_props(weight='bold', color='white')

    ax4.set_title('Summary Statistics', fontsize=12, fontweight='bold', pad=20)

    # ==================== Main Title and Layout ====================
    fig.suptitle(title, fontsize=15, fontweight='bold', y=0.98)

    # Add category explanations
    explanation = (
        "Category Definitions:\n"
        "SS: Synonymous in both frames (neutral)\n"
        "SN: Synonymous in canonical, Nonsynonymous in alternative\n"
        "NS: Nonsynonymous in canonical, Synonymous in alternative\n"
        "NN: Nonsynonymous in both frames\n"
        "STOP_F1: Introduces stop codon in alternative frame"
    )
    fig.text(0.02, 0.02, explanation, fontsize=8,
            bbox=dict(boxstyle='round', facecolor='lightgray', alpha=0.5))

    plt.tight_layout(rect=[0, 0.08, 1, 0.96])

    if output_file:
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        print(f"Mutation frequency distribution plot saved to {output_file}")

    return fig, (ax1, ax2, ax3, ax4)


def plot_mutation_frequency_by_region(mutations_df: pd.DataFrame,
                                     frequency_threshold: float = 1.0,
                                     output_file: str = None,
                                     figsize: tuple = (14, 10),
                                     title: str = 'Mutation Frequencies Across Regions'):
    """
    Compare mutation frequency distributions across different overlapping regions.

    Parameters:
    -----------
    mutations_df : pd.DataFrame
        DataFrame with columns: Region, Mutation_Category, Mutation, Count, Frequency_Pct
    frequency_threshold : float
        Threshold (in %) to highlight high-frequency mutations (default: 1.0%)
    output_file : str, optional
        Path to save figure
    figsize : tuple
        Figure size (width, height)
    title : str
        Main title for the figure

    Returns:
    --------
    fig, axes : matplotlib figure and axes objects
    """
    if mutations_df.empty:
        print("Warning: Empty mutations DataFrame provided")
        return None, None

    regions = mutations_df['Region'].unique()
    n_regions = len(regions)

    if n_regions == 0:
        print("Warning: No regions found in mutations DataFrame")
        return None, None

    categories = ['SS', 'SN', 'NS', 'NN', 'STOP_F1']
    colors = {
        'SS': '#2ecc71',
        'SN': '#e74c3c',
        'NS': '#3498db',
        'NN': '#9b59b6',
        'STOP_F1': '#e67e22'
    }

    # Create figure
    fig, axes = plt.subplots(2, 2, figsize=figsize)
    axes = axes.flatten()

    # ==================== Panel 1: Mutation Counts by Region and Category ====================
    ax = axes[0]

    # Prepare data for grouped bar chart
    count_data = mutations_df.groupby(['Region', 'Mutation_Category']).size().unstack(fill_value=0)
    # Reindex to ensure all regions and categories are present
    count_data = count_data.reindex(index=regions, columns=categories, fill_value=0)

    x = np.arange(n_regions)
    width = 0.15
    multiplier = 0

    for i, cat in enumerate(categories):
        if cat in count_data.columns:
            offset = width * multiplier
            ax.bar(x + offset, count_data[cat], width, label=cat,
                  color=colors[cat], alpha=0.8)
            multiplier += 1

    ax.set_xlabel('Region', fontsize=11, fontweight='bold')
    ax.set_ylabel('Number of Mutations', fontsize=11, fontweight='bold')
    ax.set_title('Total Mutations by Category', fontsize=12, fontweight='bold')
    ax.set_xticks(x + width * 2)
    ax.set_xticklabels(regions, rotation=45, ha='right')
    ax.legend(fontsize=8, ncol=2, loc='upper right')
    ax.grid(True, alpha=0.3, axis='y')

    # ==================== Panel 2: High-Frequency Mutations (≥threshold) ====================
    ax = axes[1]

    high_freq_data = mutations_df[mutations_df['Frequency_Pct'] >= frequency_threshold]
    high_count_data = high_freq_data.groupby(['Region', 'Mutation_Category']).size().unstack(fill_value=0)
    # Reindex to ensure all regions and categories are present
    high_count_data = high_count_data.reindex(index=regions, columns=categories, fill_value=0)

    multiplier = 0
    for i, cat in enumerate(categories):
        if cat in high_count_data.columns:
            offset = width * multiplier
            ax.bar(x + offset, high_count_data[cat], width, label=cat,
                  color=colors[cat], alpha=0.8, edgecolor='black', linewidth=1.5)
            multiplier += 1

    ax.set_xlabel('Region', fontsize=11, fontweight='bold')
    ax.set_ylabel(f'Mutations ≥{frequency_threshold}%', fontsize=11, fontweight='bold')
    ax.set_title(f'High-Frequency Mutations (≥{frequency_threshold}%)', fontsize=12, fontweight='bold')
    ax.set_xticks(x + width * 2)
    ax.set_xticklabels(regions, rotation=45, ha='right')
    ax.legend(fontsize=8, ncol=2, loc='upper right')
    ax.grid(True, alpha=0.3, axis='y')

    # ==================== Panel 3: Mean Frequency by Category ====================
    ax = axes[2]

    mean_freq_data = mutations_df.groupby(['Region', 'Mutation_Category'])['Frequency_Pct'].mean().unstack(fill_value=0)
    # Reindex to ensure all regions and categories are present
    mean_freq_data = mean_freq_data.reindex(index=regions, columns=categories, fill_value=0)

    multiplier = 0
    for i, cat in enumerate(categories):
        if cat in mean_freq_data.columns:
            offset = width * multiplier
            ax.bar(x + offset, mean_freq_data[cat], width, label=cat,
                  color=colors[cat], alpha=0.7)
            multiplier += 1

    ax.axhline(y=frequency_threshold, color='red', linestyle='--', linewidth=2, alpha=0.5,
              label=f'{frequency_threshold}% threshold')

    ax.set_xlabel('Region', fontsize=11, fontweight='bold')
    ax.set_ylabel('Mean Frequency (%)', fontsize=11, fontweight='bold')
    ax.set_title('Mean Mutation Frequency by Category', fontsize=12, fontweight='bold')
    ax.set_xticks(x + width * 2)
    ax.set_xticklabels(regions, rotation=45, ha='right')
    ax.legend(fontsize=8, ncol=2, loc='upper right')
    ax.grid(True, alpha=0.3, axis='y')

    # ==================== Panel 4: Proportion of High-Frequency Mutations ====================
    ax = axes[3]

    # Calculate proportion of high-frequency mutations for each region and category
    prop_data = []
    for region in regions:
        region_data = mutations_df[mutations_df['Region'] == region]
        props = []
        for cat in categories:
            cat_data = region_data[region_data['Mutation_Category'] == cat]
            if len(cat_data) > 0:
                prop = (cat_data['Frequency_Pct'] >= frequency_threshold).sum() / len(cat_data) * 100
            else:
                prop = 0
            props.append(prop)
        prop_data.append(props)

    prop_df = pd.DataFrame(prop_data, index=regions, columns=categories)

    multiplier = 0
    for i, cat in enumerate(categories):
        offset = width * multiplier
        ax.bar(x + offset, prop_df[cat], width, label=cat,
              color=colors[cat], alpha=0.7)
        multiplier += 1

    ax.set_xlabel('Region', fontsize=11, fontweight='bold')
    ax.set_ylabel(f'% Mutations ≥{frequency_threshold}%', fontsize=11, fontweight='bold')
    ax.set_title(f'Proportion of High-Frequency Mutations', fontsize=12, fontweight='bold')
    ax.set_xticks(x + width * 2)
    ax.set_xticklabels(regions, rotation=45, ha='right')
    ax.legend(fontsize=8, ncol=2, loc='upper right')
    ax.grid(True, alpha=0.3, axis='y')

    # ==================== Main Title ====================
    fig.suptitle(title, fontsize=14, fontweight='bold', y=0.98)

    plt.tight_layout(rect=[0, 0, 1, 0.96])

    if output_file:
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        print(f"Region comparison plot saved to {output_file}")

    return fig, axes


if __name__ == "__main__":
    print("Visualization module for overlapping frame asymmetry analysis")
    print("Functions available:")
    print("  - plot_asymmetry_scatter()")
    print("  - plot_omega_comparison()")
    print("  - plot_asymmetry_ratios()")
    print("  - plot_substitution_counts()")
    print("  - plot_beta_stop()")
    print("  - plot_raw_substitution_rates()")
    print("  - plot_mcnemar_asymmetry()")
    print("  - create_comprehensive_asymmetry_figure()")
    print("  - plot_rrna_mdp_titv()")
    print("  - plot_rrna_mdp_diversity()")
    print("  - plot_rrna_mdp_titv_difference()")
    print("  - plot_mutation_frequency_distributions()")
    print("  - plot_mutation_frequency_by_region()")
