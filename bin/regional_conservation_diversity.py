"""
Regional conservation and nucleotide diversity analysis for gene sequences.

This module provides functions to analyze per-position conservation and nucleotide
diversity within gene regions, with special focus on comparing overlapping MDP
regions to non-overlapping regions.
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import os
from collections import defaultdict
from scipy import stats
from statsmodels.stats.multitest import multipletests
from statsmodels.tsa.stattools import acf
from statannotations.Annotator import Annotator
from adjustText import adjust_text


def calculate_position_conservation(sequences, position):
    """
    Calculate conservation score for a specific position across sequences.

    Args:
        sequences (list): List of aligned sequences
        position (int): Position to analyze

    Returns:
        float: Conservation score [0-1] based on Shannon entropy
    """
    # Get all bases at this position
    bases = [seq[position] for seq in sequences if position < len(seq)]

    # Count occurrences of each base
    counts = defaultdict(int)
    for base in bases:
        if base != '-' and base != 'N':  # Ignore gaps and N's
            counts[base] += 1

    # Calculate Shannon entropy: H = -∑(pi * log2(pi))
    total = sum(counts.values())
    if total == 0:
        return np.nan

    entropy = 0
    for count in counts.values():
        p = count / total
        entropy -= p * np.log2(p)

    # Convert entropy to conservation score [0-1]
    max_entropy = np.log2(4)  # Maximum possible entropy for 4 DNA bases
    conservation = 1 - (entropy / max_entropy)

    return conservation


def calculate_position_diversity(sequences, position):
    """
    Calculate nucleotide diversity (π) for a specific position.

    Nucleotide diversity is the average number of nucleotide differences
    per site between any two sequences.

    Args:
        sequences (list): List of aligned sequences
        position (int): Position to analyze

    Returns:
        float: Nucleotide diversity [0-1]
    """
    # Get all bases at this position
    bases = [seq[position] for seq in sequences if position < len(seq)]

    # Filter out gaps and N's
    bases = [b for b in bases if b not in ['-', 'N']]

    if len(bases) < 2:
        return np.nan

    # Count occurrences of each base
    counts = defaultdict(int)
    for base in bases:
        counts[base] += 1

    # Calculate nucleotide diversity using the formula:
    # π = n/(n-1) * (1 - Σ(pi^2))
    n = len(bases)
    sum_squared_freq = sum((count/n)**2 for count in counts.values())
    diversity = (n/(n-1)) * (1 - sum_squared_freq)

    return diversity


def calculate_per_position_conservation(df, gene_col, overlapping_regions=None):
    """
    Calculate per-position conservation scores across all sequences in a gene.

    Args:
        df (pd.DataFrame): DataFrame with sequences (one row per individual)
        gene_col (str): Column name containing the aligned gene sequences
        overlapping_regions (list, optional): List of [name, start, end] for overlapping regions

    Returns:
        pd.DataFrame: DataFrame with columns: position, conservation, region_type, region_name
    """
    # Filter valid sequences
    valid_seqs = df[df[gene_col].notna() & (df[gene_col] != '')][gene_col].tolist()

    if len(valid_seqs) == 0:
        raise ValueError(f"No valid sequences found in column '{gene_col}'")

    # Get alignment length
    seq_lengths = [len(str(seq)) for seq in valid_seqs]
    max_length = max(seq_lengths)

    # Calculate conservation for each position
    results = []
    for pos in range(max_length):
        conservation = calculate_position_conservation(valid_seqs, pos)

        # Determine region type and name
        region_type = 'non-overlapping'
        region_name = 'non-overlapping'

        if overlapping_regions:
            for region in overlapping_regions:
                region_name_temp, start, end = region
                # Position is 0-indexed, regions are 1-indexed, so adjust
                if start - 1 <= pos < end:
                    region_type = 'overlapping'
                    region_name = region_name_temp
                    break

        results.append({
            'position': pos + 1,  # 1-indexed for output
            'conservation': conservation,
            'region_type': region_type,
            'region_name': region_name
        })

    return pd.DataFrame(results)


def calculate_per_position_diversity(df, gene_col, overlapping_regions=None):
    """
    Calculate per-position nucleotide diversity across all sequences in a gene.

    Args:
        df (pd.DataFrame): DataFrame with sequences (one row per individual)
        gene_col (str): Column name containing the aligned gene sequences
        overlapping_regions (list, optional): List of [name, start, end] for overlapping regions

    Returns:
        pd.DataFrame: DataFrame with columns: position, diversity, region_type, region_name
    """
    # Filter valid sequences
    valid_seqs = df[df[gene_col].notna() & (df[gene_col] != '')][gene_col].tolist()

    if len(valid_seqs) == 0:
        raise ValueError(f"No valid sequences found in column '{gene_col}'")

    # Get alignment length
    seq_lengths = [len(str(seq)) for seq in valid_seqs]
    max_length = max(seq_lengths)

    # Calculate diversity for each position
    results = []
    for pos in range(max_length):
        diversity = calculate_position_diversity(valid_seqs, pos)

        # Determine region type and name
        region_type = 'non-overlapping'
        region_name = 'non-overlapping'

        if overlapping_regions:
            for region in overlapping_regions:
                region_name_temp, start, end = region
                # Position is 0-indexed, regions are 1-indexed, so adjust
                if start - 1 <= pos < end:
                    region_type = 'overlapping'
                    region_name = region_name_temp
                    break

        results.append({
            'position': pos + 1,  # 1-indexed for output
            'diversity': diversity,
            'region_type': region_type,
            'region_name': region_name
        })

    return pd.DataFrame(results)


def calculate_position_kl_divergence(sequences, position, background_freq):
    """
    Calculate KL divergence for a specific position compared to background.

    KL divergence measures how much the nucleotide distribution at a position
    differs from the background distribution. High KL = position is "picky" about
    which bases it allows (e.g., strong purine bias, specific transitions allowed).

    Args:
        sequences (list): List of aligned sequences
        position (int): Position to analyze
        background_freq (dict): Background nucleotide frequencies {base: freq}

    Returns:
        float: KL divergence value (0 = matches background, higher = more constrained)
    """
    # Get all bases at this position
    bases = [seq[position] for seq in sequences if position < len(seq)]

    # Filter out gaps and N's
    bases = [b for b in bases if b in 'ACGT']

    if len(bases) == 0:
        return np.nan

    # Calculate observed frequencies at this position
    obs_freq = {}
    for base in 'ACGT':
        obs_freq[base] = bases.count(base) / len(bases)

    # Calculate KL divergence: D_KL(P || Q) = sum(P(i) * log(P(i) / Q(i)))
    # P = observed frequencies, Q = background frequencies
    kl_div = 0
    for base in 'ACGT':
        p = obs_freq[base]
        q = background_freq[base]

        # Add small pseudocount to avoid log(0)
        if p > 0 and q > 0:
            kl_div += p * np.log2(p / q)
        elif p > 0 and q == 0:
            # If background has 0 but observed has non-zero, use large value
            kl_div += p * 10  # Arbitrary large penalty

    return kl_div


def calculate_per_position_kl_divergence(df, gene_col, overlapping_regions=None):
    """
    Calculate per-position KL divergence across all sequences in a gene.

    Each position's nucleotide distribution is compared to the background
    distribution across the entire gene. This highlights positions with
    non-random nucleotide preferences (e.g., purine-only sites, specific
    transition biases).

    Args:
        df (pd.DataFrame): DataFrame with sequences (one row per individual)
        gene_col (str): Column name containing the aligned gene sequences
        overlapping_regions (list, optional): List of [name, start, end] for overlapping regions

    Returns:
        pd.DataFrame: DataFrame with columns: position, kl_divergence, region_type, region_name,
                     plus A_freq, C_freq, G_freq, T_freq for the position
    """
    # Filter valid sequences
    valid_seqs = df[df[gene_col].notna() & (df[gene_col] != '')][gene_col].tolist()

    if len(valid_seqs) == 0:
        raise ValueError(f"No valid sequences found in column '{gene_col}'")

    # Calculate background nucleotide frequencies across entire gene
    all_bases = []
    for seq in valid_seqs:
        all_bases.extend([b for b in str(seq) if b in 'ACGT'])

    background_freq = {}
    total_bases = len(all_bases)
    for base in 'ACGT':
        background_freq[base] = all_bases.count(base) / total_bases if total_bases > 0 else 0.25

    # Get alignment length
    seq_lengths = [len(str(seq)) for seq in valid_seqs]
    max_length = max(seq_lengths)

    # Calculate KL divergence for each position
    results = []
    for pos in range(max_length):
        # Get bases at this position
        bases = [seq[pos] for seq in valid_seqs if pos < len(seq)]
        bases_clean = [b for b in bases if b in 'ACGT']

        # Calculate position-specific frequencies
        pos_freq = {}
        if len(bases_clean) > 0:
            for base in 'ACGT':
                pos_freq[base] = bases_clean.count(base) / len(bases_clean)
        else:
            pos_freq = {base: np.nan for base in 'ACGT'}

        # Calculate KL divergence
        kl_div = calculate_position_kl_divergence(valid_seqs, pos, background_freq)

        # Determine region type and name
        region_type = 'non-overlapping'
        region_name = 'non-overlapping'

        if overlapping_regions:
            for region in overlapping_regions:
                region_name_temp, start, end = region
                # Position is 0-indexed, regions are 1-indexed, so adjust
                if start - 1 <= pos < end:
                    region_type = 'overlapping'
                    region_name = region_name_temp
                    break

        results.append({
            'position': pos + 1,  # 1-indexed for output
            'kl_divergence': kl_div,
            'region_type': region_type,
            'region_name': region_name,
            'A_freq': pos_freq['A'],
            'C_freq': pos_freq['C'],
            'G_freq': pos_freq['G'],
            'T_freq': pos_freq['T']
        })

    return pd.DataFrame(results)


def calculate_nucleotide_composition_enrichment(df, gene_col, overlapping_regions=None):
    """
    Calculate nucleotide composition enrichment at each position relative to background.

    This shows which bases are over- or under-represented at each position compared
    to the overall gene composition. Useful for identifying specific nucleotide biases
    in overlapping regions (e.g., purine-enrichment, pyrimidine-depletion).

    Args:
        df (pd.DataFrame): DataFrame with sequences (one row per individual)
        gene_col (str): Column name containing the aligned gene sequences
        overlapping_regions (list, optional): List of [name, start, end] for overlapping regions

    Returns:
        pd.DataFrame: DataFrame with position, region info, and enrichment scores for each base
                     Enrichment = log2(observed_freq / background_freq)
                     Positive = enriched, Negative = depleted, 0 = matches background
    """
    # Filter valid sequences
    valid_seqs = df[df[gene_col].notna() & (df[gene_col] != '')][gene_col].tolist()

    if len(valid_seqs) == 0:
        raise ValueError(f"No valid sequences found in column '{gene_col}'")

    # Calculate background nucleotide frequencies across entire gene
    all_bases = []
    for seq in valid_seqs:
        all_bases.extend([b for b in str(seq) if b in 'ACGT'])

    background_freq = {}
    total_bases = len(all_bases)
    for base in 'ACGT':
        count = all_bases.count(base)
        background_freq[base] = count / total_bases if total_bases > 0 else 0.25

    # Get alignment length
    seq_lengths = [len(str(seq)) for seq in valid_seqs]
    max_length = max(seq_lengths)

    # Calculate enrichment for each position
    results = []
    for pos in range(max_length):
        # Get bases at this position
        bases = [seq[pos] for seq in valid_seqs if pos < len(seq)]
        bases_clean = [b for b in bases if b in 'ACGT']

        # Calculate position-specific frequencies
        enrichment = {}
        obs_freq = {}
        if len(bases_clean) > 0:
            for base in 'ACGT':
                obs_freq[base] = bases_clean.count(base) / len(bases_clean)

                # Calculate enrichment as log2(observed / background)
                # Positive = enriched, Negative = depleted
                if obs_freq[base] > 0 and background_freq[base] > 0:
                    enrichment[base] = np.log2(obs_freq[base] / background_freq[base])
                elif obs_freq[base] == 0:
                    enrichment[base] = -10  # Very depleted
                else:
                    enrichment[base] = 10  # Very enriched
        else:
            obs_freq = {base: np.nan for base in 'ACGT'}
            enrichment = {base: np.nan for base in 'ACGT'}

        # Determine region type and name
        region_type = 'non-overlapping'
        region_name = 'non-overlapping'

        if overlapping_regions:
            for region in overlapping_regions:
                region_name_temp, start, end = region
                if start - 1 <= pos < end:
                    region_type = 'overlapping'
                    region_name = region_name_temp
                    break

        results.append({
            'position': pos + 1,
            'region_type': region_type,
            'region_name': region_name,
            'A_freq': obs_freq['A'],
            'C_freq': obs_freq['C'],
            'G_freq': obs_freq['G'],
            'T_freq': obs_freq['T'],
            'A_enrichment': enrichment['A'],
            'C_enrichment': enrichment['C'],
            'G_enrichment': enrichment['G'],
            'T_enrichment': enrichment['T'],
            'purine_freq': obs_freq['A'] + obs_freq['G'],
            'pyrimidine_freq': obs_freq['C'] + obs_freq['T']
        })

    return pd.DataFrame(results)


def calculate_autocorrelation(results_df, metric_col='conservation', max_lag=50,
                              by_region=True, by_region_name=False, alpha=0.05):
    """
    Calculate autocorrelation of a metric to detect periodic patterns.

    Autocorrelation measures the correlation between a signal and a lagged version
    of itself. Peaks at specific lags indicate periodic patterns (e.g., lag 3 for
    codon periodicity in overlapping ORFs).

    Args:
        results_df (pd.DataFrame): Output from calculate_per_position_* functions
        metric_col (str): Column to analyze ('conservation', 'diversity', 'kl_divergence')
        max_lag (int): Maximum lag to calculate (default 50). Will be automatically
                       adjusted to min(max_lag, n_positions // 2) for short regions.
        by_region (bool): Calculate separate autocorrelations for overlapping vs non-overlapping
                         (ignored if by_region_name=True)
        by_region_name (bool): Calculate separate autocorrelations for each individual region name
                              (e.g., SHLP1, SHLP6, humanin). Overrides by_region if True.
        alpha (float): Significance level for confidence intervals (default 0.05)

    Returns:
        dict: Autocorrelation results with lags, acf_values, confidence intervals, and region info
    """
    results = {}

    def _calculate_acf_for_subset(subset, region_label, max_lag):
        """Helper function to calculate ACF for a subset of data."""
        subset_clean = subset.dropna(subset=[metric_col])

        if len(subset_clean) == 0:
            return None

        # Adjust max_lag for short regions (need at least 2*max_lag points)
        effective_max_lag = min(max_lag, len(subset_clean) // 2)

        # Need at least effective_max_lag + 1 points
        if len(subset_clean) <= effective_max_lag:
            return None

        values = subset_clean[metric_col].values

        # Calculate autocorrelation
        acf_values = acf(values, nlags=effective_max_lag, alpha=alpha, fft=True)

        # acf returns (acf, confint) when alpha is specified
        if isinstance(acf_values, tuple):
            acf_vals, confint = acf_values
            result = {
                'lags': np.arange(effective_max_lag + 1),
                'acf': acf_vals,
                'confint_lower': confint[:, 0] - acf_vals,
                'confint_upper': confint[:, 1] - acf_vals,
                'n_positions': len(values),
                'metric': metric_col,
                'max_lag_used': effective_max_lag
            }
        else:
            result = {
                'lags': np.arange(effective_max_lag + 1),
                'acf': acf_values,
                'confint_lower': None,
                'confint_upper': None,
                'n_positions': len(values),
                'metric': metric_col,
                'max_lag_used': effective_max_lag
            }

        # Find significant peaks (excluding lag 0)
        significant_lags = []
        if result['confint_lower'] is not None:
            for i in range(1, len(acf_vals)):
                # Check if ACF is significantly different from 0
                if acf_vals[i] > 1.96/np.sqrt(len(values)):
                    significant_lags.append(i)

        result['significant_lags'] = significant_lags

        # Detect period-3 periodicity (codon-based)
        period_3_lags = [i for i in significant_lags if i % 3 == 0 and i > 0]
        result['has_period_3'] = len(period_3_lags) > 0
        result['period_3_lags'] = period_3_lags

        return result

    if by_region_name:
        # Calculate separately for each individual region name
        region_names = sorted(results_df['region_name'].unique())

        for region_name in region_names:
            subset = results_df[results_df['region_name'] == region_name].copy()
            result = _calculate_acf_for_subset(subset, region_name, max_lag)
            if result is not None:
                results[region_name] = result

    elif by_region:
        # Calculate separately for overlapping and non-overlapping
        for region_type in ['overlapping', 'non-overlapping']:
            subset = results_df[results_df['region_type'] == region_type].copy()
            result = _calculate_acf_for_subset(subset, region_type, max_lag)
            if result is not None:
                results[region_type] = result

    else:
        # Calculate for entire dataset
        subset = results_df.copy()
        result = _calculate_acf_for_subset(subset, 'all', max_lag)
        if result is not None:
            results['all'] = result

    return results


def plot_autocorrelation(acf_results, figsize=None, title=None,
                         mark_period_3=True, output_path=None, show_title=True,
                         col_wrap=3):
    """
    Plot autocorrelation results to visualize periodic patterns.

    Args:
        acf_results (dict): Output from calculate_autocorrelation
        figsize (tuple, optional): Figure size. If None, automatically calculated based on
                                  number of subplots (default: None for auto-sizing)
        title (str, optional): Plot title (if None, auto-generates from metric name)
        mark_period_3 (bool): Whether to mark multiples of 3 (codon periodicity)
        output_path (str, optional): Path to save figure
        show_title (bool): Whether to display the title (default: True)
        col_wrap (int): Number of columns before wrapping to next row (default: 3)

    Returns:
        matplotlib.figure.Figure: The figure object
    """
    n_subplots = len(acf_results)

    # Calculate grid dimensions
    if n_subplots == 1:
        nrows, ncols = 1, 1
    else:
        ncols = min(n_subplots, col_wrap)
        nrows = int(np.ceil(n_subplots / ncols))

    # Auto-calculate figsize if not provided
    if figsize is None:
        width_per_plot = 5
        height_per_plot = 4
        figsize = (width_per_plot * ncols, height_per_plot * nrows)

    # Create subplots
    fig, axes = plt.subplots(nrows, ncols, figsize=figsize, sharey=True,
                             squeeze=False)
    axes = axes.flatten()

    for idx, (region_name, result) in enumerate(acf_results.items()):
        ax = axes[idx]

        lags = result['lags']
        acf_vals = result['acf']

        # Plot ACF as bars
        ax.bar(lags, acf_vals, width=0.8, color='steelblue', alpha=0.7)

        # Plot confidence interval
        if result['confint_lower'] is not None:
            ax.fill_between(lags,
                           result['confint_lower'],
                           result['confint_upper'],
                           alpha=0.2, color='gray', label='95% CI')

        # Add horizontal line at 0
        ax.axhline(y=0, color='black', linestyle='-', linewidth=0.8)

        # Mark period-3 lags if requested
        if mark_period_3:
            period_3_lags = [i for i in lags if i % 3 == 0 and i > 0]
            for lag in period_3_lags:
                ax.axvline(x=lag, color='red', linestyle='--', alpha=0.3, linewidth=1)

        # Formatting
        subtitle = f"{region_name.capitalize()} (n={result['n_positions']})"
        if result.get('has_period_3', False):
            subtitle += f"\nPeriod-3 detected at lags: {result['period_3_lags'][:5]}"  # Show first 5

        ax.set_title(subtitle, fontsize=11)
        ax.set_xlabel('Lag', fontsize=10)

        # Add y-label to leftmost plots in each row
        if idx % ncols == 0:
            ax.set_ylabel('Autocorrelation', fontsize=10)

        ax.grid(True, alpha=0.3, axis='y')

    # Hide unused subplots
    for idx in range(n_subplots, len(axes)):
        axes[idx].set_visible(False)

    if show_title:
        if title:
            fig.suptitle(title, fontsize=14, y=1.00)
        else:
            metric = list(acf_results.values())[0]['metric']
            fig.suptitle(f'Autocorrelation of {metric.capitalize()}', fontsize=14, y=1.00)

    sns.despine()
    plt.tight_layout()

    if output_path:
        plt.savefig(output_path, dpi=300, bbox_inches='tight')

    return fig


def plot_nucleotide_enrichment(composition_df, overlapping_regions=None,
                               figsize=(14, 6), title='Nucleotide Enrichment (log2)',
                               mark_regions=True, output_path=None, window_size=1, show_title=True):
    """
    Plot nucleotide enrichment scores showing which bases deviate from background.

    Args:
        composition_df (pd.DataFrame): Output from calculate_nucleotide_composition_enrichment
        overlapping_regions (list, optional): List of [name, start, end] for marking regions
        figsize (tuple): Figure size
        title (str): Plot title
        mark_regions (bool): Whether to mark overlapping regions with shading
        output_path (str, optional): Path to save figure
        window_size (int): Rolling window size for smoothing (default 1 = no smoothing)
        show_title (bool): Whether to display the title (default: True)

    Returns:
        matplotlib.figure.Figure: The figure object
    """
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=figsize, sharex=True)

    # Apply smoothing if requested
    if window_size > 1:
        for col in ['A_enrichment', 'C_enrichment', 'G_enrichment', 'T_enrichment']:
            composition_df[col] = composition_df[col].rolling(window=window_size, center=True).mean()

    # Plot individual base enrichments
    ax1.plot(composition_df['position'], composition_df['A_enrichment'],
            label='A', color='#1f77b4', linewidth=1.5, alpha=0.8)
    ax1.plot(composition_df['position'], composition_df['G_enrichment'],
            label='G', color='#ff7f0e', linewidth=1.5, alpha=0.8)
    ax1.plot(composition_df['position'], composition_df['C_enrichment'],
            label='C', color='#2ca02c', linewidth=1.5, alpha=0.8)
    ax1.plot(composition_df['position'], composition_df['T_enrichment'],
            label='T', color='#d62728', linewidth=1.5, alpha=0.8)

    ax1.axhline(y=0, color='black', linestyle='--', alpha=0.5, linewidth=1)
    ax1.set_ylabel('Enrichment\n(log2 obs/bg)', fontsize=11)
    if show_title:
        ax1.set_title(title, fontsize=13)
    ax1.legend(loc='upper right', ncol=4, fontsize=9)
    ax1.grid(True, alpha=0.3)

    # Plot purine vs pyrimidine
    ax2.plot(composition_df['position'], composition_df['purine_freq'],
            label='Purine (A+G)', color='steelblue', linewidth=2, alpha=0.8)
    ax2.plot(composition_df['position'], composition_df['pyrimidine_freq'],
            label='Pyrimidine (C+T)', color='coral', linewidth=2, alpha=0.8)

    ax2.axhline(y=0.5, color='gray', linestyle='--', alpha=0.5, linewidth=1)
    ax2.set_xlabel('Position', fontsize=11)
    ax2.set_ylabel('Frequency', fontsize=11)
    ax2.legend(loc='upper right', fontsize=9)
    ax2.grid(True, alpha=0.3)

    # Mark overlapping regions on both subplots
    if mark_regions and overlapping_regions:
        colors = plt.cm.Set3(np.linspace(0, 1, len(overlapping_regions)))
        texts = []
        for i, region in enumerate(overlapping_regions):
            region_name, start, end = region
            # Add to first subplot with label (will appear in legend)
            ax1.axvspan(start, end, alpha=0.15, color=colors[i], label=region_name)
            # Add to second subplot without label (to avoid duplicate legend entries)
            ax2.axvspan(start, end, alpha=0.15, color=colors[i])

            # Add text annotation at the midpoint
            mid_pos = (start + end) / 2
            # Get y-position for text (top of the plot)
            y_max = ax1.get_ylim()[1]
            text_obj = ax1.text(mid_pos, y_max, region_name, ha='center', va='bottom',
                               fontsize=9, color='black', fontweight='bold')
            texts.append(text_obj)

        # Adjust text positions to prevent overlap
        if texts:
            adjust_text(texts, ax=ax1, only_move={'points': 'y', 'text': 'y'})

    sns.despine()
    plt.tight_layout()

    if output_path:
        plt.savefig(output_path, dpi=300, bbox_inches='tight')

    return fig


def plot_per_position_conservation(conservation_df, overlapping_regions=None,
    figsize=(14, 4), color='tab:blue',
    title='Per-Position Conservation',
    xlabel='Position', ylabel='Conservation Score',
    mark_regions=True, output_path=None, window_size=1, show_title=True):
    """
    Plot per-position conservation scores.

    Args:
        conservation_df (pd.DataFrame): Output from calculate_per_position_conservation
        overlapping_regions (list, optional): List of [name, start, end] for marking regions
        figsize (tuple): Figure size
        color (str): Line color
        title (str): Plot title
        xlabel (str): X-axis label
        ylabel (str): Y-axis label
        mark_regions (bool): Whether to mark overlapping regions with shading
        output_path (str, optional): Path to save figure
        window_size (int): Rolling window size for smoothing (default 1 = no smoothing)
        show_title (bool): Whether to display the title (default: True)

    Returns:
        matplotlib.figure.Figure: The figure object
    """
    fig, ax = plt.subplots(figsize=figsize)

    # Apply smoothing if requested
    conservation_values = conservation_df['conservation'].copy()
    if window_size > 1:
        conservation_values = conservation_values.rolling(window=window_size, center=True).mean()

    # Plot conservation scores
    ax.plot(conservation_df['position'], conservation_values,
            color=color, linewidth=1.5, alpha=0.8)

    # Mark overlapping regions
    texts = []
    if mark_regions and overlapping_regions:
        colors = plt.cm.Set3(np.linspace(0, 1, len(overlapping_regions)))
        for i, region in enumerate(overlapping_regions):
            region_name, start, end = region
            ax.axvspan(start, end, alpha=0.2, color=colors[i], label=region_name)

            # Add text annotation at the midpoint
            mid_pos = (start + end) / 2
            y_max = ax.get_ylim()[1]
            text_obj = ax.text(mid_pos, y_max, region_name, ha='center', va='bottom',
                              fontsize=9, color='black', fontweight='bold')
            texts.append(text_obj)

    ax.set_xlabel(xlabel, fontsize=12)
    ax.set_ylabel(ylabel, fontsize=12)
    if show_title:
        ax.set_title(title, fontsize=14)
    ax.grid(True, alpha=0.3, linestyle='--')

    if mark_regions and overlapping_regions:
        ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left', fontsize=9)

    # Adjust text positions to prevent overlap
    if texts:
        adjust_text(texts, ax=ax, only_move={'points': 'y', 'text': 'y'})

    sns.despine()
    plt.tight_layout()

    if output_path:
        plt.savefig(output_path, dpi=300, bbox_inches='tight')

    return fig


def plot_per_position_diversity(diversity_df, overlapping_regions=None,
                                figsize=(14, 4), color='tab:orange',
                                title='Per-Position Nucleotide Diversity',
                                xlabel='Position', ylabel='Nucleotide Diversity (π)',
                                mark_regions=True, output_path=None, window_size=1, show_title=True):
    """
    Plot per-position nucleotide diversity.

    Args:
        diversity_df (pd.DataFrame): Output from calculate_per_position_diversity
        overlapping_regions (list, optional): List of [name, start, end] for marking regions
        figsize (tuple): Figure size
        color (str): Line color
        title (str): Plot title
        xlabel (str): X-axis label
        ylabel (str): Y-axis label
        mark_regions (bool): Whether to mark overlapping regions with shading
        output_path (str, optional): Path to save figure
        window_size (int): Rolling window size for smoothing (default 1 = no smoothing)
        show_title (bool): Whether to display the title (default: True)

    Returns:
        matplotlib.figure.Figure: The figure object
    """
    fig, ax = plt.subplots(figsize=figsize)

    # Apply smoothing if requested
    diversity_values = diversity_df['diversity'].copy()
    if window_size > 1:
        diversity_values = diversity_values.rolling(window=window_size, center=True).mean()

    # Plot diversity scores
    ax.plot(diversity_df['position'], diversity_values,
            color=color, linewidth=1.5, alpha=0.8)

    # Mark overlapping regions
    texts = []
    if mark_regions and overlapping_regions:
        colors = plt.cm.Set3(np.linspace(0, 1, len(overlapping_regions)))
        for i, region in enumerate(overlapping_regions):
            region_name, start, end = region
            ax.axvspan(start, end, alpha=0.2, color=colors[i], label=region_name)

            # Add text annotation at the midpoint
            mid_pos = (start + end) / 2
            y_max = ax.get_ylim()[1]
            text_obj = ax.text(mid_pos, y_max, region_name, ha='center', va='bottom',
                              fontsize=9, color='black', fontweight='bold')
            texts.append(text_obj)

    ax.set_xlabel(xlabel, fontsize=12)
    ax.set_ylabel(ylabel, fontsize=12)
    if show_title:
        ax.set_title(title, fontsize=14)
    ax.grid(True, alpha=0.3, linestyle='--')

    if mark_regions and overlapping_regions:
        ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left', fontsize=9)

    # Adjust text positions to prevent overlap
    if texts:
        adjust_text(texts, ax=ax, only_move={'points': 'y', 'text': 'y'})

    sns.despine()
    plt.tight_layout()

    if output_path:
        plt.savefig(output_path, dpi=300, bbox_inches='tight')

    return fig


def plot_per_position_kl_divergence(kl_df, overlapping_regions=None,
                                    figsize=(14, 4), color='tab:green',
                                    title='Per-Position KL Divergence from Background',
                                    xlabel='Position', ylabel='KL Divergence (bits)',
                                    mark_regions=True, output_path=None, window_size=1, show_title=True):
    """
    Plot per-position KL divergence scores.

    Args:
        kl_df (pd.DataFrame): Output from calculate_per_position_kl_divergence
        overlapping_regions (list, optional): List of [name, start, end] for marking regions
        figsize (tuple): Figure size
        color (str): Line color
        title (str): Plot title
        xlabel (str): X-axis label
        ylabel (str): Y-axis label
        mark_regions (bool): Whether to mark overlapping regions with shading
        output_path (str, optional): Path to save figure
        window_size (int): Rolling window size for smoothing (default 1 = no smoothing)
        show_title (bool): Whether to display the title (default: True)

    Returns:
        matplotlib.figure.Figure: The figure object
    """
    fig, ax = plt.subplots(figsize=figsize)

    # Apply smoothing if requested
    kl_values = kl_df['kl_divergence'].copy()
    if window_size > 1:
        kl_values = kl_values.rolling(window=window_size, center=True).mean()

    # Plot KL divergence scores
    ax.plot(kl_df['position'], kl_values,
            color=color, linewidth=1.5, alpha=0.8)

    # Mark overlapping regions
    texts = []
    if mark_regions and overlapping_regions:
        colors = plt.cm.Set3(np.linspace(0, 1, len(overlapping_regions)))
        for i, region in enumerate(overlapping_regions):
            region_name, start, end = region
            ax.axvspan(start, end, alpha=0.2, color=colors[i], label=region_name)

            # Add text annotation at the midpoint
            mid_pos = (start + end) / 2
            y_max = ax.get_ylim()[1]
            text_obj = ax.text(mid_pos, y_max, region_name, ha='center', va='bottom',
                              fontsize=9, color='black', fontweight='bold')
            texts.append(text_obj)

    ax.set_xlabel(xlabel, fontsize=12)
    ax.set_ylabel(ylabel, fontsize=12)
    if show_title:
        ax.set_title(title, fontsize=14)
    ax.grid(True, alpha=0.3, linestyle='--')

    # Add horizontal line at KL=0 (matches background perfectly)
    ax.axhline(y=0, color='gray', linestyle='--', alpha=0.5, linewidth=1)

    if mark_regions and overlapping_regions:
        ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left', fontsize=9)

    # Adjust text positions to prevent overlap
    if texts:
        adjust_text(texts, ax=ax, only_move={'points': 'y', 'text': 'y'})

    sns.despine()
    plt.tight_layout()

    if output_path:
        plt.savefig(output_path, dpi=300, bbox_inches='tight')

    return fig


def compare_regional_distributions(results_df, metric_col='conservation',
    test='mann-whitney', figsize=(10, 6),
    palette='Set2', title=None,
    ylabel=None, output_path=None,
    show_pvalues=True, correction_method='fdr_bh',
    ylim='auto', plot_type='violin', show_points=True, show_title=True):
    """
    Compare distributions of conservation or diversity between overlapping and
    non-overlapping regions, with statistical annotations and multiple testing correction.

    Args:
        results_df (pd.DataFrame): Output from calculate_per_position_conservation
                                   or calculate_per_position_diversity
        metric_col (str): Column name to compare ('conservation' or 'diversity')
        test (str): Statistical test to use ('mann-whitney' or 't-test')
        figsize (tuple): Figure size
        palette (str): Seaborn color palette
        title (str, optional): Plot title
        ylabel (str, optional): Y-axis label
        output_path (str, optional): Path to save figure
        show_pvalues (bool): Whether to show p-values as asterisks
        correction_method (str): Multiple testing correction method. Options:
                                'fdr_bh' (Benjamini-Hochberg FDR, default),
                                'bonferroni', 'fdr_by' (Benjamini-Yekutieli),
                                'holm', 'sidak', or None for no correction
        ylim (tuple or str): Y-axis limits. Use 'auto' to automatically zoom to data range,
                            tuple (min, max) for manual limits, or None for default (0, 1).
        plot_type (str): Type of plot - 'box', 'violin', or 'both' (default: 'violin')
        show_points (bool): Whether to overlay individual data points (default: True)
        show_title (bool): Whether to display the title (default: True)

    Returns:
        tuple: (fig, stats_results) - Figure object and dictionary of statistical results
    """
    # Remove NaN values
    df_clean = results_df.dropna(subset=[metric_col])

    # Set default title and ylabel if not provided
    if title is None:
        title = f'Distribution of {metric_col.capitalize()} by Region'
    if ylabel is None:
        ylabel = metric_col.capitalize()

    # Create figure
    fig, ax = plt.subplots(figsize=figsize)

    # Create ordering for regions
    order = ['non-overlapping'] + sorted([name for name in df_clean['region_name'].unique()
                                          if name != 'non-overlapping'])

    # Create the requested plot type(s)
    if plot_type == 'box':
        sns.boxplot(data=df_clean, x='region_name', y=metric_col, hue='region_name',
                   order=order, palette=palette, legend=False, ax=ax)
    elif plot_type == 'violin':
        sns.violinplot(data=df_clean, x='region_name', y=metric_col, hue='region_name',
                      order=order, palette=palette, legend=False, ax=ax, inner='quartile')
    elif plot_type == 'both':
        # Violin plot with boxplot overlay
        sns.violinplot(data=df_clean, x='region_name', y=metric_col, hue='region_name',
                      order=order, palette=palette, legend=False, ax=ax, inner=None, alpha=0.6)
        sns.boxplot(data=df_clean, x='region_name', y=metric_col,
                   order=order, ax=ax, width=0.3, boxprops=dict(alpha=0.7),
                   showcaps=False, whiskerprops=dict(alpha=0.7),
                   flierprops=dict(alpha=0.0))  # Hide outliers from boxplot
    else:
        raise ValueError(f"Unknown plot_type: {plot_type}. Use 'box', 'violin', or 'both'")

    # Overlay individual data points if requested
    if show_points:
        sns.stripplot(data=df_clean, x='region_name', y=metric_col,
                     order=order, color='black', alpha=0.3, size=2, ax=ax, jitter=True)

    ax.set_xlabel('Region', fontsize=12)
    ax.set_ylabel(ylabel, fontsize=12)
    if show_title:
        ax.set_title(title, fontsize=14)
    plt.xticks(rotation=45, ha='right')

    # Apply y-axis limits
    if ylim == 'auto':
        # Automatically calculate smart limits with some padding
        data_min = df_clean[metric_col].min()
        data_max = df_clean[metric_col].max()
        data_range = data_max - data_min
        padding = data_range * 0.1 if data_range > 0 else 0.05
        ax.set_ylim(data_min - padding, data_max + padding)
    elif ylim is not None:
        ax.set_ylim(ylim)

    # Add statistical annotations if requested
    stats_results = {}
    if show_pvalues and len(order) > 1:
        # Prepare pairs for comparison (each overlapping region vs non-overlapping)
        pairs = []
        for region in order[1:]:  # Skip 'non-overlapping'
            pairs.append(('non-overlapping', region))

        # Perform statistical tests
        if test == 'mann-whitney':
            test_func = stats.mannwhitneyu
            test_name = 'Mann-Whitney U'
        elif test == 't-test':
            test_func = stats.ttest_ind
            test_name = 'Independent t-test'
        else:
            raise ValueError(f"Unknown test: {test}. Use 'mann-whitney' or 't-test'")

        # Calculate p-values
        formatted_pairs = []
        pvalues = []
        for pair in pairs:
            group1 = df_clean[df_clean['region_name'] == pair[0]][metric_col]
            group2 = df_clean[df_clean['region_name'] == pair[1]][metric_col]

            if len(group1) > 0 and len(group2) > 0:
                if test == 'mann-whitney':
                    stat, pval = test_func(group1, group2, alternative='two-sided')
                else:
                    stat, pval = test_func(group1, group2)

                # Store initial results with raw p-values
                stats_results[f"{pair[1]}_vs_non-overlapping"] = {
                    'test': test_name,
                    'statistic': stat,
                    'pvalue_raw': pval,
                    'n1': len(group1),
                    'n2': len(group2),
                    'mean1': group1.mean(),
                    'mean2': group2.mean()
                }

                formatted_pairs.append(pair)
                pvalues.append(pval)

        # Apply multiple testing correction
        if len(pvalues) > 0 and correction_method is not None:
            reject, pvals_corrected, _, _ = multipletests(pvalues, method=correction_method)

            # Update stats_results with corrected p-values
            for i, pair in enumerate(formatted_pairs):
                comparison_key = f"{pair[1]}_vs_non-overlapping"
                stats_results[comparison_key]['pvalue_corrected'] = pvals_corrected[i]
                stats_results[comparison_key]['significant'] = reject[i]
                stats_results[comparison_key]['correction_method'] = correction_method

            # Use corrected p-values for annotations
            pvalues_to_plot = pvals_corrected
        else:
            # No correction - use raw p-values
            pvalues_to_plot = pvalues
            for pair in formatted_pairs:
                comparison_key = f"{pair[1]}_vs_non-overlapping"
                stats_results[comparison_key]['pvalue_corrected'] = stats_results[comparison_key]['pvalue_raw']
                stats_results[comparison_key]['significant'] = stats_results[comparison_key]['pvalue_raw'] < 0.05
                stats_results[comparison_key]['correction_method'] = 'none'

        # Add annotations using Annotator
        if len(formatted_pairs) > 0:
            annotator = Annotator(ax, formatted_pairs, data=df_clean,
                                 x='region_name', y=metric_col, order=order)
            annotator.configure(test=None, text_format='star',
                              loc='outside', verbose=False)
            annotator.set_pvalues(pvalues_to_plot)
            annotator.annotate()

    sns.despine()
    plt.tight_layout()

    if output_path:
        plt.savefig(output_path, dpi=300, bbox_inches='tight')

    return fig, stats_results


def plot_split_view_comparison(results_df, metric_col='conservation',
                               figsize=(14, 6), palette='Set2',
                               invariant_threshold=None,
                               test='mann-whitney', show_pvalues=True,
                               correction_method='fdr_bh',
                               plot_type='violin', show_points=True,
                               show_title=True, title=None, output_path=None):
    """
    Create a split-view comparison plot that separates invariant sites from variable sites.

    This addresses distributions that are heavily centered (conservation near 1.0, diversity near 0.0)
    by showing:
    1. Left panel: % of invariant sites (frequency of perfect conservation/no diversity)
    2. Right panel: Distribution of variable sites only (severity when mutations occur)

    Args:
        results_df (pd.DataFrame): Output from calculate_per_position_* functions
        metric_col (str): Column name to analyze ('conservation' or 'diversity')
        figsize (tuple): Figure size (default: (14, 6))
        palette (str): Seaborn color palette
        invariant_threshold (float, optional): Threshold for defining invariant sites.
                                              If None, auto-detects based on metric:
                                              - conservation: >= 0.999
                                              - diversity: <= 0.001
        test (str): Statistical test for variable sites comparison ('mann-whitney' or 't-test')
        show_pvalues (bool): Whether to show p-values on variable sites plot
        correction_method (str): Multiple testing correction method
        plot_type (str): Plot type for variable sites ('box', 'violin', or 'both')
        show_points (bool): Whether to overlay individual points on variable sites plot
        show_title (bool): Whether to display titles
        title (str, optional): Overall figure title
        output_path (str, optional): Path to save figure

    Returns:
        tuple: (fig, stats_dict) where stats_dict contains:
               - 'invariant_percentages': DataFrame with % invariant per region
               - 'variable_stats': Statistical comparison results for variable sites
               - 'n_variable': Number of variable sites per region
    """
    # Auto-detect invariant threshold based on metric
    if invariant_threshold is None:
        if metric_col == 'conservation':
            invariant_threshold = 0.999  # Sites with conservation >= 0.999
            is_high_centered = True  # Centered near 1.0
        elif metric_col == 'diversity':
            invariant_threshold = 0.001  # Sites with diversity <= 0.001
            is_high_centered = False  # Centered near 0.0
        else:
            raise ValueError(f"metric_col must be 'conservation' or 'diversity', got '{metric_col}'")
    else:
        # If user provides threshold, infer centering from metric name
        is_high_centered = (metric_col == 'conservation')

    # Remove NaN values
    df_clean = results_df.dropna(subset=[metric_col])

    # Create ordering for regions
    order = ['non-overlapping'] + sorted([name for name in df_clean['region_name'].unique()
                                          if name != 'non-overlapping'])

    # Identify invariant vs variable sites
    if is_high_centered:
        # For conservation: invariant sites have value >= threshold
        df_clean['site_type'] = df_clean[metric_col].apply(
            lambda x: 'invariant' if x >= invariant_threshold else 'variable'
        )
    else:
        # For diversity: invariant sites have value <= threshold
        df_clean['site_type'] = df_clean[metric_col].apply(
            lambda x: 'invariant' if x <= invariant_threshold else 'variable'
        )

    # Calculate % invariant sites per region
    invariant_percentages = []
    n_variable_dict = {}
    for region in order:
        region_data = df_clean[df_clean['region_name'] == region]
        n_total = len(region_data)
        n_invariant = len(region_data[region_data['site_type'] == 'invariant'])
        n_variable = len(region_data[region_data['site_type'] == 'variable'])
        pct_invariant = (n_invariant / n_total * 100) if n_total > 0 else 0

        invariant_percentages.append({
            'region_name': region,
            'n_total': n_total,
            'n_invariant': n_invariant,
            'n_variable': n_variable,
            'pct_invariant': pct_invariant
        })
        n_variable_dict[region] = n_variable

    invariant_df = pd.DataFrame(invariant_percentages)

    # Create figure with two subplots
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=figsize)

    # Get colors from palette
    colors = sns.color_palette(palette, n_colors=len(order))
    color_dict = {region: colors[i] for i, region in enumerate(order)}

    # ===== LEFT PANEL: Bar chart of % invariant sites =====
    bars = ax1.bar(range(len(order)), invariant_df['pct_invariant'],
                   color=[color_dict[r] for r in invariant_df['region_name']],
                   alpha=0.7, edgecolor='black', linewidth=1)

    ax1.set_xticks(range(len(order)))
    ax1.set_xticklabels(order, rotation=45, ha='right')
    ax1.set_xlabel('Region', fontsize=12)
    ax1.set_ylabel('% Invariant Sites', fontsize=12)

    if show_title:
        if is_high_centered:
            ax1.set_title(f'Frequency: % Perfectly Conserved\n(≥{invariant_threshold:.3f})', fontsize=11)
        else:
            ax1.set_title(f'Frequency: % Perfectly Conserved\n(≤{invariant_threshold:.3f})', fontsize=11)

    # Add percentage labels on bars
    for i, bar in enumerate(bars):
        height = bar.get_height()
        ax1.text(bar.get_x() + bar.get_width()/2., height,
                f'{height:.1f}%',
                ha='center', va='bottom', fontsize=9)

    ax1.set_ylim(0, 100)
    ax1.grid(True, alpha=0.3, axis='y')

    # ===== RIGHT PANEL: Distribution of variable sites only =====
    variable_sites = df_clean[df_clean['site_type'] == 'variable'].copy()

    if len(variable_sites) > 0:
        # Create the requested plot type
        if plot_type == 'box':
            sns.boxplot(data=variable_sites, x='region_name', y=metric_col, hue='region_name',
                       order=order, palette=palette, legend=False, ax=ax2)
        elif plot_type == 'violin':
            sns.violinplot(data=variable_sites, x='region_name', y=metric_col, hue='region_name',
                          order=order, palette=palette, legend=False, ax=ax2, inner='quartile')
        elif plot_type == 'both':
            sns.violinplot(data=variable_sites, x='region_name', y=metric_col, hue='region_name',
                          order=order, palette=palette, legend=False, ax=ax2, inner=None, alpha=0.6)
            sns.boxplot(data=variable_sites, x='region_name', y=metric_col,
                       order=order, ax=ax2, width=0.3, boxprops=dict(alpha=0.7),
                       showcaps=False, whiskerprops=dict(alpha=0.7),
                       flierprops=dict(alpha=0.0))

        # Overlay individual data points if requested
        if show_points:
            sns.stripplot(data=variable_sites, x='region_name', y=metric_col,
                         order=order, color='black', alpha=0.3, size=2, ax=ax2, jitter=True)

        ax2.set_xlabel('Region', fontsize=12)
        ylabel = 'Conservation Score' if metric_col == 'conservation' else 'Nucleotide Diversity (π)'
        ax2.set_ylabel(ylabel, fontsize=12)

        if show_title:
            ax2.set_title(f'Severity: Variable Sites Only\n(n={len(variable_sites)} positions)', fontsize=11)

        ax2.set_xticks(range(len(order)))
        ax2.set_xticklabels(order, rotation=45, ha='right')
        ax2.grid(True, alpha=0.3, axis='y')

        # Add statistical annotations for variable sites
        variable_stats = {}
        if show_pvalues and len(order) > 1:
            # Filter order to only include regions present in variable_sites
            regions_with_variable_sites = variable_sites['region_name'].unique().tolist()
            filtered_order = [r for r in order if r in regions_with_variable_sites]

            # Prepare pairs for comparison (only regions with variable sites)
            pairs = []
            if 'non-overlapping' in filtered_order:
                for region in filtered_order:
                    if region != 'non-overlapping':
                        pairs.append(('non-overlapping', region))

            # Perform statistical tests
            if test == 'mann-whitney':
                test_func = stats.mannwhitneyu
                test_name = 'Mann-Whitney U'
            else:
                test_func = stats.ttest_ind
                test_name = 'Independent t-test'

            formatted_pairs = []
            pvalues = []

            for pair in pairs:
                group1 = variable_sites[variable_sites['region_name'] == pair[0]][metric_col]
                group2 = variable_sites[variable_sites['region_name'] == pair[1]][metric_col]

                if len(group1) > 0 and len(group2) > 0:
                    if test == 'mann-whitney':
                        stat, pval = test_func(group1, group2, alternative='two-sided')
                    else:
                        stat, pval = test_func(group1, group2)

                    variable_stats[f"{pair[1]}_vs_non-overlapping"] = {
                        'test': test_name,
                        'statistic': stat,
                        'pvalue_raw': pval,
                        'n1': len(group1),
                        'n2': len(group2),
                        'mean1': group1.mean(),
                        'mean2': group2.mean()
                    }

                    formatted_pairs.append(pair)
                    pvalues.append(pval)

            # Apply multiple testing correction
            if len(pvalues) > 0 and correction_method is not None:
                reject, pvals_corrected, _, _ = multipletests(pvalues, method=correction_method)

                for i, pair in enumerate(formatted_pairs):
                    comparison_key = f"{pair[1]}_vs_non-overlapping"
                    variable_stats[comparison_key]['pvalue_corrected'] = pvals_corrected[i]
                    variable_stats[comparison_key]['significant'] = reject[i]
                    variable_stats[comparison_key]['correction_method'] = correction_method

                pvalues_to_plot = pvals_corrected
            else:
                pvalues_to_plot = pvalues
                for pair in formatted_pairs:
                    comparison_key = f"{pair[1]}_vs_non-overlapping"
                    variable_stats[comparison_key]['pvalue_corrected'] = variable_stats[comparison_key]['pvalue_raw']
                    variable_stats[comparison_key]['significant'] = variable_stats[comparison_key]['pvalue_raw'] < 0.05
                    variable_stats[comparison_key]['correction_method'] = 'none'

            # Add annotations (using filtered order)
            if len(formatted_pairs) > 0:
                annotator = Annotator(ax2, formatted_pairs, data=variable_sites,
                                     x='region_name', y=metric_col, order=filtered_order)
                annotator.configure(test=None, text_format='star',
                                  loc='outside', verbose=False)
                annotator.set_pvalues(pvalues_to_plot)
                annotator.annotate()
    else:
        ax2.text(0.5, 0.5, 'No variable sites found',
                ha='center', va='center', transform=ax2.transAxes, fontsize=14)
        variable_stats = {}

    # Overall title
    if show_title and title:
        fig.suptitle(title, fontsize=14, y=1.02)

    sns.despine()
    plt.tight_layout()

    if output_path:
        plt.savefig(output_path, dpi=300, bbox_inches='tight')

    # Return statistics
    stats_dict = {
        'invariant_percentages': invariant_df,
        'variable_stats': variable_stats,
        'n_variable': n_variable_dict,
        'threshold': invariant_threshold
    }

    return fig, stats_dict


def plot_regional_density_comparison(results_df, metric_col='conservation',
                                     figsize=(12, 5), palette='Set2',
                                     title=None, xlabel=None,
                                     output_path=None, plot_type='kde',
                                     xlim='auto', bins=50, show_title=True):
    """
    Create overlaid density or histogram plots to compare distributions between regions.
    This is particularly useful when distributions are very similar and tightly clustered.

    Args:
        results_df (pd.DataFrame): Output from calculate_per_position_* functions
        metric_col (str): Column name to plot ('conservation', 'diversity', etc.)
        figsize (tuple): Figure size
        palette (str): Seaborn color palette
        title (str, optional): Plot title
        xlabel (str, optional): X-axis label
        output_path (str, optional): Path to save figure
        plot_type (str): 'kde' for kernel density, 'hist' for histogram, 'both' for both
        xlim (tuple or str): X-axis limits. Use 'auto' to zoom to data range, or tuple (min, max)
        bins (int): Number of bins for histogram (default 50)
        show_title (bool): Whether to display the title (default: True)

    Returns:
        matplotlib.figure.Figure: The figure object
    """
    # Remove NaN values
    df_clean = results_df.dropna(subset=[metric_col])

    # Set default title and xlabel if not provided
    if title is None:
        title = f'{metric_col.capitalize()} Distribution by Region'
    if xlabel is None:
        xlabel = metric_col.capitalize()

    # Create ordering for regions
    order = ['non-overlapping'] + sorted([name for name in df_clean['region_name'].unique()
                                          if name != 'non-overlapping'])

    # Create figure
    if plot_type == 'both':
        fig, (ax1, ax2) = plt.subplots(2, 1, figsize=figsize, sharex=True)
    else:
        fig, ax1 = plt.subplots(figsize=figsize)
        ax2 = None

    # Get colors from palette
    colors = sns.color_palette(palette, n_colors=len(order))

    # Plot KDE (kernel density estimate)
    if plot_type in ['kde', 'both']:
        for idx, region in enumerate(order):
            subset = df_clean[df_clean['region_name'] == region][metric_col]
            if len(subset) > 1:
                subset.plot(kind='density', ax=ax1, label=region, color=colors[idx],
                           linewidth=2, alpha=0.7)

        ax1.set_ylabel('Density', fontsize=11)
        ax1.set_xlabel(xlabel if plot_type == 'kde' else '', fontsize=11)
        if show_title:
            ax1.set_title(f'{title} - Kernel Density' if plot_type == 'both' else title,
                         fontsize=13)
        ax1.legend(loc='best', fontsize=9)
        ax1.grid(True, alpha=0.3, axis='y')

    # Plot histogram
    if plot_type in ['hist', 'both']:
        ax_hist = ax2 if ax2 is not None else ax1

        for idx, region in enumerate(order):
            subset = df_clean[df_clean['region_name'] == region][metric_col]
            ax_hist.hist(subset, bins=bins, alpha=0.5, label=region,
                        color=colors[idx], edgecolor='black', linewidth=0.5)

        ax_hist.set_ylabel('Count', fontsize=11)
        ax_hist.set_xlabel(xlabel, fontsize=11)
        if show_title:
            if plot_type == 'both':
                ax_hist.set_title(f'{title} - Histogram', fontsize=13)
            elif plot_type == 'hist':
                ax_hist.set_title(title, fontsize=13)
        ax_hist.legend(loc='best', fontsize=9)
        ax_hist.grid(True, alpha=0.3, axis='y')

    # Apply x-axis limits
    if xlim == 'auto':
        data_min = df_clean[metric_col].min()
        data_max = df_clean[metric_col].max()
        data_range = data_max - data_min
        padding = data_range * 0.05 if data_range > 0 else 0.02
        xlim_vals = (data_min - padding, data_max + padding)
        ax1.set_xlim(xlim_vals)
        if ax2 is not None:
            ax2.set_xlim(xlim_vals)
    elif xlim is not None:
        ax1.set_xlim(xlim)
        if ax2 is not None:
            ax2.set_xlim(xlim)

    sns.despine()
    plt.tight_layout()

    if output_path:
        plt.savefig(output_path, dpi=300, bbox_inches='tight')

    return fig


def analyze_gene_regions(df, gene_col, overlapping_regions,
                         calculate_conservation=True, calculate_diversity=True,
                         calculate_kl_divergence=True, calculate_composition=True,
                         run_autocorrelation=True,
                         plot_results=True, compare_distributions=True,
                         plot_density_comparison=False,
                         plot_split_view=False,
                         density_plot_type='kde',
                         autocorr_max_lag=50,
                         autocorr_by_region_name=False,
                         output_dir=None,
                         show_points = True,
                         show_title = True,
                        **plot_kwargs):
    """
    Complete analysis pipeline for gene regions: calculate conservation, diversity,
    KL divergence, nucleotide composition, autocorrelation, plot results, and compare distributions.

    Args:
        df (pd.DataFrame): DataFrame with sequences (one row per individual)
        gene_col (str): Column name containing the aligned gene sequences
        overlapping_regions (list): List of [name, start, end] for overlapping regions
        calculate_conservation (bool): Whether to calculate conservation
        calculate_diversity (bool): Whether to calculate diversity
        calculate_kl_divergence (bool): Whether to calculate KL divergence
        calculate_composition (bool): Whether to calculate nucleotide composition enrichment
        run_autocorrelation (bool): Whether to calculate autocorrelation
        plot_results (bool): Whether to create plots
        compare_distributions (bool): Whether to compare region distributions (violin/box plots)
        plot_density_comparison (bool): Whether to create density/histogram comparison plots (default: False)
        plot_split_view (bool): Whether to create split-view plots separating invariant from variable sites
                                (recommended for conservation/diversity). Only applies to conservation and diversity.
        density_plot_type (str): Type of density plot - 'kde', 'hist', or 'both' (default: 'kde')
        autocorr_max_lag (int): Maximum lag for autocorrelation (default 50). Automatically adjusted
                               for short regions to min(max_lag, region_length // 2).
        autocorr_by_region_name (bool): Whether to calculate autocorrelation separately for each
                                       individual overlapping region (default: False, which groups
                                       all overlapping regions together vs non-overlapping)
        output_dir (str, optional): Directory to save output files
        show_points (bool): Whether to show individual data points on distribution plots
        show_title (bool): Whether to show titles on plots
        **plot_kwargs: Additional plotting arguments

    Returns:
        dict: Dictionary containing all results (DataFrames, figures, statistics)
    """
    results = {}
    # Create dir path if needed
    if output_dir:
        os.makedirs(output_dir, exist_ok=True)
    # Calculate conservation
    if calculate_conservation:
        print(f"Calculating per-position conservation for {gene_col}...")
        conservation_df = calculate_per_position_conservation(df, gene_col, overlapping_regions)
        results['conservation_df'] = conservation_df

        if plot_results:
            print("Plotting conservation...")
            output_path = f"{output_dir}/{gene_col}_conservation.png" if output_dir else None
            fig = plot_per_position_conservation(conservation_df, overlapping_regions,
                                                 output_path=output_path, show_title=show_title, **plot_kwargs)
            results['conservation_fig'] = fig

        if compare_distributions:
            print("Comparing conservation distributions...")
            output_path = f"{output_dir}/{gene_col}_conservation_comparison.png" if output_dir else None
            fig, stats = compare_regional_distributions(conservation_df,
                                                       metric_col='conservation',
                                                       output_path=output_path, show_title=show_title, show_points=show_points)
            results['conservation_comparison_fig'] = fig
            results['conservation_stats'] = stats

        if plot_density_comparison:
            print("Creating conservation density comparison...")
            output_path = f"{output_dir}/{gene_col}_conservation_density.png" if output_dir else None
            fig = plot_regional_density_comparison(conservation_df,
                                                   metric_col='conservation',
                                                   plot_type=density_plot_type,
                                                   output_path=output_path, show_title=show_title)
            results['conservation_density_fig'] = fig

        if plot_split_view:
            print("Creating conservation split-view comparison...")
            output_path = f"{output_dir}/{gene_col}_conservation_splitview.png" if output_dir else None
            fig, stats = plot_split_view_comparison(conservation_df,
                                                    metric_col='conservation',
                                                    show_pvalues=True,
                                                    show_points=show_points,
                                                    show_title=show_title,
                                                    output_path=output_path)
            results['conservation_splitview_fig'] = fig
            results['conservation_splitview_stats'] = stats

    # Calculate diversity
    if calculate_diversity:
        print(f"Calculating per-position diversity for {gene_col}...")
        diversity_df = calculate_per_position_diversity(df, gene_col, overlapping_regions)
        results['diversity_df'] = diversity_df

        if plot_results:
            print("Plotting diversity...")
            output_path = f"{output_dir}/{gene_col}_diversity.png" if output_dir else None
            fig = plot_per_position_diversity(diversity_df, overlapping_regions,
                                             output_path=output_path, show_title=show_title, **plot_kwargs)
            results['diversity_fig'] = fig

        if compare_distributions:
            print("Comparing diversity distributions...")
            output_path = f"{output_dir}/{gene_col}_diversity_comparison.png" if output_dir else None
            fig, stats = compare_regional_distributions(diversity_df,
                                                       metric_col='diversity',
                                                       output_path=output_path,
                                                       ylabel='Nucleotide Diversity (π)', show_title=show_title, show_points=show_points)
            results['diversity_comparison_fig'] = fig
            results['diversity_stats'] = stats

        if plot_density_comparison:
            print("Creating diversity density comparison...")
            output_path = f"{output_dir}/{gene_col}_diversity_density.png" if output_dir else None
            fig = plot_regional_density_comparison(diversity_df,
                                                   metric_col='diversity',
                                                   plot_type=density_plot_type,
                                                   xlabel='Nucleotide Diversity (π)',
                                                   output_path=output_path, show_title=show_title)
            results['diversity_density_fig'] = fig

        if plot_split_view:
            print("Creating diversity split-view comparison...")
            output_path = f"{output_dir}/{gene_col}_diversity_splitview.png" if output_dir else None
            fig, stats = plot_split_view_comparison(diversity_df,
                                                    metric_col='diversity',
                                                    show_pvalues=True,
                                                    show_points=show_points,
                                                    show_title=show_title,
                                                    output_path=output_path)
            results['diversity_splitview_fig'] = fig
            results['diversity_splitview_stats'] = stats

    # Calculate KL divergence
    if calculate_kl_divergence:
        print(f"Calculating per-position KL divergence for {gene_col}...")
        kl_df = calculate_per_position_kl_divergence(df, gene_col, overlapping_regions)
        results['kl_divergence_df'] = kl_df

        if plot_results:
            print("Plotting KL divergence...")
            output_path = f"{output_dir}/{gene_col}_kl_divergence.png" if output_dir else None
            fig = plot_per_position_kl_divergence(kl_df, overlapping_regions,
                                                  output_path=output_path, show_title=show_title, **plot_kwargs)
            results['kl_divergence_fig'] = fig

        if compare_distributions:
            print("Comparing KL divergence distributions...")
            output_path = f"{output_dir}/{gene_col}_kl_divergence_comparison.png" if output_dir else None
            fig, stats = compare_regional_distributions(kl_df,
                                                       metric_col='kl_divergence',
                                                       output_path=output_path,
                                                       ylabel='KL Divergence (bits)', show_title=show_title, show_points=show_points)
            results['kl_divergence_comparison_fig'] = fig
            results['kl_divergence_stats'] = stats

        if plot_density_comparison:
            print("Creating KL divergence density comparison...")
            output_path = f"{output_dir}/{gene_col}_kl_divergence_density.png" if output_dir else None
            fig = plot_regional_density_comparison(kl_df,
                                                   metric_col='kl_divergence',
                                                   plot_type=density_plot_type,
                                                   xlabel='KL Divergence (bits)',
                                                   output_path=output_path, show_title=show_title)
            results['kl_divergence_density_fig'] = fig

    # Calculate nucleotide composition enrichment
    if calculate_composition:
        print(f"Calculating nucleotide composition enrichment for {gene_col}...")
        composition_df = calculate_nucleotide_composition_enrichment(df, gene_col, overlapping_regions)
        results['composition_df'] = composition_df

        if plot_results:
            print("Plotting nucleotide enrichment...")
            output_path = f"{output_dir}/{gene_col}_composition_enrichment.png" if output_dir else None
            fig = plot_nucleotide_enrichment(composition_df, overlapping_regions,
                                             output_path=output_path, **plot_kwargs)
            results['composition_fig'] = fig

    # Calculate autocorrelation to detect periodic patterns
    if run_autocorrelation:
        print("Calculating autocorrelation for periodic pattern detection...")

        # Conservation autocorrelation
        if calculate_conservation and 'conservation_df' in results:
            print("  - Conservation autocorrelation...")
            acf_conservation = calculate_autocorrelation(
                results['conservation_df'],
                metric_col='conservation',
                max_lag=autocorr_max_lag,
                by_region=not autocorr_by_region_name,
                by_region_name=autocorr_by_region_name
            )
            results['autocorr_conservation'] = acf_conservation

            if plot_results:
                output_path = f"{output_dir}/{gene_col}_autocorr_conservation.png" if output_dir else None
                fig = plot_autocorrelation(acf_conservation,
                                          title=f'Autocorrelation of Conservation - {gene_col}',
                                          output_path=output_path, show_title=show_title)
                results['autocorr_conservation_fig'] = fig

        # Diversity autocorrelation
        if calculate_diversity and 'diversity_df' in results:
            print("  - Diversity autocorrelation...")
            acf_diversity = calculate_autocorrelation(
                results['diversity_df'],
                metric_col='diversity',
                max_lag=autocorr_max_lag,
                by_region=not autocorr_by_region_name,
                by_region_name=autocorr_by_region_name
            )
            results['autocorr_diversity'] = acf_diversity

            if plot_results:
                output_path = f"{output_dir}/{gene_col}_autocorr_diversity.png" if output_dir else None
                fig = plot_autocorrelation(acf_diversity,
                                          title=f'Autocorrelation of Diversity - {gene_col}',
                                          output_path=output_path, show_title=show_title)
                results['autocorr_diversity_fig'] = fig

        # KL divergence autocorrelation
        if calculate_kl_divergence and 'kl_divergence_df' in results:
            print("  - KL divergence autocorrelation...")
            acf_kl = calculate_autocorrelation(
                results['kl_divergence_df'],
                metric_col='kl_divergence',
                max_lag=autocorr_max_lag,
                by_region=not autocorr_by_region_name,
                by_region_name=autocorr_by_region_name
            )
            results['autocorr_kl_divergence'] = acf_kl

            if plot_results:
                output_path = f"{output_dir}/{gene_col}_autocorr_kl_divergence.png" if output_dir else None
                fig = plot_autocorrelation(acf_kl,
                                          title=f'Autocorrelation of KL Divergence - {gene_col}',
                                          output_path=output_path, show_title=show_title)
                results['autocorr_kl_divergence_fig'] = fig

    print("Analysis complete!")
    return results


def plot_mdp_invariant_sites_summary(stats_list, gene_names,
                                     figsize=(12, 6), palette='Set2',
                                     title='MDP Invariant Sites Comparison',
                                     significance_threshold=0.05,
                                     correction_method='fdr_bh',
                                     output_path=None, show_title=True,
                                     ylabel='% Difference in Invariant Sites\n(Overlapping - Non-overlapping)'):
    """
    Create a summary plot comparing invariant site percentages across multiple MDPs.

    Shows the percent difference in invariant sites between non-overlapping regions
    and each overlapping region, with statistical significance from Mann-Whitney tests
    comparing conservation score distributions (from the variable sites analysis).

    Args:
        stats_list (list): List of conservation_splitview_stats dictionaries from
                          plot_split_view_comparison, one per gene/MDP
        gene_names (list): List of gene names corresponding to each stats dict
        figsize (tuple): Figure size (default: (12, 6))
        palette (str): Seaborn color palette for gene colors
        title (str): Plot title
        significance_threshold (float): P-value threshold for marking significance (default: 0.05)
        correction_method (str): Multiple testing correction method. Options:
                                'fdr_bh' (Benjamini-Hochberg FDR, default),
                                'bonferroni', 'fdr_by' (Benjamini-Yekutieli),
                                'holm', 'sidak', or None for no correction
        output_path (str, optional): Path to save figure
        show_title (bool): Whether to display the title
        ylabel (str): Y-axis label

    Returns:
        tuple: (fig, summary_df) where summary_df contains:
               - gene: Gene/MDP name
               - region: Overlapping region name
               - non_overlap_pct: % invariant sites in non-overlapping region
               - overlap_pct: % invariant sites in overlapping region
               - pct_difference: Absolute difference (overlap - non_overlap)
               - mann_whitney_statistic: Mann-Whitney U test statistic (from variable_stats)
               - pvalue_raw: Raw p-value from Mann-Whitney test
               - pvalue_corrected: Multiple testing corrected p-value
               - significant: Boolean indicating significance after correction

    Example:
        >>> # After running analyze_gene_regions for multiple genes
        >>> stats_list = [results_MT_ND5['conservation_splitview_stats'],
        ...               results_MT_RNR2['conservation_splitview_stats']]
        >>> gene_names = ['MT-ND5', 'MT-RNR2']
        >>> fig, summary_df = plot_mdp_invariant_sites_summary(stats_list, gene_names)
    """
    # Collect data from all genes
    summary_data = []

    for gene_name, stats in zip(gene_names, stats_list):
        invariant_df = stats['invariant_percentages']
        variable_stats = stats.get('variable_stats', {})

        # Get non-overlapping data
        non_overlap_row = invariant_df[invariant_df['region_name'] == 'non-overlapping']
        if len(non_overlap_row) == 0:
            print(f"Warning: No non-overlapping region found for {gene_name}, skipping")
            continue

        non_overlap_pct = non_overlap_row['pct_invariant'].values[0]

        # Compare each overlapping region
        for _, row in invariant_df.iterrows():
            region_name = row['region_name']
            if region_name == 'non-overlapping':
                continue

            overlap_pct = row['pct_invariant']
            pct_difference = overlap_pct - non_overlap_pct

            # Get Mann-Whitney test results from variable_stats
            comparison_key = f"{region_name}_vs_non-overlapping"
            if comparison_key in variable_stats:
                test_results = variable_stats[comparison_key]
                mann_whitney_stat = test_results.get('statistic', np.nan)
                # Get the already-corrected p-value from variable_stats
                pval_raw = test_results.get('pvalue_raw', np.nan)
                # Note: variable_stats already has correction applied within each gene
                # We'll re-apply correction across all genes below
            else:
                mann_whitney_stat = np.nan
                pval_raw = np.nan

            summary_data.append({
                'gene': gene_name,
                'region': region_name,
                'non_overlap_pct': non_overlap_pct,
                'overlap_pct': overlap_pct,
                'pct_difference': pct_difference,
                'mann_whitney_statistic': mann_whitney_stat,
                'pvalue_raw': pval_raw
            })

    summary_df = pd.DataFrame(summary_data)

    if len(summary_df) == 0:
        print("No data to plot - all genes were skipped")
        return None, summary_df

    # Apply multiple testing correction across ALL comparisons (all genes and regions)
    if correction_method is not None and len(summary_df) > 0:
        valid_pvals = summary_df['pvalue_raw'].notna()
        if valid_pvals.sum() > 0:
            pvals = summary_df.loc[valid_pvals, 'pvalue_raw'].values
            reject, pvals_corrected, _, _ = multipletests(pvals, method=correction_method)

            summary_df['pvalue_corrected'] = np.nan
            summary_df['significant'] = False
            summary_df.loc[valid_pvals, 'pvalue_corrected'] = pvals_corrected
            summary_df.loc[valid_pvals, 'significant'] = reject
        else:
            summary_df['pvalue_corrected'] = np.nan
            summary_df['significant'] = False
    else:
        # No correction
        summary_df['pvalue_corrected'] = summary_df['pvalue_raw']
        summary_df['significant'] = summary_df['pvalue_raw'] < significance_threshold

    # Create plot with regions on X-axis and genes as hue
    fig, ax = plt.subplots(figsize=figsize)

    # Maintain order of regions as they appear in the data (grouped by gene)
    # Create a custom order: for each gene, get its regions in order
    region_order = []
    for gene in gene_names:
        gene_regions = summary_df[summary_df['gene'] == gene]['region'].unique()
        for region in gene_regions:
            if region not in region_order:
                region_order.append(region)

    # Use seaborn barplot with dodge parameter for better spacing control
    sns.barplot(data=summary_df, x='region', y='pct_difference', hue='gene',
                palette=palette, alpha=0.8, edgecolor='black', linewidth=0.5,
                dodge=False, ax=ax, hue_order=gene_names, order=region_order)

    # Add significance markers
    # With dodge=False, bars are positioned at integer indices (0, 1, 2, ...)
    n_regions = len(region_order)

    # Create a mapping of (gene, region) -> bar position
    bar_positions = {}
    for idx, region in enumerate(region_order):
        region_data = summary_df[summary_df['region'] == region]
        for _, row in region_data.iterrows():
            bar_positions[(row['gene'], region)] = idx

    # Iterate through data and add significance markers
    for _, row in summary_df.iterrows():
        if row['significant']:
            gene = row['gene']
            region = row['region']
            pval = row['pvalue_corrected']
            y_val = row['pct_difference']

            # Get bar position
            bar_x = bar_positions.get((gene, region), 0)
            y_pos = y_val + (1.5 if y_val > 0 else -1.5)

            # Determine number of asterisks based on p-value
            if pval < 0.001:
                marker = '***'
            elif pval < 0.01:
                marker = '**'
            else:
                marker = '*'

            ax.text(bar_x, y_pos, marker,
                   ha='center', va='bottom' if y_val > 0 else 'top',
                   fontsize=14, fontweight='bold', color='red')

    # Formatting
    ax.set_xlabel('Overlapping Region', fontsize=10)
    ax.set_ylabel(ylabel, fontsize=7)
    if show_title:
        ax.set_title(title, fontsize=14)
    ax.set_xticks(range(n_regions))
    ax.set_xticklabels(region_order, rotation=45, ha='center')
    ax.axhline(y=0, color='black', linestyle='-', linewidth=1)
    ax.legend(title='Encompassing Gene', bbox_to_anchor=(1.05, 1), loc='upper left', fontsize=9)
    ax.grid(True, alpha=0.3, axis='y')

    # Add note about significance levels
    if show_title:
        note_text = f'* p < 0.05, ** p < 0.01, *** p < 0.001\n(Mann-Whitney test, Correction: {correction_method if correction_method else "none"})'
        ax.text(0.02, 0.98, note_text, transform=ax.transAxes,
               fontsize=8, verticalalignment='top',
               bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.3))

    sns.despine()
    plt.tight_layout()

    if output_path:
        plt.savefig(output_path, dpi=300, bbox_inches='tight')

    return fig, summary_df
