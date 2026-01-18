"""
Mitochondrial Ribosome Profiling Initiation Analysis Extension

This module extends the base mito_riboseq_analysis module with additional
functions for detecting and analyzing translation initiation sites in
mitochondrial genes, particularly in overlapping reading frames and
alternative ORFs.

Key features:
- Peak detection for identifying initiation sites
- Analysis of bicistronic/overlapping regions (e.g., MT-ATP8/MT-ATP6)
- Start codon enrichment testing
- Condition comparison for differential initiation
- Metagene analysis around start/stop codons

Author: Extended for mitochondrial microprotein initiation research
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
from scipy import signal, stats
from typing import Dict, List, Tuple, Optional, Union
import warnings
warnings.filterwarnings('ignore')

# Import base functions
from mito_riboseq_analysis import load_and_prepare_data


def detect_initiation_peaks(df: pd.DataFrame,
                            gene_id: str,
                            prominence: float = 2.0,
                            min_height: float = 1.5,
                            min_distance: int = 50,
                            window_size: int = 100) -> pd.DataFrame:
    """
    Detect potential initiation sites based on ribosome occupancy peaks.

    Logic:
    1. Filter data for the specified gene
    2. Extract RPM values (ribosome density per codon)
    3. Normalize by mean to get relative occupancy
    4. Use scipy's find_peaks to detect significant peaks
    5. Return peak positions and properties

    Parameters:
    -----------
    df : pd.DataFrame
        Data from load_and_prepare_data()
    gene_id : str
        Gene name to analyze (e.g., 'MT-ATP8', 'gene_MOTSc')
    prominence : float, default=2.0
        Minimum peak prominence for detection (relative to baseline)
    min_height : float, default=1.5
        Minimum peak height (normalized occupancy)
    min_distance : int, default=50
        Minimum distance between peaks (codons)
    window_size : int, default=100
        Number of codons from start to analyze

    Returns:
    --------
    pd.DataFrame
        Detected peaks with columns:
        - codon_position: Codon index of peak
        - codon_seq: Codon sequence at peak
        - rpm: RPM value at peak
        - normalized_rpm: Normalized RPM at peak
        - prominence: Peak prominence
        - peak_height: Peak height above baseline

    Notes:
    ------
    - Uses RPM (ribosome density) not RPM_cumsum_normgene
    - Averages across samples for each codon position
    - Focuses on initiation window (first ~100 codons)
    """
    # Filter for the specific gene
    gene_data = df[df['gene_id'] == gene_id].copy()

    if len(gene_data) == 0:
        print(f"No data found for {gene_id}")
        return pd.DataFrame()

    # Limit to initiation window
    gene_data = gene_data[gene_data['codon_index'] <= window_size]

    # Average RPM across samples for each codon position
    # Group by codon_index and aggregate
    occupancy_data = gene_data.groupby('codon_index').agg({
        'RPM_normgene': 'mean',
        'codon_seq': 'first'
    }).reset_index()

    occupancy_data = occupancy_data.sort_values('codon_index')

    occupancy = occupancy_data['RPM_normgene'].values

    # Normalize by mean
    if occupancy.mean() > 0:
        norm_occupancy = occupancy / occupancy.mean()
    else:
        norm_occupancy = occupancy

    # Detect peaks using scipy
    peaks, properties = signal.find_peaks(
        norm_occupancy,
        prominence=prominence,
        height=min_height,
        distance=min_distance
    )

    # Create results dataframe
    results = pd.DataFrame({
        'codon_position': occupancy_data.iloc[peaks]['codon_index'].values,
        'codon_seq': occupancy_data.iloc[peaks]['codon_seq'].values,
        'rpm': occupancy[peaks],
        'normalized_rpm': norm_occupancy[peaks],
        'prominence': properties['prominences'],
        'peak_height': properties['peak_heights']
    })

    return results


def test_start_codon_enrichment(df: pd.DataFrame,
                                gene_id: str,
                                start_positions: List[int],
                                window: int = 10,
                                background_window: int = 50) -> pd.DataFrame:
    """
    Test for significant enrichment at specific start codon positions.

    Logic:
    1. For each candidate start position:
       - Calculate ribosome density in peak window
       - Calculate background density from flanking regions
       - Compute enrichment ratio and statistical significance
    2. Use Mann-Whitney U test to test for significant enrichment
    3. Return enrichment statistics for each position

    Parameters:
    -----------
    df : pd.DataFrame
        Data from load_and_prepare_data()
    gene_id : str
        Gene name to analyze
    start_positions : list of int
        Codon positions to test (0-indexed)
    window : int, default=10
        Window around start codon for calculating peak
    background_window : int, default=50
        Window for background calculation (excludes peak window)

    Returns:
    --------
    pd.DataFrame
        Enrichment statistics for each position:
        - position: Codon index tested
        - codon: Codon sequence at position
        - rpm: RPM value at position
        - peak_mean: Mean RPM in peak window
        - background_mean: Mean RPM in background
        - enrichment: Fold enrichment (peak/background)
        - z_score: Z-score of enrichment
        - pvalue: Mann-Whitney U test p-value
        - significant: Whether p < 0.05

    Notes:
    ------
    - Useful for validating predicted internal ORF starts
    - Background excludes the peak region to avoid bias
    - Statistical test accounts for variability in ribosome density
    """
    # Filter for the specific gene
    gene_data = df[df['gene_id'] == gene_id].copy()

    if len(gene_data) == 0:
        print(f"No data found for {gene_id}")
        return pd.DataFrame()

    # Average RPM across samples for each codon position
    occupancy_data = gene_data.groupby('codon_index').agg({
        'RPM_normgene': 'mean',
        'codon_seq': 'first'
    }).reset_index()

    occupancy_data = occupancy_data.sort_values('codon_index')
    occupancy = occupancy_data['RPM_normgene'].values

    results = []
    for pos in start_positions:
        if pos < background_window or pos >= len(occupancy) - background_window:
            continue

        # Get peak value and local mean
        peak_start = max(0, pos - 2)
        peak_end = min(len(occupancy), pos + window)
        peak_region = occupancy[peak_start:peak_end]
        peak_value = occupancy[pos]
        peak_mean = peak_region.mean()

        # Get background (excluding the peak region)
        bg_start1 = max(0, pos - background_window)
        bg_end1 = max(0, pos - window)
        bg_start2 = min(len(occupancy), pos + window)
        bg_end2 = min(len(occupancy), pos + background_window)

        background = np.concatenate([
            occupancy[bg_start1:bg_end1],
            occupancy[bg_start2:bg_end2]
        ])

        if len(background) > 0:
            bg_mean = background.mean()
            bg_std = background.std()

            # Calculate enrichment
            enrichment = peak_value / (bg_mean + 1e-6)
            z_score = (peak_value - bg_mean) / (bg_std + 1e-6) if bg_std > 0 else 0

            # Statistical test
            try:
                _, pval = stats.mannwhitneyu(
                    peak_region,
                    background,
                    alternative='greater'
                )
            except:
                pval = 1.0
        else:
            enrichment = peak_value
            z_score = 0
            pval = 1.0
            bg_mean = 0
            peak_mean = peak_value

        # Get codon at this position
        codon_row = occupancy_data[occupancy_data['codon_index'] == pos]
        codon = codon_row['codon_seq'].values[0] if len(codon_row) > 0 else 'NNN'

        results.append({
            'position': pos,
            'codon': codon,
            'rpm': peak_value,
            'peak_mean': peak_mean,
            'background_mean': bg_mean,
            'enrichment': enrichment,
            'z_score': z_score,
            'pvalue': pval,
            'significant': pval < 0.05
        })

    return pd.DataFrame(results)


def analyze_bicistronic_region(df: pd.DataFrame,
                               gene_id: str,
                               upstream_start: int,
                               downstream_start: int,
                               window: int = 5) -> Dict:
    """
    Analyze initiation pattern in bicistronic/overlapping region.

    Logic:
    1. Extract ribosome occupancy for the gene
    2. Analyze occupancy at both upstream and downstream start positions
    3. Calculate initiation ratio (downstream/upstream)
    4. Test for enrichment at both positions
    5. Return comprehensive analysis results

    Following the approach in Wakigawa et al. Figure 4 for analyzing
    overlapping ORFs like MT-ATP8/MT-ATP6.

    Parameters:
    -----------
    df : pd.DataFrame
        Data from load_and_prepare_data()
    gene_id : str
        Gene name (e.g., 'MT-ATP8' for ATP8/ATP6 region)
    upstream_start : int
        Start codon position of upstream ORF (codon index)
    downstream_start : int
        Start codon position of downstream/internal ORF (codon index)
    window : int, default=5
        Window size around start codons for analysis

    Returns:
    --------
    dict
        Analysis results containing:
        - gene_id: Gene analyzed
        - upstream_start, downstream_start: Start positions
        - upstream_initiation, downstream_initiation: RPM at start codons
        - upstream_mean, downstream_mean: Mean RPM in windows
        - upstream_codon, downstream_codon: Codon sequences
        - initiation_ratio: Downstream/upstream initiation strength
        - occupancy_profile: Full normalized occupancy profile
        - enrichment_stats: Statistical enrichment test results

    Notes:
    ------
    - Identifies differential initiation between overlapping ORFs
    - Initiation ratio > 0.5 suggests significant internal initiation
    - Used to validate alternative reading frame translation
    """
    # Filter for the specific gene
    gene_data = df[df['gene_id'] == gene_id].copy()

    if len(gene_data) == 0:
        print(f"No data found for {gene_id}")
        return {}

    # Average RPM across samples for each codon position
    occupancy_data = gene_data.groupby('codon_index').agg({
        'RPM_normgene': 'mean',
        'codon_seq': 'first'
    }).reset_index()

    occupancy_data = occupancy_data.sort_values('codon_index')
    occupancy = occupancy_data['RPM_normgene'].values

    # Normalize
    if occupancy.mean() > 0:
        norm_occupancy = occupancy / occupancy.mean()
    else:
        norm_occupancy = occupancy

    results = {
        'gene_id': gene_id,
        'upstream_start': upstream_start,
        'downstream_start': downstream_start,
        'total_codons': len(occupancy)
    }

    # Upstream ORF initiation
    if upstream_start < len(occupancy):
        upstream_region_start = max(0, upstream_start)
        upstream_region_end = min(len(occupancy), upstream_start + window)
        upstream_region = occupancy[upstream_region_start:upstream_region_end]

        results['upstream_initiation'] = occupancy[upstream_start]
        results['upstream_mean'] = upstream_region.mean()
        results['upstream_max'] = upstream_region.max()

        # Get codon
        codon_row = occupancy_data[occupancy_data['codon_index'] == upstream_start]
        if len(codon_row) > 0:
            results['upstream_codon'] = codon_row['codon_seq'].values[0]

    # Downstream/internal ORF initiation
    if downstream_start < len(occupancy):
        downstream_region_start = max(0, downstream_start)
        downstream_region_end = min(len(occupancy), downstream_start + window)
        downstream_region = occupancy[downstream_region_start:downstream_region_end]

        results['downstream_initiation'] = occupancy[downstream_start]
        results['downstream_mean'] = downstream_region.mean()
        results['downstream_max'] = downstream_region.max()

        # Get codon
        codon_row = occupancy_data[occupancy_data['codon_index'] == downstream_start]
        if len(codon_row) > 0:
            results['downstream_codon'] = codon_row['codon_seq'].values[0]

    # Calculate initiation ratio
    if 'upstream_initiation' in results and results['upstream_initiation'] > 0:
        results['initiation_ratio'] = (
            results.get('downstream_initiation', 0) / results['upstream_initiation']
        )
    else:
        results['initiation_ratio'] = 0

    # Store full profile
    results['occupancy_profile'] = norm_occupancy
    results['raw_counts'] = occupancy

    # Test for enrichment at both positions
    enrichment_test = test_start_codon_enrichment(
        df,
        gene_id,
        [upstream_start, downstream_start]
    )
    results['enrichment_stats'] = enrichment_test

    return results


def compare_initiation_conditions(df1: pd.DataFrame,
                                  df2: pd.DataFrame,
                                  gene_id: str,
                                  start_positions: List[int],
                                  label1: str = "Condition 1",
                                  label2: str = "Condition 2") -> pd.DataFrame:
    """
    Compare initiation between two conditions.

    Logic:
    1. Extract occupancy data for the same gene from both conditions
    2. Normalize both profiles
    3. Calculate fold change and statistical differences
    4. Return comparison metrics for each start position

    Useful for detecting condition-dependent initiation (e.g., control vs
    retapamulin treatment, or mtIF3 dependency).

    Parameters:
    -----------
    df1 : pd.DataFrame
        Data from first condition (e.g., treatment)
    df2 : pd.DataFrame
        Data from second condition (e.g., control)
    gene_id : str
        Gene or region name to compare
    start_positions : list of int
        Start codon positions to compare
    label1 : str, default="Condition 1"
        Label for first condition
    label2 : str, default="Condition 2"
        Label for second condition

    Returns:
    --------
    pd.DataFrame
        Comparison statistics:
        - position: Codon index
        - codon: Codon sequence
        - condition1, condition2: Normalized RPM values
        - fold_change: condition1 / condition2
        - log2_fold_change: Log2 fold change
        - difference: Absolute difference

    Notes:
    ------
    - Both conditions are normalized to mean=1 for fair comparison
    - Log2 fold change is calculated with pseudocount (+1) to handle zeros
    - Useful for identifying condition-specific initiation events
    """
    # Get data from both conditions
    gene1 = df1[df1['gene_id'] == gene_id].copy()
    gene2 = df2[df2['gene_id'] == gene_id].copy()

    if len(gene1) == 0 or len(gene2) == 0:
        print(f"Data not available in one or both conditions for {gene_id}")
        return pd.DataFrame()

    # Average across samples for each condition
    occ1_data = gene1.groupby('codon_index').agg({
        'RPM_normgene': 'mean',
        'codon_seq': 'first'
    }).reset_index()

    occ2_data = gene2.groupby('codon_index').agg({
        'RPM_normgene': 'mean',
        'codon_seq': 'first'
    }).reset_index()

    occ1_data = occ1_data.sort_values('codon_index')
    occ2_data = occ2_data.sort_values('codon_index')

    occ1 = occ1_data['RPM_normgene'].values
    occ2 = occ2_data['RPM_normgene'].values

    # Normalize both
    if occ1.mean() > 0:
        occ1 = occ1 / occ1.mean()
    if occ2.mean() > 0:
        occ2 = occ2 / occ2.mean()

    results = []
    for pos in start_positions:
        if pos >= min(len(occ1), len(occ2)):
            continue

        # Get codon
        codon_row = occ1_data[occ1_data['codon_index'] == pos]
        codon = codon_row['codon_seq'].values[0] if len(codon_row) > 0 else 'NNN'

        result = {
            'position': pos,
            'codon': codon,
            label1: occ1[pos],
            label2: occ2[pos],
            'fold_change': occ1[pos] / (occ2[pos] + 1e-6),
            'log2_fold_change': np.log2(occ1[pos] + 1) - np.log2(occ2[pos] + 1),
            'difference': occ1[pos] - occ2[pos]
        }
        results.append(result)

    return pd.DataFrame(results)


def plot_initiation_profile(df: pd.DataFrame,
                            gene_id: str,
                            peaks: Optional[pd.DataFrame] = None,
                            start_positions: Optional[List[int]] = None,
                            window_size: int = 100,
                            figsize: Tuple[int, int] = (14, 6),
                            save_path: Optional[Path] = None,
                            title: Optional[bool] = False) -> plt.Figure:
    """
    Plot ribosome occupancy profile with detected initiation sites.

    Logic:
    1. Extract ribosome occupancy for the gene
    2. Plot normalized RPM profile
    3. Mark detected peaks and known start positions
    4. Add annotations for key positions

    Parameters:
    -----------
    df : pd.DataFrame
        Data from load_and_prepare_data()
    gene_id : str
        Gene name to plot
    peaks : pd.DataFrame, optional
        Detected peaks from detect_initiation_peaks()
    start_positions : list of int, optional
        Known/predicted start positions to mark
    window_size : int, default=100
        Number of codons from start to show
    figsize : tuple, default=(14, 6)
        Figure size (width, height)
    save_path : Path, optional
        Path to save figure

    Returns:
    --------
    matplotlib.figure.Figure
        The generated figure

    Notes:
    ------
    - Y-axis shows normalized RPM (mean = 1.0)
    - Red dashed lines mark known start positions
    - Orange triangles mark automatically detected peaks
    - Useful for visual validation of initiation sites
    """
    gene_data = df[df['gene_id'] == gene_id].copy()

    if len(gene_data) == 0:
        print(f"No data for {gene_id}")
        return None

    # Limit to window
    gene_data = gene_data[gene_data['codon_index'] <= window_size]

    # Average across samples
    occupancy_data = gene_data.groupby('codon_index').agg({
        'RPM_normgene': 'mean',
        'codon_seq': 'first'
    }).reset_index()

    occupancy_data = occupancy_data.sort_values('codon_index')

    # Convert codon_index to int
    occupancy_data['codon_index'] = occupancy_data['codon_index'].astype(int)

    occupancy = occupancy_data['RPM_normgene'].values
    if occupancy.mean() > 0:
        norm_occ = occupancy / occupancy.mean()
    else:
        norm_occ = occupancy

    codons = occupancy_data['codon_index'].values

    # Create figure
    fig, ax = plt.subplots(figsize=figsize)

    # Plot main occupancy line
    ax.plot(codons, norm_occ, linewidth=2, alpha=0.8,
           color='#2E86AB', label='Ribosome occupancy')
    ax.fill_between(codons, norm_occ, alpha=0.2, color='#2E86AB')

    # Mark known start positions
    if start_positions:
        for i, pos in enumerate(start_positions):
            if pos <= window_size:
                label = 'Known start' if i == 0 else None
                #ax.axvline(pos, color='red', linestyle='--',
                #         linewidth=2.5, alpha=0.7, label=label, zorder=0)

                # Add position label
                pos_idx = occupancy_data[occupancy_data['codon_index'] == pos].index
                if len(pos_idx) > 0:
                    codon = occupancy_data.loc[pos_idx[0], 'codon_seq']
                    #ax.text(pos, ax.get_ylim()[1] * 0.95,
                    #      f'{pos}\n{codon}',
                    #       ha='center', va='top', fontsize=9,
                    #       bbox=dict(boxstyle='round', facecolor='white',
                    #               alpha=0.8, edgecolor='red'))

    # Mark detected peaks
    if peaks is not None and len(peaks) > 0:
        peak_pos = peaks['codon_position'].values
        peak_heights = []
        for pos in peak_pos:
            pos_idx = occupancy_data[occupancy_data['codon_index'] == pos].index
            if len(pos_idx) > 0:
                idx = pos_idx[0]
                peak_heights.append(norm_occ[idx])

        if peak_heights:
            ax.scatter(peak_pos[:len(peak_heights)], peak_heights,
                      s=150, marker='^', color='orange',
                      edgecolors='black', linewidths=2,
                      zorder=5, label='Detected initiation')

    ax.set_xlabel('Codon position', fontsize=12, fontweight='bold')
    ax.set_ylabel('Ribosome occupancy\n(normalized)', fontsize=12, fontweight='bold')
    if title:
        ax.set_title(f'{gene_id.replace("gene_", "")} - Initiation Site Analysis',
                fontsize=14, fontweight='bold', pad=15)
    ax.legend(loc='best', framealpha=0.95, fontsize=10)
    ax.grid(alpha=0.3, linestyle='--')
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    plt.tight_layout()

    if save_path:
        fig.savefig(save_path, dpi=300, bbox_inches='tight')
        print(f"Saved to {save_path}")

    return fig


def plot_bicistronic_initiation(df: pd.DataFrame,
                                gene_id: str,
                                upstream_start: int,
                                downstream_start: int,
                                df_control: Optional[pd.DataFrame] = None,
                                label_main: str = "Sample",
                                label_control: str = "Control",
                                figsize: Tuple[int, int] = (14, 7),
                                save_path: Optional[Path] = None) -> plt.Figure:
    """
    Plot initiation pattern in bicistronic region.

    Logic:
    1. Analyze the bicistronic region
    2. Plot normalized occupancy profile
    3. Mark upstream and downstream start codons
    4. Optionally overlay control condition for comparison
    5. Annotate initiation strengths

    Following Wakigawa et al. Figure 4 style for overlapping ORF analysis.

    Parameters:
    -----------
    df : pd.DataFrame
        Data from load_and_prepare_data() for main condition
    gene_id : str
        Gene name (e.g., 'MT-ATP8' for ATP8/ATP6 region)
    upstream_start : int
        Upstream ORF start position (codons)
    downstream_start : int
        Downstream ORF start position (codons)
    df_control : pd.DataFrame, optional
        Data for control condition
    label_main : str, default="Sample"
        Label for main condition
    label_control : str, default="Control"
        Label for control condition
    figsize : tuple, default=(14, 7)
        Figure size
    save_path : Path, optional
        Path to save figure

    Returns:
    --------
    matplotlib.figure.Figure
        The generated figure

    Notes:
    ------
    - Gray shaded region shows overlap between ORFs
    - Red dashed line: upstream start codon
    - Orange dashed line: downstream/internal start codon
    - Annotations show initiation strength at each position
    - Initiation ratio displayed in title
    """
    # Analyze this condition
    results = analyze_bicistronic_region(
        df, gene_id, upstream_start, downstream_start
    )

    if not results:
        return None

    fig, ax = plt.subplots(figsize=figsize)

    occupancy = results['occupancy_profile']
    codons = np.arange(len(occupancy))

    # Plot this condition
    ax.plot(codons, occupancy, linewidth=2.5,
           label=label_main,
           color='#2E86AB', alpha=0.9)

    # Plot control if provided
    if df_control is not None:
        control_results = analyze_bicistronic_region(
            df_control, gene_id, upstream_start, downstream_start
        )
        if control_results and 'occupancy_profile' in control_results:
            control_occ = control_results['occupancy_profile']
            ax.plot(codons[:len(control_occ)], control_occ,
                   linewidth=2.5,
                   label=label_control,
                   color='#A23B72', alpha=0.9)

    # Mark start codons
    upstream_codon = results.get('upstream_codon', 'AUG')
    downstream_codon = results.get('downstream_codon', 'AUG')

    ax.axvline(upstream_start, color='#D62828', linestyle='--',
              linewidth=3, alpha=0.8,
              label=f'Upstream start ({upstream_codon})')
    ax.axvline(downstream_start, color='#F77F00', linestyle='--',
              linewidth=3, alpha=0.8,
              label=f'Internal start ({downstream_codon})')

    # Highlight overlap region
    if downstream_start > upstream_start:
        ax.axvspan(upstream_start, downstream_start,
                  alpha=0.15, color='gray', label='Overlap region', zorder=0)

    # Add annotations for initiation strength
    y_max = ax.get_ylim()[1]
    if 'upstream_initiation' in results and upstream_start < len(occupancy):
        ax.annotate(f'Init: {results["upstream_initiation"]:.2f}',
                   xy=(upstream_start, occupancy[upstream_start]),
                   xytext=(upstream_start - 20, y_max * 0.85),
                   fontsize=9, ha='center',
                   bbox=dict(boxstyle='round', facecolor='white', alpha=0.8),
                   arrowprops=dict(arrowstyle='->', color='red', lw=1.5))

    if 'downstream_initiation' in results and downstream_start < len(occupancy):
        ax.annotate(f'Init: {results["downstream_initiation"]:.2f}',
                   xy=(downstream_start, occupancy[downstream_start]),
                   xytext=(downstream_start + 20, y_max * 0.85),
                   fontsize=9, ha='center',
                   bbox=dict(boxstyle='round', facecolor='white', alpha=0.8),
                   arrowprops=dict(arrowstyle='->', color='orange', lw=1.5))

    ax.set_xlabel('Codon position', fontsize=12, fontweight='bold')
    ax.set_ylabel('Mitoribosome occupancy', fontsize=12, fontweight='bold')

    title = f'{gene_id} - Bicistronic Initiation Pattern'
    if 'initiation_ratio' in results:
        title += f'\nInternal/Upstream ratio: {results["initiation_ratio"]:.2f}'
    ax.set_title(title, fontsize=14, fontweight='bold', pad=15)

    ax.legend(loc='best', framealpha=0.95, fontsize=10)
    ax.grid(alpha=0.3, linestyle='--')
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    plt.tight_layout()

    if save_path:
        fig.savefig(save_path, dpi=300, bbox_inches='tight')
        print(f"Saved to {save_path}")

    return fig


def plot_metagene_initiation(df: pd.DataFrame,
                             gene_ids: List[str],
                             window: int = 50,
                             align_to: str = 'start',
                             figsize: Tuple[int, int] = (12, 6),
                             save_path: Optional[Path] = None) -> plt.Figure:
    """
    Create metagene plot around start or stop codons.

    Logic:
    1. For each gene, extract occupancy around alignment position
       - For 'start': extracts from start codon (position 0) downstream
       - For 'stop': extracts window centered on stop codon
    2. Normalize each gene's profile
    3. Calculate mean and SEM across all genes
    4. Plot composite initiation/termination signature

    Similar to Wakigawa et al. Figure 2F-G metagene analysis.

    Parameters:
    -----------
    df : pd.DataFrame
        Data from load_and_prepare_data()
    gene_ids : list of str
        Genes to include in metagene (e.g., all canonical genes)
    window : int, default=50
        Window size around alignment position (codons)
    align_to : str, default='start'
        'start' or 'stop' - codon to align to
    figsize : tuple, default=(12, 6)
        Figure size
    save_path : Path, optional
        Path to save figure

    Returns:
    --------
    matplotlib.figure.Figure
        The generated metagene figure

    Notes:
    ------
    - Position 0 marks the alignment codon (start or stop)
    - For 'start' mode: shows downstream region (0 to +window)
    - For 'stop' mode: shows region centered on stop (-window to +window)
    - Shaded region shows SEM across genes
    - Reveals consensus initiation/termination signature
    - Differences between canonical and overlapping genes can be detected
    """
    fig, ax = plt.subplots(figsize=figsize)

    all_profiles = []

    for gene_id in gene_ids:
        gene_data = df[df['gene_id'] == gene_id].copy()

        if len(gene_data) == 0:
            continue

        # Average across samples
        occupancy_data = gene_data.groupby('codon_index').agg({
            'RPM_normgene': 'mean'
        }).reset_index()

        occupancy_data = occupancy_data.sort_values('codon_index')
        occupancy = occupancy_data['RPM_normgene'].values

        # Normalize
        if occupancy.mean() > 0:
            norm_occ = occupancy / occupancy.mean()
        else:
            norm_occ = occupancy

        # Extract window around start (position 0) or stop
        if align_to == 'start':
            align_pos = 0
            # For start codon, extract from position 0 onwards (downstream)
            start_idx = 0
            end_idx = min(len(norm_occ), window)
        else:  # stop
            # For stop codon, extract window centered on the stop
            align_pos = len(norm_occ) - 1
            start_idx = max(0, align_pos - window)
            end_idx = min(len(norm_occ), align_pos + window)

        if end_idx - start_idx > 0:
            profile = norm_occ[start_idx:end_idx]
            all_profiles.append(profile)

    if not all_profiles:
        print("No data to plot")
        return None

    # Calculate mean and SEM
    max_len = max(len(p) for p in all_profiles)
    padded = []
    for p in all_profiles:
        if len(p) < max_len:
            padded.append(np.pad(p, (0, max_len - len(p)), constant_values=np.nan))
        else:
            padded.append(p)

    profiles_array = np.array(padded)
    mean_profile = np.nanmean(profiles_array, axis=0)
    sem_profile = np.nanstd(profiles_array, axis=0) / np.sqrt(np.sum(~np.isnan(profiles_array), axis=0))

    # Set positions based on alignment mode
    if align_to == 'start':
        # For start codon, positions are 0 onwards (downstream region)
        positions = np.arange(0, len(mean_profile))
    else:
        # For stop codon, center on position 0
        positions = np.arange(-window, -window + len(mean_profile))

    # Plot mean with error band
    ax.plot(positions, mean_profile, linewidth=3,
           color='#2E86AB', label=f'Mean (n={len(all_profiles)})')
    ax.fill_between(positions,
                   mean_profile - sem_profile,
                   mean_profile + sem_profile,
                   alpha=0.3, color='#2E86AB')

    # Mark alignment position
    ax.axvline(0, color='red', linestyle='--', linewidth=3,
              alpha=0.7, label=f'{align_to.capitalize()} codon', zorder=0)

    # Highlight region around alignment position
    if align_to == 'start':
        # For start codon, highlight the initiation region (0 to +5)
        ax.axvspan(0, 5, alpha=0.1, color='red', zorder=0)
    else:
        # For stop codon, highlight region around it (-5 to +5)
        ax.axvspan(-5, 5, alpha=0.1, color='red', zorder=0)

    ax.set_xlabel(f'Distance from {align_to} codon (codons)',
                 fontsize=12, fontweight='bold')
    ax.set_ylabel('Mean ribosome occupancy', fontsize=12, fontweight='bold')
    ax.set_title(f'Metagene Profile - {align_to.capitalize()} Codon',
                fontsize=14, fontweight='bold', pad=15)
    ax.legend(loc='best', framealpha=0.95, fontsize=11)
    ax.grid(alpha=0.3, linestyle='--')
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    plt.tight_layout()

    if save_path:
        fig.savefig(save_path, dpi=300, bbox_inches='tight')
        print(f"Saved to {save_path}")

    return fig


# Convenience functions for quick analysis

def quick_initiation_analysis(table_path: Union[str, Path],
                              gene_id: str,
                              known_starts: Optional[List[int]] = None,
                              output_dir: Optional[Path] = None,
                              sep: str = ',',
                              title: Optional[bool] = False,
                              figsize: Optional[Tuple[int, int]] = (10, 6)) -> Tuple[pd.DataFrame, plt.Figure]:
    """
    Quick initiation analysis for a single gene.

    Logic:
    1. Load the data
    2. Detect peaks automatically
    3. Test known start positions if provided
    4. Generate visualization

    Parameters:
    -----------
    table_path : str or Path
        Path to cumsum table CSV file
    gene_id : str
        Gene to analyze
    known_starts : list of int, optional
        Known start positions to test
    output_dir : Path, optional
        Directory to save results
    sep : str, default=','
        CSV separator
    title : bool, default=False
        Whether to add title to plot
    figsize : tuple, default=(10, 6)
        Figure size (width, height) in inches

    Returns:
    --------
    tuple : (peaks_df, figure)
        Detected peaks and visualization

    Notes:
    ------
    - One-line function for quick analysis
    - Automatically detects peaks and generates plot
    - Optionally validates known start positions
    """
    # Load data
    df = load_and_prepare_data(table_path, sep=sep)

    # Detect peaks
    peaks = detect_initiation_peaks(df, gene_id)

    print(f"\nDetected {len(peaks)} peak(s) in {gene_id}:")
    print(peaks)

    # Test known starts if provided
    if not known_starts:
        print(f'No known start positions provided for {gene_id}. Using detected peaks only.')
        known_starts = peaks['codon_position'].tolist()

    enrichment = test_start_codon_enrichment(df, gene_id, known_starts)
    print("\nEnrichment at known starts:")
    print(enrichment)


    # Plot
    fig = plot_initiation_profile(df, gene_id, peaks=peaks, start_positions=known_starts, title=title, figsize = figsize)

    if output_dir:
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)
        fig.savefig(output_dir / f"{gene_id}_initiation.pdf",
                   dpi=300, bbox_inches='tight')

    return peaks, fig


def batch_initiation_analysis(df: pd.DataFrame,
                               gene_ids: List[str],
                               known_starts: Optional[Dict[str, List[int]]] = None,
                               window_size: int = 100,
                               output_dir: Optional[Path] = None,
                               ncols: int = 3,
                               figsize_per_plot: Tuple[float, float] = (5, 4)) -> Tuple[Dict[str, pd.DataFrame], plt.Figure]:
    """
    Perform initiation analysis for multiple genes and create a facet grid.

    Logic:
    1. For each gene, detect peaks and calculate occupancy
    2. Create a multi-panel figure with one subplot per gene
    3. Plot initiation profile for each gene in its own panel
    4. Return all peaks and the combined figure

    Parameters:
    -----------
    df : pd.DataFrame
        Data from load_and_prepare_data()
    gene_ids : list of str
        List of genes to analyze
    known_starts : dict, optional
        Dictionary mapping gene_id to list of known start positions
        Example: {'MT-ATP8': [0, 43], 'gene_MOTSc': [0]}
    window_size : int, default=100
        Number of codons from start to show
    output_dir : Path, optional
        Directory to save results
    ncols : int, default=3
        Number of columns in the facet grid
    figsize_per_plot : tuple, default=(5, 4)
        Size of each individual subplot (width, height)

    Returns:
    --------
    tuple : (peaks_dict, figure)
        - peaks_dict: Dictionary mapping gene_id to detected peaks DataFrame
        - figure: Combined facet grid figure

    Notes:
    ------
    - Automatically arranges plots in a grid
    - Each panel shows one gene's initiation profile
    - Useful for comparing initiation patterns across multiple genes
    - Color-codes overlapping genes differently from canonical genes
    """
    # If gene_id is empty, use all
    if gene_ids == []:
        gene_ids = df['gene_id'].unique().tolist()
    # Calculate grid dimensions
    nrows = int(np.ceil(len(gene_ids) / ncols))
    fig_width = ncols * figsize_per_plot[0]
    fig_height = nrows * figsize_per_plot[1]

    # Create figure with subplots
    fig, axes = plt.subplots(nrows, ncols, figsize=(fig_width, fig_height))

    # Flatten axes array for easier iteration
    if nrows == 1 and ncols == 1:
        axes = np.array([axes])
    axes = axes.flatten() if isinstance(axes, np.ndarray) else [axes]

    all_peaks = {}

    for idx, gene_id in enumerate(gene_ids):
        ax = axes[idx]

        # Filter for the specific gene
        gene_data = df[df['gene_id'] == gene_id].copy()

        if len(gene_data) == 0:
            ax.text(0.5, 0.5, f'No data for\n{gene_id}',
                   ha='center', va='center', transform=ax.transAxes,
                   fontsize=10, color='red')
            ax.set_xticks([])
            ax.set_yticks([])
            continue

        # Limit to window
        gene_data = gene_data[gene_data['codon_index'] <= window_size]

        # Average across samples
        occupancy_data = gene_data.groupby('codon_index').agg({
            'RPM_normgene': 'mean',
            'codon_seq': 'first',
            'transcript_class': 'first'
        }).reset_index()

        occupancy_data = occupancy_data.sort_values('codon_index')
        occupancy_data['codon_index'] = occupancy_data['codon_index'].astype(int)

        occupancy = occupancy_data['RPM_normgene'].values
        if occupancy.mean() > 0:
            norm_occ = occupancy / occupancy.mean()
        else:
            norm_occ = occupancy

        codons = occupancy_data['codon_index'].values

        # Determine color based on transcript class
        trans_class = occupancy_data['transcript_class'].iloc[0]
        color = '#E74C3C' if trans_class == 'Overlapping' else '#3498DB'

        # Plot occupancy profile
        ax.plot(codons, norm_occ, linewidth=2, alpha=0.8,
               color=color, label='Occupancy')
        ax.fill_between(codons, norm_occ, alpha=0.2, color=color)

        # Detect peaks
        peaks = detect_initiation_peaks(df, gene_id, window_size=window_size)
        all_peaks[gene_id] = peaks

        # Mark detected peaks
        if peaks is not None and len(peaks) > 0:
            peak_pos = peaks['codon_position'].values
            peak_heights = []
            for pos in peak_pos:
                pos_idx = occupancy_data[occupancy_data['codon_index'] == pos].index
                if len(pos_idx) > 0:
                    idx = pos_idx[0]
                    peak_heights.append(norm_occ[idx])

            if peak_heights:
                ax.scatter(peak_pos[:len(peak_heights)], peak_heights,
                          s=100, marker='^', color='orange',
                          edgecolors='black', linewidths=1.5,
                          zorder=5)

        # Mark known starts if provided
        if known_starts and gene_id in known_starts:
            for pos in known_starts[gene_id]:
                if pos <= window_size:
                    ax.axvline(pos, color='red', linestyle='--',
                             linewidth=1.5, alpha=0.5, zorder=0)

        # Format subplot
        gene_display = gene_id.replace('gene_', '')
        ax.set_title(gene_display, fontsize=11, fontweight='bold')
        ax.set_xlabel('Codon position', fontsize=9)
        ax.set_ylabel('Norm. occupancy', fontsize=9)
        ax.tick_params(labelsize=8)
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.grid(alpha=0.2, linestyle='--')

    # Hide unused subplots
    for idx in range(len(gene_ids), len(axes)):
        axes[idx].axis('off')

    plt.tight_layout()

    if output_dir:
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)
        fig.savefig(output_dir / "batch_initiation_analysis.jpg",
                   dpi=300, bbox_inches='tight')

        # Also save peaks summary
        peaks_summary = []
        for gene_id, peaks_df in all_peaks.items():
            for _, row in peaks_df.iterrows():
                peaks_summary.append({
                    'gene_id': gene_id,
                    **row.to_dict()
                })
        if peaks_summary:
            pd.DataFrame(peaks_summary).to_csv(
                output_dir / "batch_peaks_summary.csv",
                index=False
            )

    return all_peaks, fig


def compare_to_canonical_metagene(df: pd.DataFrame,
                                  alternative_gene_id: str,
                                  canonical_gene_ids: Optional[List[str]] = None,
                                  window: int = 50,
                                  method: str = 'correlation',
                                  plot: bool = True,
                                  plot_metrics: bool = True,
                                  figsize: Tuple[int, int] = (12, 5),
                                  save_path: Optional[Path] = None) -> Dict:
    """
    Compare an alternative gene's initiation profile to the canonical metagene.

    Logic:
    1. Calculate the metagene profile for canonical genes
    2. Extract the alternative gene's initiation profile
    3. Align both profiles (normalize and truncate to same length)
    4. Calculate similarity metrics:
       - Pearson correlation
       - Spearman correlation
       - Peak position similarity
       - Shape similarity (cosine similarity)
    5. Optionally plot comparison
    6. Return comprehensive comparison metrics

    This helps determine if an alternative reading frame shows a credible
    translation initiation signature similar to known canonical genes.

    Parameters:
    -----------
    df : pd.DataFrame
        Data from load_and_prepare_data()
    alternative_gene_id : str
        Alternative gene to test (e.g., 'gene_MOTSc', 'gene_SHLP2')
    canonical_gene_ids : list of str, optional
        List of canonical genes for metagene. If None, uses all canonical genes
    window : int, default=50
        Window size for comparison (codons from start)
    method : str, default='correlation'
        Primary comparison method: 'correlation', 'cosine', or 'combined'
    plot : bool, default=True
        Whether to generate comparison plot
    plot_metrics : bool, default=True
        Whether to show metrics panel alongside coverage plot (only if plot=True)
    figsize : tuple, default=(12, 5)
        Figure size if plotting
    save_path : Path, optional
        Path to save figure

    Returns:
    --------
    dict
        Comparison results containing:
        - pearson_r: Pearson correlation coefficient
        - pearson_p: Pearson p-value
        - spearman_r: Spearman correlation coefficient
        - spearman_p: Spearman p-value
        - cosine_similarity: Cosine similarity score
        - peak_position_alt: Peak position in alternative gene
        - peak_position_canonical: Peak position in canonical metagene
        - peak_distance: Absolute difference in peak positions
        - peak_distance_from_start_canonical: Distance of canonical peak from start codon
        - peak_distance_from_start_alternative: Distance of alternative peak from start codon
        - start_proximity_score: Score based on peak proximity to start (0-1)
        - coverage_alternative: Mean RPM coverage of alternative gene
        - coverage_canonical_mean: Mean RPM coverage across canonical genes
        - coverage_ratio: Ratio of alternative to canonical coverage
        - coverage_adequacy_score: Score based on coverage level (0-1)
        - alternative_profile: Normalized profile of alternative gene
        - canonical_metagene: Normalized metagene profile
        - similarity_score: Combined similarity score (0-1, includes start proximity and coverage)
        - interpretation: Text interpretation of results

    Notes:
    ------
    - High correlation (>0.5) suggests similar initiation pattern
    - Similar peak positions indicate comparable ribosome recruitment
    - Peak proximity to start codon (0-5 codons) is expected for genuine translation
    - Adequate coverage (RPM) indicates active translation vs. noise
    - Cosine similarity measures shape similarity independent of magnitude
    - Combined score integrates multiple metrics including start proximity and coverage
    - Useful for validating novel microprotein predictions
    """
    # Get canonical genes if not provided
    if canonical_gene_ids is None:
        canonical_data = df[df['transcript_class'] == 'Canonical']
        canonical_gene_ids = canonical_data['gene_id'].unique().tolist()

    if len(canonical_gene_ids) == 0:
        raise ValueError("No canonical genes found or provided")

    # Calculate canonical metagene profile and collect coverage data
    canonical_profiles = []
    canonical_coverages = []  # Track mean RPM for each canonical gene
    for gene_id in canonical_gene_ids:
        gene_data = df[df['gene_id'] == gene_id].copy()

        if len(gene_data) == 0:
            continue

        # Average across samples
        occupancy_data = gene_data.groupby('codon_index').agg({
            'RPM': 'mean'
        }).reset_index()

        occupancy_data = occupancy_data.sort_values('codon_index')
        occupancy = occupancy_data['RPM'].values

        # Store coverage (mean RPM before normalization)
        canonical_coverages.append(occupancy.mean())

        # Normalize
        if occupancy.mean() > 0:
            norm_occ = occupancy / occupancy.mean()
        else:
            norm_occ = occupancy

        # Extract window from start
        start_idx = 0
        end_idx = min(len(norm_occ), window)

        if end_idx - start_idx > 0:
            profile = norm_occ[start_idx:end_idx]
            canonical_profiles.append(profile)

    if len(canonical_profiles) == 0:
        raise ValueError("Could not generate canonical metagene profile")

    # Calculate mean canonical profile
    max_len = max(len(p) for p in canonical_profiles)
    padded = []
    for p in canonical_profiles:
        # Pad shorter profiles with NaN for averaging
        if len(p) < max_len:
            padded.append(np.pad(p, (0, max_len - len(p)), constant_values=np.nan))
        else:
            padded.append(p)

    profiles_array = np.array(padded)
    canonical_metagene = np.nanmean(profiles_array, axis=0) # Mean profile
    canonical_sem = np.nanstd(profiles_array, axis=0) / np.sqrt(np.sum(~np.isnan(profiles_array), axis=0)) # SEM

    # Get alternative gene profile
    alt_gene_data = df[df['gene_id'] == alternative_gene_id].copy() # Alternative gene data

    if len(alt_gene_data) == 0:
        raise ValueError(f"No data found for alternative gene {alternative_gene_id}")

    # Average across samples
    alt_occupancy_data = alt_gene_data.groupby('codon_index').agg({
        'RPM': 'mean'
    }).reset_index()

    alt_occupancy_data = alt_occupancy_data.sort_values('codon_index')
    alt_occupancy = alt_occupancy_data['RPM'].values

    # Store alternative gene coverage (mean RPM before normalization)
    alternative_coverage = alt_occupancy.mean()

    # Normalize
    if alt_occupancy.mean() > 0:
        alt_norm_occ = alt_occupancy / alt_occupancy.mean()
    else:
        alt_norm_occ = alt_occupancy

    # Extract window
    alt_start_idx = 0
    alt_end_idx = min(len(alt_norm_occ), window)
    alternative_profile = alt_norm_occ[alt_start_idx:alt_end_idx]

    # Align profiles to same length for comparison
    min_length = min(len(canonical_metagene), len(alternative_profile)) # Minimum length
    canonical_aligned = canonical_metagene[:min_length] # Aligned canonical profile
    alternative_aligned = alternative_profile[:min_length] # Aligned alternative profile

    # Calculate similarity metrics
    results = {
        'alternative_gene': alternative_gene_id,
        'n_canonical_genes': len(canonical_profiles),
        'comparison_length': min_length,
        'alternative_profile': alternative_profile,
        'canonical_metagene': canonical_metagene,
        'canonical_sem': canonical_sem
    }

    # Pearson correlation
    if len(canonical_aligned) > 1:
        pearson_r, pearson_p = stats.pearsonr(canonical_aligned, alternative_aligned)
        results['pearson_r'] = pearson_r
        results['pearson_p'] = pearson_p
    else:
        results['pearson_r'] = np.nan
        results['pearson_p'] = np.nan

    # Spearman correlation
    if len(canonical_aligned) > 1:
        spearman_r, spearman_p = stats.spearmanr(canonical_aligned, alternative_aligned)
        results['spearman_r'] = spearman_r
        results['spearman_p'] = spearman_p
    else:
        results['spearman_r'] = np.nan
        results['spearman_p'] = np.nan

    # Cosine similarity
    dot_product = np.dot(canonical_aligned, alternative_aligned)
    norm_canonical = np.linalg.norm(canonical_aligned)
    norm_alternative = np.linalg.norm(alternative_aligned)

    if norm_canonical > 0 and norm_alternative > 0:
        cosine_sim = dot_product / (norm_canonical * norm_alternative)
        results['cosine_similarity'] = cosine_sim
    else:
        results['cosine_similarity'] = np.nan

    # Peak positions
    canonical_peak_pos = np.argmax(canonical_aligned)
    alternative_peak_pos = np.argmax(alternative_aligned)

    results['peak_position_canonical'] = canonical_peak_pos
    results['peak_position_alternative'] = alternative_peak_pos
    results['peak_distance'] = abs(canonical_peak_pos - alternative_peak_pos)

    # Peak height comparison
    results['peak_height_canonical'] = canonical_aligned[canonical_peak_pos]
    results['peak_height_alternative'] = alternative_aligned[alternative_peak_pos]
    results['peak_height_ratio'] = alternative_aligned[alternative_peak_pos] / (canonical_aligned[canonical_peak_pos] + 1e-6)

    # Peak distance from start codon (position 0)
    # Real translation initiation should show peaks very close to start (0-5 codons)
    results['peak_distance_from_start_canonical'] = canonical_peak_pos
    results['peak_distance_from_start_alternative'] = alternative_peak_pos

    # Calculate start proximity score (1.0 for peak at start, decreasing with distance)
    # Define expected ranges:
    # - 0-3 codons: Excellent (score = 1.0)
    # - 4-6 codons: Good (score = 0.8-0.9)
    # - 7-10 codons: Moderate (score = 0.5-0.7)
    # - >10 codons: Poor (score < 0.5, decreasing)
    def calculate_start_proximity_score(peak_position):
        if peak_position <= 3:
            return 1.0
        elif peak_position <= 6:
            # Linear decay from 1.0 to 0.7
            return 1.0 - (peak_position - 3) * 0.1
        elif peak_position <= 10:
            # Linear decay from 0.7 to 0.3
            return 0.7 - (peak_position - 6) * 0.1
        else:
            # Exponential decay for positions > 10
            return max(0.0, 0.3 * np.exp(-(peak_position - 10) / 10))

    canonical_start_proximity = calculate_start_proximity_score(canonical_peak_pos)
    alternative_start_proximity = calculate_start_proximity_score(alternative_peak_pos)

    # Use the average of both (alternative gene should be similar to canonical pattern)
    results['start_proximity_score_canonical'] = canonical_start_proximity
    results['start_proximity_score_alternative'] = alternative_start_proximity
    results['start_proximity_score'] = np.mean([canonical_start_proximity, alternative_start_proximity])

    # Coverage analysis
    canonical_mean_coverage = np.mean(canonical_coverages)
    results['coverage_alternative'] = alternative_coverage
    results['coverage_canonical_mean'] = canonical_mean_coverage
    results['coverage_ratio'] = alternative_coverage / (canonical_mean_coverage + 1e-6)

    # Calculate coverage adequacy score
    # Real translation should show reasonable coverage, not just noise
    # Score based on:
    # 1. Absolute coverage level (RPM)
    # 2. Relative coverage compared to canonical genes
    def calculate_coverage_adequacy_score(alt_coverage, canonical_mean, ratio):
        # Component 1: Absolute coverage score
        # Very low coverage (<5 RPM) is likely noise
        # Moderate coverage (5-20 RPM) is acceptable
        # High coverage (>20 RPM) is strong evidence
        if alt_coverage < 5:
            absolute_score = alt_coverage / 5  # Linear from 0 to 1
        elif alt_coverage < 20:
            absolute_score = 0.7 + (alt_coverage - 5) / 15 * 0.2  # 0.7 to 0.9
        else:
            absolute_score = min(1.0, 0.9 + (alt_coverage - 20) / 80 * 0.1)  # 0.9 to 1.0

        # Component 2: Relative coverage score
        # Alternative gene should have comparable coverage to canonical genes
        # Ratio of 0.1-2.0 is acceptable (alternative genes can be less expressed)
        if ratio < 0.05:
            relative_score = 0.0  # Too low, likely noise
        elif ratio < 0.1:
            relative_score = ratio / 0.1 * 0.5  # 0 to 0.5
        elif ratio < 2.0:
            relative_score = 0.5 + (min(ratio, 1.0) - 0.1) / 0.9 * 0.5  # 0.5 to 1.0
        else:
            relative_score = max(0.5, 1.0 - (ratio - 2.0) / 3.0 * 0.5)  # Decreasing for very high ratios

        # Weight absolute score more heavily (70/30)
        return 0.7 * absolute_score + 0.3 * relative_score

    results['coverage_adequacy_score'] = calculate_coverage_adequacy_score(
        alternative_coverage, canonical_mean_coverage, results['coverage_ratio']
    )

    # Calculate combined similarity score (0-1 scale)
    # Weight different metrics
    score_components = []

    # Correlation component (0-1, where 1 is perfect positive correlation)
    if not np.isnan(results['pearson_r']):
        score_components.append((results['pearson_r'] + 1) / 2)  # Scale from [-1,1] to [0,1]

    # Cosine similarity component (already 0-1 for positive values)
    if not np.isnan(results['cosine_similarity']):
        score_components.append(max(0, results['cosine_similarity']))

    # Peak position similarity (1 when peaks align, decreases with distance)
    if min_length > 0:
        peak_similarity = 1 - (results['peak_distance'] / min_length)
        score_components.append(max(0, peak_similarity))

    # Start proximity component (critical for validating real translation initiation)
    # Weight this heavily as peak distance from start is a strong biological indicator
    score_components.append(results['start_proximity_score'])
    score_components.append(results['start_proximity_score'])  # Add twice for 2x weight

    # Coverage adequacy component (real translation requires substantial coverage)
    # Weight this moderately - it's important but not as critical as start proximity
    score_components.append(results['coverage_adequacy_score'])

    if score_components:
        results['similarity_score'] = np.mean(score_components)
    else:
        results['similarity_score'] = np.nan

    # Interpretation
    # Consider both overall similarity AND start proximity
    base_interpretation = ""
    if results['similarity_score'] >= 0.7:
        base_interpretation = "STRONG evidence of canonical-like initiation signature"
    elif results['similarity_score'] >= 0.5:
        base_interpretation = "MODERATE evidence of canonical-like initiation signature"
    elif results['similarity_score'] >= 0.3:
        base_interpretation = "WEAK evidence of canonical-like initiation signature"
    else:
        base_interpretation = "MINIMAL evidence of canonical-like initiation signature"

    # Add qualifiers based on start proximity and coverage
    qualifiers = []

    # Start proximity qualifier
    if alternative_peak_pos > 10:
        qualifiers.append(f"Peak at position {alternative_peak_pos}, far from start codon")
    elif alternative_peak_pos > 6:
        qualifiers.append(f"Peak at position {alternative_peak_pos}, moderately distant from start")
    else:
        qualifiers.append(f"Peak at position {alternative_peak_pos}, near start codon")

    # Coverage qualifier
    if alternative_coverage < 5:
        qualifiers.append(f"Very low coverage ({alternative_coverage:.1f} RPM), possible noise")
    elif alternative_coverage < 10:
        qualifiers.append(f"Low coverage ({alternative_coverage:.1f} RPM)")
    elif alternative_coverage >= 20:
        qualifiers.append(f"Strong coverage ({alternative_coverage:.1f} RPM)")
    else:
        qualifiers.append(f"Moderate coverage ({alternative_coverage:.1f} RPM)")

    # Combine interpretation
    results['interpretation'] = base_interpretation + " (" + "; ".join(qualifiers) + ")"

    # Create comparison plot
    if plot:
        # Create figure with 1 or 2 panels based on plot_metrics
        if plot_metrics:
            fig, (ax1, ax2) = plt.subplots(1, 2, figsize=figsize)
        else:
            fig, ax1 = plt.subplots(1, 1, figsize=figsize)

        positions = np.arange(len(canonical_metagene))
        alt_positions = np.arange(len(alternative_profile))

        # Panel 1: Overlaid profiles
        ax1.plot(positions, canonical_metagene, linewidth=3,
                color='#3498DB', label=f'Canonical metagene (n={len(canonical_profiles)})',
                alpha=0.8)
        ax1.fill_between(positions,
                        canonical_metagene - canonical_sem,
                        canonical_metagene + canonical_sem,
                        alpha=0.1, color='#3498DB')

        ax1.plot(alt_positions, alternative_profile, linewidth=2.5,
                color='#E74C3C', label=f'{alternative_gene_id.replace("gene_", "")}',
                alpha=0.8, linestyle='--')

        # Mark expected initiation region (0-5 codons from start)
        ax1.axvspan(-0.5, 5.5, alpha=0.05, color='green',
                   label='Expected initiation region (0-5 codons)')

        # Mark start codon position
        ax1.axvline(0, color='black', linestyle='-', linewidth=2,
                   alpha=0.5, label='Start codon (position 0)')

        # Mark peaks
        ax1.axvline(canonical_peak_pos, color='#3498DB', linestyle=':',
                   linewidth=2, alpha=0.5, label=f'Canonical peak ({canonical_peak_pos})')
        ax1.axvline(alternative_peak_pos, color='#E74C3C', linestyle=':',
                   linewidth=2, alpha=0.5, label=f'Alt peak ({alternative_peak_pos})')

        ax1.set_xlabel('Codon position from start', fontsize=8, fontweight='bold')
        ax1.set_ylabel('Normalized ribosome occupancy', fontsize=8, fontweight='bold')
        ax1.set_title('', fontsize=12, fontweight='bold')
        ax1.legend(loc='best', fontsize=9, framealpha=0.9)
        ax1.grid(alpha=0.3, linestyle='--')
        ax1.spines['top'].set_visible(False)
        ax1.spines['right'].set_visible(False)

        # Panel 2: Metrics summary (optional)
        if plot_metrics:
            ax2.axis('off')

            metrics_text = f"""
Similarity Metrics:

Pearson r: {results['pearson_r']:.3f} (p={results['pearson_p']:.3e})
Spearman r: {results['spearman_r']:.3f} (p={results['spearman_p']:.3e})
Cosine similarity: {results['cosine_similarity']:.3f}

Peak Analysis:
  Canonical peak: position {canonical_peak_pos}
  Alternative peak: position {alternative_peak_pos}
  Peak distance: {results['peak_distance']} codons
  Peak height ratio: {results['peak_height_ratio']:.2f}

Start Proximity Analysis:
  Canonical distance from start: {results['peak_distance_from_start_canonical']} codons
  Alternative distance from start: {results['peak_distance_from_start_alternative']} codons
  Start proximity score: {results['start_proximity_score']:.3f}

Coverage Analysis:
  Alternative coverage: {results['coverage_alternative']:.1f} RPM
  Canonical mean coverage: {results['coverage_canonical_mean']:.1f} RPM
  Coverage ratio: {results['coverage_ratio']:.2f}
  Coverage adequacy score: {results['coverage_adequacy_score']:.3f}

Overall Similarity Score: {results['similarity_score']:.3f}

{results['interpretation']}
            """

            ax2.text(0.1, 0.5, metrics_text,
                    fontsize=10, family='monospace',
                    verticalalignment='center',
                    bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.3))

        plt.tight_layout()

        if save_path:
            fig.savefig(save_path, dpi=300, bbox_inches='tight')
            print(f"Saved comparison to {save_path}")

        results['figure'] = fig

    return results


def batch_compare_to_canonical(df: pd.DataFrame,
                               alternative_gene_ids: Optional[List[str]] = None,
                               canonical_gene_ids: Optional[List[str]] = None,
                               window: int = 50,
                               plot_individual: bool = False,
                               plot_metrics: bool = True,
                               indiv_figsize: Tuple[int, int] = (12, 5),
                               plot_summary: bool = True,
                               output_dir: Optional[Path] = None,
                               figsize: Tuple[int, int] = (14, 10)) -> Tuple[pd.DataFrame, Optional[plt.Figure]]:
    """
    Batch compare multiple alternative genes to the canonical metagene.

    Logic:
    1. If alternative_gene_ids not provided, use all overlapping genes
    2. For each alternative gene, call compare_to_canonical_metagene()
    3. Collect all similarity metrics into a summary DataFrame
    4. Optionally create summary visualization showing all comparisons
    5. Sort results by similarity score (highest to lowest)
    6. Return results table and optional summary figure

    Parameters:
    -----------
    df : pd.DataFrame
        Data from load_and_prepare_data()
    alternative_gene_ids : list of str, optional
        List of alternative genes to validate. If None, uses all overlapping genes
    canonical_gene_ids : list of str, optional
        List of canonical genes for metagene. If None, uses all canonical genes
    window : int, default=50
        Window size for comparison (codons from start)
    plot_individual : bool, default=False
        Whether to save individual comparison plots for each gene
    plot_metrics : bool, default=True
        Whether to show metrics panel in individual plots (only if plot_individual=True)
    indiv_figsize : tuple, default=(12, 5)
        Figure size for individual gene comparison plots
    plot_summary : bool, default=True
        Whether to create a summary plot showing all comparisons
    output_dir : Path, optional
        Directory to save results and plots
    figsize : tuple, default=(14, 10)
        Figure size for summary plot (larger to accommodate 4 panels)

    Returns:
    --------
    tuple : (results_df, summary_figure)
        - results_df: DataFrame with all comparison metrics, sorted by similarity score
          Columns include: gene_id, gene_name, similarity_score, correlation metrics,
          peak positions, peak_distance_from_start (canonical/alternative),
          start_proximity_score, coverage metrics (alternative/canonical/ratio),
          coverage_adequacy_score, interpretation, n_canonical_genes
        - summary_figure: Optional matplotlib figure with summary visualization

    Notes:
    ------
    - Results sorted by similarity_score (descending)
    - Useful for prioritizing which alternative genes are most likely real
    - Summary plot shows 4 panels: similarity scores, correlations, peak distances,
      and start proximity vs coverage with biological thresholds
    - Individual plots saved if plot_individual=True and output_dir provided
    - New metrics include start proximity and coverage analysis for robust validation
    - Genes in upper-left quadrant of panel 4 (near start + high coverage) are strongest candidates
    """
    # Get alternative genes if not provided
    if alternative_gene_ids is None:
        alt_data = df[df['transcript_class'] == 'Overlapping']
        alternative_gene_ids = alt_data['gene_id'].unique().tolist()

    if len(alternative_gene_ids) == 0:
        raise ValueError("No alternative genes found or provided")

    print(f"Comparing {len(alternative_gene_ids)} alternative genes to canonical metagene...")

    # Collect results for all genes
    all_results = []

    for alt_gene in alternative_gene_ids:
        print(f"  Analyzing {alt_gene}...", end=" ")

        try:
            # Determine output path for individual plot
            if plot_individual and output_dir:
                output_dir_path = Path(output_dir)
                output_dir_path.mkdir(parents=True, exist_ok=True)
                gene_safe = alt_gene.replace('gene_', '').replace('/', '_')
                save_path = output_dir_path / f"{gene_safe}_vs_canonical.jpg"
            else:
                save_path = None

            # Compare to canonical
            result = compare_to_canonical_metagene(
                df,
                alt_gene,
                canonical_gene_ids=canonical_gene_ids,
                window=window,
                plot=plot_individual,
                plot_metrics=plot_metrics,
                save_path=save_path,
                figsize=indiv_figsize
            )

            # Extract key metrics
            all_results.append({
                'gene_id': alt_gene,
                'gene_name': alt_gene.replace('gene_', ''),
                'similarity_score': result['similarity_score'],
                'pearson_r': result['pearson_r'],
                'pearson_p': result['pearson_p'],
                'spearman_r': result['spearman_r'],
                'spearman_p': result['spearman_p'],
                'cosine_similarity': result['cosine_similarity'],
                'peak_position_canonical': result['peak_position_canonical'],
                'peak_position_alternative': result['peak_position_alternative'],
                'peak_distance': result['peak_distance'],
                'peak_height_ratio': result['peak_height_ratio'],
                'peak_distance_from_start_canonical': result['peak_distance_from_start_canonical'],
                'peak_distance_from_start_alternative': result['peak_distance_from_start_alternative'],
                'start_proximity_score': result['start_proximity_score'],
                'coverage_alternative': result['coverage_alternative'],
                'coverage_canonical_mean': result['coverage_canonical_mean'],
                'coverage_ratio': result['coverage_ratio'],
                'coverage_adequacy_score': result['coverage_adequacy_score'],
                'interpretation': result['interpretation'],
                'n_canonical_genes': result['n_canonical_genes']
            })

            print(f"Score: {result['similarity_score']:.3f}")

        except Exception as e:
            print(f"ERROR: {str(e)}")
            all_results.append({
                'gene_id': alt_gene,
                'gene_name': alt_gene.replace('gene_', ''),
                'similarity_score': np.nan,
                'pearson_r': np.nan,
                'pearson_p': np.nan,
                'spearman_r': np.nan,
                'spearman_p': np.nan,
                'cosine_similarity': np.nan,
                'peak_position_canonical': np.nan,
                'peak_position_alternative': np.nan,
                'peak_distance': np.nan,
                'peak_height_ratio': np.nan,
                'peak_distance_from_start_canonical': np.nan,
                'peak_distance_from_start_alternative': np.nan,
                'start_proximity_score': np.nan,
                'coverage_alternative': np.nan,
                'coverage_canonical_mean': np.nan,
                'coverage_ratio': np.nan,
                'coverage_adequacy_score': np.nan,
                'interpretation': f'ERROR: {str(e)}',
                'n_canonical_genes': np.nan
            })

    # Create results DataFrame
    results_df = pd.DataFrame(all_results)

    # Sort by similarity score (descending)
    results_df = results_df.sort_values('similarity_score', ascending=False, na_position='last')
    results_df = results_df.reset_index(drop=True)

    # Save results table
    if output_dir:
        output_dir_path = Path(output_dir)
        output_dir_path.mkdir(parents=True, exist_ok=True)
        results_df.to_csv(output_dir_path / "alternative_genes_validation_summary.csv", index=False)
        print(f"\nSaved results to {output_dir_path / 'alternative_genes_validation_summary.csv'}")

    # Create summary plot
    summary_fig = None
    if plot_summary and len(results_df) > 0:
        summary_fig = _create_validation_summary_plot(results_df, figsize=figsize)

        if output_dir:
            summary_fig.savefig(output_dir_path / "alternative_genes_validation_summary.jpg",
                              dpi=300, bbox_inches='tight')
            print(f"Saved summary plot to {output_dir_path / 'alternative_genes_validation_summary.jpg'}")

    return results_df, summary_fig


def _create_validation_summary_plot(results_df: pd.DataFrame,
                                    figsize: Tuple[int, int] = (14, 10)) -> plt.Figure:
    """
    Create summary visualization for batch validation results.

    Creates 4-panel plot showing:
    1. Similarity scores (bar chart)
    2. Correlation metrics (scatter plot)
    3. Peak position similarity (bar chart)
    4. Start proximity vs coverage (scatter plot with biological thresholds)

    Internal helper function for batch_compare_to_canonical.
    """
    # Remove rows with NaN similarity scores
    plot_df = results_df[~results_df['similarity_score'].isna()].copy()

    if len(plot_df) == 0:
        print("No valid results to plot")
        return None

    # Create figure with 4 panels
    fig = plt.figure(figsize=figsize)
    gs = fig.add_gridspec(3, 2, hspace=0.4, wspace=0.3)
    ax1 = fig.add_subplot(gs[0, :])  # Top: full width
    ax2 = fig.add_subplot(gs[1, 0])   # Middle left
    ax3 = fig.add_subplot(gs[1, 1])   # Middle right
    ax4 = fig.add_subplot(gs[2, :])   # Bottom: full width

    # Color code by interpretation
    color_map = {
        'STRONG': '#27AE60',      # Green
        'MODERATE': '#F39C12',    # Orange
        'WEAK': '#E67E22',        # Dark orange
        'MINIMAL': '#E74C3C'      # Red
    }

    def get_color(interpretation):
        for key in color_map:
            if key in interpretation:
                return color_map[key]
        return '#95A5A6'  # Gray for unknown

    plot_df['color'] = plot_df['interpretation'].apply(get_color)

    # Panel 1: Similarity scores bar chart
    y_pos = np.arange(len(plot_df))
    ax1.barh(y_pos, plot_df['similarity_score'].values,
            color=plot_df['color'].values, alpha=0.8, edgecolor='black', linewidth=1)

    ax1.set_yticks(y_pos)
    ax1.set_yticklabels(plot_df['gene_name'].values, fontsize=9)
    ax1.set_xlabel('Similarity Score', fontsize=11, fontweight='bold')
    ax1.set_title('Alternative Genes Validation - Similarity to Canonical Metagene',
                 fontsize=13, fontweight='bold', pad=10)
    ax1.axvline(0.7, color='green', linestyle='--', linewidth=1.5, alpha=0.5, label='Strong threshold')
    ax1.axvline(0.5, color='orange', linestyle='--', linewidth=1.5, alpha=0.5, label='Moderate threshold')
    ax1.axvline(0.3, color='red', linestyle='--', linewidth=1.5, alpha=0.5, label='Weak threshold')
    ax1.set_xlim(0, 1.0)
    ax1.legend(loc='lower right', fontsize=8)
    ax1.grid(alpha=0.3, axis='x', linestyle='--')
    ax1.spines['top'].set_visible(False)
    ax1.spines['right'].set_visible(False)

    # Panel 2: Correlation scatter plot
    ax2.scatter(plot_df['pearson_r'].values, plot_df['cosine_similarity'].values,
               s=150, c=plot_df['color'].values, alpha=0.7,
               edgecolors='black', linewidths=1.5)

    # Add gene labels
    for idx, row in plot_df.iterrows():
        ax2.annotate(row['gene_name'],
                    (row['pearson_r'], row['cosine_similarity']),
                    xytext=(5, 5), textcoords='offset points',
                    fontsize=8, alpha=0.8)

    ax2.set_xlabel('Pearson Correlation', fontsize=10, fontweight='bold')
    ax2.set_ylabel('Cosine Similarity', fontsize=10, fontweight='bold')
    ax2.set_title('Correlation Metrics', fontsize=11, fontweight='bold')
    ax2.axhline(0.5, color='gray', linestyle='--', linewidth=1, alpha=0.3)
    ax2.axvline(0.5, color='gray', linestyle='--', linewidth=1, alpha=0.3)
    ax2.grid(alpha=0.3, linestyle='--')
    ax2.spines['top'].set_visible(False)
    ax2.spines['right'].set_visible(False)

    # Panel 3: Peak distance
    y_pos = np.arange(len(plot_df))
    bars = ax3.barh(y_pos, plot_df['peak_distance'].values,
                   color=plot_df['color'].values, alpha=0.8,
                   edgecolor='black', linewidth=1)

    ax3.set_yticks(y_pos)
    ax3.set_yticklabels(plot_df['gene_name'].values, fontsize=9)
    ax3.set_xlabel('Peak Distance (codons)', fontsize=10, fontweight='bold')
    ax3.set_title('Peak Position Similarity', fontsize=11, fontweight='bold')
    ax3.axvline(10, color='orange', linestyle='--', linewidth=1.5, alpha=0.5,
               label='10 codon threshold')
    ax3.legend(loc='lower right', fontsize=8)
    ax3.grid(alpha=0.3, axis='x', linestyle='--')
    ax3.spines['top'].set_visible(False)
    ax3.spines['right'].set_visible(False)

    # Panel 4: Start proximity and coverage scatter plot
    scatter = ax4.scatter(plot_df['peak_distance_from_start_alternative'].values,
                         plot_df['coverage_alternative'].values,
                         s=200, c=plot_df['color'].values, alpha=0.7,
                         edgecolors='black', linewidths=1.5)

    # Add gene labels
    for idx, row in plot_df.iterrows():
        ax4.annotate(row['gene_name'],
                    (row['peak_distance_from_start_alternative'], row['coverage_alternative']),
                    xytext=(5, 5), textcoords='offset points',
                    fontsize=8, alpha=0.8)

    # Mark expected initiation region
    ax4.axvspan(-0.5, 5.5, alpha=0.15, color='green', label='Expected initiation region (0-5 codons)')

    # Mark low coverage threshold
    ax4.axhline(5, color='orange', linestyle='--', linewidth=1.5, alpha=0.5,
               label='Low coverage threshold (5 RPM)')

    ax4.set_xlabel('Peak Distance from Start Codon (codons)', fontsize=10, fontweight='bold')
    ax4.set_ylabel('Mean Coverage (RPM, log scale)', fontsize=10, fontweight='bold')
    ax4.set_title('Start Proximity vs Coverage - Biological Validation Metrics',
                 fontsize=11, fontweight='bold')
    ax4.legend(loc='best', fontsize=8)
    ax4.grid(alpha=0.3, linestyle='--', which='both')  # Show grid for both major and minor ticks
    ax4.spines['top'].set_visible(False)
    ax4.spines['right'].set_visible(False)

    # Set logarithmic scale for y-axis
    ax4.set_yscale('log')

    # Set reasonable axis limits
    ax4.set_xlim(-1, max(15, plot_df['peak_distance_from_start_alternative'].max() + 2))
    # For log scale, set lower limit to a small positive value instead of 0
    y_min = max(0.1, plot_df['coverage_alternative'].min() * 0.5)
    y_max = plot_df['coverage_alternative'].max() * 2
    ax4.set_ylim(y_min, y_max)

    plt.tight_layout()

    return fig


def analyze_all_overlapping_regions(table_path: Union[str, Path],
                                   output_dir: Optional[Path] = None,
                                   control_table_path: Optional[Union[str, Path]] = None,
                                   sep: str = ',') -> Dict:
    """
    Analyze all known overlapping regions.

    Logic:
    1. Load the main dataset
    2. Optionally load control dataset
    3. Analyze each known overlapping region:
       - MT-ATP8/MT-ATP6
       - MT-ND4L/MT-ND4
    4. Generate plots and export results

    Parameters:
    -----------
    table_path : str or Path
        Path to cumsum table CSV file
    output_dir : Path, optional
        Directory to save results
    control_table_path : str or Path, optional
        Path to control condition cumsum table
    sep : str, default=','
        CSV separator

    Returns:
    --------
    dict
        Analysis results for all regions

    Notes:
    ------
    - Analyzes all known bicistronic regions in human mitochondria
    - Start positions are approximate and may need adjustment
    - Compares to control if provided
    """
    # Load data
    df = load_and_prepare_data(table_path, sep=sep)

    # Load control if provided
    df_control = None
    if control_table_path:
        df_control = load_and_prepare_data(control_table_path, sep=sep)

    results = {}

    # Known overlapping regions and their approximate start positions
    # Note: These positions may need adjustment based on your specific annotation
    region_configs = {
        'MT-ATP8': {
            'upstream_start': 0,      # MT-ATP8 start
            'downstream_start': 43,   # MT-ATP6 start (overlaps by ~46 nt = ~15 codons)
            'description': 'MT-ATP8/MT-ATP6 overlap'
        },
        'MT-ND4L': {
            'upstream_start': 0,      # MT-ND4L start
            'downstream_start': 95,   # MT-ND4 start (overlaps by ~7 nt = ~2 codons)
            'description': 'MT-ND4L/MT-ND4 overlap'
        }
    }

    for gene_id, config in region_configs.items():
        print(f"\nAnalyzing {config['description']}...")

        # Analyze
        region_results = analyze_bicistronic_region(
            df,
            gene_id,
            config['upstream_start'],
            config['downstream_start']
        )

        if region_results:
            results[gene_id] = region_results

            print(f"  Upstream initiation: {region_results.get('upstream_initiation', 'N/A'):.2f}")
            print(f"  Downstream initiation: {region_results.get('downstream_initiation', 'N/A'):.2f}")
            print(f"  Initiation ratio: {region_results.get('initiation_ratio', 'N/A'):.2f}")

            # Plot
            fig = plot_bicistronic_initiation(
                df,
                gene_id,
                config['upstream_start'],
                config['downstream_start'],
                df_control=df_control
            )

            if output_dir and fig:
                output_dir = Path(output_dir)
                output_dir.mkdir(parents=True, exist_ok=True)
                fig.savefig(
                    output_dir / f"{gene_id}_bicistronic.pdf",
                    dpi=300, bbox_inches='tight'
                )
                plt.close(fig)

            # Export enrichment stats
            if output_dir and 'enrichment_stats' in region_results:
                enrichment_df = region_results['enrichment_stats']
                enrichment_df.to_csv(
                    output_dir / f"{gene_id}_enrichment.csv",
                    index=False
                )

    return results


if __name__ == "__main__":
    print("Mitochondrial Initiation Analysis Extension")
    print("=" * 60)
    print("\nThis module extends mito_riboseq_analysis with initiation-specific functions.")
    print("\nKey functions:")
    print("  - detect_initiation_peaks(): Detect initiation sites via peak finding")
    print("  - test_start_codon_enrichment(): Test specific positions for enrichment")
    print("  - analyze_bicistronic_region(): Analyze overlapping ORFs")
    print("  - compare_initiation_conditions(): Compare two conditions")
    print("  - plot_initiation_profile(): Visualize initiation patterns")
    print("  - plot_bicistronic_initiation(): Plot overlapping ORF initiation")
    print("  - plot_metagene_initiation(): Create metagene plots")
    print("\nConvenience functions:")
    print("  - quick_initiation_analysis(): One-line analysis for a gene")
    print("  - batch_initiation_analysis(): Analyze multiple genes in facet grid")
    print("  - compare_to_canonical_metagene(): Validate alternative genes vs canonical")
    print("  - batch_compare_to_canonical(): Batch validate all alternative genes")
    print("  - analyze_all_overlapping_regions(): Batch analyze all overlapping regions")
    print("\nExample usage:")
    print("  # Single gene analysis")
    print("  from mito_riboseq_analysis_extension import quick_initiation_analysis")
    print("  peaks, fig = quick_initiation_analysis('data.csv', 'MT-ATP8', known_starts=[0, 43])")
    print("\n  # Batch analysis with facet grid")
    print("  from mito_riboseq_analysis_extension import batch_initiation_analysis, load_and_prepare_data")
    print("  df = load_and_prepare_data('data.csv')")
    print("  genes = ['MT-ND1', 'MT-ND2', 'gene_MOTSc', 'gene_SHLP2']")
    print("  peaks_dict, fig = batch_initiation_analysis(df, genes)")
    print("\n  # Validate single alternative gene")
    print("  from mito_riboseq_analysis_extension import compare_to_canonical_metagene")
    print("  results = compare_to_canonical_metagene(df, 'gene_MOTSc')")
    print("  print(f'Similarity score: {results[\"similarity_score\"]:.3f}')")
    print("  print(results['interpretation'])")
    print("\n  # Batch validate all alternative genes")
    print("  from mito_riboseq_analysis_extension import batch_compare_to_canonical")
    print("  results_df, fig = batch_compare_to_canonical(df, output_dir='validation_results/')")
    print("  print(results_df[['gene_name', 'similarity_score', 'interpretation']])")
