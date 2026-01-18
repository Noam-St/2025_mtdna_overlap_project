"""
Mitochondrial Ribosome Profiling Analysis Module

This script provides functions to analyze translation initiation signatures
from mitochondrial ribosome profiling data, comparing overlapping reading 
frames (microproteins) to canonical mitochondrial genes.

Key concepts:
- RPM_cumsum_normgene: Cumulative sum of RPM-normalized ribosome footprint reads
- Initiation signature: The pattern of ribosome density at the 5' end of transcripts
- Overlapping genes: Alternative reading frames (gene_*), e.g., MOTSc, SHLP family
- Canonical genes: Standard mitochondrial protein-coding genes (ENST*)

Author: Generated for mitochondrial microprotein research
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path


def load_and_prepare_data(table_path, sep = ','):
    """
    Load the ribosome profiling cumulative sum data and categorize transcripts.
    
    Logic:
    1. Read CSV with ribosome profiling data
    2. Classify transcripts as 'Overlapping' (gene_*) or 'Canonical' (ENST*)
    3. Return DataFrame with additional 'transcript_class' column
    
    Parameters:
    -----------
    table_path : str or Path
        Path to the tabular file containing mito_cumsum_table data
        
    Returns:
    --------
    pd.DataFrame
        Original data with added 'transcript_class' column
        
    Notes:
    ------
    - Overlapping genes start with 'gene_' prefix
    - Canonical genes use Ensembl transcript IDs (ENST)
    - This classification is critical for downstream comparative analysis
    """
    df = pd.read_csv(table_path, sep=sep)
    
    # Classify transcripts based on their naming convention
    # Overlapping reading frames/microproteins: start with 'gene_'
    # Canonical genes: Ensembl IDs starting with 'ENST'
    df['transcript_class'] = df['transcript_id'].apply(
        lambda x: 'Overlapping' if x.startswith('gene_') else 'Canonical'
    )
    
    return df


def calculate_initiation_metrics(df, window_size=30):
    """
    Calculate metrics for translation initiation signatures.
    
    Logic:
    1. For each transcript-sample combination, extract the first N codons
    2. Calculate multiple metrics:
       - peak_position: Codon with HIGHEST individual RPM (ribosome density)
       - peak_rpm: The RPM value at that peak position
       - max_cumsum_change: Position with largest increase in cumsum (steepest slope)
       - initiation_slope: Overall rate of ribosome accumulation
    3. Return summary statistics for comparative analysis
    
    IMPORTANT: We use RPM (not RPM_cumsum_normgene) for peak_position because:
    - RPM = ribosome density at each codon
    - RPM_cumsum_normgene = cumulative sum (always increasing, max at end)
    - Peak initiation signature = codon with highest individual density
    
    Parameters:
    -----------
    df : pd.DataFrame
        Data from load_and_prepare_data()
    window_size : int, default=30
        Number of codons from start to analyze for initiation signature
        (Based on typical elongation phase analysis starting at +30bp)
        
    Returns:
    --------
    pd.DataFrame
        Summary metrics per transcript and sample:
        - peak_position: Codon index with maximum RPM (highest density)
        - peak_rpm: The RPM value at the peak position
        - max_cumsum_change_position: Codon with largest increase in cumsum
        - max_cumsum_change: Size of that increase
        - final_cumsum: Total cumulative RPM in window
        - initiation_slope: Linear fit slope of cumsum
        
    Notes:
    ------
    - Window size of 30 codons (~90 nt) captures the initiation/early elongation region
    - For mt-mRNAs, peaks often appear downstream (+30bp) rather than at position 0
    - Peak position uses RPM to find where ribosome density is highest
    - Cumsum change identifies where ribosomes accumulate most rapidly
    """
    metrics = []
    
    for (transcript, sample, trans_class), group in df.groupby(
        ['transcript_id', 'sample', 'transcript_class']
    ):
        # Extract initiation window (first N codons)
        window_data = group[group['codon_index'] <= window_size].copy()
        
        if len(window_data) == 0:
            continue
        
        # Sort by codon_index to ensure proper order
        window_data = window_data.sort_values('codon_index')
            
        # CORRECTED: Find peak position using RPM (individual density), not cumsum
        peak_idx = window_data.loc[window_data['RPM_normgene'].idxmax(), 'codon_index']
        peak_rpm_val = window_data['RPM_normgene'].max()
        
        # Calculate position with largest change in cumsum (steepest increase)
        # This identifies where ribosomes are accumulating most rapidly
        if len(window_data) > 1:
            window_data['cumsum_diff'] = window_data['RPM_cumsum_normgene'].diff()
            max_change_idx = window_data.loc[window_data['cumsum_diff'].idxmax(), 'codon_index']
            max_change_val = window_data['cumsum_diff'].max()
        else:
            max_change_idx = np.nan
            max_change_val = np.nan
        
        # Final cumulative value in window
        final_cumsum = window_data['RPM_cumsum_normgene'].max()
        
        # Calculate slope of cumulative sum (overall rate of ribosome accumulation)
        # Linear fit: RPM_cumsum_normgene ~ codon_index
        if len(window_data) > 1:
            slope = np.polyfit(window_data['codon_index'], 
                             window_data['RPM_cumsum_normgene'], 1)[0]
        else:
            slope = np.nan
            
        metrics.append({
            'transcript_id': transcript,
            'sample': sample,
            'transcript_class': trans_class,
            'peak_position': peak_idx,
            'peak_rpm': peak_rpm_val,
            'max_cumsum_change_position': max_change_idx,
            'max_cumsum_change': max_change_val,
            'final_cumsum': final_cumsum,
            'initiation_slope': slope
        })
    
    return pd.DataFrame(metrics)


def plot_per_codon_profile(df, transcript_id, 
                           figsize=(10, 6),
                           colors=None,
                           save_path=None,
                           show_window=50):
    """
    Plot per-codon cumulative RPM profile for a specific transcript.
    
    Logic:
    1. Filter data for specified transcript
    2. Plot RPM_cumsum_normgene vs codon_index for each sample (replicate)
    3. Visualize the ribosome density accumulation along the transcript
    
    Parameters:
    -----------
    df : pd.DataFrame
        Data from load_and_prepare_data()
    transcript_id : str
        Specific transcript to plot (e.g., 'gene_MOTSc')
    figsize : tuple, default=(10, 6)
        Figure size (width, height) in inches
    colors : list or None, default=None
        Colors for each sample. If None, uses default palette
    save_path : str or None, default=None
        If provided, save figure to this path
    show_window : int or None, default=50
        If provided, only show first N codons (useful for initiation region)
        If None, show entire transcript
        
    Returns:
    --------
    matplotlib.figure.Figure
        The generated figure object
        
    Notes:
    ------
    - RPM_cumsum_normgene increases monotonically (cumulative sum property)
    - Steep increases indicate regions of high ribosome density
    - For initiation analysis, focus on first 30-50 codons
    """
    fig, ax = plt.subplots(figsize=figsize)
    
    # Filter for specific transcript
    transcript_data = df[df['transcript_id'] == transcript_id].copy()
    transcript_name = transcript_data['gene_id'].iloc[0] if 'gene_id' in transcript_data.columns else transcript_id
    
    if len(transcript_data) == 0:
        print(f"Warning: No data found for transcript {transcript_id}")
        return fig
    
    # Apply window if specified
    if show_window is not None:
        transcript_data = transcript_data[transcript_data['codon_index'] <= show_window]
    
    # Get unique samples
    samples = transcript_data['sample'].unique()
    
    # Set colors
    if colors is None:
        colors = sns.color_palette("husl", len(samples))
    
    # Plot each sample
    for idx, sample in enumerate(samples):
        sample_data = transcript_data[transcript_data['sample'] == sample]
        ax.plot(sample_data['codon_index'], 
               sample_data['RPM'],
               marker='o',
               markersize=3,
               label=sample,
               color=colors[idx],
               linewidth=2,
               alpha=0.8)
    
    # Get transcript class for title
    trans_class = transcript_data['transcript_class'].iloc[0]
    
    # Add more xaxis ticks if window is small
    if show_window is not None and show_window <= 100:
        ax.set_xticks(np.arange(0, show_window + 1, 5))

    # Formatting
    ax.set_xlabel('Codon Index (Position)', fontsize=12, fontweight='bold')
    ax.set_ylabel('RPM', fontsize=12, fontweight='bold')
    ax.set_title(f'{transcript_name} ({trans_class})\nPer-Codon Ribosome Profiling', 
                fontsize=14, fontweight='bold')
    ax.legend(title='Sample', frameon=True, fancybox=True)
    ax.grid(True, alpha=0.3, linestyle='--')
    
    # Despine and remove grid lines for cleaner look
    sns.despine(ax=ax)
    ax.grid(False)

    plt.tight_layout()
    
    # Save if path provided
    if save_path:
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
        print(f"Saved figure to {save_path}")
    
    return fig


def plot_all_transcripts(df, output_dir=None, 
                         figsize=(10, 6),
                         colors=None,
                         show_window=50):
    """
    Generate per-codon profiles for all transcripts in the dataset.
    
    Logic:
    1. Iterate through all unique transcripts
    2. Call plot_per_codon_profile() for each
    3. Optionally save all figures to output directory
    
    Parameters:
    -----------
    df : pd.DataFrame
        Data from load_and_prepare_data()
    output_dir : str or None, default=None
        Directory to save all figures. If None, figures are only displayed
    figsize : tuple, default=(10, 6)
        Figure size for each plot
    colors : list or None, default=None
        Colors for samples
    show_window : int or None, default=50
        Window size for each plot
        
    Returns:
    --------
    dict
        Dictionary mapping transcript_id to figure object
        
    Notes:
    ------
    - Useful for batch generation of all transcript profiles
    - Creates one figure per transcript
    - Organize output by creating subdirectories for overlapping vs canonical
    """
    figures = {}
    
    # Create output directory if needed
    if output_dir:
        output_path = Path(output_dir)
        output_path.mkdir(parents=True, exist_ok=True)
        
        # Create subdirectories for each class
        (output_path / 'overlapping').mkdir(exist_ok=True)
        (output_path / 'canonical').mkdir(exist_ok=True)
    
    # Generate plots for each transcript
    for transcript in df['transcript_id'].unique():
        trans_class = df[df['transcript_id'] == transcript]['transcript_class'].iloc[0]
        gene_id = df[df['transcript_id'] == transcript]['gene_id'].iloc[0] if 'gene_id' in df.columns else transcript
        # Determine save path
        if output_dir:
            subdir = 'overlapping' if trans_class == 'Overlapping' else 'canonical'
            save_path = output_path / subdir / f'{gene_id}_profile.png'
        else:
            save_path = None
        
        # Generate plot
        fig = plot_per_codon_profile(
            df, 
            transcript,
            figsize=figsize,
            colors=colors,
            save_path=save_path,
            show_window=show_window
        )
        
        figures[transcript] = fig
        plt.close(fig)  # Close to save memory
    
    print(f"Generated {len(figures)} transcript profiles")
    return figures


def plot_initiation_summary(df, window_size=30,
                            figsize=(14, 10),
                            colors=None,
                            save_path=None):
    """
    Create summary plots comparing initiation signatures between transcript classes.
    
    Logic:
    1. Extract initiation window for all transcripts
    2. Calculate average profiles for Overlapping vs Canonical
    3. Create multi-panel figure:
       - Panel A: Mean cumsum profiles with confidence intervals
       - Panel B: Peak position distribution (where maximum density occurs)
       - Panel C: Initiation slope comparison (rate of density increase)
       - Panel D: Heatmap of individual transcript patterns
    
    Parameters:
    -----------
    df : pd.DataFrame
        Data from load_and_prepare_data()
    window_size : int, default=30
        Number of codons to include in initiation analysis
    figsize : tuple, default=(14, 10)
        Figure size (width, height) for the multi-panel plot
    colors : dict or None, default=None
        Dictionary mapping 'Overlapping' and 'Canonical' to colors
        Example: {'Overlapping': '#E74C3C', 'Canonical': '#3498DB'}
    save_path : str or None, default=None
        If provided, save the summary figure to this path
        
    Returns:
    --------
    matplotlib.figure.Figure
        Multi-panel summary figure
        
    Notes:
    ------
    - This is the key comparative analysis between overlapping and canonical genes
    - Differences in peak position reveal distinct translation initiation patterns
    - For mt-mRNAs, overlapping genes may show different initiation dynamics
    - Statistical differences can indicate functional translation of alternative ORFs
    """
    # Set default colors
    if colors is None:
        colors = {
            'Overlapping': '#E74C3C',  # Red
            'Canonical': '#3498DB'      # Blue
        }
    
    # Create figure with subplots
    fig = plt.figure(figsize=figsize)
    gs = fig.add_gridspec(2, 2, hspace=0.3, wspace=0.3)
    ax1 = fig.add_subplot(gs[0, :])  # Top: full width
    ax2 = fig.add_subplot(gs[1, 0])
    ax3 = fig.add_subplot(gs[1, 1])
    
    # === Panel A: Mean Relative RPM profiles ===
    for trans_class in ['Overlapping', 'Canonical']:
        class_data = df[df['transcript_class'] == trans_class].copy()
        
        # Filter to initiation window
        init_data = class_data[class_data['codon_index'] <= window_size]
        
        # Calculate mean and SEM across transcripts for each codon position
        # Group by codon_index and sample, then aggregate
        summary = init_data.groupby(['codon_index', 'sample']).agg({
            'RPM_normgene': ['mean', 'sem']
        }).reset_index()
        summary.columns = ['codon_index', 'sample', 'mean', 'sem']
        
        # Average across samples
        final_summary = summary.groupby('codon_index').agg({
            'mean': 'mean',
            'sem': lambda x: np.sqrt(np.sum(x**2)) / len(x)  # Propagate error
        }).reset_index()
        
        # Plot mean with confidence interval
        ax1.plot(final_summary['codon_index'], 
                final_summary['mean'],
                label=trans_class,
                color=colors[trans_class],
                linewidth=2.5)
        
        ax1.fill_between(
            final_summary['codon_index'],
            final_summary['mean'] - final_summary['sem'],
            final_summary['mean'] + final_summary['sem'],
            alpha=0.2,
            color=colors[trans_class]
        )
    
    ax1.set_xlabel('Codon Index', fontsize=12, fontweight='bold')
    ax1.set_ylabel('Mean Relative RPM', fontsize=12, fontweight='bold')
    ax1.set_title('Initiation Signature Comparison', fontsize=14, fontweight='bold')
    ax1.legend(frameon=True, fancybox=True)
    ax1.grid(True, alpha=0.3, linestyle='--')
    
    # Remove top and right spines
    sns.despine(ax=ax1)
    
    # === Panel B: Peak position distribution ===
    metrics = calculate_initiation_metrics(df, window_size)
    
    for trans_class in ['Overlapping', 'Canonical']:
        class_metrics = metrics[metrics['transcript_class'] == trans_class]
        ax2.hist(class_metrics['peak_position'], 
                bins=np.arange(0.5, window_size + 1.5, 1),
                alpha=0.6,
                label=trans_class,
                color=colors[trans_class],
                edgecolor='black')
    
    ax2.set_xlabel('Peak Position (Codon Index)', fontsize=11, fontweight='bold')
    ax2.set_ylabel('Frequency', fontsize=11, fontweight='bold')
    ax2.set_title('Distribution of Peak Positions', fontsize=12, fontweight='bold')
    ax2.legend()
    ax2.grid(True, alpha=0.3, axis='y', linestyle='--')
    
    # === Panel C: Initiation slope comparison ===
    plot_data = []
    for trans_class in ['Overlapping', 'Canonical']:
        class_metrics = metrics[metrics['transcript_class'] == trans_class]
        for slope in class_metrics['initiation_slope'].dropna():
            plot_data.append({
                'Class': trans_class,
                'Slope': slope
            })
    
    plot_df = pd.DataFrame(plot_data)
    
    # Violin plot for slope comparison
    parts = ax3.violinplot(
        [plot_df[plot_df['Class'] == 'Overlapping']['Slope'].values,
         plot_df[plot_df['Class'] == 'Canonical']['Slope'].values],
        positions=[0, 1],
        showmeans=True,
        showmedians=True
    )
    
    # Color the violin plots
    for idx, pc in enumerate(parts['bodies']):
        class_name = ['Overlapping', 'Canonical'][idx]
        pc.set_facecolor(colors[class_name])
        pc.set_alpha(0.6)
    
    ax3.set_xticks([0, 1])
    ax3.set_xticklabels(['Overlapping', 'Canonical'])
    ax3.set_ylabel('Initiation Slope (RPM/codon)', fontsize=11, fontweight='bold')
    ax3.set_title('Rate of Ribosome Density Increase', fontsize=12, fontweight='bold')
    ax3.grid(True, alpha=0.3, axis='y', linestyle='--')
    
    plt.tight_layout()
    
    # Save if path provided
    if save_path:
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
        print(f"Saved summary figure to {save_path}")
    
    return fig


def plot_initiation_ridgeplot(df, window_size=30,
                            figsize=(12, 10),
                            cmap='YlOrRd',
                            save_path=None,
                            normalize_per_transcript=True):
    """
    Create ridgeplot showing initiation patterns for all transcripts.
    
    Logic:
    1. Extract initiation window (first N codons) for all transcripts
    2. Create overlapping density curves for each transcript
    3. Optionally normalize each transcript to [0, 1] for pattern visualization
    4. Separate ridgeplots for Overlapping vs Canonical genes
    
    Parameters:
    -----------
    df : pd.DataFrame
        Data from load_and_prepare_data()
    window_size : int, default=30
        Number of codons to include
    figsize : tuple, default=(12, 10)
        Figure size
    cmap : str, default='YlOrRd'
        Colormap for ridgeplot (Yellow-Orange-Red shows density well)
    save_path : str or None, default=None
        Save path for figure
    normalize_per_transcript : bool, default=True
        If True, normalize each transcript's RPM to [0, 1]
        Allows comparison of initiation patterns independent of absolute density
        
    Returns:
    --------
    matplotlib.figure.Figure
        Ridgeplot figure
        
    Notes:
    ------
    - Normalization reveals initiation pattern shapes
    - Overlapping curves allow visual comparison of pattern similarity
    - Distinct patterns between overlapping and canonical may indicate
      different translation initiation mechanisms
    """
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=figsize)
    
    for idx, (trans_class, ax) in enumerate([('Overlapping', ax1), 
                                              ('Canonical', ax2)]):
        # Filter data
        class_data = df[df['transcript_class'] == trans_class].copy()
        init_data = class_data[class_data['codon_index'] <= window_size]
        
        # Get unique transcripts
        transcripts = init_data['gene_id'].unique()
        n_transcripts = len(transcripts)
        
        # Get colormap
        colors = plt.cm.get_cmap(cmap)(np.linspace(0.3, 0.9, n_transcripts))
        
        # Calculate vertical spacing
        vertical_spacing = 1.0
        
        # Plot each transcript as a ridge
        for i, transcript in enumerate(transcripts):
            transcript_data = init_data[init_data['gene_id'] == transcript]
            
            # Average across samples for each codon position
            profile = transcript_data.groupby('codon_index')['RPM_normgene'].mean()
            
            # Normalize if requested
            if normalize_per_transcript and profile.max() > 0:
                profile = profile / profile.max()
            
            # Create x-axis (codon indices)
            x = profile.index.values
            y = profile.values
            
            # Offset vertically for ridge effect
            baseline = i * vertical_spacing
            y_offset = y + baseline
            
            # Fill area under curve
            ax.fill_between(x, baseline, y_offset, 
                           alpha=0.7, 
                           color=colors[i],
                           edgecolor='black',
                           linewidth=0.5)
            
            # Add transcript label on the left
            if len(y) > 0:
                ax.text(-2, baseline + 0.3, transcript.replace('gene_', '').replace('MT-', ''), 
                       fontsize=12, 
                       ha='right', 
                       va='center')
        
        # Set axis properties
        ax.set_xlim(-0.5, window_size + 0.5)
        ax.set_ylim(-0.5, n_transcripts * vertical_spacing + 0.5)
        ax.set_xlabel('Codon Index', fontsize=11, fontweight='bold')
        ax.set_title(f'{trans_class} Genes\nInitiation Patterns', 
                    fontsize=12, fontweight='bold')
        
        # Despine
        sns.despine(ax=ax, left=True)
        # Remove y-axis ticks (transcript names are on the plot)
        ax.set_yticks([])
        # Add more xaxis ticks if window is small
        if window_size <= 100:
            ax.set_xticks(np.arange(0, window_size + 1, 5))
        
        # Add grid for x-axis only
        ax.grid(True, alpha=0.3, axis='x', linestyle='--')
        ax.set_axisbelow(True)
    
    plt.tight_layout()
    
    if save_path:
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
        print(f"Saved ridgeplot to {save_path}")
    
    return fig


def get_transcript_summary_stats(df):
    """
    Calculate summary statistics for all transcripts.
    
    Logic:
    1. For each transcript, calculate: length, total RPM, mean per-codon RPM
    2. Provide overview statistics by transcript class
    
    Parameters:
    -----------
    df : pd.DataFrame
        Data from load_and_prepare_data()
        
    Returns:
    --------
    pd.DataFrame
        Summary statistics per transcript
        
    Notes:
    ------
    - Useful for quality control and understanding dataset composition
    - Can identify outliers or problematic transcripts
    """
    stats = []
    
    for (transcript, trans_class), group in df.groupby(['gene_id', 'transcript_class']):
        # Calculate per-transcript metrics
        length = group['codon_index'].max()
        
        # Average across samples
        sample_stats = []
        for sample in group['sample'].unique():
            sample_data = group[group['sample'] == sample]
            total_rpm = sample_data['RPM'].max()
            mean_rpm = sample_data['RPM'].mean()
            sample_stats.append({
                'total_rpm': total_rpm,
                'mean_rpm': mean_rpm
            })
        
        # Average the sample statistics
        avg_total_rpm = np.mean([s['total_rpm'] for s in sample_stats])
        avg_mean_rpm = np.mean([s['mean_rpm'] for s in sample_stats])
        
        stats.append({
            'transcript_id': transcript,
            'transcript_class': trans_class,
            'length_codons': length,
            'avg_total_rpm': avg_total_rpm,
            'avg_mean_rpm_per_codon': avg_mean_rpm
        })
    
    return pd.DataFrame(stats)

def occupancy_with_canonical_bg(df, ov_transcript_id, can_transcript_id, n_codons=50, metric = 'RPM', figsize = (12,4), savefig_path=None):
    """
    Plots the codon occupancy of an overlapping gene against the background
    occupancy of the canonical gene it overlaps.

    Parameters:
    -----------
    df : pd.DataFrame
        The dataframe containing columns: 'transcript_id', 'codon_index', 
        'codon_seq', 'RPM' (or 'codon_count_sum').
    ov_transcript_id : str
        The transcript_id of the overlapping gene (e.g., 'gene_MOTSc').
    can_transcript_id : str
        The transcript_id of the canonical gene (e.g., 'MT-ND1').
    n_codons : int
        Number of codons from the start of the overlapping gene to plot.

    Returns:
    --------
    fig : matplotlib.figure.Figure
        The generated plot figure.
    """
    
    # 1. Prepare Data: Calculate mean RPM per codon position across samples
    # We assume 'codon_count_sum' or 'RPM' is the metric. Using RPM based on your context.
    if metric not in df.columns and 'codon_count_sum' in df.columns:
        metric = 'codon_count_sum'
        
    # Aggregate data to get one value per codon index per transcript
    df_agg = df.groupby(['gene_id', 'codon_index', 'codon_seq'])[metric].mean().reset_index()
    
    # Extract specific transcript data
    ov_data = df_agg[df_agg['gene_id'] == ov_transcript_id].sort_values('codon_index')
    can_data = df_agg[df_agg['gene_id'] == can_transcript_id].sort_values('codon_index')

    if ov_data.empty:
        raise ValueError(f"Overlapping transcript {ov_transcript_id} not found in data.")
    if can_data.empty:
        raise ValueError(f"Canonical transcript {can_transcript_id} not found in data.")

    # 2. Reconstruct Sequences to find alignment
    # We create a string of the sequence to find the offset
    ov_seq = "".join(ov_data['codon_seq'].tolist())
    can_seq = "".join(can_data['codon_seq'].tolist())

    # Find the starting nucleotide position of the overlapping gene within the canonical gene
    # We assume the overlapping gene is a substring of the canonical gene (on the same strand)
    start_nucleotide_idx = can_seq.find(ov_seq)

    if start_nucleotide_idx == -1:
        print(f"Warning: Exact sequence of {ov_transcript_id} not found in {can_transcript_id}.")
        print("Attempting to align based on first 15 nucleotides...")
        start_nucleotide_idx = can_seq.find(ov_seq[:15])
        if start_nucleotide_idx == -1:
            raise ValueError("Could not align overlapping gene to canonical gene sequence.")

    # 3. Map Canonical Codons to Overlapping Frame
    # The overlapping gene starts at 'start_nucleotide_idx' (0-based) in the canonical sequence.
    # Canonical codon indices are 1-based in the data.
    # One codon = 3 nucleotides.
    
    # We want to extract the canonical background for the window requested
    aligned_can_rpms = []
    aligned_ov_rpms = []
    
    # Limit to available data length
    max_codons = min(n_codons, len(ov_data))
    
    for i in range(max_codons):
        # Get Overlapping RPM
        ov_val = ov_data.iloc[i][metric]
        aligned_ov_rpms.append(ov_val)
        
        # Calculate which canonical codon corresponds to this overlapping codon
        # Current overlapping nucleotide index relative to canonical start:
        current_nuc_index = start_nucleotide_idx + (i * 3)
        
        # Convert nucleotide index to Canonical Codon Index (1-based)
        # Example: Nuc index 0,1,2 -> Codon 1. Nuc index 3,4,5 -> Codon 2.
        can_codon_idx = (current_nuc_index // 3) + 1
        
        # Fetch the RPM for this canonical codon
        # We interpret "background" as the canonical codon that covers the start of the overlapping codon
        can_val_row = can_data[can_data['codon_index'] == can_codon_idx]
        
        if not can_val_row.empty:
            aligned_can_rpms.append(can_val_row.iloc[0][metric])
        else:
            aligned_can_rpms.append(0) # Fill with 0 if we run off the end of canonical gene

    # 4. Plotting
    sns.set_style("whitegrid")
    fig, ax = plt.subplots(figsize=figsize)
    
    x_axis = range(1, max_codons + 1)
    
    # Plot Canonical Background (Area)
    ax.fill_between(x_axis, aligned_can_rpms, color='gray', alpha=0.3, 
                    label=f'Background: {can_transcript_id} (Canonical)')
    ax.plot(x_axis, aligned_can_rpms, color='gray', alpha=0.5, linestyle='--')

    # Plot Overlapping Signal (Bars/Line)
    ax.plot(x_axis, aligned_ov_rpms, color='#d62728', linewidth=2.5, marker='o', markersize=4,
            label=f'Signal: {ov_transcript_id.replace("gene_","")} (Overlapping)')
    

    # Aesthetics
    ax.set_title(f'Ribosome Occupancy: {ov_transcript_id} vs {can_transcript_id} Background', fontsize=14)
    ax.set_xlabel(f'Codon Position ({ov_transcript_id})', fontsize=12)
    ax.set_ylabel(f'Mean {metric}', fontsize=12)

    # Move legend to the center
    ax.legend(loc='upper center',  ncol=2, frameon=True, fancybox=True, framealpha=0.3)

    # Add more xticks if n_codons is small
    if n_codons <= 100:
        ax.set_xticks(np.arange(1, max_codons + 1, 3))
    # Add text regarding alignment
    frame_shift = start_nucleotide_idx % 3
    info_text = (f"Alignment Offset: {start_nucleotide_idx} nt\n"
                 f"Reading Frame Shift: +{frame_shift}")
    plt.text(0.02, 0.95, info_text, transform=ax.transAxes, 
             bbox=dict(facecolor='white', alpha=0.8, edgecolor='none'))

    # Despine and remove grid lines for cleaner look
    sns.despine(ax=ax)
    ax.grid(False)

    plt.tight_layout()

    if savefig_path:
        plt.savefig(savefig_path, dpi=300, bbox_inches='tight')
        print(f"Saved figure to {savefig_path}")


def plot_frame_distribution(df, gene_id, gene_label=None):
    """
    Plots the reading frame distribution (Triplet Periodicity) for a specific transcript
    based on pre-computed P-site codon counts.

    Parameters:
    -----------
    df : pd.DataFrame
        Must contain columns: 'transcript_id', 'position_1_count', 
        'position_2_count', 'position_3_count', 'codon_index'.
    transcript_id : str
        The ID of the gene to analyze (e.g., 'gene_MOTSc').
    gene_label : str (optional)
        A label for the plot title (e.g., 'MOTSc'). If None, uses transcript_id.

    Returns:
    --------
    fig : matplotlib.figure.Figure
        A figure containing two subplots:
        1. Bar chart of total reads per frame.
        2. Line chart of frame percentages across the gene body.
    """
    
    # 1. Filter data for the specific transcript
    subset = df[df['gene_id'] == gene_id].copy()
    
    if subset.empty:
        print(f"Error: Transcript {gene_id} not found in dataframe.")
        return None
        
    # Ensure data is sorted by codon index for the line plot
    subset = subset.sort_values('codon_index')
    
    # 2. Calculate Global Totals (The Bar Chart)
    f0_total = subset['position_1_count'].sum() # Frame 0 (In-frame)
    f1_total = subset['position_2_count'].sum() # Frame 1 (+1 shift)
    f2_total = subset['position_3_count'].sum() # Frame 2 (+2 shift)
    
    total_counts = f0_total + f1_total + f2_total
    
    if total_counts == 0:
        print(f"No counts found for {gene_id}.")
        return None

    # Calculate percentages
    f0_pct = (f0_total / total_counts) * 100
    f1_pct = (f1_total / total_counts) * 100
    f2_pct = (f2_total / total_counts) * 100
    
    # 3. Prepare Line Plot Data (Rolling Average for smoothness)
    # We want to see if the frame dominance is consistent across the gene
    # Using a rolling window to smooth out noise (e.g., window of 10 codons)
    window = 10
    subset['F0_roll'] = subset['position_1_count'].rolling(window).sum()
    subset['F1_roll'] = subset['position_2_count'].rolling(window).sum()
    subset['F2_roll'] = subset['position_3_count'].rolling(window).sum()
    
    # Normalize rows to 100% to see relative frame usage per window
    row_sums = subset[['F0_roll', 'F1_roll', 'F2_roll']].sum(axis=1)
    subset['F0_norm'] = (subset['F0_roll'] / row_sums) * 100
    subset['F1_norm'] = (subset['F1_roll'] / row_sums) * 100
    subset['F2_norm'] = (subset['F2_roll'] / row_sums) * 100

    # 4. Plotting
    fig = plt.figure(figsize=(12, 5))
    
    # --- Subplot 1: Global Frame Distribution (Bar Chart) ---
    ax1 = plt.subplot(1, 2, 1)
    frames = ['Frame 0\n(Pos 1)', 'Frame 1\n(Pos 2)', 'Frame 2\n(Pos 3)']
    counts = [f0_total, f1_total, f2_total]
    pcts = [f0_pct, f1_pct, f2_pct]
    colors = ['#1f77b4', '#ff7f0e', '#2ca02c'] # Blue, Orange, Green
    
    bars = ax1.bar(frames, counts, color=colors, alpha=0.8, edgecolor='black')
    
    # Add percentage labels on top of bars
    max_height = max(counts)
    for bar, pct in zip(bars, pcts):
        height = bar.get_height()
        ax1.text(bar.get_x() + bar.get_width()/2., height + (max_height*0.02),
                 f'{pct:.1f}%', ha='center', va='bottom', fontsize=11, fontweight='bold')
        
    ax1.set_ylabel('Total Ribosome Footprints (Corrected)', fontsize=11)
    ax1.set_title(f'Global Frame Distribution: {gene_label or transcript_id}', fontsize=13)
    ax1.grid(axis='y', linestyle='--', alpha=0.5)

    # --- Subplot 2: Frame Consistency Across Gene (Line Plot) ---
    ax2 = plt.subplot(1, 2, 2)
    
    # Plot the smoothed percentage lines
    ax2.plot(subset['codon_index'], subset['F0_norm'], label='Frame 0', color='#1f77b4', linewidth=2)
    ax2.plot(subset['codon_index'], subset['F1_norm'], label='Frame 1', color='#ff7f0e', linewidth=1.5, alpha=0.7)
    ax2.plot(subset['codon_index'], subset['F2_norm'], label='Frame 2', color='#2ca02c', linewidth=1.5, alpha=0.7)
    
    ax2.set_xlabel('Codon Index')
    ax2.set_ylabel('% of Reads (Rolling Window)')
    ax2.set_title('Frame Periodicity Across Gene Body', fontsize=13)
    ax2.legend(loc='upper right')
    ax2.set_ylim(0, 100)
    ax2.grid(True, linestyle='--', alpha=0.3)
    
    plt.tight_layout()
    return fig

