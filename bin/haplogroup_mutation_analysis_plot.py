import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from collections import defaultdict
import os
import re
from typing import Dict, List, Tuple, Optional
import warnings
warnings.filterwarnings('ignore')

from Bio.Seq import Seq
from Bio.Data import CodonTable
from adjustText import adjust_text

def get_amino_acid_change(position: int, ref_allele: str, alt_allele: str,
                          reference_seq: str, orf_start: int = 1,
                          codon_table_id: int = 2) -> Tuple[str, str, int, str]:
    """
    Calculate amino acid change for a nucleotide substitution using Biopython.

    Parameters:
    -----------
    position : int
        1-based position of mutation in the sequence
    ref_allele : str
        Reference nucleotide
    alt_allele : str
        Alternative nucleotide
    reference_seq : str
        Full reference sequence
    orf_start : int
        1-based position where the ORF starts (default 1)
    codon_table_id : int
        NCBI codon table ID (default 2 = Vertebrate Mitochondrial)
        1 = Standard, 2 = Vertebrate Mitochondrial, etc.

    Returns:
    --------
    Tuple[str, str, int, str]
        (ref_aa, alt_aa, codon_position, change_type)
        change_type is 'synonymous', 'nonsynonymous', 'nonsense', or 'unknown'
    """
    # Get codon table
    try:
        table = CodonTable.unambiguous_dna_by_id[codon_table_id]
    except KeyError:
        table = CodonTable.unambiguous_dna_by_id[2]  # Default to vertebrate mitochondrial

    # Convert to 0-based indexing
    pos_0based = position - 1
    orf_start_0based = orf_start - 1

    # Check if position is within the sequence
    if pos_0based < 0 or pos_0based >= len(reference_seq):
        return ('?', '?', 0, 'unknown')

    # Calculate position relative to ORF start
    rel_pos = pos_0based - orf_start_0based

    if rel_pos < 0:
        return ('?', '?', 0, 'upstream')

    # Determine codon number (0-based) and position within codon (0, 1, or 2)
    codon_num = rel_pos // 3
    codon_pos = rel_pos % 3  # 0 = first, 1 = second, 2 = third position

    # Get codon start position in sequence
    codon_start = orf_start_0based + (codon_num * 3)

    # Check if we can extract the full codon
    if codon_start < 0 or codon_start + 3 > len(reference_seq):
        return ('?', '?', codon_pos + 1, 'incomplete')

    # Extract reference codon
    ref_codon = reference_seq[codon_start:codon_start + 3].upper()

    # Build mutant codon
    codon_list = list(ref_codon)
    codon_list[codon_pos] = alt_allele.upper()
    alt_codon = ''.join(codon_list)

    # Translate codons using Biopython
    try:
        ref_aa = str(Seq(ref_codon).translate(table=codon_table_id))
        alt_aa = str(Seq(alt_codon).translate(table=codon_table_id))
    except Exception:
        return ('?', '?', codon_pos + 1, 'unknown')

    # Determine change type
    if ref_aa == alt_aa:
        change_type = 'synonymous'
    elif alt_aa == '*':
        change_type = 'nonsense'
    elif ref_aa == '*':
        change_type = 'readthrough'
    else:
        change_type = 'nonsynonymous'

    # Return amino acid number (1-based)
    aa_num = codon_num + 1
    return (ref_aa, alt_aa, aa_num, change_type)

def format_aa_change(ref_aa: str, alt_aa: str, aa_num: int, change_type: str) -> str:
    """Format amino acid change as standard notation (e.g., M15V or M15= for synonymous)."""
    if ref_aa == '?' or alt_aa == '?':
        return ''
    if change_type == 'synonymous':
        return f"{ref_aa}{aa_num}="
    return f"{ref_aa}{aa_num}{alt_aa}"

def load_haplogroup_associations(excel_file: str) -> pd.DataFrame:
    """
    Load haplogroup mutation association data from Excel file.
    
    Parameters:
    -----------
    excel_file : str
        Path to Excel file containing haplogroup mutation data
        
    Returns:
    --------
    pd.DataFrame
        DataFrame with mutation association data
    """
    try:
        # Read the Associations sheet
        df = pd.read_excel(excel_file, sheet_name='Associations')
        print(f"Loaded associations data with shape: {df.shape}")
        print(f"Columns: {list(df.columns)}")
        return df
    except Exception as e:
        print(f"Error loading Excel file: {e}")
        return pd.DataFrame()

def parse_mutation_string(mutation: str) -> Tuple[Optional[int], Optional[str], Optional[str]]:
    """
    Parse mutation string to extract position and allele changes.
    
    Expected formats: A123G, 123A>G, A123T, etc.
    
    Parameters:
    -----------
    mutation : str
        Mutation string
        
    Returns:
    --------
    Tuple[Optional[int], Optional[str], Optional[str]]
        (position, reference_allele, alternative_allele)
    """
    if pd.isna(mutation):
        return None, None, None
    
    mutation = str(mutation).strip()
    
    # Pattern 1: A123G (reference + position + alternative)
    pattern1 = re.match(r'^([ATGC])(\d+)([ATGC])$', mutation, re.IGNORECASE)
    if pattern1:
        ref_allele = pattern1.group(1).upper()
        position = int(pattern1.group(2))
        alt_allele = pattern1.group(3).upper()
        return position, ref_allele, alt_allele
    
    # Pattern 2: 123A>G (position + ref > alt)
    pattern2 = re.match(r'^(\d+)([ATGC])>([ATGC])$', mutation, re.IGNORECASE)
    if pattern2:
        position = int(pattern2.group(1))
        ref_allele = pattern2.group(2).upper()
        alt_allele = pattern2.group(3).upper()
        return position, ref_allele, alt_allele
    
    # Pattern 3: Try to extract just position if it's a simple number
    pattern3 = re.match(r'^(\d+)$', mutation)
    if pattern3:
        position = int(pattern3.group(1))
        return position, None, None
    
    # Pattern 4: More complex patterns - try to extract position
    position_match = re.search(r'(\d+)', mutation)
    if position_match:
        position = int(position_match.group(1))
        return position, None, None
    
    print(f"Could not parse mutation: {mutation}")
    return None, None, None

def reshape_haplogroup_data(df: pd.DataFrame) -> pd.DataFrame:
    """
    Reshape data from wide format (separate frequency columns) to long format.
    
    Parameters:
    -----------
    df : pd.DataFrame
        Input dataframe with Freq_X columns
        
    Returns:
    --------
    pd.DataFrame
        Reshaped dataframe with one row per gene-mutation-haplogroup combination
    """
    # Get all frequency columns
    freq_columns = [col for col in df.columns if col.startswith('Freq_')]
    
    if not freq_columns:
        print("No frequency columns found!")
        return pd.DataFrame()
    
    # Extract haplogroup names from column names
    haplogroup_mapping = {}
    for col in freq_columns:
        haplogroup = col.replace('Freq_', '')
        haplogroup_mapping[col] = haplogroup
    
    print(f"Found {len(freq_columns)} haplogroups: {list(haplogroup_mapping.values())}")
    
    # Parse mutation information
    mutation_info = df['Mutation'].apply(parse_mutation_string)
    df['Position'] = [info[0] for info in mutation_info]
    df['Ref_Allele'] = [info[1] for info in mutation_info]
    df['Alt_Allele'] = [info[2] for info in mutation_info]
    
    # Reshape to long format
    id_vars = ['Gene', 'Mutation', 'Position', 'Ref_Allele', 'Alt_Allele', 
               'N_Carriers', 'N_Total', 'P_Value', 'P_Corrected', 
               'Significant_Uncorrected', 'Significant_Corrected']
    
    # Only include id_vars that actually exist in the dataframe
    id_vars = [col for col in id_vars if col in df.columns]
    
    df_long = pd.melt(df, id_vars=id_vars, 
                      value_vars=freq_columns,
                      var_name='Frequency_Column', 
                      value_name='Frequency')
    
    # Map frequency column names to haplogroup names
    df_long['Haplogroup'] = df_long['Frequency_Column'].map(haplogroup_mapping)
    
    # Remove rows with missing or zero frequencies
    df_long = df_long.dropna(subset=['Frequency'])
    df_long = df_long[df_long['Frequency'] > 0]
    
    # Remove rows with missing positions
    df_long = df_long.dropna(subset=['Position'])
    
    print(f"Reshaped data: {len(df_long)} gene-mutation-haplogroup combinations")
    
    return df_long

def get_gene_sequences(hs_pop_df: pd.DataFrame, gene_list: List[str]) -> Dict[str, str]:
    """
    Extract reference sequences for each gene from the population dataframe.
    
    Parameters:
    -----------
    hs_pop_df : pd.DataFrame
        Population dataframe with gene sequences
    gene_list : List[str]
        List of gene names
        
    Returns:
    --------
    Dict[str, str]
        Dictionary mapping gene names to their reference sequences
    """
    gene_sequences = {}
    
    for gene in gene_list:
        seq_col = f"{gene}"
        if seq_col in hs_pop_df.columns:
            # Use the most common sequence as reference
            sequences = hs_pop_df[seq_col].dropna()
            if not sequences.empty:
                reference_seq = sequences.mode().iloc[0] if not sequences.mode().empty else sequences.iloc[0]
                gene_sequences[gene] = reference_seq
                print(f"Found reference sequence for {gene}: {len(reference_seq)} bp")
            else:
                print(f"Warning: No sequences found for {gene}")
        else:
            print(f"Warning: Column {seq_col} not found in dataframe")
    
    return gene_sequences

def create_substitution_type_map():
    """Create mapping for different types of substitutions."""
    return {
        # Transitions
        ('A', 'G'): 'Transition', ('G', 'A'): 'Transition',
        ('C', 'T'): 'Transition', ('T', 'C'): 'Transition',
        # Transversions
        ('A', 'C'): 'Transversion', ('C', 'A'): 'Transversion',
        ('A', 'T'): 'Transversion', ('T', 'A'): 'Transversion',
        ('G', 'C'): 'Transversion', ('C', 'G'): 'Transversion',
        ('G', 'T'): 'Transversion', ('T', 'G'): 'Transversion',
    }

def plot_gene_mutations(gene_data: pd.DataFrame, gene_name: str, reference_seq: str,
                       output_dir: str = "plots", figsize: Tuple[int, int] = (16, 10),
                       gene_start_position: int = None,
                       show_absolute_positions: bool = True,
                       positions_are_absolute: bool = False,
                       orf_start: int = 1,
                       codon_table_id: int = 2,
                       single_panel: bool = False,
                       show_aa_change: bool = True,
                       annotate_mutations: List[str] = None) -> None:
    """
    Create a comprehensive mutation plot for a single gene.

    Parameters:
    -----------
    gene_data : pd.DataFrame
        Mutation data for the specific gene
    gene_name : str
        Name of the gene
    reference_seq : str
        Reference sequence for the gene
    output_dir : str
        Directory to save plots
    figsize : Tuple[int, int]
        Figure size (width, height)
    gene_start_position : int, optional
        Start position of the gene in mtDNA (1-based, Cambridge reference).
        If provided, absolute positions will be displayed alongside relative positions.
        Common values: RNR1=1671, RNR2=3230, ND4L-ND4=10470, CYTB=14747
    show_absolute_positions : bool
        Whether to show absolute mtDNA positions (requires gene_start_position)
    positions_are_absolute : bool
        If True, the 'Position' column in gene_data contains absolute mtDNA positions
        and will be converted to gene-relative positions for plotting.
        Requires gene_start_position to be set.
    orf_start : int
        1-based position where the ORF/reading frame starts within the gene sequence.
        Default is 1 (sequence starts at first codon). Used for amino acid translation.
    codon_table_id : int
        NCBI codon table ID for translation (default 2 = Vertebrate Mitochondrial).
        1 = Standard, 2 = Vertebrate Mitochondrial, 3 = Yeast Mitochondrial, etc.
    single_panel : bool
        If True, only generate the frequency scatter plot (ax2) as a single panel.
        Default is False (generates full 4-panel figure).
    show_aa_change : bool
        If True, include amino acid change (e.g., '*25W') in annotation boxes.
        Default is True. Set to False for non-coding regions or to simplify annotations.
    annotate_mutations : List[str], optional
        List of specific mutations to annotate regardless of frequency threshold.
        Matches against both relative (e.g., 'C119T') and absolute (e.g., 'C11674T') mutation names.
        These will be added to the high-frequency mutations for annotation.
    """
    if gene_data.empty:
        print(f"No mutation data found for {gene_name}")
        return

    # Create output directory if it doesn't exist
    os.makedirs(output_dir, exist_ok=True)

    # Handle absolute position input - convert to relative positions
    if positions_are_absolute:
        if gene_start_position is None:
            raise ValueError("gene_start_position must be provided when positions_are_absolute=True")
        # Make a copy to avoid modifying original data
        gene_data = gene_data.copy()
        # Store original absolute position and convert to relative
        gene_data['Absolute_Position_Original'] = gene_data['Position']
        gene_data['Position'] = gene_data['Position'] - gene_start_position + 1
        print(f"Converted absolute positions to gene-relative (subtracted {gene_start_position - 1})")
        print(f"  Position range: {gene_data['Position'].min()} - {gene_data['Position'].max()}")

    # Determine if we can show absolute positions
    use_absolute = show_absolute_positions and gene_start_position is not None

    # Set up the plot based on single_panel mode
    if single_panel:
        # Single panel: only frequency scatter plot
        fig, ax2 = plt.subplots(1, 1, figsize=figsize)
        ax1 = ax3 = ax4 = None  # Not used in single panel mode
    else:
        # Full 4-panel figure
        fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=figsize, height_ratios=[1, 2])

    
    # Get substitution type mapping
    sub_type_map = create_substitution_type_map()
    
    # Add substitution type to data
    if 'Ref_Allele' in gene_data.columns and 'Alt_Allele' in gene_data.columns:
        gene_data = gene_data.copy()
        gene_data['Substitution_Type'] = gene_data.apply(
            lambda row: sub_type_map.get((row['Ref_Allele'], row['Alt_Allele']), 'Unknown') 
            if pd.notna(row['Ref_Allele']) and pd.notna(row['Alt_Allele']) else 'Unknown', 
            axis=1
        )
    else:
        gene_data['Substitution_Type'] = 'Unknown'
    
    # Find the top haplogroup for each position and get total population frequency
    position_data = []
    for pos in gene_data['Position'].unique():
        pos_data = gene_data[gene_data['Position'] == pos]
        if not pos_data.empty:
            # Get the haplogroup with highest frequency at this position
            top_row = pos_data.loc[pos_data['Frequency'].idxmax()]
            # Get total population frequency from N_Carriers/N_Total
            total_frequency = top_row['N_Carriers'] / top_row['N_Total'] if top_row['N_Total'] > 0 else 0

            # Calculate absolute position if gene_start_position is provided
            abs_pos = gene_start_position + pos - 1 if use_absolute else None

            # Build mutation string with absolute position
            ref_allele = top_row.get('Ref_Allele', 'N')
            alt_allele = top_row.get('Alt_Allele', 'N')
            mutation_absolute = f"{ref_allele}{abs_pos}{alt_allele}" if use_absolute else None

            # Calculate amino acid change
            ref_aa, alt_aa, aa_num, change_type = get_amino_acid_change(
                position=pos, ref_allele=ref_allele, alt_allele=alt_allele,
                reference_seq=reference_seq, orf_start=orf_start,
                codon_table_id=codon_table_id
            )
            aa_change = format_aa_change(ref_aa, alt_aa, aa_num, change_type)

            position_data.append({
                'position': pos,
                'absolute_position': abs_pos,
                'frequency': top_row['Frequency'],
                'total_frequency': total_frequency,
                'n_carriers': top_row['N_Carriers'],
                'n_total': top_row['N_Total'],
                'haplogroup': top_row['Haplogroup'],
                'ref_allele': ref_allele,
                'alt_allele': alt_allele,
                'substitution_type': top_row['Substitution_Type'],
                'mutation': top_row['Mutation'],
                'aa_change': aa_change,
                'aa_change_type': change_type,
                'mutation_absolute': mutation_absolute,
                'p_value': top_row.get('P_Value', 1.0),
                'significant': top_row.get('Significant_Corrected', False)
            })
    
    position_df = pd.DataFrame(position_data)
    
    if position_df.empty:
        print(f"No position data available for {gene_name}")
        return
    
    # Get sequence length (needed for multiple plots)
    seq_length = len(reference_seq)

    # Plot 1: Reference sequence visualization (first 300 bp) - skip in single panel mode
    if not single_panel:
        display_length = min(300, seq_length)

        # Color code nucleotides
        nucleotide_colors = {'A': 'red', 'T': 'blue', 'G': 'green', 'C': 'orange'}

        for i, nucleotide in enumerate(reference_seq[:display_length]):
            color = nucleotide_colors.get(nucleotide.upper(), 'gray')
            ax1.barh(0, 1, left=i, color=color, alpha=0.7)

        # Mark mutation positions on sequence
        mutation_positions = position_df['position'].values
        for pos in mutation_positions:
            if pos <= display_length:
                ax1.axvline(x=pos-1, color='black', linestyle='-', alpha=0.8, linewidth=2)

        ax1.set_xlim(0, display_length)
        ax1.set_ylim(-0.5, 0.5)

        # Build xlabel with absolute positions if available
        if use_absolute:
            abs_start = gene_start_position
            abs_end_display = gene_start_position + display_length - 1
            abs_end_total = gene_start_position + seq_length - 1
            ax1.set_xlabel(f'Gene pos 1-{display_length} of {seq_length} bp\n'
                          f'(mtDNA pos {abs_start}-{abs_end_display} of {abs_start}-{abs_end_total})')
        else:
            ax1.set_xlabel(f'Position (showing 1-{display_length} of {seq_length} bp)')
        ax1.set_yticks([])

        # Add nucleotide legend
        legend_elements = [plt.Rectangle((0,0),1,1, facecolor=color, alpha=0.7, label=nuc)
                          for nuc, color in nucleotide_colors.items()]
        ax1.legend(handles=legend_elements, loc='upper right', bbox_to_anchor=(1, 1))
    
    # Plot 2: Mutation frequencies with top haplogroup per position

    
    # Create scatter plot with different symbols for substitution types
    substitution_markers = {'Transition': 'o', 'Transversion': 's', 'Unknown': '^'}
    substitution_colors = {'Transition': 'blue', 'Transversion': 'red', 'Unknown': 'gray'}
    
    # Plot top haplogroup frequencies
    for sub_type in position_df['substitution_type'].unique():
        mask = position_df['substitution_type'] == sub_type
        subset = position_df[mask]
        
        # Use different sizes for significant vs non-significant
        sizes = [150 if sig else 100 for sig in subset['significant']]
        
        scatter = ax2.scatter(subset['position'], subset['frequency'], 
                            marker=substitution_markers.get(sub_type, 'o'),
                            c=substitution_colors.get(sub_type, 'gray'),
                            s=sizes, alpha=0.7, label=f'{sub_type} (Top Haplogroup)', 
                            edgecolors='black')
    
    # Add total frequency line plot
    ax2.scatter(position_df['position'], position_df['total_frequency'], 
            alpha=0.8, label='Total Frequency (All Haplogroups)', marker= 'x', c ='black', s=50, edgecolors='white')
    
    # Add small circles for total frequencies
    ax2.scatter(position_df['position'], position_df['total_frequency'], 
               c='black', s=50, alpha=0.6, marker='o', edgecolors='white')
    
    # Add frequency threshold line
    ax2.axhline(y=0.01, color='gray', linestyle='--', alpha=0.5, label='1% threshold')

    # Set xlabel based on position mode
    if use_absolute:
        ax2.set_xlabel(f'Position in Gene (mtDNA: {gene_start_position}-{gene_start_position + seq_length - 1})')
    else:
        ax2.set_xlabel('Position in Gene')
    ax2.set_ylabel('Frequency')

    # Add padding to x-axis to prevent markers from being clipped at edges
    x_padding = seq_length * 0.03  # 3% padding on each side
    ax2.set_xlim(-x_padding, seq_length + x_padding)
    ax2.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
    ax2.grid(True, alpha=0.3)

    # Add secondary x-axis for absolute positions if available
    if use_absolute:
        ax2_abs = ax2.twiny()
        ax2_abs.set_xlim(gene_start_position - x_padding, gene_start_position + seq_length + x_padding)
        ax2_abs.set_xlabel('mtDNA Position (Cambridge Reference)')

    # Annotate high-frequency mutations using adjustText for automatic positioning
    high_freq_threshold = max(0.05, position_df['frequency'].quantile(0.8))  # At least 5% or top 20%
    high_freq_mutations = position_df[position_df['frequency'] >= high_freq_threshold]

    # Add explicitly requested mutations regardless of frequency
    if annotate_mutations:
        # Match against both relative mutation name and absolute mutation name
        explicit_mask = (
            position_df['mutation'].isin(annotate_mutations) |
            position_df['mutation_absolute'].isin(annotate_mutations)
        )
        explicit_mutations = position_df[explicit_mask]
        # Combine with high frequency mutations, avoiding duplicates
        high_freq_mutations = pd.concat([high_freq_mutations, explicit_mutations]).drop_duplicates(subset=['position'])

    # Get max frequency to determine y-axis bounds
    max_freq = position_df['frequency'].max() if not position_df.empty else 1.0

    # Collect text objects and their target points for adjust_text
    texts = []
    x_points = []
    y_points = []

    for _, row in high_freq_mutations.iterrows():
        sig_marker = "*" if row['significant'] else ""
        # Build annotation text with nucleotide mutation and optionally amino acid change
        aa_line = f" ({row['aa_change']})" if show_aa_change and row.get('aa_change') else ""

        # Only show mtDNA absolute position if positions were NOT already absolute in input
        if use_absolute and row['mutation_absolute'] and not positions_are_absolute:
            annotation_text = (f"{row['mutation']}{sig_marker}{aa_line}\n"
                              f"mtDNA: {row['mutation_absolute']}\n"
                              f"{row['haplogroup']}: {row['frequency']:.1%}\n"
                              f"Total: {row['total_frequency']:.1%}")
        else:
            annotation_text = (f"{row['mutation']}{sig_marker}{aa_line}\n"
                              f"{row['haplogroup']}: {row['frequency']:.1%}\n"
                              f"Total: {row['total_frequency']:.1%}")

        # Store point coordinates
        x_points.append(row['position'])
        y_points.append(row['frequency'])

        # Determine box color based on significance (with transparency)
        box_color = 'yellow' if row['significant'] else 'lightblue'

        # Calculate initial offset for text (away from BOTH markers)
        # Consider both the top haplogroup frequency and total frequency positions
        freq_gap = abs(row['frequency'] - row['total_frequency'])

        # Determine if there's enough gap between the two markers
        if freq_gap > 0.15:
            # Markers are far apart - place text outside the gap (above top or below total)
            if row['frequency'] > row['total_frequency']:
                # Top marker is higher - place text above top marker
                initial_y_offset = 0.08 * max_freq
                va = 'bottom'
            else:
                # Total is higher (rare) - place text below top marker
                initial_y_offset = -0.25 * max_freq
                va = 'top'
        else:
            # Markers are close together - place text clearly outside both
            # Place above if there's room, otherwise below
            if row['frequency'] < max_freq * 0.7:
                # Room above - place text above both markers
                initial_y_offset = 0.12 * max_freq
                va = 'bottom'
            else:
                # Near top - place text well below both markers
                initial_y_offset = -0.30 * max_freq
                va = 'top'

        # Slight x offset to avoid sitting directly on the marker
        initial_x_offset = seq_length * 0.02

        # Create text annotation with initial offset from the data point
        txt = ax2.text(row['position'] + initial_x_offset,
                      row['frequency'] + initial_y_offset,
                      annotation_text,
                      fontsize=8, ha='left', va=va,
                      bbox=dict(boxstyle='round,pad=0.3',
                               facecolor=box_color,
                               alpha=0.4,  # Partial transparency
                               edgecolor='gray',
                               linewidth=0.5))
        texts.append(txt)

    # Use adjust_text to automatically position annotations avoiding overlaps
    if texts:
        # Get ALL scatter points to avoid (both top haplogroup and total frequency points)
        all_x = list(position_df['position'].values) * 2  # Both frequency types at same x
        all_y = list(position_df['frequency'].values) + list(position_df['total_frequency'].values)

        adjust_text(texts,
                   x=all_x,  # Pass ALL points to avoid, not just annotated ones
                   y=all_y,
                   ax=ax2,
                   expand_points=(2.0, 2.0),  # Increased expansion around points
                   expand_text=(1.5, 1.5),    # Increased expansion around text
                   force_points=(1.0, 1.0),   # Increased repulsion from points
                   force_text=(0.8, 0.8),     # Increased repulsion between texts
                   lim=1000,  # More iterations for better placement
                   only_move={'points': 'y', 'text': 'xy'},
                   avoid_self=True)
    
    # Calculate haplogroup counts (needed for summary stats even in single panel mode)
    haplogroup_counts = position_df['haplogroup'].value_counts()

    # Plot 3: Haplogroup distribution - skip in single panel mode
    if not single_panel:
        colors = plt.cm.Set3(np.linspace(0, 1, len(haplogroup_counts)))

        bars = ax3.bar(range(len(haplogroup_counts)), haplogroup_counts.values, color=colors)
        ax3.set_xticks(range(len(haplogroup_counts)))
        ax3.set_xticklabels(haplogroup_counts.index, rotation=45)
        ax3.set_ylabel('Number of Positions')
        ax3.grid(True, alpha=0.3)

        # Add value labels on bars
        for bar, count in zip(bars, haplogroup_counts.values):
            ax3.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.1,
                    str(count), ha='center', va='bottom', fontweight='bold')

    # Plot 4: All haplogroup frequencies for all positions (heatmap style) - skip in single panel mode
    if not single_panel:
        # Create a pivot table for heatmap
        pivot_data = gene_data.pivot_table(index='Position', columns='Haplogroup',
                                          values='Frequency', fill_value=0)

        if not pivot_data.empty:
            # Only show positions and haplogroups with some signal
            pivot_data = pivot_data.loc[(pivot_data > 0.01).any(axis=1)]  # Positions with >1% in any haplogroup
            pivot_data = pivot_data.loc[:, (pivot_data > 0.01).any(axis=0)]  # Haplogroups with >1% somewhere

            if not pivot_data.empty and pivot_data.shape[0] > 1 and pivot_data.shape[1] > 1:
                im = ax4.imshow(pivot_data.T, aspect='auto', cmap='YlOrRd', interpolation='nearest')

                # Set ticks and labels - include absolute positions if available
                ax4.set_xticks(range(len(pivot_data.index)))
                if use_absolute:
                    # Show both relative and absolute positions
                    abs_labels = [f"{pos}\n({gene_start_position + pos - 1})"
                                 for pos in pivot_data.index]
                    ax4.set_xticklabels(abs_labels, rotation=45, fontsize=7)
                else:
                    ax4.set_xticklabels(pivot_data.index, rotation=45)
                ax4.set_yticks(range(len(pivot_data.columns)))
                ax4.set_yticklabels(pivot_data.columns)

                if use_absolute:
                    ax4.set_xlabel('Gene Position (mtDNA Position)')
                else:
                    ax4.set_xlabel('Position')
                ax4.set_ylabel('Haplogroup')

                # Add colorbar
                cbar = plt.colorbar(im, ax=ax4, shrink=0.8)
                cbar.set_label('Frequency')
            else:
                ax4.text(0.5, 0.5, 'Insufficient data for heatmap',
                        ha='center', va='center', transform=ax4.transAxes)
        else:
            ax4.text(0.5, 0.5, 'No data for heatmap',
                    ha='center', va='center', transform=ax4.transAxes)

    # Add figure title with gene location info
    if use_absolute:
        abs_end = gene_start_position + seq_length - 1
        title_text = (f'{gene_name} Haplogroup Mutations\n'
                     f'(mtDNA positions {gene_start_position}-{abs_end}, {seq_length} bp)')
    else:
        title_text = f'{gene_name} Haplogroup Mutations'

    if single_panel:
        ax2.set_title(title_text, fontsize=14, fontweight='bold')
    else:
        fig.suptitle(title_text, fontsize=14, fontweight='bold', y=1.02)

    sns.despine(fig=fig, trim=False)
    plt.tight_layout()

    # Save plot with different filename for single panel
    if single_panel:
        output_file = os.path.join(output_dir, f"{gene_name}_haplogroup_mutations_frequency.png")
    else:
        output_file = os.path.join(output_dir, f"{gene_name}_haplogroup_mutations.png")
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    plt.show()

    print(f"Plot saved: {output_file}")

    # Print summary statistics
    print(f"\nSummary for {gene_name}:")
    print(f"  Gene length: {seq_length} bp")
    if use_absolute:
        print(f"  mtDNA position: {gene_start_position}-{gene_start_position + seq_length - 1}")
    print(f"  Number of mutation positions: {len(position_df)}")
    print(f"  Total mutations across haplogroups: {len(gene_data)}")
    print(f"  Top haplogroup frequency range: {position_df['frequency'].min():.1%} - {position_df['frequency'].max():.1%}")
    print(f"  Total frequency range: {position_df['total_frequency'].min():.1%} - {position_df['total_frequency'].max():.1%}")
    print(f"  Most common haplogroup: {haplogroup_counts.index[0]} ({haplogroup_counts.iloc[0]} positions)")
    print(f"  Substitution types: {dict(position_df['substitution_type'].value_counts())}")
    sig_count = position_df['significant'].sum() if 'significant' in position_df.columns else 0
    print(f"  Significant mutations (corrected): {sig_count}/{len(position_df)}")

    # Additional statistics about total vs top haplogroup frequencies
    mean_ratio = (position_df['frequency'] / position_df['total_frequency']).mean()
    print(f"  Mean ratio (top haplogroup / total frequency): {mean_ratio:.2f}")
    print(f"  Positions with >10% total frequency: {(position_df['total_frequency'] > 0.1).sum()}")

def analyze_haplogroup_mutations(excel_file: str, hs_pop_df: pd.DataFrame,
                                ignore_trna: bool = True,
                                output_dir: str = "haplogroup_plots",
                                min_frequency: float = 0.01,
                                figsize: Tuple[int, int] = (16, 10),
                                gene_start_positions: Dict[str, int] = None,
                                show_absolute_positions: bool = True,
                                positions_are_absolute: bool = False,
                                orf_starts: Dict[str, int] = None,
                                codon_table_id: int = 2,
                                single_panel: bool = False,
                                show_aa_change: bool = True,
                                annotate_mutations: Dict[str, List[str]] = None) -> Dict[str, pd.DataFrame]:
    """
    Main function to analyze haplogroup-specific mutations across genes.

    Parameters:
    -----------
    excel_file : str
        Path to Excel file with haplogroup association data (Associations sheet)
    hs_pop_df : pd.DataFrame
        Population dataframe with gene sequences (columns: {gene_name} or {gene_name}_seq)
    ignore_trna : bool
        Whether to ignore tRNA genes (genes starting with 'TRN')
    output_dir : str
        Directory to save output plots
    min_frequency : float
        Minimum frequency threshold for plotting
    figsize : Tuple[int, int]
        Figure size for individual gene plots
    gene_start_positions : Dict[str, int], optional
        Dictionary mapping gene names to their mtDNA start positions (1-based).
        Example: {'RNR1': 1671, 'RNR2': 3230, 'ND4L-ND4': 10470, 'CYTB': 14747}
        If provided, absolute mtDNA positions will be shown in plots.
    show_absolute_positions : bool
        Whether to show absolute mtDNA positions (requires gene_start_positions)
    positions_are_absolute : bool
        If True, the 'Position' column in the data contains absolute mtDNA positions
        and will be converted to gene-relative positions for plotting.
        Requires gene_start_positions to be set.
    orf_starts : Dict[str, int], optional
        Dictionary mapping gene names to their ORF start positions (1-based).
        Example: {'RNR1': 1, 'CYTB': 1}. Default is 1 for all genes.
        Used for amino acid translation in annotations.
    codon_table_id : int
        NCBI codon table ID for translation (default 2 = Vertebrate Mitochondrial).
        1 = Standard, 2 = Vertebrate Mitochondrial, etc.
    single_panel : bool
        If True, only generate the frequency scatter plot as a single panel.
        Default is False (generates full 4-panel figure).
    show_aa_change : bool
        If True, include amino acid change in annotation boxes.
        Default is True. Set to False for non-coding regions or to simplify annotations.
    annotate_mutations : Dict[str, List[str]], optional
        Dictionary mapping gene names to lists of mutations to annotate regardless of frequency.
        Example: {'ND4_alt_seq': ['C11674T'], 'CYTB': ['G15043A', 'A15326G']}
        Matches against both relative and absolute mutation names.

    Returns:
    --------
    Dict[str, pd.DataFrame]
        Dictionary mapping gene names to their mutation data
    """
    print("Starting haplogroup mutation analysis...")
    
    # Load and prepare data
    raw_data = load_haplogroup_associations(excel_file)
    if raw_data.empty:
        print("Failed to load data. Exiting.")
        return {}
    
    # Reshape data from wide to long format
    mutation_data = reshape_haplogroup_data(raw_data)
    if mutation_data.empty:
        print("Failed to reshape mutation data. Exiting.")
        return {}
    
    # Filter by frequency threshold
    mutation_data = mutation_data[mutation_data['Frequency'] >= min_frequency]
    print(f"Filtered to {len(mutation_data)} gene-haplogroup combinations with frequency >= {min_frequency:.1%}")
    
    # Get unique genes
    all_genes = mutation_data['Gene'].unique()
    print(f"Found {len(all_genes)} genes with mutations")
    
    # Filter out tRNA genes if requested
    if ignore_trna:
        genes_to_analyze = [gene for gene in all_genes if not str(gene).upper().startswith('TRN')]
        excluded_trna = [gene for gene in all_genes if str(gene).upper().startswith('TRN')]
        if excluded_trna:
            print(f"Ignoring {len(excluded_trna)} tRNA genes: {excluded_trna}")
    else:
        genes_to_analyze = all_genes
    
    print(f"Analyzing {len(genes_to_analyze)} genes: {list(genes_to_analyze)}")
    
    # Get gene sequences
    gene_sequences = get_gene_sequences(hs_pop_df, genes_to_analyze)
    
    # Create output directory
    os.makedirs(output_dir, exist_ok=True)
    
    # Analyze each gene
    gene_results = {}
    
    for gene in genes_to_analyze:
        print(f"\nAnalyzing {gene}...")
        
        # Get mutation data for this gene
        gene_mutation_data = mutation_data[mutation_data['Gene'] == gene].copy()
        
        if gene_mutation_data.empty:
            print(f"No mutation data for {gene}")
            continue
        
        # Get reference sequence
        reference_seq = gene_sequences.get(gene, "")
        if not reference_seq:
            print(f"No reference sequence found for {gene}, skipping...")
            continue
        
        # Get gene start position if available
        gene_start = None
        if gene_start_positions:
            gene_start = gene_start_positions.get(gene)
            if gene_start:
                print(f"  Using mtDNA start position: {gene_start}")

        # Get ORF start position for this gene (default to 1)
        orf_start = 1
        if orf_starts:
            orf_start = orf_starts.get(gene, 1)

        # Get gene-specific mutations to annotate
        gene_annotate = None
        if annotate_mutations and gene in annotate_mutations:
            gene_annotate = annotate_mutations[gene]

        # Create plot for this gene
        try:
            plot_gene_mutations(gene_mutation_data, gene, reference_seq, output_dir, figsize,
                               gene_start_position=gene_start,
                               show_absolute_positions=show_absolute_positions,
                               positions_are_absolute=positions_are_absolute,
                               orf_start=orf_start,
                               codon_table_id=codon_table_id,
                               single_panel=single_panel,
                               show_aa_change=show_aa_change,
                               annotate_mutations=gene_annotate)
            gene_results[gene] = gene_mutation_data
        except Exception as e:
            print(f"Error plotting {gene}: {e}")
            continue
    
    # Create summary statistics
    print(f"\n{'='*60}")
    print("ANALYSIS SUMMARY")
    print(f"{'='*60}")
    print(f"Total genes analyzed: {len(gene_results)}")
    print(f"Total gene-haplogroup combinations: {sum(len(df) for df in gene_results.values())}")
    print(f"Output directory: {output_dir}")
    
    # Summary by gene
    print(f"\nMutations per gene:")
    for gene, data in gene_results.items():
        n_positions = len(data['Position'].unique())
        n_haplogroups = len(data['Haplogroup'].unique())
        max_freq = data['Frequency'].max()
        sig_count = data.get('Significant_Corrected', pd.Series([False]*len(data))).sum()
        print(f"  {gene}: {len(data)} combinations, {n_positions} positions, "
              f"{n_haplogroups} haplogroups, max freq: {max_freq:.1%}, significant: {sig_count}")
    
    return gene_results

def create_summary_plot(gene_results: Dict[str, pd.DataFrame], output_dir: str = "haplogroup_plots"):
    """
    Create a summary plot showing mutation patterns across all genes.
    
    Parameters:
    -----------
    gene_results : Dict[str, pd.DataFrame]
        Results from analyze_haplogroup_mutations
    output_dir : str
        Directory to save the summary plot
    """
    if not gene_results:
        print("No results to summarize")
        return
    
    # Compile summary statistics
    summary_data = []
    
    for gene, data in gene_results.items():
        # Calculate statistics per position (taking max frequency across haplogroups)
        position_stats = data.groupby('Position').agg({
            'Frequency': ['max', 'sum'],  # max for top haplogroup, sum for total
            'Haplogroup': 'count',
            'Significant_Corrected': 'any' if 'Significant_Corrected' in data.columns else lambda x: False
        }).reset_index()
        
        # Flatten column names
        position_stats.columns = ['Position', 'Max_Frequency', 'Total_Frequency', 'Haplogroup_Count', 'Any_Significant']
        
        # Handle substitution type counts safely
        n_transitions = 0
        n_transversions = 0
        if 'Substitution_Type' in data.columns:
            n_transitions = len(data[data['Substitution_Type'] == 'Transition'])
            n_transversions = len(data[data['Substitution_Type'] == 'Transversion'])
        
        summary_data.append({
            'gene': gene,
            'n_combinations': len(data),
            'n_positions': len(data['Position'].unique()),
            'n_haplogroups': len(data['Haplogroup'].unique()),
            'max_frequency': data['Frequency'].max(),
            'mean_frequency': data['Frequency'].mean(),
            'max_total_frequency': position_stats['Total_Frequency'].max(),
            'mean_total_frequency': position_stats['Total_Frequency'].mean(),
            'n_significant': data['Significant_Corrected'].sum() if 'Significant_Corrected' in data.columns else 0,
            'n_transitions': n_transitions,
            'n_transversions': n_transversions
        })
    
    summary_df = pd.DataFrame(summary_data)
    
    # Create summary plots
    fig, axes = plt.subplots(2, 3, figsize=(18, 12))
    
    # Plot 1: Number of positions per gene
    ax1 = axes[0, 0]
    bars1 = ax1.bar(summary_df['gene'], summary_df['n_positions'])
    ax1.set_xlabel('Gene')
    ax1.set_ylabel('Number of Positions')
    ax1.tick_params(axis='x', rotation=45)
    ax1.grid(True, alpha=0.3)
    
    # Add value labels
    for bar, value in zip(bars1, summary_df['n_positions']):
        ax1.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.1,
                str(value), ha='center', va='bottom')
    
    # Plot 2: Maximum frequency per gene
    ax2 = axes[0, 1]
    bars2 = ax2.bar(summary_df['gene'], summary_df['max_frequency'])
    ax2.set_xlabel('Gene')
    ax2.set_ylabel('Maximum Frequency')
    ax2.tick_params(axis='x', rotation=45)
    ax2.grid(True, alpha=0.3)
    
    # Plot 3: Number of haplogroups per gene
    ax3 = axes[0, 2]
    bars3 = ax3.bar(summary_df['gene'], summary_df['n_haplogroups'])
    ax3.set_xlabel('Gene')
    ax3.set_ylabel('Number of Haplogroups')
    ax3.tick_params(axis='x', rotation=45)
    ax3.grid(True, alpha=0.3)
    
    # Plot 4: Transition vs Transversion counts
    ax4 = axes[1, 0]
    width = 0.35
    x = np.arange(len(summary_df))
    
    bars4a = ax4.bar(x - width/2, summary_df['n_transitions'], width, 
                     label='Transitions', alpha=0.7, color='blue')
    bars4b = ax4.bar(x + width/2, summary_df['n_transversions'], width, 
                     label='Transversions', alpha=0.7, color='red')
    
    ax4.set_xlabel('Gene')
    ax4.set_ylabel('Count')
    ax4.set_xticks(x)
    ax4.set_xticklabels(summary_df['gene'], rotation=45)
    ax4.legend()
    ax4.grid(True, alpha=0.3)
    
    # Plot 5: Significant mutations per gene
    ax5 = axes[1, 1]
    bars5 = ax5.bar(summary_df['gene'], summary_df['n_significant'])
    ax5.set_xlabel('Gene')
    ax5.set_ylabel('Number of Significant Mutations')
    ax5.tick_params(axis='x', rotation=45)
    ax5.grid(True, alpha=0.3)
    
    # Plot 6: Total mutation frequency per gene
    ax6 = axes[1, 2]
    bars6 = ax6.bar(summary_df['gene'], summary_df['max_total_frequency'])
    ax6.set_xlabel('Gene')
    ax6.set_ylabel('Maximum Total Frequency')
    ax6.tick_params(axis='x', rotation=45)
    ax6.grid(True, alpha=0.3)
    
    plt.tight_layout()
    
    # Save summary plot
    summary_file = os.path.join(output_dir, "haplogroup_mutation_summary.png")
    plt.savefig(summary_file, dpi=300, bbox_inches='tight')
    plt.show()
    
    print(f"Summary plot saved: {summary_file}")
    
    # Save summary statistics
    summary_csv = os.path.join(output_dir, "haplogroup_mutation_summary.csv")
    summary_df.to_csv(summary_csv, index=False)
    print(f"Summary statistics saved: {summary_csv}")

# Example usage
if __name__ == "__main__":
    print("Haplogroup mutation analysis script loaded successfully!")
    print("\nTo use this script:")
    print("1. Ensure your Excel file has an 'Associations' sheet with the specified columns")
    print("2. Ensure your population dataframe has columns named '{gene_name}' or '{gene_name}_seq'")
    print("3. Call analyze_haplogroup_mutations() with your data")
    print("4. Optionally call create_summary_plot() for an overview")

    print("\n" + "="*60)
    print("EXAMPLE USAGE WITH ABSOLUTE POSITIONS")
    print("="*60)
    print("""
# Define gene start positions (mtDNA Cambridge reference, 1-based)
gene_start_positions = {
    'RNR1': 1671,       # 12S rRNA
    'RNR2': 3230,       # 16S rRNA
    'ND1': 3307,        # NADH dehydrogenase subunit 1
    'ND2': 4470,        # NADH dehydrogenase subunit 2
    'COX1': 5904,       # Cytochrome c oxidase subunit I
    'COX2': 7586,       # Cytochrome c oxidase subunit II
    'ATP8': 8366,       # ATP synthase F0 subunit 8
    'ATP6': 8527,       # ATP synthase F0 subunit 6
    'COX3': 9207,       # Cytochrome c oxidase subunit III
    'ND3': 10059,       # NADH dehydrogenase subunit 3
    'ND4L': 10470,      # NADH dehydrogenase subunit 4L
    'ND4': 10760,       # NADH dehydrogenase subunit 4
    'ND4L-ND4': 10470,  # Combined ND4L-ND4 region
    'ND5': 12337,       # NADH dehydrogenase subunit 5
    'ND6': 14149,       # NADH dehydrogenase subunit 6
    'CYTB': 14747,      # Cytochrome b
}

# Run analysis with absolute position display
results = analyze_haplogroup_mutations(
    excel_file="haplogroup_associations.xlsx",
    hs_pop_df=hs_pop_df,
    ignore_trna=True,
    output_dir="haplogroup_mutation_plots",
    min_frequency=0.01,
    figsize=(16, 10),
    gene_start_positions=gene_start_positions,  # Enable absolute positions
    show_absolute_positions=True
)

# Create summary plot
create_summary_plot(results, output_dir="haplogroup_mutation_plots")
""")

    print("\n" + "="*60)
    print("SINGLE GENE PLOT WITH ABSOLUTE POSITIONS")
    print("="*60)
    print("""
# Plot a single gene with absolute positions
plot_gene_mutations(
    gene_data=gene_mutation_data,
    gene_name='RNR1',
    reference_seq=reference_sequence,
    output_dir='plots',
    figsize=(16, 10),
    gene_start_position=1671,  # RNR1 starts at mtDNA position 1671
    show_absolute_positions=True
)
""")

    print("\n" + "="*60)
    print("WITHOUT ABSOLUTE POSITIONS (relative only)")
    print("="*60)
    print("""
# Run analysis without absolute positions (original behavior)
results = analyze_haplogroup_mutations(
    excel_file="haplogroup_associations.xlsx",
    hs_pop_df=hs_pop_df,
    ignore_trna=True,
    output_dir="haplogroup_mutation_plots",
    min_frequency=0.01,
    figsize=(16, 10),
    gene_start_positions=None,  # No absolute positions
    show_absolute_positions=False
)
""")
