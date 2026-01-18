import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from collections import Counter
import logging
from typing import List, Tuple, Dict, Optional, Union
import os
from Bio.Seq import reverse_complement
from scipy.cluster.hierarchy import dendrogram, linkage
from scipy.spatial.distance import pdist
from scipy import stats
import matplotlib.patches as mpatches
import random
from multiprocessing import Pool, cpu_count
from functools import partial

# Try to import logomaker, provide fallback if not available
try:
    import logomaker
    LOGOMAKER_AVAILABLE = True
except ImportError:
    LOGOMAKER_AVAILABLE = False
    print("Warning: logomaker not available. Install with: pip install logomaker")

# Try to import tqdm for progress bars
try:
    from tqdm import tqdm
    TQDM_AVAILABLE = True
except ImportError:
    TQDM_AVAILABLE = False
    # Provide a dummy tqdm that just returns the iterable
    def tqdm(iterable, **kwargs):
        return iterable

# Set up logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

# Mitochondrial genetic code - differs from nuclear genetic code
MITOCHONDRIAL_GENETIC_CODE = {
    'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L',
    'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S',
    'TAT': 'Y', 'TAC': 'Y', 'TAA': '*', 'TAG': '*',
    'TGT': 'C', 'TGC': 'C', 'TGA': 'W', 'TGG': 'W',
    'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L',
    'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
    'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
    'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
    'ATT': 'I', 'ATC': 'I', 'ATA': 'M', 'ATG': 'M',  # ATA is Met in mitochondria
    'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
    'AAT': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K',
    'AGT': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',  # AGA/AGG are stop in mitochondria
    'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V',
    'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
    'GAT': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
    'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G'
}

# Nuclear/Standard genetic code
NUCLEAR_GENETIC_CODE = {
    'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L',
    'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S',
    'TAT': 'Y', 'TAC': 'Y', 'TAA': '*', 'TAG': '*',
    'TGT': 'C', 'TGC': 'C', 'TGA': '*', 'TGG': 'W',  # TGA is stop in nuclear code
    'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L',
    'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
    'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
    'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
    'ATT': 'I', 'ATC': 'I', 'ATA': 'I', 'ATG': 'M',  # ATA is Ile in nuclear code
    'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
    'AAT': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K',
    'AGT': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',  # AGA/AGG are Arg in nuclear code
    'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V',
    'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
    'GAT': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
    'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G'
}

def translate_sequence(sequence: str, genetic_code: str = 'mitochondrial') -> str:
    """Translate DNA sequence using specified genetic code."""
    if len(sequence) % 3 != 0:
        logger.warning(f"Sequence length {len(sequence)} is not divisible by 3")
        # Trim to nearest multiple of 3
        sequence = sequence[:len(sequence) - (len(sequence) % 3)]

    # Select genetic code
    if genetic_code.lower() == 'mitochondrial':
        codon_table = MITOCHONDRIAL_GENETIC_CODE
    elif genetic_code.lower() == 'nuclear':
        codon_table = NUCLEAR_GENETIC_CODE
    else:
        raise ValueError("genetic_code must be 'mitochondrial' or 'nuclear'")

    protein = []
    for i in range(0, len(sequence), 3):
        codon = sequence[i:i+3].upper()
        if len(codon) == 3:
            amino_acid = codon_table.get(codon, 'X')
            protein.append(amino_acid)
        else:
            break
    return ''.join(protein)

def build_synonymous_codon_table(genetic_code: str = 'mitochondrial') -> Dict[str, List[str]]:
    """
    Build a reverse codon table mapping each amino acid to its synonymous codons.

    Parameters:
    -----------
    genetic_code : str
        'mitochondrial' or 'nuclear' genetic code

    Returns:
    --------
    Dictionary mapping amino acids to list of synonymous codons
    """
    # Select genetic code
    if genetic_code.lower() == 'mitochondrial':
        codon_table = MITOCHONDRIAL_GENETIC_CODE
    elif genetic_code.lower() == 'nuclear':
        codon_table = NUCLEAR_GENETIC_CODE
    else:
        raise ValueError("genetic_code must be 'mitochondrial' or 'nuclear'")

    # Build reverse table
    synonymous_codons = {}
    for codon, aa in codon_table.items():
        if aa not in synonymous_codons:
            synonymous_codons[aa] = []
        synonymous_codons[aa].append(codon)

    return synonymous_codons

def shuffle_synonymous_codons(dna_sequence: str, genetic_code: str = 'mitochondrial') -> str:
    """
    Shuffle synonymous codons in a DNA sequence while preserving the amino acid sequence.

    Parameters:
    -----------
    dna_sequence : str
        DNA sequence to shuffle
    genetic_code : str
        Genetic code to use for translation

    Returns:
    --------
    DNA sequence with shuffled synonymous codons
    """
    # Ensure sequence length is divisible by 3
    if len(dna_sequence) % 3 != 0:
        dna_sequence = dna_sequence[:len(dna_sequence) - (len(dna_sequence) % 3)]

    # Translate original sequence
    aa_sequence = translate_sequence(dna_sequence, genetic_code)

    # Get synonymous codon table
    synonymous_codons = build_synonymous_codon_table(genetic_code)

    # Build shuffled DNA sequence
    shuffled_dna = []
    for aa in aa_sequence:
        # Get synonymous codons for this amino acid
        possible_codons = synonymous_codons.get(aa, [])
        if possible_codons:
            # Randomly select a synonymous codon
            selected_codon = random.choice(possible_codons)
            shuffled_dna.append(selected_codon)
        else:
            # If no synonymous codons found, keep original
            #logger.warning(f"No synonymous codons found for amino acid {aa}")
            pass
    return ''.join(shuffled_dna)

def calculate_ghost_peptide_conservation(
    dna_sequences: pd.Series,
    genetic_code: str = 'mitochondrial',
    frame_shift: int = 1
) -> List[float]:
    """
    Calculate conservation scores for "ghost peptides" - the same DNA sequences
    translated in a different reading frame. This serves as a background control.

    Parameters:
    -----------
    dna_sequences : pd.Series
        DNA sequences to translate in alternate frame
    genetic_code : str
        Genetic code to use for translation
    frame_shift : int
        Reading frame shift (+1 or +2)

    Returns:
    --------
    List of conservation scores for ghost peptides
    """
    # Shift sequences by the specified frame
    shifted_sequences = dna_sequences.apply(lambda x: x[frame_shift:] if len(x) > frame_shift else None).dropna()

    # Translate shifted sequences
    ghost_aa_sequences = []
    for seq in shifted_sequences:
        if pd.isna(seq) or len(seq) < 3:
            continue

        # Ensure sequence length is divisible by 3
        seq_len = len(seq)
        if seq_len % 3 != 0:
            # Trim to nearest multiple of 3
            seq = seq[:seq_len - (seq_len % 3)]
            logger.debug(f"Trimmed shifted sequence from {seq_len} to {len(seq)} bp (divisible by 3)")

        if len(seq) < 3:
            continue

        try:
            aa_seq = translate_sequence(seq, genetic_code)
            if aa_seq:
                ghost_aa_sequences.append(aa_seq)
        except Exception as e:
            logger.warning(f"Could not translate shifted sequence: {e}")
            continue

    if not ghost_aa_sequences:
        logger.warning("No valid ghost peptide sequences generated")
        return []

    # Ensure all sequences have the same length
    seq_lengths = [len(seq) for seq in ghost_aa_sequences]
    min_length = min(seq_lengths)
    ghost_aa_sequences = [seq[:min_length] for seq in ghost_aa_sequences]

    # Calculate conservation scores for ghost peptides
    ghost_conservation_scores = []

    for pos in range(min_length):
        # Get amino acids at this position
        amino_acids = [seq[pos] for seq in ghost_aa_sequences if pos < len(seq)]

        # Count amino acids
        aa_counts = Counter(amino_acids)
        total_seqs = len(amino_acids)

        # Calculate frequencies
        aa_frequencies = {aa: count/total_seqs for aa, count in aa_counts.items()}

        # Calculate Shannon entropy
        entropy = 0
        for freq in aa_frequencies.values():
            if freq > 0:
                entropy -= freq * np.log2(freq)

        # Calculate conservation (information content)
        max_entropy = np.log2(min(20.0, len(aa_frequencies)))
        conservation = max_entropy - entropy if max_entropy > 0 else 0

        ghost_conservation_scores.append(conservation)

    logger.info(f"Ghost peptide: {len(ghost_aa_sequences)} sequences, "
                f"{min_length} positions, mean conservation: {np.mean(ghost_conservation_scores):.3f} bits")

    return ghost_conservation_scores

def _run_single_permutation(perm_idx: int, dna_sequences: List[str], genetic_code: str) -> Optional[float]:
    """
    Helper function to run a single permutation. Designed to be called in parallel.

    Parameters:
    -----------
    perm_idx : int
        Permutation index (for logging/debugging)
    dna_sequences : List[str]
        List of DNA sequences to shuffle
    genetic_code : str
        Genetic code to use

    Returns:
    --------
    Mean conservation score for this permutation, or None if failed
    """
    try:
        # Shuffle synonymous codons for each DNA sequence
        shuffled_sequences = []
        for seq in dna_sequences:
            if pd.isna(seq) or len(seq) == 0:
                continue
            try:
                shuffled_seq = shuffle_synonymous_codons(seq, genetic_code)
                shuffled_sequences.append(shuffled_seq)
            except Exception:
                continue

        # Calculate conservation for Frame +1 of shuffled sequences
        if shuffled_sequences:
            shuffled_series = pd.Series(shuffled_sequences)
            permuted_conservation = calculate_ghost_peptide_conservation(
                shuffled_series, genetic_code, frame_shift=1
            )

            if permuted_conservation:
                # Return mean conservation for this permutation
                return np.mean(permuted_conservation)
    except Exception as e:
        logger.warning(f"Permutation {perm_idx} failed: {e}")
        return None

    return None

def create_information_matrix(amino_acid_sequences: List[str]) -> pd.DataFrame:
    """
    Create information content matrix for logomaker.
    
    Parameters:
    -----------
    amino_acid_sequences : List[str]
        List of amino acid sequences
        
    Returns:
    --------
    pd.DataFrame with positions as index and amino acids as columns
    """
    if not amino_acid_sequences:
        return pd.DataFrame()
    
    # Get sequence length
    seq_length = len(amino_acid_sequences[0])
    
    # Standard amino acids plus stop (as X)
    amino_acids = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 
                   'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y', 'X']
    
    # Initialize information matrix
    info_matrix = pd.DataFrame(0.0, index=range(1, seq_length + 1), columns=amino_acids)
    
    # Calculate information content for each position
    for pos in range(seq_length):
        # Get amino acids at this position
        amino_acids_at_pos = [seq[pos] if pos < len(seq) else 'X' for seq in amino_acid_sequences]
        # Convert stop codons to X
        amino_acids_at_pos = ['X' if aa == '*' else aa for aa in amino_acids_at_pos]
        
        # Count amino acids
        aa_counts = Counter(amino_acids_at_pos)
        total_seqs = len(amino_acids_at_pos)
        
        # Calculate frequencies
        aa_frequencies = {aa: count/total_seqs for aa, count in aa_counts.items()}
        
        # Calculate conservation (Shannon entropy)
        entropy = 0
        for freq in aa_frequencies.values():
            if freq > 0:
                entropy -= freq * np.log2(freq)
        # Convert to conservation score (bits)
 #       if len(aa_frequencies) == 1:
        max_entropy = np.log2(20.0)  # Assume max entropy for 20 amino acids
#        else:
 #           max_entropy =np.log2(min(20.0, len(aa_frequencies)))  # Limit to 20 amino acids
        conservation = max_entropy - entropy if max_entropy > 0 else 0
        
        # Fill information matrix (frequency * conservation for each amino acid)
        for aa, freq in aa_frequencies.items():
            if aa in info_matrix.columns:
                info_matrix.loc[pos + 1, aa] = freq * conservation
    
    return info_matrix

def create_microprotein_sequence_logo(
    sequences: pd.Series,
    genetic_code: str = 'mitochondrial',
    strand: str = '+',
    title: Optional[str] = None,
    figsize: Tuple[int, int] = (12, 6),
    color_scheme: str = 'chemistry',
    show_conservation: bool = True,
    min_frequency: float = 0.01,
    output_file: Optional[str] = None,
    dpi: int = 300,
    font_scale: float = 1.0,
    show_numbers: bool = True,
    background_color: str = 'white',
    information_box: bool = True,
    statistical_test: bool = False,
    n_permutations: int = 100,
    use_parallel: bool = True,
    n_jobs: Optional[int] = None
) -> Dict:
    """
    Create a sequence logo for amino acid sequences from a micro-protein region using logomaker.

    Parameters:
    -----------
    sequences : pd.Series
        Series containing DNA sequences for the micro-protein region.
    genetic_code : str
        'mitochondrial' or 'nuclear' genetic code to use for translation.
    strand : str
        '+' for forward strand, '-' for reverse strand.
    title : str, optional
        Title for the plot.
    figsize : Tuple[int, int]
        Figure size for the plot.
    color_scheme : str
        Color scheme for amino acids.
    show_conservation : bool
        Whether to show conservation scores above the logo.
    min_frequency : float
        Minimum frequency for an amino acid to be displayed.
    output_file : str, optional
        File to save the plot.
    dpi : int
        DPI for saving the plot.
    font_scale : float
        Scaling factor for fonts.
    show_numbers : bool
        Whether to show position numbers.
    background_color : str
        Background color.
    information_box : bool
        Whether to show the information box with statistics.
    statistical_test : bool
        If True, perform a permutation test comparing the conservation scores
        to a null distribution generated by shuffling synonymous codons and translating
        in +1 reading frame. This distinguishes conservation due to functional constraints
        from background conservation of the underlying DNA region.
    n_permutations : int
        Number of permutations to generate for the null distribution (default: 100).
        Only used if statistical_test=True.
    use_parallel : bool
        If True (default), use multiprocessing to run permutations in parallel.
        This can significantly speed up the statistical test.
    n_jobs : int, optional
        Number of parallel jobs to run. If None (default), uses all available CPU cores.
        Only used if use_parallel=True and statistical_test=True.

    Returns:
    --------
    Dict with analysis results.
    """
    
    if not LOGOMAKER_AVAILABLE:
        raise ImportError("logomaker is required. Install with: pip install logomaker")
    
    microprotein_name = sequences.name if sequences.name else "MicroProtein"
    logger.info(f"Creating sequence logo for {microprotein_name} with {len(sequences)} sequences")
    
    # Handle reverse strand
    if strand == '-':
        sequences = sequences.apply(reverse_complement)
        logger.info("Applied reverse complement for negative strand")

    # Keep DNA sequences for statistical test (after reverse complement if applicable)
    dna_sequences = sequences.copy()

    # Translate sequences to amino acids
    amino_acid_sequences = []
    for seq in sequences:
        if pd.isna(seq) or len(seq) == 0:
            continue
        try:
            aa_seq = translate_sequence(seq, genetic_code)
            if aa_seq:
                amino_acid_sequences.append(aa_seq)
        except Exception as e:
            logger.warning(f"Could not translate sequence {seq[:20]}...: {e}")
            continue
    
    if not amino_acid_sequences:
        logger.error("No valid amino acid sequences generated")
        return {}
    
    logger.info(f"Successfully translated {len(amino_acid_sequences)} sequences")
    
    # Ensure all sequences have the same length
    seq_lengths = [len(seq) for seq in amino_acid_sequences]
    max_length = max(seq_lengths)
    min_length = min(seq_lengths)
    
    if max_length != min_length:
        logger.warning(f"Sequences have variable lengths ({min_length}-{max_length}). Using minimum length.")
        amino_acid_sequences = [seq[:min_length] for seq in amino_acid_sequences]
        max_length = min_length
    
    # Create information matrix for logomaker
    info_matrix = create_information_matrix(amino_acid_sequences)
    
    # Calculate position-wise statistics for results
    position_frequencies = []
    conservation_scores = []
    
    for pos in range(max_length):
        amino_acids_at_pos = [seq[pos] for seq in amino_acid_sequences if pos < len(seq)]
        amino_acids_at_pos = ['X' if aa == '*' else aa for aa in amino_acids_at_pos]
        
        aa_counts = Counter(amino_acids_at_pos)
        total_seqs = len(amino_acids_at_pos)
        aa_frequencies = {aa: count/total_seqs for aa, count in aa_counts.items()}
        
        # Calculate conservation
        entropy = 0
        for freq in aa_frequencies.values():
            if freq > 0:
                entropy -= freq * np.log2(freq)
         # Convert to conservation score (bits)
     #   if len(aa_frequencies) == 1:
        max_entropy = np.log2(20.0)  # Assume max entropy for 20 amino acids
    #    else:
     #       max_entropy =np.log2(min(20.0, len(aa_frequencies)))  # Limit to 20 amino acids
        conservation = max_entropy - entropy if max_entropy > 0 else 0
        
        position_frequencies.append(aa_frequencies)
        conservation_scores.append(conservation)

    # Perform statistical test if requested
    statistical_test_results = None
    if statistical_test:
        parallel_str = f" (parallel on {n_jobs if n_jobs else cpu_count()} cores)" if use_parallel else " (sequential)"
        progress_str = " with progress bar" if TQDM_AVAILABLE else ""
        logger.info(f"Performing permutation test with {n_permutations} permutations{parallel_str}{progress_str}...")

        # Calculate mean conservation for actual Frame 0
        mean_actual = np.mean(conservation_scores)
        std_actual = np.std(conservation_scores)

        # Convert dna_sequences to list for multiprocessing
        dna_list = dna_sequences.dropna().tolist()

        # Generate null distribution by permutation
        if use_parallel:
            # Use multiprocessing for parallel execution
            n_processes = n_jobs if n_jobs else cpu_count()
            logger.info(f"  Using {n_processes} parallel processes...")

            # Create partial function with fixed arguments
            permutation_func = partial(_run_single_permutation,
                                      dna_sequences=dna_list,
                                      genetic_code=genetic_code)

            # Run permutations in parallel with progress bar
            with Pool(processes=n_processes) as pool:
                # Use imap for better progress tracking
                results = list(tqdm(
                    pool.imap(permutation_func, range(n_permutations)),
                    total=n_permutations,
                    desc="  Permutations",
                    unit="perm",
                    disable=not TQDM_AVAILABLE
                ))

            # Filter out None values
            null_distribution = [r for r in results if r is not None]

        else:
            # Sequential execution with progress bar
            null_distribution = []
            for perm_idx in tqdm(range(n_permutations),
                                desc="  Permutations",
                                unit="perm",
                                disable=not TQDM_AVAILABLE):
                result = _run_single_permutation(perm_idx, dna_list, genetic_code)
                if result is not None:
                    null_distribution.append(result)

        if null_distribution:
            # Calculate p-value: proportion of permutations with conservation >= actual
            null_distribution = np.array(null_distribution)
            mean_background = np.mean(null_distribution)
            std_background = np.std(null_distribution)

            # Two-tailed test: count how many permutations are as extreme as actual
            # For conservation, we expect actual > background if functionally constrained
            n_greater_equal = np.sum(null_distribution >= mean_actual)
            n_less_equal = np.sum(null_distribution <= mean_actual)

            # One-tailed p-value (testing if actual > background)
            p_value_one_tailed = n_greater_equal / len(null_distribution)

            # Two-tailed p-value (testing if actual != background)
            p_value_two_tailed = 2 * min(p_value_one_tailed, 1 - p_value_one_tailed)

            # Use one-tailed test since we expect functional constraint to increase conservation
            p_value = p_value_one_tailed

            # Calculate effect size (Cohen's d)
            pooled_std = np.sqrt((std_actual**2 + std_background**2) / 2)
            cohens_d = (mean_actual - mean_background) / pooled_std if pooled_std > 0 else 0

            # Determine significance and interpretation
            is_significant = p_value < 0.05
            if is_significant:
                if mean_actual > mean_background:
                    interpretation = 'Actual conservation is significantly HIGHER than permuted background (functional constraint)'
                else:
                    interpretation = 'Actual conservation is significantly LOWER than permuted background'
            else:
                interpretation = 'No significant difference from permuted background (may not be functionally constrained)'

            statistical_test_results = {
                'test_type': 'Permutation test',
                'n_permutations': len(null_distribution),
                'null_distribution': null_distribution.tolist(),
                'mean_background': mean_background,
                'std_background': std_background,
                'mean_actual': mean_actual,
                'std_actual': std_actual,
                'p_value': p_value,
                'p_value_two_tailed': p_value_two_tailed,
                'cohens_d': cohens_d,
                'significant': is_significant,
                'interpretation': interpretation,
                'n_greater_equal': int(n_greater_equal),
                'n_less_equal': int(n_less_equal)
            }

            logger.info(f"Permutation test results:")
            logger.info(f"  Test: Permutation test ({len(null_distribution)} permutations)")
            logger.info(f"  Mean actual conservation: {mean_actual:.3f} ± {std_actual:.3f} bits")
            logger.info(f"  Mean background (permuted): {mean_background:.3f} ± {std_background:.3f} bits")
            logger.info(f"  p-value (one-tailed): {p_value:.4f}")
            logger.info(f"  p-value (two-tailed): {p_value_two_tailed:.4f}")
            logger.info(f"  Cohen's d: {cohens_d:.3f}")
            logger.info(f"  Permutations >= actual: {n_greater_equal}/{len(null_distribution)}")
            logger.info(f"  Interpretation: {interpretation}")
        else:
            logger.warning("Could not generate null distribution for permutation test")

    # Define color schemes for logomaker
    color_schemes = {
        'chemistry': {
            'A': '#C8C8C8', 'I': '#C8C8C8', 'L': '#C8C8C8', 'M': '#C8C8C8', 'F': '#C8C8C8',
            'W': '#C8C8C8', 'Y': '#C8C8C8', 'V': '#C8C8C8', 'P': '#C8C8C8',
            'S': '#00FF00', 'T': '#00FF00', 'Q': '#00FF00', 'N': '#00FF00', 'C': '#00FF00',
            'K': '#0000FF', 'R': '#0000FF', 'H': '#0000FF',
            'D': '#FF0000', 'E': '#FF0000',
            'G': '#FF8C00', 'X': '#000000'
        },
        'hydrophobicity': {
            'A': '#1f77b4', 'I': '#1f77b4', 'L': '#1f77b4', 'M': '#1f77b4', 'F': '#1f77b4',
            'W': '#1f77b4', 'Y': '#1f77b4', 'V': '#1f77b4', 'P': '#1f77b4',
            'S': '#ff7f0e', 'T': '#ff7f0e', 'Q': '#ff7f0e', 'N': '#ff7f0e', 'C': '#ff7f0e',
            'K': '#ff7f0e', 'R': '#ff7f0e', 'H': '#ff7f0e', 'D': '#ff7f0e', 'E': '#ff7f0e',
            'G': '#2ca02c', 'X': '#000000'
        },
        'charge': {
            'K': '#0000FF', 'R': '#0000FF', 'H': '#ADD8E6',
            'D': '#FF0000', 'E': '#FF0000',
            'A': '#808080', 'I': '#808080', 'L': '#808080', 'M': '#808080', 'F': '#808080',
            'W': '#808080', 'Y': '#808080', 'V': '#808080', 'P': '#808080', 'S': '#808080',
            'T': '#808080', 'Q': '#808080', 'N': '#808080', 'C': '#808080', 'G': '#808080',
            'X': '#000000'
        }
    }
    
    # Set up the plot
    plt.style.use('default')
    fig = plt.figure(figsize=figsize, facecolor=background_color)
    
    if show_conservation:
        # Create subplots with conservation plot on top
        gs = fig.add_gridspec(2, 1, height_ratios=[1, 4], hspace=0.1)
        ax_cons = fig.add_subplot(gs[0])
        ax_logo = fig.add_subplot(gs[1])
        
        # Plot conservation scores
        positions = list(range(1, max_length + 1))
        ax_cons.plot(positions, conservation_scores, 'b-', linewidth=2, alpha=0.8)
        ax_cons.fill_between(positions, conservation_scores, alpha=0.3, color='lightblue')
        ax_cons.set_ylabel('Conservation\n(bits)', fontsize=10, fontweight='bold')
        ax_cons.set_xlim(0.5, max_length + 0.5)
        ax_cons.grid(True, alpha=0.3)
        ax_cons.spines['top'].set_visible(False)
        ax_cons.spines['right'].set_visible(False)
        ax_cons.set_xticklabels([])
        
        if title is None:
            title = f'Sequence Logo for {microprotein_name}'
        ax_cons.set_title(title, fontsize=14, fontweight='bold', pad=15)
    else:
        ax_logo = fig.add_subplot(1, 1, 1)
        if title is None:
            title = f'Sequence Logo for {microprotein_name}'
        ax_logo.set_title(title, fontsize=14, fontweight='bold', pad=15)
    
    # Create the sequence logo using logomaker
    if color_scheme in color_schemes:
        colors = color_schemes[color_scheme]
    else:
        colors = color_schemes['chemistry']
        logger.warning(f"Unknown color scheme '{color_scheme}', using 'chemistry'")
    # Create logo
    logo = logomaker.Logo(info_matrix, 
                         ax=ax_logo,
                         color_scheme=colors,
                         font_name='monospace',
                         alpha=0.8,
                         vpad=0.1,
                         show_spines = False)
    
    # Customize the logo
    logo.style_spines(visible=False)
    logo.style_spines(spines=['left', 'bottom'], visible=True)
    
    # Set axis labels and formatting
    ax_logo.set_ylabel('Information Content (bits)', fontsize=12, fontweight='bold')
    ax_logo.set_xlabel('Position', fontsize=12, fontweight='bold')
    
    if show_numbers:
        n_aa = len(range(1, max_length + 1))
        sep = 1
        if n_aa > 50:
            sep = 3
        elif n_aa > 100:
            sep = 7

        ax_logo.set_xticks(range(1, max_length + 1, sep))
        ax_logo.set_xticklabels(range(1, max_length + 1, sep))
    
    
    # Add grid
    ax_logo.grid(True, alpha=0.2, axis='both', linestyle='-', linewidth=0.5)
    
    # Add color legend based on scheme
    legend_elements = []
    if color_scheme == 'chemistry':
        legend_elements = [
            plt.Rectangle((0,0),1,1, facecolor='#C8C8C8', label='Hydrophobic'),
            plt.Rectangle((0,0),1,1, facecolor='#00FF00', label='Polar'),
            plt.Rectangle((0,0),1,1, facecolor='#0000FF', label='Basic'),
            plt.Rectangle((0,0),1,1, facecolor='#FF0000', label='Acidic'),
            plt.Rectangle((0,0),1,1, facecolor='#FF8C00', label='Glycine'),
            plt.Rectangle((0,0),1,1, facecolor='#000000', label='Stop (X)')
        ]
    elif color_scheme == 'hydrophobicity':
        legend_elements = [
            plt.Rectangle((0,0),1,1, facecolor='#1f77b4', label='Hydrophobic'),
            plt.Rectangle((0,0),1,1, facecolor='#ff7f0e', label='Hydrophilic'),
            plt.Rectangle((0,0),1,1, facecolor='#2ca02c', label='Glycine'),
            plt.Rectangle((0,0),1,1, facecolor='#000000', label='Stop (X)')
        ]
    elif color_scheme == 'charge':
        legend_elements = [
            plt.Rectangle((0,0),1,1, facecolor='#0000FF', label='Positive'),
            plt.Rectangle((0,0),1,1, facecolor='#FF0000', label='Negative'),
            plt.Rectangle((0,0),1,1, facecolor='#808080', label='Neutral'),
            plt.Rectangle((0,0),1,1, facecolor='#000000', label='Stop (X)')
        ]
    
    if legend_elements:
        ax_logo.legend(handles=legend_elements, loc='upper right', fontsize=9, 
                      framealpha=0.8, fancybox=True)
    
    # Add sequence information box
    info_text = f"Sequences: {len(amino_acid_sequences)}\n"
    info_text += f"Length: {max_length} AA\n"
    info_text += f"Genetic code: {genetic_code}\n"
    info_text += f"Mean conservation: {np.mean(conservation_scores):.2f} bits\n"
    info_text += f"Max conservation: {max(conservation_scores):.2f} bits"
    
    # Calculate the most conserved position
    if conservation_scores:
        most_conserved_pos = np.argmax(conservation_scores) + 1
        info_text += f"\nMost conserved: pos {most_conserved_pos}"
    
    # Count invariant positions
    invariant_positions = sum(1 for freq_dict in position_frequencies
                             if len(freq_dict) > 0 and max(freq_dict.values()) > 0.95)
    if invariant_positions > 0:
        info_text += f"\nInvariant positions: {invariant_positions} (*)"

    # Add statistical test results if performed
    if statistical_test_results:
        info_text += f"\n\n--- Statistical Test ---"
        info_text += f"\n{statistical_test_results['test_type']}"
        info_text += f"\nBackground: {statistical_test_results['mean_background']:.2f} bits"
        info_text += f"\np-value: {statistical_test_results['p_value']:.4f}"
        if statistical_test_results['significant']:
            info_text += f"\n*** SIGNIFICANT ***" if statistical_test_results['mean_actual'] > statistical_test_results['mean_background'] else f"\n* Significant *"

    if information_box:
        ax_logo.text(0.02, 0.98, info_text, transform=ax_logo.transAxes,
                fontsize=8 if statistical_test_results else 9, verticalalignment='top',
                bbox=dict(boxstyle='round,pad=0.5',
                         facecolor='lightgreen' if (statistical_test_results and statistical_test_results['significant'] and statistical_test_results['mean_actual'] > statistical_test_results['mean_background']) else 'lightgray',
                         alpha=0.5))
    
    plt.tight_layout()
    
    if output_file:
        plt.savefig(output_file, dpi=dpi, bbox_inches='tight', facecolor=background_color)
        logger.info(f"Sequence logo saved to {output_file}")
    else:
        plt.show()
    
    # Calculate variant counts exceeding min_frequency threshold
    variants_per_position = []
    for freq_dict in position_frequencies:
        # Count how many amino acids exceed the min_frequency threshold
        n_variants = sum(1 for freq in freq_dict.values() if freq >= min_frequency)
        variants_per_position.append(n_variants)

    # Calculate summary statistics for variants
    total_variants = sum(variants_per_position) if variants_per_position else 0
    variants_summary = {
        'total_variants': total_variants,  # Total number of variants exceeding threshold across all positions
        'mean_variants_per_position': np.mean(variants_per_position) if variants_per_position else 0,
        'median_variants_per_position': np.median(variants_per_position) if variants_per_position else 0,
        'max_variants_per_position': max(variants_per_position) if variants_per_position else 0,
        'min_variants_per_position': min(variants_per_position) if variants_per_position else 0,
        'total_positions': len(variants_per_position),
        'positions_with_multiple_variants': sum(1 for v in variants_per_position if v > 1)
    }

    # Prepare results
    results = {
        'microprotein_name': microprotein_name,
        'n_sequences': len(amino_acid_sequences),
        'sequence_length': max_length,
        'genetic_code': genetic_code,
        'strand': strand,
        'position_frequencies': position_frequencies,
        'conservation_scores': conservation_scores,
        'mean_conservation': np.mean(conservation_scores),
        'max_conservation': max(conservation_scores),
        'most_conserved_position': np.argmax(conservation_scores) + 1 if conservation_scores else None,
        'amino_acid_sequences': amino_acid_sequences,
        'color_scheme': color_scheme,
        'invariant_positions': invariant_positions,
        'information_matrix': info_matrix,
        'statistical_test': statistical_test_results,
        'min_frequency_threshold': min_frequency,
        'variants_per_position': variants_per_position,
        'variants_summary': variants_summary
    }

    return results

def create_multiple_sequence_logos(
    sequences_dict: Dict[str, pd.Series],
    genetic_codes: Dict[str, str],
    output_dir: Optional[str] = None,
    **kwargs
) -> Tuple[Dict[str, Dict], pd.DataFrame]:
    """
    Create sequence logos for multiple micro-protein regions.

    Parameters:
    -----------
    sequences_dict : Dict[str, pd.Series]
        Dictionary mapping micro-protein names to their sequence Series.
    genetic_codes : Dict[str, str]
        Dictionary mapping micro-protein names to their genetic codes.
    output_dir : str, optional
        Directory to save individual logo files.
    **kwargs
        Additional arguments to pass to create_microprotein_sequence_logo.

    Returns:
    --------
    Tuple containing:
        - Dict mapping micro-protein names to their logo analysis results
        - DataFrame with summary statistics including variant counts for each region
    """
    results = {}

    if output_dir and not os.path.exists(output_dir):
        os.makedirs(output_dir)
        logger.info(f"Created output directory: {output_dir}")

    for mp_name, sequences in sequences_dict.items():
        logger.info(f"Creating sequence logo for {mp_name}")

        # Determine genetic code
        genetic_code = genetic_codes.get(mp_name, 'mitochondrial')

        # Set output file if directory provided
        output_file = None
        if output_dir:
            safe_name = mp_name.replace('/', '_').replace('\\', '_')
            output_file = os.path.join(output_dir, f"{safe_name}_sequence_logo.png")

        # Create the logo
        try:
            mp_results = create_microprotein_sequence_logo(
                sequences=sequences,
                genetic_code=genetic_code,
                output_file=output_file,
                **kwargs
            )
            results[mp_name] = mp_results
        except Exception as e:
            logger.error(f"Error creating logo for {mp_name}: {e}")
            continue

    logger.info(f"Created sequence logos for {len(results)} micro-proteins")

    # Create summary DataFrame with variant information
    summary_data = []
    for mp_name, mp_results in results.items():
        variants_summary = mp_results.get('variants_summary', {})
        summary_row = {
            'microprotein_name': mp_name,
            'n_sequences': mp_results['n_sequences'],
            'sequence_length_aa': mp_results['sequence_length'],
            'genetic_code': mp_results['genetic_code'],
            'mean_conservation_bits': mp_results['mean_conservation'],
            'max_conservation_bits': mp_results['max_conservation'],
            'invariant_positions': mp_results.get('invariant_positions', 0),
            'min_frequency_threshold': mp_results.get('min_frequency_threshold', 0.01),
            'total_variants_above_threshold': variants_summary.get('total_variants', 0),
            'mean_variants_per_position': variants_summary.get('mean_variants_per_position', 0),
            'median_variants_per_position': variants_summary.get('median_variants_per_position', 0),
            'max_variants_per_position': variants_summary.get('max_variants_per_position', 0),
            'min_variants_per_position': variants_summary.get('min_variants_per_position', 0),
            'positions_with_multiple_variants': variants_summary.get('positions_with_multiple_variants', 0),
            'percent_positions_with_multiple_variants': (
                100 * variants_summary.get('positions_with_multiple_variants', 0) /
                variants_summary.get('total_positions', 1) if variants_summary.get('total_positions', 0) > 0 else 0
            )
        }

        # Add statistical test information if available
        if mp_results.get('statistical_test'):
            stat_test = mp_results['statistical_test']
            summary_row['stat_test_pvalue'] = stat_test.get('p_value', np.nan)
            summary_row['stat_test_significant'] = stat_test.get('significant', False)
            summary_row['mean_background_conservation'] = stat_test.get('mean_background', np.nan)

        summary_data.append(summary_row)

    summary_df = pd.DataFrame(summary_data)

    # Log summary
    logger.info(f"\nVariant Summary (>{summary_df['min_frequency_threshold'].iloc[0]*100:.0f}% frequency threshold):")
    for _, row in summary_df.iterrows():
        logger.info(f"  {row['microprotein_name']}: "
                   f"Total={int(row['total_variants_above_threshold'])} variants, "
                   f"Mean={row['mean_variants_per_position']:.2f} variants/position, "
                   f"{row['percent_positions_with_multiple_variants']:.1f}% positions have multiple variants")

    return results, summary_df

def compare_conservation_across_microproteins(
    results_dict: Dict[str, Dict],
    output_file: Optional[str] = None,
    figsize: Tuple[int, int] = (14, 8)
) -> None:
    """
    Create a comparison plot showing conservation across multiple micro-proteins.
    
    Parameters:
    -----------
    results_dict : Dict[str, Dict]
        Dictionary mapping micro-protein names to their sequence logo results.
    output_file : str, optional
        File to save the comparison plot.
    figsize : Tuple[int, int]
        Figure size for the plot.
    """
    if not results_dict:
        logger.error("No results provided for comparison")
        return
    
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=figsize)
    
    # Plot 1: Mean conservation comparison
    mp_names = list(results_dict.keys())
    mean_conservations = [results['mean_conservation'] for results in results_dict.values()]
    genetic_codes = [results['genetic_code'] for results in results_dict.values()]
    invariant_counts = [results.get('invariant_positions', 0) for results in results_dict.values()]
    
    colors = ['blue' if gc == 'mitochondrial' else 'orange' for gc in genetic_codes]
    
    bars = ax1.bar(range(len(mp_names)), mean_conservations, color=colors, alpha=0.7)
    ax1.set_xticks(range(len(mp_names)))
    ax1.set_xticklabels(mp_names, rotation=45, ha='right')
    ax1.set_ylabel('Mean Conservation (bits)')
    ax1.set_title('Conservation Comparison Across Micro-Proteins')
    ax1.grid(True, alpha=0.3)
    
    # Add invariant position count on top of bars
    for i, (bar, inv_count) in enumerate(zip(bars, invariant_counts)):
        if inv_count > 0:
            ax1.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.02, 
                    f'{inv_count}*', ha='center', va='bottom', fontsize=9, 
                    fontweight='bold', color='red')
    
    # Add genetic code legend
    from matplotlib.patches import Patch
    legend_elements = [
        Patch(facecolor='blue', alpha=0.7, label='Mitochondrial'),
        Patch(facecolor='orange', alpha=0.7, label='Nuclear')
    ]
    ax1.legend(handles=legend_elements, loc='upper right')
    ax1.text(0.98, 0.85, '* = Invariant positions', transform=ax1.transAxes, 
            ha='right', va='top', fontsize=9, color='red', fontweight='bold')
    
    # Plot 2: Conservation profiles
    for mp_name, results in results_dict.items():
        scores = results['conservation_scores']
        positions = list(range(1, len(scores) + 1))
        gc = results['genetic_code']
        color = 'blue' if gc == 'mitochondrial' else 'orange'
        ax2.plot(positions, scores, label=f"{mp_name} ({gc})", 
                color=color, alpha=0.8, linewidth=2)
        
        # Mark invariant positions with asterisks
        freq_data = results.get('position_frequencies', [])
        for pos, freq_dict in enumerate(freq_data):
            if len(freq_dict) > 0 and max(freq_dict.values()) > 0.95:
                if pos < len(scores):
                    ax2.scatter(pos + 1, scores[pos], marker='*', s=100, 
                              color=color, edgecolor='red', linewidth=1)
    
    ax2.set_xlabel('Position')
    ax2.set_ylabel('Conservation (bits)')
    ax2.set_title('Conservation Profiles')
    ax2.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
    ax2.grid(True, alpha=0.3)
    
    plt.tight_layout()
    
    if output_file:
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        logger.info(f"Comparison plot saved to {output_file}")
    else:
        plt.show()

def extract_microprotein_sequences(
    dataframe: pd.DataFrame,
    sequence_column: str,
    microprotein_regions: List[List],
    microprotein_names: Optional[List[str]] = None
) -> Dict[str, pd.Series]:
    """
    Extract micro-protein sequences from a dataframe based on coordinates.
    
    Parameters:
    -----------
    dataframe : pd.DataFrame
        DataFrame containing the sequences.
    sequence_column : str
        Name of the column containing sequences.
    microprotein_regions : List[List]
        List of [start, end] coordinates (1-based) for each micro-protein.
    microprotein_names : List[str], optional
        Names for each micro-protein. If None, will use generic names.
    
    Returns:
    --------
    Dict mapping micro-protein names to their sequence Series.
    """
    sequences_dict = {}
    
    for i, coords in enumerate(microprotein_regions):
        if len(coords) >= 3:
            start, end, name, strand = coords[0], coords[1], coords[2], coords[3] if len(coords) > 3 else '+'
        else:
            start, end = coords[0], coords[1]
            name = microprotein_names[i] if microprotein_names and i < len(microprotein_names) else f"MP_{i+1}"
            strand = '+'  # Default strand
        
        # Extract sequences (convert to 0-based indexing)
        mp_sequences = dataframe[sequence_column].apply(
            lambda x: x[start:end] if isinstance(x, str) and len(x) >= end else None
        ).dropna()
        if strand == '-':
            mp_sequences = mp_sequences.apply(reverse_complement)
            logger.info(f"Applied reverse complement for {name} ({start}-{end})")
        mp_sequences.name = name
        sequences_dict[name] = mp_sequences
        
        logger.info(f"Extracted {len(mp_sequences)} sequences for {name} ({start}-{end}) ")
    
    return sequences_dict

def create_mdp_heatmap_summary(mdp_results_dict, figsize=(16, 10), output_file=None):
    """
    Create a heatmap-based summary of MDP conservation patterns with statistical analysis.
    
    Parameters:
    -----------
    mdp_results_dict : dict
        Dictionary mapping MDP names to their analysis results
    figsize : tuple
        Figure size
    output_file : str, optional
        Path to save the plot
    """
    
    # Prepare data for heatmap
    mdp_names = list(mdp_results_dict.keys())
    max_length = max(len(results['conservation_scores']) for results in mdp_results_dict.values())
    
    # Create conservation matrix (pad shorter sequences with NaN)
    conservation_matrix = []
    for mdp_name in mdp_names:
        scores = mdp_results_dict[mdp_name]['conservation_scores']
        # Pad with NaN if shorter than max length
        padded_scores = scores + [np.nan] * (max_length - len(scores))
        conservation_matrix.append(padded_scores)
    
    conservation_df = pd.DataFrame(conservation_matrix, 
                                 index=mdp_names,
                                 columns=[f'Pos_{i+1}' for i in range(max_length)])
    
    # Create figure with subplots
    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=figsize,
                                                gridspec_kw={'width_ratios': [3, 1],
                                                            'height_ratios': [1, 2]})
    
    # 1. Main heatmap of conservation scores
    mask = conservation_df.isna()
    sns.heatmap(conservation_df, mask=mask, annot=False, cmap='viridis', 
                cbar_kws={'label': 'Conservation (bits)'}, ax=ax1,
                xticklabels=5, yticklabels=True)
    ax1.set_xlabel('Position', fontsize=12)
    ax1.set_ylabel('MDP', fontsize=12)
    
    # 2. Summary statistics
    summary_stats = []
    for mdp_name, results in mdp_results_dict.items():
        # Extract statistical test results if available
        stat_test = results.get('statistical_test', None)

        stats = {
            'MDP': mdp_name,
            'Length': results['sequence_length'],
            'Mean_Conservation': results['mean_conservation'],
            'Max_Conservation': results['max_conservation'],
            'Min_Conservation': min(results['conservation_scores']),
            'Conservation_CV': np.std(results['conservation_scores']) / np.mean(results['conservation_scores']),
            'N_Sequences': results['n_sequences'],
            'Genetic_Code': results.get('genetic_code', 'unknown'),
            'Stat_Test_Performed': stat_test is not None,
            'Stat_Test_Significant': stat_test['significant'] if stat_test else False,
            'Stat_Test_PValue': stat_test['p_value'] if stat_test else np.nan,
            'Mean_Background': stat_test['mean_background'] if stat_test else np.nan,
            'Test_Type': stat_test['test_type'] if stat_test else None
        }
        summary_stats.append(stats)

    stats_df = pd.DataFrame(summary_stats)

    # Plot mean conservation comparison with statistical significance indication
    # Color bars based on statistical test results
    bar_colors = []
    for _, row in stats_df.iterrows():
        if row['Stat_Test_Performed']:
            if row['Stat_Test_Significant']:
                # Green for significantly higher than background
                if row['Mean_Conservation'] > row['Mean_Background']:
                    bar_colors.append('green')
                else:
                    bar_colors.append('red')
            else:
                # Gray for not significant
                bar_colors.append('gray')
        else:
            # Orange for no test performed
            bar_colors.append('orange')

    bars = ax2.barh(stats_df['MDP'], stats_df['Mean_Conservation'],
                   color=bar_colors, alpha=0.7)

    # Add background conservation as vertical lines if statistical test was performed
    for i, (_, row) in enumerate(stats_df.iterrows()):
        if row['Stat_Test_Performed'] and not pd.isna(row['Mean_Background']):
            ax2.plot([row['Mean_Background'], row['Mean_Background']],
                    [i - 0.4, i + 0.4], 'k--', linewidth=2, alpha=0.5)

    # Add standard deviation error bars
    ax2.errorbar(stats_df['Mean_Conservation'], stats_df['MDP'],
                xerr=stats_df['Conservation_CV'], fmt='none',
                ecolor='black', elinewidth=1, capsize=3, alpha=0.8)

    ax2.set_xlabel('Mean Conservation (bits)', fontsize=12)
    ax2.grid(True, alpha=0.3, axis='x')

    # Add legend for statistical test
    from matplotlib.patches import Patch
    legend_elements = [
        plt.Line2D([0], [0], color='black', linestyle='--', label='Background Mean')
    ]
    ax2.legend(handles=legend_elements, fontsize=7, loc='lower right', framealpha=0.3, fancybox=True)

    # Add values on bars with significance indicator
    for i, (bar, row) in enumerate(zip(bars, stats_df.iterrows())):
        _, row_data = row
        text = f'{row_data["Mean_Conservation"]:.2f}'
        if row_data['Stat_Test_Performed'] and row_data['Stat_Test_Significant']:
            text += ' ***'
        ax2.text(2.4, bar.get_y() + bar.get_height()/2, text,
                va='center', fontsize=9)
    
    # 3. Conservation variability analysis
    ax3.scatter(stats_df['Mean_Conservation'], stats_df['Conservation_CV'], 
               s=stats_df['Length']*10, alpha=0.7,
               c=['blue' if gc == 'mitochondrial' else 'orange' 
                  for gc in stats_df['Genetic_Code']])
    
    for i, mdp in enumerate(stats_df['MDP']):
        ax3.annotate(mdp, (stats_df.iloc[i]['Mean_Conservation'], 
                          stats_df.iloc[i]['Conservation_CV']),
                    xytext=(5, 5), textcoords='offset points', fontsize=9)
    
    ax3.set_xlabel('Mean Conservation (bits)', fontsize=12)
    ax3.set_ylabel('Conservation CV', fontsize=12)
    ax3.grid(True, alpha=0.3)
    
    # 4. Hierarchical clustering of conservation patterns
    # Calculate distance matrix (only for positions that exist in all MDPs)
    min_length = min(len(results['conservation_scores']) for results in mdp_results_dict.values())
    truncated_matrix = []
    for mdp_name in mdp_names:
        scores = mdp_results_dict[mdp_name]['conservation_scores'][:min_length]
        truncated_matrix.append(scores)
    
    if len(truncated_matrix) > 1 and min_length > 0:
        # Perform hierarchical clustering
        distance_matrix = pdist(truncated_matrix, metric='euclidean')
        linkage_matrix = linkage(distance_matrix, method='ward')
        
        # Create dendrogram
        dendro = dendrogram(linkage_matrix, labels=mdp_names, ax=ax4, orientation='left')
        ax4.set_xlabel('Distance', fontsize=12)
    else:
        ax4.text(0.5, 0.5, 'Insufficient data\nfor clustering', 
                ha='center', va='center', transform=ax4.transAxes, fontsize=12)
        ax4.set_title('Conservation Pattern\nClustering', fontsize=12, fontweight='bold')
    
    # Add legend for genetic codes
    from matplotlib.patches import Patch
    legend_elements = [Patch(facecolor='blue', label='Mitochondrial'),
                      Patch(facecolor='orange', label='Nuclear')]
    #fig.legend(handles=legend_elements, loc='upper right', bbox_to_anchor=(0.98, 0.95),
    #          title='Genetic Code', title_fontsize=10, fontsize=10)
    

    plt.tight_layout()
    
    if output_file:
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        print(f"Heatmap summary saved to {output_file}")
    
    plt.show()
    
    return stats_df, conservation_df

def create_mdp_selection_summary(mdp_results_dict, dnds_results_dict=None, 
                                figsize=(18, 12), output_file=None):
    """
    Create a comprehensive selection analysis summary for MDPs combining conservation
    and dN/dS analysis results.
    
    Parameters:
    -----------
    mdp_results_dict : dict
        Dictionary with conservation analysis results
    dnds_results_dict : dict, optional
        Dictionary with dN/dS analysis results
    figsize : tuple
        Figure size
    output_file : str, optional
        Path to save the plot
    """
    
    # Create figure with complex layout
    fig = plt.figure(figsize=figsize)
    
    # Define grid layout
    gs = fig.add_gridspec(4, 4, height_ratios=[1, 1, 1, 0.3], width_ratios=[2, 1, 1, 1],
                         hspace=0.4, wspace=0.3)
    
    # Color schemes
    mito_color = '#2E86AB'  # Blue for mitochondrial
    nuclear_color = '#F24236'  # Red for nuclear
    
    # Prepare summary data
    summary_data = []
    for mdp_name, results in mdp_results_dict.items():
        data = {
            'MDP': mdp_name,
            'Length': results['sequence_length'],
            'Mean_Conservation': results['mean_conservation'],
            'Max_Conservation': results['max_conservation'],
            'Conservation_Range': results['max_conservation'] - min(results['conservation_scores']),
            'Conservation_CV': np.std(results['conservation_scores']) / np.mean(results['conservation_scores']),
            'N_Sequences': results['n_sequences'],
            'Genetic_Code': results.get('genetic_code', 'unknown'),
            'Most_Conserved_Pos': np.argmax(results['conservation_scores']) + 1
        }
        
        # Add dN/dS data if available
        if dnds_results_dict and mdp_name in dnds_results_dict:
            dnds_data = dnds_results_dict[mdp_name]
            data.update({
                'dN_dS_Ratio': dnds_data.get('dnds_ratio', np.nan),
                'dN_Rate': dnds_data.get('nonsynonymous_rate', np.nan),
                'dS_Rate': dnds_data.get('synonymous_rate', np.nan),
                'P_Value': dnds_data.get('p_value', np.nan),
                'Is_Significant': dnds_data.get('is_significant', False)
            })
        else:
            data.update({
                'dN_dS_Ratio': np.nan,
                'dN_Rate': np.nan,
                'dS_Rate': np.nan,
                'P_Value': np.nan,
                'Is_Significant': False
            })
        
        summary_data.append(data)
    
    df = pd.DataFrame(summary_data)
    
    # 1. Conservation profiles comparison (top left)
    ax1 = fig.add_subplot(gs[0, :3])
    colors = [mito_color if gc == 'mitochondrial' else nuclear_color 
              for gc in df['Genetic_Code']]
    
    for i, (mdp_name, results) in enumerate(mdp_results_dict.items()):
        positions = range(1, len(results['conservation_scores']) + 1)
        ax1.plot(positions, results['conservation_scores'], 
                label=f"{mdp_name} ({results.get('genetic_code', 'unknown')})",
                linewidth=2, alpha=0.8, color=colors[i])
        
        # Mark most conserved position
        max_pos = np.argmax(results['conservation_scores'])
        max_val = results['conservation_scores'][max_pos]
        ax1.scatter(max_pos + 1, max_val, color=colors[i], s=100, 
                   marker='*', edgecolor='black', linewidth=1, zorder=5)
    
    ax1.set_xlabel('Position', fontsize=12)
    ax1.set_ylabel('Conservation (bits)', fontsize=12)
    ax1.legend(bbox_to_anchor=(1.05, 1), loc='upper left', fontsize=10)
    ax1.grid(True, alpha=0.3)
    
    # 2. dN/dS ratios (top right)
    ax2 = fig.add_subplot(gs[0, 3])
    if not df['dN_dS_Ratio'].isna().all():
        valid_df = df.dropna(subset=['dN_dS_Ratio'])
        colors_dnds = [mito_color if gc == 'mitochondrial' else nuclear_color 
                      for gc in valid_df['Genetic_Code']]
        
        bars = ax2.barh(range(len(valid_df)), valid_df['dN_dS_Ratio'], 
                       color=colors_dnds, alpha=0.7)
        ax2.axvline(x=1, color='black', linestyle='--', alpha=0.7, label='Neutral')
        ax2.set_yticks(range(len(valid_df)))
        ax2.set_yticklabels(valid_df['MDP'], fontsize=10)
        ax2.set_xlabel('dN/dS Ratio', fontsize=12)
        ax2.set_title('Selection Pressure', fontsize=12, fontweight='bold')
        
        # Mark significant results
        for i, (idx, row) in enumerate(valid_df.iterrows()):
            if row['Is_Significant']:
                ax2.text(row['dN_dS_Ratio'] + 0.05, i, '*', 
                        fontsize=16, fontweight='bold', va='center')
        
        ax2.legend(fontsize=9)
    else:
        ax2.text(0.5, 0.5, 'No dN/dS data\navailable', ha='center', va='center',
                transform=ax2.transAxes, fontsize=12)
        ax2.set_title('Selection Pressure', fontsize=12, fontweight='bold')
    
    # 3. Conservation vs Selection correlation (middle left)
    ax3 = fig.add_subplot(gs[1, :2])
    if not df['dN_dS_Ratio'].isna().all():
        valid_df = df.dropna(subset=['dN_dS_Ratio'])
        colors_scatter = [mito_color if gc == 'mitochondrial' else nuclear_color 
                         for gc in valid_df['Genetic_Code']]
        
        scatter = ax3.scatter(valid_df['Mean_Conservation'], valid_df['dN_dS_Ratio'],
                             s=valid_df['Length']*15, c=colors_scatter, alpha=0.7,
                             edgecolors='black', linewidth=1)
        
        # Add correlation line if enough data points
        if len(valid_df) > 2:
            z = np.polyfit(valid_df['Mean_Conservation'], valid_df['dN_dS_Ratio'], 1)
            p = np.poly1d(z)
            ax3.plot(valid_df['Mean_Conservation'], p(valid_df['Mean_Conservation']),
                    "r--", alpha=0.8, linewidth=2)
            
            # Calculate correlation
            r, p_val = stats.pearsonr(valid_df['Mean_Conservation'], valid_df['dN_dS_Ratio'])
            ax3.text(0.05, 0.95, f'r = {r:.3f}\np = {p_val:.3f}', 
                    transform=ax3.transAxes, fontsize=10,
                    bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
        
        # Add labels
        for i, row in valid_df.iterrows():
            ax3.annotate(row['MDP'], (row['Mean_Conservation'], row['dN_dS_Ratio']),
                        xytext=(5, 5), textcoords='offset points', fontsize=9)
        
        ax3.axhline(y=1, color='black', linestyle='--', alpha=0.5)
        ax3.set_xlabel('Mean Conservation (bits)', fontsize=12)
        ax3.set_ylabel('dN/dS Ratio', fontsize=12)
        ax3.set_title('Conservation vs Selection Pressure\n(bubble size = length)', 
                     fontsize=12, fontweight='bold')
        ax3.grid(True, alpha=0.3)
    else:
        ax3.text(0.5, 0.5, 'No dN/dS data available\nfor correlation analysis', 
                ha='center', va='center', transform=ax3.transAxes, fontsize=12)
        ax3.set_title('Conservation vs Selection', fontsize=12, fontweight='bold')
    
    # 4. Conservation distribution by genetic code (middle right)
    ax4 = fig.add_subplot(gs[1, 2:])
    
    mito_conservation = df[df['Genetic_Code'] == 'mitochondrial']['Mean_Conservation']
    nuclear_conservation = df[df['Genetic_Code'] == 'nuclear']['Mean_Conservation']
    
    data_to_plot = []
    labels = []
    colors_box = []
    
    if len(mito_conservation) > 0:
        data_to_plot.append(mito_conservation)
        labels.append('Mitochondrial')
        colors_box.append(mito_color)
    
    if len(nuclear_conservation) > 0:
        data_to_plot.append(nuclear_conservation)
        labels.append('Nuclear')
        colors_box.append(nuclear_color)
    
    if data_to_plot:
        bp = ax4.boxplot(data_to_plot, labels=labels, patch_artist=True)
        for patch, color in zip(bp['boxes'], colors_box):
            patch.set_facecolor(color)
            patch.set_alpha(0.7)
        
        # Add statistical test if both groups present
        if len(data_to_plot) == 2 and len(mito_conservation) > 0 and len(nuclear_conservation) > 0:
            stat, p_val = stats.mannwhitneyu(mito_conservation, nuclear_conservation)
            ax4.text(0.5, 0.95, f'Mann-Whitney U\np = {p_val:.3f}', 
                    transform=ax4.transAxes, ha='center', fontsize=10,
                    bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
    
    ax4.set_ylabel('Mean Conservation (bits)', fontsize=12)
    ax4.set_title('Conservation by\nGenetic Code', fontsize=12, fontweight='bold')
    ax4.grid(True, alpha=0.3, axis='y')
    
    # 5. Summary statistics table (bottom left)
    ax5 = fig.add_subplot(gs[2, :2])
    ax5.axis('off')
    
    # Prepare table data
    table_data = []
    for _, row in df.iterrows():
        table_row = [
            row['MDP'],
            f"{row['Length']} AA",
            f"{row['Mean_Conservation']:.2f}",
            f"{row['Conservation_CV']:.3f}",
            row['Genetic_Code'][:4],  # Abbreviated
            f"{row['dN_dS_Ratio']:.2f}" if not pd.isna(row['dN_dS_Ratio']) else 'N/A',
            '*' if row['Is_Significant'] else ''
        ]
        table_data.append(table_row)
    
    headers = ['MDP', 'Length', 'Mean\nCons', 'Cons\nCV', 'Code', 'dN/dS', 'Sig']
    
    table = ax5.table(cellText=table_data, colLabels=headers,
                     cellLoc='center', loc='center',
                     colWidths=[0.15, 0.12, 0.12, 0.12, 0.12, 0.12, 0.08])
    table.auto_set_font_size(False)
    table.set_fontsize(9)
    table.scale(1, 1.5)
    
    # Color code genetic codes in table
    for i, row in enumerate(table_data):
        if row[4] == 'mito':  # mitochondrial
            for j in range(len(headers)):
                table[(i+1, j)].set_facecolor('#E3F2FD')  # Light blue
        else:  # nuclear
            for j in range(len(headers)):
                table[(i+1, j)].set_facecolor('#FFEBEE')  # Light red
    
    ax5.set_title('Summary Statistics', fontsize=12, fontweight='bold', pad=20)
    
    # 6. Functional classification plot (bottom right)
    ax6 = fig.add_subplot(gs[2, 2:])
    
    # Classify MDPs by conservation and selection
    classifications = []
    for _, row in df.iterrows():
        if pd.isna(row['dN_dS_Ratio']):
            if row['Mean_Conservation'] > df['Mean_Conservation'].median():
                classification = 'Highly Conserved'
            else:
                classification = 'Moderately Conserved'
        else:
            if row['dN_dS_Ratio'] > 1 and row['Is_Significant']:
                classification = 'Positive Selection'
            elif row['dN_dS_Ratio'] < 1 and row['Is_Significant']:
                classification = 'Negative Selection'
            elif row['Mean_Conservation'] > df['Mean_Conservation'].median():
                classification = 'Highly Conserved'
            else:
                classification = 'Moderately Conserved'
        classifications.append(classification)
    
    df['Classification'] = classifications
    
    # Count classifications
    class_counts = df['Classification'].value_counts()
    colors_pie = ['#FF6B6B', '#4ECDC4', '#45B7D1', '#96CEB4', '#FFEAA7']
    
    wedges, texts, autotexts = ax6.pie(class_counts.values, labels=class_counts.index,
                                      autopct='%1.0f%%', colors=colors_pie[:len(class_counts)],
                                      startangle=90)
    ax6.set_title('Functional Classification', fontsize=12, fontweight='bold')
    
    # 7. Legend and summary (bottom)
    ax7 = fig.add_subplot(gs[3, :])
    ax7.axis('off')
    
    # Create comprehensive legend
    legend_elements = [
        mpatches.Patch(color=mito_color, label='Mitochondrial Genetic Code'),
        mpatches.Patch(color=nuclear_color, label='Nuclear Genetic Code'),
        plt.Line2D([0], [0], marker='*', color='w', markerfacecolor='black', 
                  markersize=10, label='Most Conserved Position'),
        plt.Line2D([0], [0], color='black', linestyle='--', label='Neutral Selection (dN/dS=1)'),
        mpatches.Patch(color='white', label='* = Statistically Significant')
    ]
    
    ax7.legend(handles=legend_elements, loc='center', ncol=5, fontsize=11,
              frameon=True, fancybox=True, shadow=True)
    
    # Add overall summary text
    n_total = len(df)
    n_mito = len(df[df['Genetic_Code'] == 'mitochondrial'])
    n_nuclear = len(df[df['Genetic_Code'] == 'nuclear'])
    n_with_dnds = len(df.dropna(subset=['dN_dS_Ratio']))
    n_significant = len(df[df['Is_Significant'] == True])
    
    summary_text = f"""Dataset Summary: {n_total} MDPs analyzed ({n_mito} mitochondrial, {n_nuclear} nuclear genetic code)
    Selection Analysis: {n_with_dnds} MDPs with dN/dS data, {n_significant} showing significant selection pressure
    Mean Conservation Range: {df['Mean_Conservation'].min():.2f} - {df['Mean_Conservation'].max():.2f} bits"""
    
    fig.text(0.5, 0.02, summary_text, ha='center', fontsize=10,
             bbox=dict(boxstyle='round', facecolor='lightgray', alpha=0.8))
    
    # Overall title
    fig.suptitle('Mitochondrial-Derived Peptides: Comprehensive Selection and Conservation Analysis', 
                fontsize=16, fontweight='bold', y=0.98)
    
    if output_file:
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        print(f"Selection summary plot saved to {output_file}")
    
    plt.show()
    
    return df
# Example usage and demonstration
def example_usage():
    """Demonstrate how to use the logomaker-based sequence logo functions."""
    
    print("Micro-Protein Sequence Logo Creator - Logomaker Version")
    print("======================================================")
    print()
    
    if not LOGOMAKER_AVAILABLE:
        print("ERROR: logomaker package is required!")
        print("Install with: pip install logomaker")
        return
    
    print("Key Features:")
    print("=============")
    print("• Professional sequence logos using logomaker package")
    print("• Proper letter stacking and scaling")
    print("• Stop codons displayed as 'X'")
    print("• Invariant positions marked with red asterisks")
    print("• Multiple color schemes available")
    print("• Conservation plots above logos")
    print()
    
    print("Installation:")
    print("=============")
    print("pip install logomaker")
    print("pip install tqdm  # Optional: for progress bars during permutation tests")
    print()
    
    print("Example Usage:")
    print("==============")
    print("""
# Basic usage
results = create_microprotein_sequence_logo(
    sequences=your_sequences,
    genetic_code='nuclear',
    title='MOTS-c Conservation',
    color_scheme='chemistry',
    show_conservation=True,
    output_file='motsc_logo.png'
)

# Multiple logos - now returns both results dict and summary DataFrame
results_dict, summary_df = create_multiple_sequence_logos(
    sequences_dict=sequences_dict,
    genetic_codes=genetic_codes,
    output_dir='logos',
    color_scheme='chemistry'
)

# View the variant summary DataFrame
print(summary_df)

# Access variant information for a specific region
print(f"Total variants exceeding 1%: {summary_df.loc[0, 'total_variants_above_threshold']}")
print(f"Mean variants per position: {summary_df.loc[0, 'mean_variants_per_position']:.2f}")

# Compare conservation
compare_conservation_across_microproteins(
    results_dict=results_dict,
    output_file='conservation_comparison.png'
)
    """)
    
    return True

if __name__ == "__main__":
    example_usage()
