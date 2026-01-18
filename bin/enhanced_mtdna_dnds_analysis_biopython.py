import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from collections import Counter, defaultdict
import logging
from typing import List, Tuple, Dict, Optional, Union
from scipy import stats
from scipy.stats import false_discovery_control
import warnings
from Bio.Seq import reverse_complement
from matplotlib.legend import Legend
import clean_amb_and_gaps as cag
from importlib import reload

reload(cag)

# NEW: Biopython imports for dN/dS calculation
from Bio.codonalign.codonseq import cal_dn_ds, CodonSeq
from Bio.Seq import Seq
from Bio.Data.CodonTable import generic_by_id


# Set up logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

# Define ambiguous nucleotide codes
AMBIGUOUS_CODES = set('RYKMSWBDHVN')

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

# Start codons for different genetic codes
MITOCHONDRIAL_START_CODONS = {'ATG', 'ATA', 'ATT'}  # ATT sometimes reported but less established
NUCLEAR_START_CODONS = {'ATG'}  # Nuclear genes typically use only ATG

# Stop codons for different genetic codes
MITOCHONDRIAL_STOP_CODONS = {'TAA', 'TAG'}
NUCLEAR_STOP_CODONS = {'TAA', 'TAG', 'TGA'}

def has_ambiguous_bases(codon):
    """Check if codon contains any ambiguous bases"""
    return any(base.upper() in AMBIGUOUS_CODES for base in codon)

def translate_sequence(sequence: str, genetic_code: str = 'mitochondrial', allow_att_start = False) -> str:
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
            if i == 0 and allow_att_start:
                if codon == 'ATT':
                    amino_acid = 'M'
            protein.append(amino_acid)
        else:
            break
    return ''.join(protein)

def validate_coding_region(sequence: str, strand: str = '+', genetic_code: str = 'mitochondrial', allow_alternative_starts: bool = True) -> Dict:
    """
    Validate that a sequence represents a valid coding region.
    
    Parameters:
    -----------
    sequence : str
        DNA sequence to validate
    strand : str
        '+' for forward strand, '-' for reverse strand
    genetic_code : str
        'mitochondrial' or 'nuclear' genetic code to use
    allow_alternative_starts : bool
        Whether to allow alternative start codons (currently ATT for mitochondrial, though less well-established)
    
    Returns:
    --------
    Dict with validation results
    """
    validation_results = {
        'is_valid': False,
        'length_divisible_by_3': False,
        'has_start_codon': False,
        'has_premature_stop': False,
        'start_codon': None,
        'stop_codons_found': [],
        'reading_frame': None,
        'genetic_code': genetic_code,
        'issues': []
    }
    
    # Select appropriate codon sets
    if genetic_code.lower() == 'mitochondrial':
        start_codons = MITOCHONDRIAL_START_CODONS
        stop_codons = MITOCHONDRIAL_STOP_CODONS
        alternative_starts = {'ATT'} if allow_alternative_starts else set()
    elif genetic_code.lower() == 'nuclear':
        start_codons = NUCLEAR_START_CODONS
        stop_codons = NUCLEAR_STOP_CODONS
        alternative_starts = set()  # Nuclear code typically only uses ATG
    else:
        raise ValueError("genetic_code must be 'mitochondrial' or 'nuclear'")
    
    # Strand handling is already done in main function, so we assume sequence is in correct orientation
    sequence = sequence.upper()
    # Check if length is divisible by 3
    if len(sequence) % 3 == 0:
        validation_results['length_divisible_by_3'] = True
    else:
        validation_results['issues'].append(f"Length {len(sequence)} not divisible by 3")
        # Try different reading frames
        for frame in [1, 2]:
            if (len(sequence) - frame) % 3 == 0:
                validation_results['reading_frame'] = frame
                sequence = sequence[frame:]
                validation_results['issues'].append(f"Using reading frame +{frame}")
                break
    
    if len(sequence) < 3:
        validation_results['issues'].append("Sequence too short")
        return validation_results
    
    # Check start codon
    start_codon = sequence[:3]
    if start_codon in start_codons:
        validation_results['has_start_codon'] = True
        validation_results['start_codon'] = start_codon
    elif start_codon in alternative_starts:
        validation_results['has_start_codon'] = True
        validation_results['start_codon'] = start_codon
        validation_results['issues'].append(f"Using potentially alternative start codon: {start_codon}")
    else:
        validation_results['issues'].append(f"Invalid start codon: {start_codon} (expected one of {start_codons})")
    
    # Check for premature stop codons (excluding the last codon)
    stop_positions = []
    for i in range(0, len(sequence) - 3, 3):
        codon = sequence[i:i+3]
        if codon in stop_codons:
            stop_positions.append(i // 3)
            validation_results['stop_codons_found'].append((i // 3, codon))
    
    if stop_positions:
        validation_results['has_premature_stop'] = True
        validation_results['issues'].append(f"Premature stop codons at positions: {stop_positions}")
    
    # Overall validation
    validation_results['is_valid'] = (
        validation_results['length_divisible_by_3'] and
        validation_results['has_start_codon'] and
        not validation_results['has_premature_stop']
    )
    
    return validation_results

def calculate_synonymous_nonsynonymous(seq1: str, seq2: str, strand: str = '+',
                                      genetic_code: str = 'mitochondrial',
                                      method: str = 'NG86') -> Dict:
    """
    Calculate synonymous and non-synonymous differences between two sequences using Biopython.

    MODIFIED: Now uses Biopython's cal_dn_ds for dN/dS calculation instead of custom implementation.
    Tracks fallback statistics (accessible via get_fallback_stats()).

    Parameters:
    -----------
    seq1, seq2 : str
        DNA sequences to compare
    strand : str
        '+' for forward strand, '-' for reverse strand
    genetic_code : str
        'mitochondrial' or 'nuclear' genetic code to use
    method : str
        Biopython method for dN/dS calculation: 'NG86' (default), 'LWL85', 'ML', 'YN00'

    Returns:
    --------
    Dict with counts of synonymous and non-synonymous differences
    """
    # Initialize counters as function attributes if they don't exist
    if not hasattr(calculate_synonymous_nonsynonymous, 'total_calls'):
        calculate_synonymous_nonsynonymous.total_calls = 0
        calculate_synonymous_nonsynonymous.fallback_calls = 0
        calculate_synonymous_nonsynonymous.fallback_reasons = {
            'internal_stops': 0,
            'too_short_after_strip': 0,
            'too_short_min_length': 0,
            'negative_rates': 0,
            'exception': 0
        }

    # Increment total calls
    calculate_synonymous_nonsynonymous.total_calls += 1

    # Map genetic code to NCBI codon table ID
    if genetic_code.lower() == 'mitochondrial':
        codon_table_id = 2  # Vertebrate mitochondrial code
        codon_table = MITOCHONDRIAL_GENETIC_CODE
        stop_codons = MITOCHONDRIAL_STOP_CODONS
    elif genetic_code.lower() == 'nuclear':
        codon_table_id = 1  # Standard genetic code
        codon_table = NUCLEAR_GENETIC_CODE
        stop_codons = NUCLEAR_STOP_CODONS
    else:
        raise ValueError("genetic_code must be 'mitochondrial' or 'nuclear'")
    
    # Prepare sequences
    seq1, seq2 = seq1.upper(), seq2.upper()
    
    # Ensure sequences are same length and divisible by 3
    min_len = min(len(seq1), len(seq2))
    min_len = min_len - (min_len % 3)
    seq1, seq2 = seq1[:min_len], seq2[:min_len]
    
    if min_len < 3:
        return {
            'synonymous_sites': 0,
            'nonsynonymous_sites': 0,
            'synonymous_differences': 0,
            'nonsynonymous_differences': 0,
            'total_sites': 0,
            'total_differences': 0
        }
    
    # Strip terminal stop codons for Biopython (it doesnt handle terminal stops well)
    seq1_for_bp = seq1
    seq2_for_bp = seq2
    stripped_stop = False

    # Check if sequence ends with a stop codon
    if len(seq1) >= 3 and len(seq2) >= 3:
        last_codon1 = seq1[-3:]
        last_codon2 = seq2[-3:]
        seq1_for_bp = seq1[:-3]
        seq2_for_bp = seq2[:-3]
        stripped_stop = True
        #logger.info(f"Stripped terminal stop codons for Biopython: {last_codon1}, {last_codon2}")

    # Count synonymous and nonsynonymous sites (keep existing site counting logic)
    synonymous_sites = 0
    nonsynonymous_sites = 0
    
    for i in range(0, min_len, 3):
        codon = seq1[i:i+3]
        if len(codon) == 3:# and not has_ambiguous_bases(codon) and '-' not in codon:
            syn_sites = count_synonymous_sites(codon, codon_table)
            synonymous_sites += syn_sites
            nonsynonymous_sites += (3 - syn_sites)
    
    # Determine if we can use Biopython based on sequence quality checks
    use_biopython = True

    # Check for internal stop codons
    has_internal_stops = False
    for i in range(0, len(seq1_for_bp) - 3, 3):
        if seq1_for_bp[i:i+3] in stop_codons or seq2_for_bp[i:i+3] in stop_codons:
            has_internal_stops = True
            break

    if has_internal_stops:
        logger.warning("Internal stop codons detected, cannot use Biopython dN/dS calculation. Falling back to custom implementation.")
        use_biopython = False
        calculate_synonymous_nonsynonymous.fallback_reasons['internal_stops'] += 1

    if use_biopython and (len(seq1_for_bp) < 3 or len(seq2_for_bp) < 3):
        logger.warning("Sequences too short after stripping terminal stop codons. Falling back to custom implementation.")
        use_biopython = False
        calculate_synonymous_nonsynonymous.fallback_reasons['too_short_after_strip'] += 1

    # Check minimum sequence length for reliable Biopython calculation
    # Require at least 12 nucleotides (4 codons) for stable dN/dS estimation
    if use_biopython and (len(seq1_for_bp) < 12 or len(seq2_for_bp) < 12):
        logger.warning(f"Sequences too short for reliable Biopython calculation (seq1: {len(seq1_for_bp)}bp, seq2: {len(seq2_for_bp)}bp). Minimum 12bp required. Falling back to custom implementation.")
        use_biopython = False
        calculate_synonymous_nonsynonymous.fallback_reasons['too_short_min_length'] += 1

    # Use Biopython to calculate dN and dS (on sequences without terminal stop codons)
    if use_biopython:
        try:
            # Create CodonSeq objects
            codon_seq1 = CodonSeq(Seq(seq1_for_bp))
            codon_seq2 = CodonSeq(Seq(seq2_for_bp))

            # Iterate codon by codon, and remove codons with ambiguous bases or gaps
            filtered_codons1 = []
            filtered_codons2 = []
            for i in range(0, len(seq1_for_bp), 3):
                codon1 = seq1_for_bp[i:i+3]
                codon2 = seq2_for_bp[i:i+3]
                if len(codon1) != 3 or len(codon2) != 3:
                    continue
                # Remove codon if it contains anything other than A,T,G,C
                #if has_ambiguous_bases(codon1) or has_ambiguous_bases(codon2):
                    continue
                #if '-' in codon1 or '-' in codon2:
                    continue
                filtered_codons1.append(codon1)
                filtered_codons2.append(codon2)
            codon_seq1 = CodonSeq(Seq(''.join(filtered_codons1)))
            codon_seq2 = CodonSeq(Seq(''.join(filtered_codons2)))
            #logger.debug(f"Filtered codons for Biopython: {len(filtered_codons1)} codons remaining after removing ambiguous/gap codons.")
            # Get codon table from Biopython
            bp_codon_table = generic_by_id[codon_table_id]


            # Calculate dN and dS using Biopython
            dn, ds = cal_dn_ds(codon_seq1, codon_seq2, method=method, codon_table=bp_codon_table)

            # Check for problematic dN/dS values (negative or non-finite)
            if ds < 0 or dn < 0 or not np.isfinite(ds) or not np.isfinite(dn):
                logger.warning(f"Biopython returned problematic values (dN={dn:.6f}, dS={ds:.6f}). Negative or non-finite rates detected. Falling back to custom implementation.")
                use_biopython = False
                calculate_synonymous_nonsynonymous.fallback_reasons['negative_rates'] += 1
            else:
                # Calculate differences from Biopython-derived rates and site counts
                # differences = rate × sites
                synonymous_differences = ds * synonymous_sites if synonymous_sites > 0 else 0
                nonsynonymous_differences = dn * nonsynonymous_sites if nonsynonymous_sites > 0 else 0

                logger.debug(f"Biopython {method}: dN={dn:.6f}, dS={ds:.6f}")

        except Exception as e:
            logger.warning(f"Biopython calculation failed: {e}. Falling back to custom implementation.")
            use_biopython = False
            calculate_synonymous_nonsynonymous.fallback_reasons['exception'] += 1

    # Fall back to custom implementation if Biopython was not used or failed
    if not use_biopython:
        calculate_synonymous_nonsynonymous.fallback_calls += 1
        synonymous_differences = 0
        nonsynonymous_differences = 0
        
        # Compare codon by codon (original logic)
        for i in range(0, min_len, 3):
            codon1 = seq1[i:i+3]
            codon2 = seq2[i:i+3]
            
            if len(codon1) != 3 or len(codon2) != 3:
                continue
                
            # Skip if either codon contains ambiguous bases
            if has_ambiguous_bases(codon1) or has_ambiguous_bases(codon2):
                continue
                
            # Skip if either codon contains gaps '-'
            if '-' in codon1 or '-' in codon2:
                continue

            aa1 = codon_table.get(codon1, 'X')
            aa2 = codon_table.get(codon2, 'X')
            
            # Skip stop codons
            if aa1 == '*' or aa2 == '*':
                continue

            # Count differences between codons
            codon_differences = sum(1 for a, b in zip(codon1, codon2) if a != b)
            
            if codon_differences == 1:
                # Single nucleotide difference
                if aa1 == aa2:
                    synonymous_differences += 1
                else:
                    nonsynonymous_differences += 1
            elif codon_differences > 1:
                # Multiple differences - use Nei-Gojobori method approximation
                syn_count, nonsyn_count = estimate_multiple_changes(codon1, codon2, codon_table)
                synonymous_differences += syn_count
                nonsynonymous_differences += nonsyn_count
    
    return {
        'synonymous_sites': synonymous_sites,
        'nonsynonymous_sites': nonsynonymous_sites,
        'synonymous_differences': synonymous_differences,
        'nonsynonymous_differences': nonsynonymous_differences,
        'total_sites': synonymous_sites + nonsynonymous_sites,
        'total_differences': synonymous_differences + nonsynonymous_differences
    }

def get_fallback_stats() -> Dict:
    """
    Get statistics on how often calculate_synonymous_nonsynonymous fell back to custom implementation.

    Returns:
    --------
    Dict with fallback statistics including:
        - total_calls: Total number of times the function was called
        - fallback_calls: Number of times custom implementation was used
        - biopython_calls: Number of times Biopython was successfully used
        - fallback_percentage: Percentage of calls that used fallback
        - fallback_reasons: Breakdown of reasons for fallback
    """
    if not hasattr(calculate_synonymous_nonsynonymous, 'total_calls'):
        return {
            'total_calls': 0,
            'fallback_calls': 0,
            'biopython_calls': 0,
            'fallback_percentage': 0.0,
            'fallback_reasons': {}
        }

    total = calculate_synonymous_nonsynonymous.total_calls
    fallback = calculate_synonymous_nonsynonymous.fallback_calls
    biopython = total - fallback
    fallback_pct = (fallback / total * 100) if total > 0 else 0.0

    return {
        'total_calls': total,
        'fallback_calls': fallback,
        'biopython_calls': biopython,
        'fallback_percentage': fallback_pct,
        'fallback_reasons': calculate_synonymous_nonsynonymous.fallback_reasons.copy()
    }

def reset_fallback_stats():
    """Reset fallback statistics counters."""
    if hasattr(calculate_synonymous_nonsynonymous, 'total_calls'):
        calculate_synonymous_nonsynonymous.total_calls = 0
        calculate_synonymous_nonsynonymous.fallback_calls = 0
        calculate_synonymous_nonsynonymous.fallback_reasons = {
            'internal_stops': 0,
            'too_short_after_strip': 0,
            'too_short_min_length': 0,
            'negative_rates': 0,
            'exception': 0
        }

def log_fallback_stats():
    """Log a summary of fallback statistics."""
    stats = get_fallback_stats()

    if stats['total_calls'] == 0:
        logger.info("No calls to calculate_synonymous_nonsynonymous yet.")
        return

    logger.info("=" * 60)
    logger.info("Biopython dN/dS Calculation Statistics")
    logger.info("=" * 60)
    logger.info(f"Total function calls: {stats['total_calls']}")
    logger.info(f"Biopython used successfully: {stats['biopython_calls']} ({100-stats['fallback_percentage']:.1f}%)")
    logger.info(f"Fallback to custom implementation: {stats['fallback_calls']} ({stats['fallback_percentage']:.1f}%)")

    if stats['fallback_calls'] > 0:
        logger.info("\nFallback reasons breakdown:")
        for reason, count in stats['fallback_reasons'].items():
            if count > 0:
                pct = (count / stats['fallback_calls'] * 100)
                logger.info(f"  - {reason.replace('_', ' ').title()}: {count} ({pct:.1f}% of fallbacks)")
    logger.info("=" * 60)

def count_synonymous_sites(codon: str, codon_table: Dict) -> float:
    """Count the number of synonymous sites in a codon."""
    if codon not in codon_table:
        return 0
    
    original_aa = codon_table[codon]
    synonymous_count = 0
    
    bases = ['A', 'T', 'G', 'C']
    
    # Check each position
    for pos in range(3):
        for base in bases:
            if base != codon[pos]:
                new_codon = codon[:pos] + base + codon[pos+1:]
                if new_codon in codon_table:
                    new_aa = codon_table[new_codon]
                    if new_aa == original_aa:
                        synonymous_count += 1
    
    return synonymous_count / 9  # Normalize by possible changes (3 positions × 3 alternatives)

def count_nonsynonymous_sites(codon: str, codon_table: Dict) -> float:
    """Count the number of non-synonymous sites in a codon."""
    return 3 - count_synonymous_sites(codon, codon_table)

def estimate_multiple_changes(codon1: str, codon2: str, codon_table: Dict) -> Tuple[float, float]:
    """Estimate synonymous and non-synonymous changes for multiple differences."""
    # Simple approximation - in practice, you might want to use more sophisticated methods
    # like the Nei-Gojobori or Yang-Nielsen methods
    
    differences = sum(1 for a, b in zip(codon1, codon2) if a != b)
    
    if differences == 0:
        return 0, 0
    
    # Estimate based on the proportion of synonymous sites
    syn_sites1 = count_synonymous_sites(codon1, codon_table)
    syn_sites2 = count_synonymous_sites(codon2, codon_table)
    avg_syn_sites = (syn_sites1 + syn_sites2) / 2
    
    # Rough approximation
    expected_syn = differences * (avg_syn_sites / 3)
    expected_nonsyn = differences - expected_syn
    
    return expected_syn, expected_nonsyn

def find_best_match_coordinates(reference_sequence: str, expected_sequence: str, 
                               original_start: int, original_end: int, 
                               search_window: int = 50, min_similarity: float = 0.8) -> Tuple[int, int, float, str]:
    """
    Search for the best match of an expected sequence within a window around original coordinates.
    
    Parameters:
    -----------
    reference_sequence : str
        The full reference sequence to search in
    expected_sequence : str
        The expected sequence to find
    original_start : int
        Original start position (1-based, absolute coordinates)
    original_end : int
        Original end position (1-based, absolute coordinates)
    search_window : int
        Number of nucleotides to search on each side of original coordinates
    min_similarity : float
        Minimum similarity threshold to consider a match
    
    Returns:
    --------
    Tuple of (best_start, best_end, best_similarity, best_match_sequence)
    """
    def calculate_similarity(seq1: str, seq2: str) -> float:
        """Calculate sequence similarity (proportion of matching nucleotides)."""
        if len(seq1) != len(seq2):
            return 0.0
        matches = sum(1 for a, b in zip(seq1.upper(), seq2.upper()) if a == b)
        return matches / len(seq1)
    
    expected_sequence = expected_sequence.upper()
    reference_sequence = reference_sequence.upper()
    expected_length = len(expected_sequence)
    
    # Define search boundaries (convert to 0-based for searching)
    search_start = max(0, original_start - 1 - search_window)
    search_end = min(len(reference_sequence), original_end + search_window)
    
    best_similarity = 0
    best_start = original_start
    best_end = original_end
    best_match = reference_sequence[original_start:original_end] if original_end <= len(reference_sequence) else ""
    
    # Search within the window
    for start_pos in range(search_start, search_end - expected_length + 1):
        end_pos = start_pos + expected_length
        if end_pos > len(reference_sequence):
            continue
            
        candidate_sequence = reference_sequence[start_pos:end_pos]
        similarity = calculate_similarity(expected_sequence, candidate_sequence)
        
        if similarity > best_similarity:
            best_similarity = similarity
            best_start = start_pos + 1  # Convert back to 1-based
            best_end = end_pos
            best_match = candidate_sequence
    
    logger.info(f"Coordinate search: original ({original_start}-{original_end}), "
                f"best match ({best_start}-{best_end}), similarity: {best_similarity:.3f}")
    
    return best_start, best_end, best_similarity, best_match

def calculate_dnds_statistical_tests(synonymous_diffs: float, nonsynonymous_diffs: float, 
                                    synonymous_sites: float, nonsynonymous_sites: float,
                                    test_method: str = 'fisher') -> Dict[str, float]:
    """
    Calculate statistical tests for dN/dS = 1 (neutrality test).
    
    Provides multiple test options, with Fisher's exact test being more appropriate
    for small counts (common in microprotein analysis).
    
    Parameters:
    -----------
    synonymous_diffs : float
        Number of synonymous differences
    nonsynonymous_diffs : float
        Number of nonsynonymous differences
    synonymous_sites : float
        Number of synonymous sites
    nonsynonymous_sites : float
        Number of nonsynonymous sites
    test_method : str
        Statistical test method: 'fisher', 'z_test', 'z_test_robust', or 'all'
    
    Returns:
    --------
    Dict with test statistics and p-values
    """
    
    # Calculate rates
    syn_rate = synonymous_diffs / synonymous_sites if synonymous_sites > 0 else 0
    nonsyn_rate = nonsynonymous_diffs / nonsynonymous_sites if nonsynonymous_sites > 0 else 0
    
    results = {
        'synonymous_rate': syn_rate,
        'nonsynonymous_rate': nonsyn_rate,
        'dnds_ratio': nonsyn_rate / syn_rate if syn_rate > 0 else float('inf')
    }
    
    # Handle edge cases
    if syn_rate == 0 or nonsyn_rate == 0:
        results.update({
            'test_method': test_method,
            'statistic': 0.0,
            'p_value': 1.0,
            'is_significant': False,
            'interpretation': 'No variation - cannot test neutrality'
        })
        return results
    
    if test_method == 'fisher' or test_method == 'all':
        # Fisher's exact test: compare observed syn/nonsyn ratio to expected under neutrality
        # Under neutrality, we expect dN/dS = 1, so dN*syn_sites = dS*nonsyn_sites
        # We can test if the ratio of differences deviates from the ratio of sites
        
        # Expected synonymous differences under neutrality
        total_diffs = synonymous_diffs + nonsynonymous_diffs
        total_sites = synonymous_sites + nonsynonymous_sites
        
        if total_sites > 0 and total_diffs > 0:
            expected_syn_diffs = total_diffs * (synonymous_sites / total_sites)
            
            # Use Fisher's exact test on the 2x2 contingency table:
            # [observed_syn, observed_nonsyn] vs [expected_syn, expected_nonsyn]
            try:
                _, p_fisher = stats.fisher_exact([
                    [int(synonymous_diffs), int(nonsynonymous_diffs)],
                    [int(expected_syn_diffs), int(total_diffs - expected_syn_diffs)]
                ])
                print([
                    [int(synonymous_diffs), int(nonsynonymous_diffs)],
                    [int(expected_syn_diffs), int(total_diffs - expected_syn_diffs)]
                ])
                results.update({
                    'fisher_p_value': p_fisher,
                    'fisher_is_significant': p_fisher < 0.05
                })

            except:
                results.update({
                    'fisher_p_value': 1.0,
                    'fisher_is_significant': False
                })
    
    if test_method == 'z_test' or test_method == 'all':
        # Original Z-test (problematic for small counts)
        dnds = results['dnds_ratio']
        if dnds != float('inf') and dnds > 0:
            ln_dnds = np.log(dnds)
            if synonymous_diffs > 0 and nonsynonymous_diffs > 0:
                se_ln_dnds = np.sqrt(1/synonymous_diffs + 1/nonsynonymous_diffs)
                z_stat = ln_dnds / se_ln_dnds
                p_z = 2 * (1 - stats.norm.cdf(abs(z_stat)))
            else:
                z_stat, p_z = 0.0, 1.0
        else:
            z_stat, p_z = 0.0, 1.0
            
        results.update({
            'z_statistic': z_stat,
            'z_p_value': p_z,
            'z_is_significant': p_z < 0.05
        })
    
    if test_method == 'z_test_robust' or test_method == 'all':
        # Robust Z-test using sites instead of differences for variance
        dnds = results['dnds_ratio']
        if dnds != float('inf') and dnds > 0:
            ln_dnds = np.log(dnds)
            # Use total sites for more stable variance estimate
            if synonymous_sites > 0 and nonsynonymous_sites > 0:
                # Add small constant to avoid division by zero and reduce extreme SEs
                pseudo_count = 0.5
                se_ln_dnds_robust = np.sqrt(
                    (1/(synonymous_sites + pseudo_count)) + 
                    (1/(nonsynonymous_sites + pseudo_count))
                )
                z_stat_robust = ln_dnds / se_ln_dnds_robust
                p_z_robust = 2 * (1 - stats.norm.cdf(abs(z_stat_robust)))
            else:
                z_stat_robust, p_z_robust = 0.0, 1.0
        else:
            z_stat_robust, p_z_robust = 0.0, 1.0
            
        results.update({
            'z_robust_statistic': z_stat_robust,
            'z_robust_p_value': p_z_robust,
            'z_robust_is_significant': p_z_robust < 0.05
        })
    
    # Set primary test results based on method
    if test_method == 'fisher':
        results.update({
            'test_method': 'Fisher exact test',
            'statistic': results.get('fisher_p_value', 1.0),  # For Fisher's test, report p-value as statistic
            'p_value': results.get('fisher_p_value', 1.0),
            'is_significant': results.get('fisher_is_significant', False)
        })
    elif test_method == 'z_test_robust':
        results.update({
            'test_method': 'Robust Z-test',
            'statistic': results.get('z_robust_statistic', 0.0),
            'p_value': results.get('z_robust_p_value', 1.0),
            'is_significant': results.get('z_robust_is_significant', False)
        })
    elif test_method == 'z_test':
        results.update({
            'test_method': 'Z-test',
            'statistic': results.get('z_statistic', 0.0),
            'p_value': results.get('z_p_value', 1.0),
            'is_significant': results.get('z_is_significant', False)
        })
    else:  # 'all'
        # Use Fisher's test as primary, but include all results
        results.update({
            'test_method': 'Multiple tests',
            'statistic': results.get('fisher_p_value', 1.0),
            'p_value': results.get('fisher_p_value', 1.0),
            'is_significant': results.get('fisher_is_significant', False)
        })
    
    return results

def analyze_protein_dnds(sequences: List[str], reference_seq: str, name: str,
                        genetic_code: str, strand: str, min_frequency: float,
                        statistical_test: str = 'fisher', allow_att_start = False,
                        dnds_method: str = 'NG86', weight_by_frequency: bool = True) -> Dict:
    """
    Analyze dN/dS for a protein region.
    
    Parameters:
    -----------
    sequences : List[str]
        List of sequences for the protein region
    reference_seq : str
        Reference sequence for the protein
    name : str
        Name of the protein
    genetic_code : str
        Genetic code to use
    strand : str
        Strand information
    min_frequency : float
        Minimum frequency for variants
    statistical_test : str
        Statistical test method: 'fisher', 'z_test', 'z_test_robust', or 'all'
    allow_att_start : bool
        Allow ATT as start codon
    dnds_method : str
        Biopython method for dN/dS calculation: 'NG86', 'LWL85', 'ML', 'YN00'
    weight_by_frequency : bool
        If True, weight differences and sites by variant frequency (default: True).
        If False, treat all variants equally regardless of frequency.
        Use False for population-level data where rare mutations should be considered equally.

    Returns:
    --------
    Dict with analysis results
    """
    
    # Calculate dN/dS for each variant vs reference
    # Store both weighted and unweighted for comparison
    synonymous_diffs_weighted = []
    nonsynonymous_diffs_weighted = []
    synonymous_sites_weighted = []
    nonsynonymous_sites_weighted = []

    synonymous_diffs_unweighted = []
    nonsynonymous_diffs_unweighted = []
    synonymous_sites_unweighted = []
    nonsynonymous_sites_unweighted = []

    variant_counts = Counter(sequences)
    total_sequences = len(sequences)

    # Track variants above frequency threshold and variants with mutations
    variants_above_frequency = 0
    variants_with_synonymous = 0
    variants_with_nonsynonymous = 0
    all_variants_with_synonymous = 0
    all_variants_with_nonsynonymous = 0

    for variant_seq, count in variant_counts.items():
        frequency = count / total_sequences

        # Count all variants with mutations (regardless of frequency threshold)
        if variant_seq != reference_seq:
            diff_results_all = calculate_synonymous_nonsynonymous(
                reference_seq, variant_seq, strand, genetic_code, method=dnds_method
            )
            if diff_results_all['synonymous_differences'] > 0:
                all_variants_with_synonymous += 1
            if diff_results_all['nonsynonymous_differences'] > 0:
                all_variants_with_nonsynonymous += 1

        if frequency < min_frequency:
            # Report rare variant skips
            print(f"Skipping rare variant with frequency {frequency:.4f}")
            continue  # Skip rare variants

        # Count variants above frequency threshold (excluding reference)
        if variant_seq != reference_seq:
            variants_above_frequency += 1

        if variant_seq == reference_seq:
            continue  # Skip reference sequence

        # Calculate differences using Biopython (use already calculated if same variant)
        diff_results = calculate_synonymous_nonsynonymous(
            reference_seq, variant_seq, strand, genetic_code, method=dnds_method
        )

        # Count variants with observed mutations (above frequency threshold)
        if diff_results['synonymous_differences'] > 0:
            variants_with_synonymous += 1
        if diff_results['nonsynonymous_differences'] > 0:
            variants_with_nonsynonymous += 1

        # Always calculate both weighted and unweighted for comparison
        # Weighted version
        synonymous_diffs_weighted.append(diff_results['synonymous_differences'] * frequency)
        nonsynonymous_diffs_weighted.append(diff_results['nonsynonymous_differences'] * frequency)
        synonymous_sites_weighted.append(diff_results['synonymous_sites'] * frequency)
        nonsynonymous_sites_weighted.append(diff_results['nonsynonymous_sites'] * frequency)

        # Unweighted version
        synonymous_diffs_unweighted.append(diff_results['synonymous_differences'])
        nonsynonymous_diffs_unweighted.append(diff_results['nonsynonymous_differences'])
        synonymous_sites_unweighted.append(diff_results['synonymous_sites'])
        nonsynonymous_sites_unweighted.append(diff_results['nonsynonymous_sites'])

    # Select which version to use for main results
    if weight_by_frequency:
        synonymous_diffs = synonymous_diffs_weighted
        nonsynonymous_diffs = nonsynonymous_diffs_weighted
        synonymous_sites = synonymous_sites_weighted
        nonsynonymous_sites = nonsynonymous_sites_weighted
    else:
        synonymous_diffs = synonymous_diffs_unweighted
        nonsynonymous_diffs = nonsynonymous_diffs_unweighted
        synonymous_sites = synonymous_sites_unweighted
        nonsynonymous_sites = nonsynonymous_sites_unweighted
    
    # Calculate summary statistics
    total_syn_diffs = sum(synonymous_diffs)
    total_nonsyn_diffs = sum(nonsynonymous_diffs)
    total_syn_sites = sum(synonymous_sites) if synonymous_sites else 1  # Avoid division by zero
    total_nonsyn_sites = sum(nonsynonymous_sites) if nonsynonymous_sites else 1
    
    # Calculate rates
    syn_rate = total_syn_diffs / total_syn_sites if total_syn_sites > 0 else 0
    nonsyn_rate = total_nonsyn_diffs / total_nonsyn_sites if total_nonsyn_sites > 0 else 0
    
    # Calculate dN/dS ratio
    dnds_ratio = nonsyn_rate / syn_rate if syn_rate > 0 else float('inf')
    
    # Calculate statistical tests for neutrality
    test_results = calculate_dnds_statistical_tests(
        total_syn_diffs, total_nonsyn_diffs, 
        total_syn_sites, total_nonsyn_sites, 
        statistical_test
    )
    
    # Translate reference sequence to amino acids
    amino_acid_sequence = translate_sequence(reference_seq, genetic_code, allow_att_start = allow_att_start)

    # Periodic logging of fallback statistics at milestones
    if hasattr(calculate_synonymous_nonsynonymous, 'total_calls'):
        total = calculate_synonymous_nonsynonymous.total_calls
        # Log at these milestones: 100, 500, 1000, 5000, 10000, etc.
        milestones = [100, 500, 1000, 5000, 10000, 50000, 100000]
        if total in milestones:
            log_fallback_stats()

    # Calculate the alternative weighting method for comparison
    if weight_by_frequency:
        # Main results are weighted, so calculate unweighted for comparison
        alt_syn_diffs = synonymous_diffs_unweighted
        alt_nonsyn_diffs = nonsynonymous_diffs_unweighted
        alt_syn_sites = synonymous_sites_unweighted
        alt_nonsyn_sites = nonsynonymous_sites_unweighted
        alt_method = 'unweighted'
    else:
        # Main results are unweighted, so calculate weighted for comparison
        alt_syn_diffs = synonymous_diffs_weighted
        alt_nonsyn_diffs = nonsynonymous_diffs_weighted
        alt_syn_sites = synonymous_sites_weighted
        alt_nonsyn_sites = nonsynonymous_sites_weighted
        alt_method = 'weighted'

    # Calculate alternative results
    alt_total_syn_diffs = sum(alt_syn_diffs)
    alt_total_nonsyn_diffs = sum(alt_nonsyn_diffs)
    alt_total_syn_sites = sum(alt_syn_sites) if alt_syn_sites else 1
    alt_total_nonsyn_sites = sum(alt_nonsyn_sites) if alt_nonsyn_sites else 1

    alt_syn_rate = alt_total_syn_diffs / alt_total_syn_sites if alt_total_syn_sites > 0 else 0
    alt_nonsyn_rate = alt_total_nonsyn_diffs / alt_total_nonsyn_sites if alt_total_nonsyn_sites > 0 else 0
    alt_dnds_ratio = alt_nonsyn_rate / alt_syn_rate if alt_syn_rate > 0 else float('inf')

    alt_test_results = calculate_dnds_statistical_tests(
        alt_total_syn_diffs, alt_total_nonsyn_diffs,
        alt_total_syn_sites, alt_total_nonsyn_sites,
        statistical_test
    )

    return {
        'name': name,
        'length': len(reference_seq),
        'reference_sequence': reference_seq,
        'amino_acid_sequence': amino_acid_sequence,
        'genetic_code': genetic_code,
        'unique_variants': len(variant_counts),
        'variants_above_frequency': variants_above_frequency,
        'variants_with_synonymous': variants_with_synonymous,
        'variants_with_nonsynonymous': variants_with_nonsynonymous,
        'all_variants_with_synonymous': all_variants_with_synonymous,
        'all_variants_with_nonsynonymous': all_variants_with_nonsynonymous,
        'total_sequences': total_sequences,
        'synonymous_differences': total_syn_diffs,
        'nonsynonymous_differences': total_nonsyn_diffs,
        'synonymous_sites': total_syn_sites,
        'nonsynonymous_sites': total_nonsyn_sites,
        'synonymous_rate': syn_rate,
        'nonsynonymous_rate': nonsyn_rate,
        'dnds_ratio': dnds_ratio,
        'statistic': test_results.get('z_statistic', test_results.get('statistic', 0.0)),
        'z_statistic': test_results.get('z_statistic', test_results.get('statistic', 0.0)),  # For backward compatibility
        'p_value': test_results['p_value'],
        'is_significant': test_results['is_significant'],
        'test_method': test_results['test_method'],
        'all_test_results': test_results,  # Store all test results for comparison
        'weighting_method': 'weighted' if weight_by_frequency else 'unweighted',
        # Alternative weighting results for comparison
        'comparison_results': {
            'weighting_method': alt_method,
            'synonymous_differences': alt_total_syn_diffs,
            'nonsynonymous_differences': alt_total_nonsyn_diffs,
            'synonymous_sites': alt_total_syn_sites,
            'nonsynonymous_sites': alt_total_nonsyn_sites,
            'synonymous_rate': alt_syn_rate,
            'nonsynonymous_rate': alt_nonsyn_rate,
            'dnds_ratio': alt_dnds_ratio,
            'statistic': alt_test_results.get('z_statistic', alt_test_results.get('statistic', 0.0)),
            'p_value': alt_test_results['p_value'],
            'is_significant': alt_test_results['is_significant'],
            'test_method': alt_test_results['test_method']
        }
    }

def analyze_microprotein_dnds(
    df: pd.DataFrame,
    microprotein_regions: List[List],
    gene_name: str,
    sequence_column: str,
    reference_sequence: Optional[str] = None,
    expected_sequences: Optional[Dict[str, str]] = None,
    genetic_codes: Optional[Dict[str, str]] = None,
    min_frequency: float = 0.01,
    coordinate_search_window: int = 50,
    min_coordinate_similarity: float = 0.8,
    plot_results: bool = True,
    output_file: Optional[str] = None,
    figsize: Tuple[int, int] = (15, 10),
    within_protein_mode: bool = False,
    gene_column : str = None,
    encompassing_genetic_code: str = 'mitochondrial',
    encompassing_gene_strand : str = '+',
    encompassing_gene_start : int = 1,
    encompassing_gene_reference_sequence: Optional[str] = None,
    statistical_test: str = 'z_test_robust',
    max_sequences : int = 50_000,
    max_seq_len : int = 10_000,
    allow_att_start : bool = False,
    dnds_method: str = 'NG86',
    threshold = 0.05,
    weight_by_frequency: bool = True
   ) -> Dict:
    """
    Analyze dN/dS ratios for micro-protein regions within DNA sequences.
    
    Parameters:
    -----------
    df : pd.DataFrame
        DataFrame with sequences. Should have columns for sequence data.
    microprotein_regions : List[List]
        List of [start, end, name, strand] for each micro-protein.
        Coordinates are absolute positions (1-based indexing) within the sequences.
    gene_name : str
        Name of the gene being analyzed.
    sequence_column : str
        Column name containing the sequences.
    reference_sequence : str, optional
        Reference sequence. If None, most common sequence is used.
        This should be the full sequence that the absolute coordinates refer to.
    expected_sequences : Dict[str, str], optional
        Dictionary mapping micro-protein names to their expected DNA sequences.
        If provided, coordinates will be adjusted to best match these sequences.
    genetic_codes : Dict[str, str], optional
        Dictionary mapping micro-protein names to genetic code type ('mitochondrial' or 'nuclear').
        If None, 'mitochondrial' is used for all.
    min_frequency : float
        Minimum allele frequency to consider for analysis.
    coordinate_search_window : int
        Search window size when looking for expected sequences.
    min_coordinate_similarity : float
        Minimum similarity threshold for coordinate adjustment.
    plot_results : bool
        Whether to generate plots.
    output_file : str, optional
        File to save results.
    figsize : Tuple[int, int]
        Figure size for plots.
    within_protein_mode : bool
        Whether to analyze the entire gene sequence as an encompassing protein 
        for all microproteins (useful for alternative reading frames within genes).
    encompassing_genetic_code : str
        Genetic code to use for the encompassing protein analysis ('mitochondrial' or 'nuclear').
        Only used when within_protein_mode=True
    encompassing_gene_strand : str
        Strand of the encompassing gene ('+' or '-'). Only used when within_protein_mode=True
    encompassing_gene_start : int
        Start position of the encompassing gene (1-based, absolute coordinates).
    max_sequences : int
        Maximum number of sequences to analyze in within_protein_mode.
        If the number of sequences exceeds this, random subsampling will be performed.
    max_seq_len : int
        Maximum length of sequences to analyze in within_protein_mode.
        If the length of the reference sequence exceeds this, it will be truncated.
    statistical_test : str
        Statistical test method: 'fisher' (Fisher's exact test - better for small counts),
        'z_test' (original Z-test), 'z_test_robust' (Z-test with pseudo-counts), or 'all'.
        Default 'fisher' is recommended for microprotein analysis.
    dnds_method : str
        Biopython method for dN/dS calculation: 'NG86' (default, Nei-Gojobori 1986),
        'LWL85' (Li et al. 1985), 'ML' (Goldman-Yang 1994), 'YN00' (Yang-Nielsen 2000)
    threshold : float
        Threshold for filtering the entire sequence
    weight_by_frequency : bool
        If True, weight differences and sites by variant frequency (default: True).
        If False, treat all variants equally regardless of frequency.
        Use False for population-level data where rare mutations should be considered equally.

    Returns:
    --------
    Dict with analysis results
    """
    logger.info(f"Starting dN/dS analysis for {gene_name}")
    logger.info(f"Using Biopython method: {dnds_method}")
    logger.info(f"Analyzing {len(microprotein_regions)} micro-protein regions")
    logger.info(f"Dataset contains {len(df)} sequences")
    logger.info(f"Within protein mode: {within_protein_mode}")
    
    # Extract sequences
    sequences = df[sequence_column].tolist()
    
    # Extract gene sequences
    gene_sequences = df[f"{gene_name}_seq"].tolist() if f"{gene_name}_seq" in df.columns else sequences
    if f"{gene_name}_seq" not in df.columns:
        logger.warning(f"Column '{gene_name}_seq' not found in DataFrame. Using {sequence_column} as gene sequences.")
    # Determine reference sequence
    if reference_sequence is None:
        # Use most common sequence
        sequence_counts = Counter(sequences)
        reference_sequence = sequence_counts.most_common(1)[0][0]
        logger.info("Using most common sequence as reference")
        # Generate gene reference
        encompassing_gene_reference_sequence = Counter(gene_sequences).most_common(1)[0][0]
        logger.info("Using most common gene sequence as reference")
    else:
        logger.info("Using provided reference sequence")
        if encompassing_gene_reference_sequence is not None: pass        
        else:
            logger.warning("No reference gene sequence provided, using reference sequence for gene analysis")
            encompassing_gene_reference_sequence = reference_sequence
            

    if within_protein_mode and (len(encompassing_gene_reference_sequence) > max_seq_len or len(gene_sequences) > max_sequences):
        logger.info(f"Large encompassing protein detected ({len(encompassing_gene_reference_sequence)} bp) with {len(gene_sequences)} samples. "
                f"Subsampling to {max_sequences} sequences for performance in within_protein_mode.")
        # Keep the reference sequence if it's in the data
        reference_in_data = reference_sequence in sequences
    
        # Randomly subsample
        np.random.seed(42)  # For reproducibility
        subsample_indices = np.random.choice(len(sequences), size=max_sequences, replace=False)
        sequences = [sequences[i] for i in subsample_indices]
        gene_sequences = [gene_sequences[i] for i in subsample_indices]        
        # Ensure reference sequence is included if it was in original data
        if reference_in_data and reference_sequence not in sequences:
            sequences[0] = reference_sequence
    # **NEW: Performance optimization for long encompassing proteins only**
    original_sequences_count = len(sequences)

    logger.info(f"Subsampled from {original_sequences_count} to {len(sequences)} sequences")
    print(f"Using {len(gene_sequences)} sequences for gene {gene_name}")
      
    results = {
        'gene_name': gene_name,
        'microprotein_results': {},
        'validation_results': {},
        'coordinate_adjustments': {},
        'summary_statistics': {},
        'encompassing_protein_results': {}  # New section for encompassing protein results
    }

    # Store raw P-values for multiple testing correction
    raw_p_values = []
    microprotein_names_order = []

    # Analyze each micro-protein region
    for region_info in microprotein_regions:
        if len(region_info) == 3:
            start, end, name = region_info
            strand = '+'  # Default to forward strand
        elif len(region_info) == 4:
            start, end, name, strand = region_info
        else:
            logger.error(f"Invalid region format: {region_info}")
            continue
            
        logger.info(f"Analyzing micro-protein: {name} ({start}-{end}, strand: {strand})")
        
        # Determine genetic code for this micro-protein
        if genetic_codes and name in genetic_codes:
            genetic_code = genetic_codes[name]
        else:
            genetic_code = 'mitochondrial'  # Default
        
        logger.info(f"Using {genetic_code} genetic code for {name}")
        
        # Use absolute coordinates directly
        # Ensure coordinates are within sequence bounds
        sequence_length = len(reference_sequence)
        if start < 1 or end > sequence_length:
            logger.warning(f"Micro-protein {name} coordinates ({start}-{end}) extend beyond sequence boundaries (1-{sequence_length}). Skipping.")
            continue
        
        # Extract reference sequence for this region using absolute coordinates (convert to 0-based indexing)
        if strand == '+':
            ref_region_seq = reference_sequence[start:end]
        else:
            ref_region_seq = reverse_complement(reference_sequence[start:end])
        logger.info(f"Extracted reference sequence for {name}: {ref_region_seq}")
        
        # Adjust coordinates if expected sequence is provided
        coordinate_adjusted = False
        original_coords = (start, end)
        
        if expected_sequences and name in expected_sequences:
            expected_seq = expected_sequences[name]
            logger.info(f"Searching for expected sequence for {name}")
            
            best_start, best_end, similarity, best_match = find_best_match_coordinates(
                reference_sequence, expected_seq, start, end, 
                coordinate_search_window, min_coordinate_similarity
            )
            
            if similarity >= min_coordinate_similarity and (best_start != start or best_end != end):
                logger.info(f"Adjusting coordinates for {name}: ({start}-{end}) -> ({best_start}-{best_end})")
                start, end = best_start, best_end
                ref_region_seq = best_match
                coordinate_adjusted = True
            else:
                logger.info(f"No significant coordinate adjustment needed for {name} (similarity: {similarity:.3f})")
            
            results['coordinate_adjustments'][name] = {
                'original_coords': original_coords,
                'adjusted_coords': (start, end),
                'similarity': similarity,
                'coordinate_adjusted': coordinate_adjusted,
                'expected_sequence': expected_seq,
                'found_sequence': best_match
            }
        
        # Validate the coding region
        validation = validate_coding_region(ref_region_seq, strand, genetic_code)
        results['validation_results'][name] = validation
        
        if not validation['is_valid']:
            logger.warning(f"Micro-protein {name} failed validation: {validation['issues']}")
            # Continue with analysis but flag the issues
        
        # Translate reference sequence to amino acids
        amino_acid_sequence = translate_sequence(ref_region_seq, genetic_code, allow_att_start=allow_att_start)
        logger.info(f"{name}: {' '.join([ref_region_seq[i:i+3] for i in range(0, len(ref_region_seq) + 1, 3)])}")
        # Align amino acid sequence with the above codon separated sequence

        logger.info(f"{name}: {'  ' + '   '.join(amino_acid_sequence[i: i + 1] for i in range(0, len(amino_acid_sequence) + 1))}")

        
        # Collect sequences for this region from all samples using absolute coordinates
        region_sequences = []
        for seq in sequences:
            if end <= len(seq):
                if strand == '+':
                    region_seq = seq[start:end]  # Convert to 0-based indexing
                else:
                    region_seq = reverse_complement(seq[start:end])
                region_sequences.append(region_seq)
        logger.info(f"Collected {len(region_sequences)} sequences for {name} from dataset")
        # Log the first 5 sequences for verification
        for i, seq in enumerate(region_sequences[:5]):
            logger.info(f"Sample sequence {i+1} for {name}: {seq}")
        if len(region_sequences) == 0:
            logger.warning(f"No valid sequences found for {name}")
            continue
        
        variant_counts = Counter(region_sequences)
        logger.info(f"Found {len(variant_counts)} unique variants for {name} before filtering")
        cleanup = cag.clean_sequences_adaptive(region_sequences, threshold = threshold)
        region_sequences = cleanup.cleaned_sequences
        # Report the number of codon positions and sequences removed during cleaning
        logger.info(f'out of {cleanup.original_seq_count} sequences, {cleanup.final_seq_count} were kept after cleaning for {name}')
        logger.info(f'out of {cleanup.original_codon_count} codon positions, {cleanup.final_codon_count} were kept after cleaning for {name}')
        # Calculate dN/dS for each variant vs reference using Biopython
        # Store both weighted and unweighted for comparison
        synonymous_diffs_weighted = []
        nonsynonymous_diffs_weighted = []
        synonymous_sites_weighted = []
        nonsynonymous_sites_weighted = []

        synonymous_diffs_unweighted = []
        nonsynonymous_diffs_unweighted = []
        synonymous_sites_unweighted = []
        nonsynonymous_sites_unweighted = []

        variant_counts = Counter(region_sequences)
        total_sequences = len(region_sequences)

        logger.info(f"Found {len(variant_counts)} unique variants for {name} after filtering")

        # Track variants above frequency threshold and variants with mutations
        variants_above_frequency = 0
        variants_with_synonymous = 0
        variants_with_nonsynonymous = 0
        all_variants_with_synonymous = 0
        all_variants_with_nonsynonymous = 0

        for variant_seq, count in variant_counts.items():
            frequency = count / total_sequences

            # Count all variants with mutations (regardless of frequency threshold)
            if variant_seq != ref_region_seq:
                diff_results_all = calculate_synonymous_nonsynonymous(
                    ref_region_seq, variant_seq, strand, genetic_code, method=dnds_method
                )
                if diff_results_all['synonymous_differences'] > 0:
                    all_variants_with_synonymous += 1
                if diff_results_all['nonsynonymous_differences'] > 0:
                    all_variants_with_nonsynonymous += 1

            if frequency < min_frequency:
                continue  # Skip rare variants

            # Count variants above frequency threshold (excluding reference)
            if variant_seq != ref_region_seq:
                variants_above_frequency += 1

            if variant_seq == ref_region_seq:
                continue  # Skip reference sequence

            # Calculate differences using Biopython (use already calculated if same variant)
            diff_results = calculate_synonymous_nonsynonymous(
                ref_region_seq, variant_seq, strand, genetic_code, method=dnds_method
            )

            # Count variants with observed mutations (above frequency threshold)
            if diff_results['synonymous_differences'] > 0:
                variants_with_synonymous += 1
            if diff_results['nonsynonymous_differences'] > 0:
                variants_with_nonsynonymous += 1

            # Always calculate both weighted and unweighted for comparison
            # Weighted version
            synonymous_diffs_weighted.append(diff_results['synonymous_differences'] * frequency)
            nonsynonymous_diffs_weighted.append(diff_results['nonsynonymous_differences'] * frequency)
            synonymous_sites_weighted.append(diff_results['synonymous_sites'] * frequency)
            nonsynonymous_sites_weighted.append(diff_results['nonsynonymous_sites'] * frequency)

            # Unweighted version
            synonymous_diffs_unweighted.append(diff_results['synonymous_differences'])
            nonsynonymous_diffs_unweighted.append(diff_results['nonsynonymous_differences'])
            synonymous_sites_unweighted.append(diff_results['synonymous_sites'])
            nonsynonymous_sites_unweighted.append(diff_results['nonsynonymous_sites'])

        # Select which version to use for main results
        if weight_by_frequency:
            synonymous_diffs = synonymous_diffs_weighted
            nonsynonymous_diffs = nonsynonymous_diffs_weighted
            synonymous_sites = synonymous_sites_weighted
            nonsynonymous_sites = nonsynonymous_sites_weighted
        else:
            synonymous_diffs = synonymous_diffs_unweighted
            nonsynonymous_diffs = nonsynonymous_diffs_unweighted
            synonymous_sites = synonymous_sites_unweighted
            nonsynonymous_sites = nonsynonymous_sites_unweighted
        
        # Calculate summary statistics
        total_syn_diffs = sum(synonymous_diffs)
        total_nonsyn_diffs = sum(nonsynonymous_diffs)
        total_syn_sites = sum(synonymous_sites) if synonymous_sites else 1  # Avoid division by zero
        total_nonsyn_sites = sum(nonsynonymous_sites) if nonsynonymous_sites else 1
        
        logger.info(f"Total synonymous sites for {name}: {total_syn_sites}")
        logger.info(f"Total non-synonymous sites for {name}: {total_nonsyn_sites}")
        # Calculate rates
        syn_rate = total_syn_diffs / total_syn_sites if total_syn_sites > 0 else 0
        nonsyn_rate = total_nonsyn_diffs / total_nonsyn_sites if total_nonsyn_sites > 0 else 0
        
        # Calculate dN/dS ratio
        dnds_ratio = nonsyn_rate / syn_rate if syn_rate > 0 else float('inf')
        
        logger.info(f"Total synonymous differences for {name}: {total_syn_diffs}")
        logger.info(f"Total non-synonymous differences for {name}: {total_nonsyn_diffs}")
        logger.info(f"Synonymous rate for {name}: {syn_rate:.6f}")
        logger.info(f"Non-synonymous rate for {name}: {nonsyn_rate:.6f}")
        
        # Calculate statistical tests for neutrality
        test_results = calculate_dnds_statistical_tests(
            total_syn_diffs, total_nonsyn_diffs, 
            total_syn_sites, total_nonsyn_sites, 
            statistical_test
        )
        
        # Calculate confidence intervals (using bootstrap or analytical methods)
        # For simplicity, using Poisson confidence intervals
        def poisson_ci(count, alpha=0.05):
            from scipy.stats import chi2
            lower = chi2.ppf(alpha/2, 2*count) / 2 if count > 0 else 0
            upper = chi2.ppf(1-alpha/2, 2*count+2) / 2
            return lower, upper
        
        syn_ci = poisson_ci(total_syn_diffs)
        nonsyn_ci = poisson_ci(total_nonsyn_diffs)

        # Calculate the alternative weighting method for comparison
        if weight_by_frequency:
            # Main results are weighted, so calculate unweighted for comparison
            alt_syn_diffs = synonymous_diffs_unweighted
            alt_nonsyn_diffs = nonsynonymous_diffs_unweighted
            alt_syn_sites = synonymous_sites_unweighted
            alt_nonsyn_sites = nonsynonymous_sites_unweighted
            alt_method = 'unweighted'
        else:
            # Main results are unweighted, so calculate weighted for comparison
            alt_syn_diffs = synonymous_diffs_weighted
            alt_nonsyn_diffs = nonsynonymous_diffs_weighted
            alt_syn_sites = synonymous_sites_weighted
            alt_nonsyn_sites = nonsynonymous_sites_weighted
            alt_method = 'weighted'

        # Calculate alternative results
        alt_total_syn_diffs = sum(alt_syn_diffs)
        alt_total_nonsyn_diffs = sum(alt_nonsyn_diffs)
        alt_total_syn_sites = sum(alt_syn_sites) if alt_syn_sites else 1
        alt_total_nonsyn_sites = sum(alt_nonsyn_sites) if alt_nonsyn_sites else 1

        alt_syn_rate = alt_total_syn_diffs / alt_total_syn_sites if alt_total_syn_sites > 0 else 0
        alt_nonsyn_rate = alt_total_nonsyn_diffs / alt_total_nonsyn_sites if alt_total_nonsyn_sites > 0 else 0
        alt_dnds_ratio = alt_nonsyn_rate / alt_syn_rate if alt_syn_rate > 0 else float('inf')

        alt_test_results = calculate_dnds_statistical_tests(
            alt_total_syn_diffs, alt_total_nonsyn_diffs,
            alt_total_syn_sites, alt_total_nonsyn_sites,
            statistical_test
        )

        # Store raw P-value for multiple testing correction
        raw_p_values.append(test_results['p_value'])
        microprotein_names_order.append(name)

        results['microprotein_results'][name] = {
            'coordinates': {
                'absolute_start': start,
                'absolute_end': end,
                'strand': strand,
                'coordinate_adjusted': coordinate_adjusted
            },
            'sequence_info': {
                'length': len(ref_region_seq),
                'reference_sequence': ref_region_seq,
                'amino_acid_sequence': amino_acid_sequence,
                'genetic_code': genetic_code,
                'unique_variants': len(variant_counts),
                'variants_above_frequency': variants_above_frequency,
                'variants_with_synonymous': variants_with_synonymous,
                'variants_with_nonsynonymous': variants_with_nonsynonymous,
                'all_variants_with_synonymous': all_variants_with_synonymous,
                'all_variants_with_nonsynonymous': all_variants_with_nonsynonymous,
                'total_sequences': total_sequences
            },
            'mutation_counts': {
                'synonymous_differences': total_syn_diffs,
                'nonsynonymous_differences': total_nonsyn_diffs,
                'synonymous_sites': total_syn_sites,
                'nonsynonymous_sites': total_nonsyn_sites
            },
            'rates': {
                'synonymous_rate': syn_rate,
                'nonsynonymous_rate': nonsyn_rate,
                'dnds_ratio': dnds_ratio
            },
            'statistical_tests': {
                'statistic': test_results.get('z_statistic', test_results.get('statistic', 0.0)),
                'p_value': test_results['p_value'],
                'is_significant': test_results['is_significant'],
                'test_method': test_results['test_method'],
                'all_test_results': test_results  # Store all test results for comparison
            },
            'confidence_intervals': {
                'synonymous_ci': syn_ci,
                'nonsynonymous_ci': nonsyn_ci
            },
            'validation': validation,
            'weighting_method': 'weighted' if weight_by_frequency else 'unweighted',
            # Alternative weighting results for comparison
            'comparison_results': {
                'weighting_method': alt_method,
                'mutation_counts': {
                    'synonymous_differences': alt_total_syn_diffs,
                    'nonsynonymous_differences': alt_total_nonsyn_diffs,
                    'synonymous_sites': alt_total_syn_sites,
                    'nonsynonymous_sites': alt_total_nonsyn_sites
                },
                'rates': {
                    'synonymous_rate': alt_syn_rate,
                    'nonsynonymous_rate': alt_nonsyn_rate,
                    'dnds_ratio': alt_dnds_ratio
                },
                'statistical_tests': {
                    'statistic': alt_test_results.get('z_statistic', alt_test_results.get('statistic', 0.0)),
                    'p_value': alt_test_results['p_value'],
                    'is_significant': alt_test_results['is_significant'],
                    'test_method': alt_test_results['test_method']
                }
            }
        }
        
        # Analyze encompassing protein if in within_protein_mode
        if within_protein_mode:
            logger.info(f"Analyzing encompassing gene {gene_name} for microprotein {name}")
            
            # Use the entire gene sequence as the encompassing protein
            # Extract full gene sequences for all samples
            # Analyze gene-level dN/dS
            if gene_sequences:
                # Clean up gene sequences
                cleanup_gene = cag.clean_sequences_adaptive(gene_sequences, threshold = threshold)
                gene_sequences = cleanup_gene.cleaned_sequences
                logger.info(f'Cleaned up gene sequences - {len(gene_sequences)} left')
                enc_results = analyze_protein_dnds(
                    gene_sequences, encompassing_gene_reference_sequence, gene_name,
                    encompassing_genetic_code, encompassing_gene_strand, min_frequency, statistical_test,
                    allow_att_start=allow_att_start, dnds_method=dnds_method,
                    weight_by_frequency=weight_by_frequency
                )
                results['encompassing_protein_results'][name] = enc_results
                logger.info(f"Encompassing gene {gene_name} dN/dS ratio: {enc_results['dnds_ratio']:.4f}")
        
        logger.info(f"  dN/dS ratio for {name}: {dnds_ratio:.4f}")
        logger.info(f"  Synonymous rate: {syn_rate:.6f}")
        logger.info(f"  Non-synonymous rate: {nonsyn_rate:.6f}")
        logger.info(f"  Statistical test ({test_results['test_method']}): statistic={test_results.get('z_statistic', test_results.get('statistic', 0.0)):.4f}, p-value: {test_results['p_value']:.4e}")
        logger.info(f"  Amino acid sequence: {amino_acid_sequence}")

    # Apply Benjamini-Hochberg multiple test correction
    if len(raw_p_values) > 0:
        logger.info(f"Applying Benjamini-Hochberg correction to {len(raw_p_values)} P-values")

        # Apply false discovery rate control using Benjamini-Hochberg procedure
        corrected_p_values = false_discovery_control(raw_p_values, method='bh')

        # Update results with corrected P-values (overwriting raw P-values)
        for i, name in enumerate(microprotein_names_order):
            if name in results['microprotein_results']:
                raw_p_val = raw_p_values[i]
                corrected_p_val = corrected_p_values[i]

                # Overwrite the P-value with the corrected one
                results['microprotein_results'][name]['statistical_tests']['p_value'] = corrected_p_val
                # Update significance based on corrected P-value
                results['microprotein_results'][name]['statistical_tests']['is_significant'] = corrected_p_val < 0.05

                logger.info(f"  {name}: Raw P-value = {raw_p_val:.4e}, Corrected P-value = {corrected_p_val:.4e}")

    # Generate summary statistics
    all_dnds = [result['rates']['dnds_ratio'] for result in results['microprotein_results'].values() 
                if not np.isinf(result['rates']['dnds_ratio'])]
    
    results['summary_statistics'] = {
        'mean_dnds': np.mean(all_dnds) if all_dnds else np.nan,
        'median_dnds': np.median(all_dnds) if all_dnds else np.nan,
        'std_dnds': np.std(all_dnds) if all_dnds else np.nan,
        'n_microproteins': len(results['microprotein_results']),
        'n_significant': sum(1 for result in results['microprotein_results'].values() 
                           if result['statistical_tests']['is_significant'])
    }
    
    # Create results DataFrame
    results_data = []
    for name, mp_results in results['microprotein_results'].items():
        row_data = {
            'MicroProtein': name,
            'Gene': gene_name,
            'AbsoluteStart': mp_results['coordinates']['absolute_start'],
            'AbsoluteEnd': mp_results['coordinates']['absolute_end'],
            'Strand': mp_results['coordinates']['strand'],
            'GeneticCode': mp_results['sequence_info']['genetic_code'],
            'CoordinateAdjusted': mp_results['coordinates']['coordinate_adjusted'],
            'Length': mp_results['sequence_info']['length'],
            'AminoAcidSequence': mp_results['sequence_info']['amino_acid_sequence'],
            'UniqueVariants': mp_results['sequence_info']['unique_variants'],
            'VariantsAboveFrequency': mp_results['sequence_info']['variants_above_frequency'],
            'VariantsWithSynonymous': mp_results['sequence_info']['variants_with_synonymous'],
            'VariantsWithNonsynonymous': mp_results['sequence_info']['variants_with_nonsynonymous'],
            'AllVariantsWithSynonymous': mp_results['sequence_info']['all_variants_with_synonymous'],
            'AllVariantsWithNonsynonymous': mp_results['sequence_info']['all_variants_with_nonsynonymous'],
            'SynonymousDifferences': mp_results['mutation_counts']['synonymous_differences'],
            'NonsynonymousDifferences': mp_results['mutation_counts']['nonsynonymous_differences'],
            'SynonymousRate': mp_results['rates']['synonymous_rate'],
            'NonsynonymousRate': mp_results['rates']['nonsynonymous_rate'],
            'dN_dS_Ratio': mp_results['rates']['dnds_ratio'],
            'Statistic': mp_results['statistical_tests']['statistic'],
            'P_Value': mp_results['statistical_tests']['p_value'],
            'Is_Significant': mp_results['statistical_tests']['is_significant'],
            'Test_Method': mp_results['statistical_tests']['test_method'],
            'IsValid': mp_results['validation']['is_valid'],
            'ValidationIssues': '; '.join(mp_results['validation']['issues'])
        }
        
        # Add encompassing protein data if available
        if name in results['encompassing_protein_results']:
            enc_results = results['encompassing_protein_results'][name]
            row_data.update({
                'EncompassingProtein': enc_results['name'],
                'EncompassingProtein_Length': enc_results['length'],
                'EncompassingProtein_dN_dS_Ratio': enc_results['dnds_ratio'],
                'EncompassingProtein_Statistic': enc_results.get('z_statistic', enc_results.get('statistic', 0.0)),
                'EncompassingProtein_P_Value': enc_results['p_value'],
                'EncompassingProtein_Is_Significant': enc_results['is_significant'],
                'EncompassingProtein_Test_Method': enc_results.get('test_method', 'Unknown'),
                'EncompassingProtein_AminoAcidSequence': enc_results['amino_acid_sequence']
            })
        
        results_data.append(row_data)
    # **NEW: Add encompassing proteins as separate entries if within_protein_mode is True**
    if within_protein_mode and results['encompassing_protein_results']:
        for name, enc_results in results['encompassing_protein_results'].items():
            # Only add one encompassing protein entry (they're all the same gene)
            if name == list(results['encompassing_protein_results'].keys())[0]:
                enc_row_data = {
                    'MicroProtein': f"{enc_results['name']}_encompassing",
                    'Gene': gene_name,
                    'AbsoluteStart': encompassing_gene_start,  # Full gene
                    'AbsoluteEnd': encompassing_gene_start + len(encompassing_gene_reference_sequence),  # Full gene
                    'Strand': '+',
                    'GeneticCode': encompassing_genetic_code,
                    'CoordinateAdjusted': False,
                    'Length': enc_results['length'],
                    'AminoAcidSequence': enc_results['amino_acid_sequence'],
                    'UniqueVariants': enc_results['unique_variants'],
                    'VariantsAboveFrequency': enc_results['variants_above_frequency'],
                    'VariantsWithSynonymous': enc_results['variants_with_synonymous'],
                    'VariantsWithNonsynonymous': enc_results['variants_with_nonsynonymous'],
                    'AllVariantsWithSynonymous': enc_results['all_variants_with_synonymous'],
                    'AllVariantsWithNonsynonymous': enc_results['all_variants_with_nonsynonymous'],
                    'SynonymousDifferences': enc_results['synonymous_differences'],
                    'NonsynonymousDifferences': enc_results['nonsynonymous_differences'],
                    'SynonymousRate': enc_results['synonymous_rate'],
                    'NonsynonymousRate': enc_results['nonsynonymous_rate'],
                    'dN_dS_Ratio': enc_results['dnds_ratio'],
                    'Statistic': enc_results.get('z_statistic', enc_results.get('statistic', 0.0)),
                    'P_Value': enc_results['p_value'],
                    'Is_Significant': enc_results['is_significant'],
                    'Test_Method': enc_results.get('test_method', 'Unknown'),
                    'IsValid': True,  # Assume encompassing proteins are valid
                    'ValidationIssues': '',
                    'Type': 'EncompassingProtein'  # **NEW: Distinguish type**
                }
                results_data.append(enc_row_data)

    results_df = pd.DataFrame(results_data)
    results['results_df'] = results_df
        
    # Plotting with significance indicators
    if plot_results and len(results_data) > 0:
        fig, axes = plt.subplots(2, 2, figsize=figsize)
        
        # **NEW: Separate microproteins and encompassing proteins**
        mp_data = [r for r in results_data if r.get('Type', 'MicroProtein') == 'MicroProtein']
        enc_data = [r for r in results_data if r.get('Type', 'MicroProtein') == 'EncompassingProtein']
        all_plot_data = mp_data + enc_data
        
        # Plot 1: dN/dS ratios with significance indicators
        ax = axes[0, 0]
        valid_dnds = [r['dN_dS_Ratio'] for r in all_plot_data if not np.isinf(r['dN_dS_Ratio'])]
        names = [r['MicroProtein'] for r in all_plot_data if not np.isinf(r['dN_dS_Ratio'])]
        genetic_codes_plot = [r['GeneticCode'] for r in all_plot_data if not np.isinf(r['dN_dS_Ratio'])]
        is_significant_plot = [r['Is_Significant'] for r in all_plot_data if not np.isinf(r['dN_dS_Ratio'])]
        types_plot = [r.get('Type', 'MicroProtein') for r in all_plot_data if not np.isinf(r['dN_dS_Ratio'])]
        
        if valid_dnds:
            bars = ax.bar(range(len(valid_dnds)), valid_dnds)
            ax.axhline(y=1, color='red', linestyle='--', alpha=0.7, label='Neutral (dN/dS=1)')
            ax.set_xticks(range(len(valid_dnds)))
            ax.set_xticklabels([f"{name}\n({gc})" for name, gc in zip(names, genetic_codes_plot)], 
                            rotation=45, fontsize=9)
            ax.set_ylabel('dN/dS Ratio')
            ax.set_title('dN/dS Ratios for Micro-Proteins and Encompassing Proteins')
            ax.legend()
            
            # **NEW: Color bars based on type, selection pressure and genetic code**
            for i, (bar, ratio, gc, is_sig, ptype) in enumerate(zip(bars, valid_dnds, genetic_codes_plot, is_significant_plot, types_plot)):
                if ptype == 'EncompassingProtein':
                    # Use different styling for encompassing proteins
                    if ratio > 1:
                        color = 'darkred' if gc == 'mitochondrial' else 'darkorange'
                    elif ratio < 1:
                        color = 'darkblue' if gc == 'mitochondrial' else 'darkslateblue'
                    else:
                        color = 'darkgray'
                    bar.set_edgecolor('black')
                    bar.set_linewidth(2)
                else:
                    # Original coloring for microproteins
                    if ratio > 1:
                        color = 'red' if gc == 'mitochondrial' else 'orange'
                    elif ratio < 1:
                        color = 'blue' if gc == 'mitochondrial' else 'lightblue'
                    else:
                        color = 'gray'
                
                bar.set_color(color)
                
                # Add asterisk for significant results
                if is_sig:
                    ax.text(i, ratio + max(valid_dnds) * 0.02, '*', ha='center', va='bottom', 
                        fontsize=16, fontweight='bold', color='black')
        
        ax.grid(False)

        # Plot 2: Synonymous vs Non-synonymous differences
        ax = axes[0, 1]
        syn_diffs = [r['SynonymousDifferences'] for r in all_plot_data]
        nonsyn_diffs = [r['NonsynonymousDifferences'] for r in all_plot_data]
        names = [r['MicroProtein'] for r in all_plot_data]
        types_plot = [r.get('Type', 'MicroProtein') for r in all_plot_data]
        
        x = np.arange(len(names))
        width = 0.35
        
        # **NEW: Different alpha for encompassing proteins**
        syn_bars = ax.bar(x - width/2, syn_diffs, width, label='Synonymous', alpha=0.7)
        nonsyn_bars = ax.bar(x + width/2, nonsyn_diffs, width, label='Non-synonymous', alpha=0.7)
        
        # Style encompassing proteins differently
        for i, ptype in enumerate(types_plot):
            if ptype == 'EncompassingProtein':
                syn_bars[i].set_edgecolor('black')
                syn_bars[i].set_linewidth(2)
                nonsyn_bars[i].set_edgecolor('black')
                nonsyn_bars[i].set_linewidth(2)
        
        ax.set_xticks(x)
        ax.set_xticklabels(names, rotation=45)
        ax.set_ylabel('Number of Differences')
        ax.set_title('Synonymous vs Non-synonymous Differences')
        ax.legend()
        ax.grid(False)
        
        # Plot 3: Mutation rates
        ax = axes[1, 0]
        syn_rates = [r['SynonymousRate'] for r in all_plot_data]
        nonsyn_rates = [r['NonsynonymousRate'] for r in all_plot_data]
        
        syn_bars = ax.bar(x - width/2, syn_rates, width, label='Synonymous rate', alpha=0.7)
        nonsyn_bars = ax.bar(x + width/2, nonsyn_rates, width, label='Non-synonymous rate', alpha=0.7)
        
        # Style encompassing proteins differently
        for i, ptype in enumerate(types_plot):
            if ptype == 'EncompassingProtein':
                syn_bars[i].set_edgecolor('black')
                syn_bars[i].set_linewidth(2)
                nonsyn_bars[i].set_edgecolor('black')
                nonsyn_bars[i].set_linewidth(2)
        
        ax.set_xticks(x)
        ax.set_xticklabels(names, rotation=45)
        ax.set_ylabel('Mutation Rate')
        ax.set_title('Synonymous vs Non-synonymous Rates')
        ax.legend()
        ax.grid(False)
        
        # Plot 4: P-values and significance
        ax = axes[1, 1]
        p_values = [r['P_Value'] for r in all_plot_data]
        
        bars = ax.bar(range(len(names)), -np.log10(p_values), alpha=0.7)
        ax.axhline(y=-np.log10(0.05), color='red', linestyle='--', alpha=0.7, label='p=0.05')
        ax.set_xticks(range(len(names)))
        ax.set_xticklabels(names, rotation=45)
        ax.set_ylabel('-log10(p-value)')
        ax.set_title('Statistical Significance of dN/dS Tests')
        ax.legend()
        
        # Color bars based on significance and type
        for i, (bar, p_val, ptype) in enumerate(zip(bars, p_values, types_plot)):
            if p_val < 0.05:
                bar.set_color('red')
            else:
                bar.set_color('gray')
            
            if ptype == 'EncompassingProtein':
                bar.set_edgecolor('black')
                bar.set_linewidth(2)
        
        sns.despine()
        ax.grid(False)
        plt.tight_layout()
        
        if output_file:
            plt.savefig(output_file, dpi=300, bbox_inches='tight')
            logger.info(f"Plot saved to {output_file}")
        else:
            plt.show()
    
    return results

def plot_combined_microprotein_analysis(
    results_dict: Dict[str, Dict],
    output_file: Optional[str] = None,
    figsize: Tuple[int, int] = (16, 8),
    colors: Optional[List[str]] = None,
    show_genetic_codes: bool = True,
    show_gene_labels: bool = True,
    title: Optional[str] = None,
    display_summary_text: bool = True,
    include_encompassing_proteins: bool = True,  # **NEW PARAMETER**
    ylim = None,
    legend_outside: bool = False,
    break_axis_threshold: float = 50.0  # **NEW PARAMETER**
) -> None:
    """
    Plot combined results from multiple analyze_microprotein_dnds runs.

    Parameters:
    -----------
    include_encompassing_proteins : bool
        Whether to include encompassing proteins in the plots (if available).
    break_axis_threshold : float
        If the difference between max and min y-values exceeds this threshold,
        use broken axes for better visualization. Default is 50. Set to None to disable.
    """
    logger.info(f"Creating combined plot for {len(results_dict)} gene analyses")

    # Try to import brokenaxes for broken axis functionality
    try:
        from brokenaxes import brokenaxes
        brokenaxes_available = True
    except ImportError:
        brokenaxes_available = False
        if break_axis_threshold is not None:
            logger.warning("brokenaxes library not available. Install with: pip install brokenaxes")
    
    # Extract data from all results
    all_data = []
    gene_colors = {}
    
    if colors is None:
        colors = plt.cm.Set1(np.linspace(0, 1, len(results_dict)))
    
    for i, (gene_name, results) in enumerate(results_dict.items()):
        gene_color = colors[i % len(colors)]
        gene_colors[gene_name] = gene_color
        
        if 'microprotein_results' not in results:
            logger.warning(f"No microprotein results found for {gene_name}")
            continue
            
        # **NEW: Process microproteins**
        for mp_name, mp_results in results['microprotein_results'].items():
            # Create label for this micro-protein
            label_parts = []
            if show_gene_labels:
                label_parts.append(gene_name)
            label_parts.append(mp_name)
            if show_genetic_codes:
                genetic_code = mp_results['sequence_info'].get('genetic_code', 'unknown')
                label_parts.append(f"({genetic_code})")
            
            full_label = mp_name
            
            data_point = {
                'gene': gene_name,
                'microprotein': mp_name,
                'full_label': full_label,
                'dnds_ratio': mp_results['rates']['dnds_ratio'],
                'syn_rate': mp_results['rates']['synonymous_rate'],
                'nonsyn_rate': mp_results['rates']['nonsynonymous_rate'],
                'genetic_code': mp_results['sequence_info'].get('genetic_code', 'unknown'),
                'is_valid': mp_results['validation']['is_valid'],
                'is_significant': mp_results['statistical_tests']['is_significant'],
                'p_value': mp_results['statistical_tests']['p_value'],
                'color': gene_color,
                'type': 'MicroProtein'  # **NEW: Mark as microprotein**
            }
            all_data.append(data_point)
        
        # **NEW: Process encompassing proteins if requested and available**
        if include_encompassing_proteins and 'encompassing_protein_results' in results:
            # Only add one encompassing protein per gene (they're all the same)
            if results['encompassing_protein_results']:
                first_mp = list(results['encompassing_protein_results'].keys())[0]
                enc_results = results['encompassing_protein_results'][first_mp]
                
                # Create label for encompassing protein
                label_parts = []
                if show_gene_labels:
                    label_parts.append(gene_name)
                label_parts.append(f"Full Gene")
                if show_genetic_codes:
                    enc_genetic_code = enc_results.get('genetic_code', 'unknown')
                    label_parts.append(f"({enc_genetic_code})")
                
                full_label = ':'.join(label_parts[:2])
                
                data_point = {
                    'gene': gene_name,
                    'microprotein': f"Full Gene",
                    'full_label': full_label,
                    'dnds_ratio': enc_results['dnds_ratio'],
                    'syn_rate': enc_results['synonymous_rate'],
                    'nonsyn_rate': enc_results['nonsynonymous_rate'],
                    'genetic_code': enc_results.get('genetic_code', 'unknown'),
                    'is_valid': True,  # Assume encompassing proteins are valid
                    'is_significant': enc_results['is_significant'],
                    'p_value': enc_results['p_value'],
                    'color': gene_color,
                    'type': 'EncompassingProtein'  # **NEW: Mark as encompassing protein**
                }
                all_data.append(data_point)
    
    if not all_data:
        logger.error("No valid data found for plotting")
        return
    
    # Filter out infinite dN/dS ratios for plotting
    valid_data = [d for d in all_data if not np.isinf(d['dnds_ratio']) and not np.isnan(d['dnds_ratio'])]
    
    if not valid_data:
        logger.error("No valid dN/dS ratios found for plotting")
        return
    
    logger.info(f"Plotting {len(valid_data)} items from {len(results_dict)} genes")

    # Extract data for plotting
    labels = [d['full_label'] for d in valid_data]
    dnds_ratios = [d['dnds_ratio'] for d in valid_data]
    syn_rates = [d['syn_rate'] for d in valid_data]
    nonsyn_rates = [d['nonsyn_rate'] for d in valid_data]
    colors_plot = [d['color'] for d in valid_data]
    genetic_codes = [d['genetic_code'] for d in valid_data]
    is_valid = [d['is_valid'] for d in valid_data]
    is_significant = [d['is_significant'] for d in valid_data]
    types = [d['type'] for d in valid_data]  # **NEW: Get types**

    x_positions = np.arange(len(labels))

    # **NEW: Check if we need broken axes**
    use_broken_axes = False
    if break_axis_threshold is not None and brokenaxes_available:
        # Check dN/dS ratios range
        min_dnds = min(dnds_ratios)
        max_dnds = max(dnds_ratios)
        dnds_range = max_dnds - min_dnds

        # Check rates range
        all_rates = syn_rates + nonsyn_rates
        min_rate = min(all_rates)
        max_rate = max(all_rates)
        rates_range = max_rate - min_rate

        if dnds_range > break_axis_threshold or rates_range > break_axis_threshold:
            use_broken_axes = True
            logger.info(f"Using broken axes (dN/dS range: {dnds_range:.2f}, rates range: {rates_range:.2f})")

            # Calculate appropriate breaks for dN/dS plot
            if dnds_range > break_axis_threshold:
                # Find natural break points by looking at gaps in the data
                sorted_dnds = sorted(dnds_ratios)
                gaps = []
                for i in range(len(sorted_dnds) - 1):
                    gap_size = sorted_dnds[i + 1] - sorted_dnds[i]
                    if gap_size > break_axis_threshold * 0.3:  # At least 30% of threshold
                        gaps.append((sorted_dnds[i], sorted_dnds[i + 1], gap_size))

                if gaps:
                    # Use the largest gap
                    largest_gap = max(gaps, key=lambda x: x[2])
                    dnds_ylims = ((0, largest_gap[0] * 1.1), (largest_gap[1] * 0.9, max_dnds * 1.1))
                else:
                    # No clear gap, divide roughly in middle
                    median_val = np.median(dnds_ratios)
                    dnds_ylims = ((0, median_val * 1.5), (median_val * 2, max_dnds * 1.1))
            else:
                dnds_ylims = None

            # Calculate appropriate breaks for rates plot
            if rates_range > break_axis_threshold:
                sorted_rates = sorted(all_rates)
                gaps = []
                for i in range(len(sorted_rates) - 1):
                    gap_size = sorted_rates[i + 1] - sorted_rates[i]
                    if gap_size > break_axis_threshold * 0.3:
                        gaps.append((sorted_rates[i], sorted_rates[i + 1], gap_size))

                if gaps:
                    largest_gap = max(gaps, key=lambda x: x[2])
                    rates_ylims = ((0, largest_gap[0] * 1.1), (largest_gap[1] * 0.9, max_rate * 1.1))
                else:
                    median_val = np.median(all_rates)
                    rates_ylims = ((0, median_val * 1.5), (median_val * 2, max_rate * 1.1))
            else:
                rates_ylims = None

    # Create the plot (with or without broken axes)
    if use_broken_axes:
        from matplotlib.gridspec import GridSpec

        fig = plt.figure(figsize=figsize)
        if title is None:
            title = f"Micro-Protein Selection Analysis ({len(results_dict)} Genes)"
        fig.suptitle(title, fontsize=16, fontweight='bold', y=0.98)

        # Create GridSpec for subplot layout
        gs = GridSpec(1, 2, figure=fig)

        # Create broken axes for left plot (dN/dS)
        if dnds_ylims is not None:
            ax1 = brokenaxes(xlims=((0, len(labels)),), ylims=dnds_ylims,
                           subplot_spec=gs[0], fig=fig, despine=False)
        else:
            ax1 = fig.add_subplot(gs[0])

        # Create broken axes for right plot (rates)
        if rates_ylims is not None:
            ax2 = brokenaxes(xlims=((0, len(labels)),), ylims=rates_ylims,
                           subplot_spec=gs[1], fig=fig, despine=False)
        else:
            ax2 = fig.add_subplot(gs[1])
    else:
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=figsize)
        if title is None:
            title = f"Micro-Protein Selection Analysis ({len(results_dict)} Genes)"
        fig.suptitle(title, fontsize=16, fontweight='bold')
    
    # Plot 1: dN/dS ratios with significance indicators
    bars1 = ax1.bar(x_positions, dnds_ratios, color=colors_plot, alpha=0.7, edgecolor='black', linewidth=0.5)

    # Add horizontal line at dN/dS = 1 (neutral selection)
    ax1.axhline(y=1, color='red', linestyle='--', alpha=0.8, linewidth=2, label='Neutral (dN/dS=1)')

    # **NEW: Extract individual bar patches (needed for broken axes compatibility)**
    from matplotlib.container import BarContainer
    if isinstance(bars1, BarContainer):
        bar_patches1 = bars1.patches
    elif isinstance(bars1, list) and len(bars1) > 0 and isinstance(bars1[0], BarContainer):
        # If bars1 is a list of BarContainers (from broken axes), flatten to get all patches
        bar_patches1 = [patch for container in bars1 for patch in container.patches]
    else:
        bar_patches1 = bars1

    # **NEW: Color bars based on selection pressure, validation status, and type**
    for i, (bar, ratio, valid, gc, is_sig, ptype) in enumerate(zip(bar_patches1, dnds_ratios, is_valid, genetic_codes, is_significant, types)):
        if not valid:
            # Add hatching for invalid sequences
            bar.set_hatch('///')
            bar.set_alpha(0.5)
        
        # **NEW: Different styling for encompassing proteins**
        if ptype == 'EncompassingProtein':
            bar.set_edgecolor('black')
            bar.set_linewidth(3)  # Thicker border for encompassing proteins
            
        # Add asterisk for significant results
        if is_sig:
            try:
                ax1.text(i, ratio + max(dnds_ratios) * 0.02, '*', ha='center', va='bottom',
                        fontsize=16, fontweight='bold', color='black')
            except ValueError:
                # Skip text if coordinate falls in broken axis gap
                pass

        # Add text annotation for extreme values
        if ratio > 2:
            try:
                ax1.text(i, ratio + 0.1, f'{ratio:.2f}', ha='center', va='bottom', fontsize=8, fontweight='bold')
            except ValueError:
                # Skip text if coordinate falls in broken axis gap
                pass
    
    ax1.set_xticks(x_positions)
    ax1.set_xticklabels(labels, rotation=45, fontsize=10)
    ax1.set_ylabel('dN/dS Ratio', fontsize=12, fontweight='bold')
    ax1.set_title('', fontsize=14, fontweight='bold')
    
    if legend_outside:
        ax1.legend(loc='upper left', bbox_to_anchor=(1, 1), fancybox=True, fontsize=10)
    else:
        ax1.legend(loc='lower right', fancybox=True, fontsize=10, framealpha=0.5, edgecolor='black')
    ax1.grid(False)

    # **NEW: Only set ylim if not using broken axes and ylim is specified**
    if not use_broken_axes:
        if not ylim:
            ax1.set_ylim(0, max(dnds_ratios) * 1.1)
        else:
            ax1.set_ylim(ylim)
    
    # Add background coloring for different selection types (skip for broken axes)
    if not use_broken_axes:
        y_max = ax1.get_ylim()[1]
        ax1.axhspan(0, 1, alpha=0.1, color='blue', label='Purifying selection')
        ax1.axhspan(1, y_max, alpha=0.1, color='red', label='Positive selection')
    
    # Plot 2: Synonymous and Non-synonymous rates
    width = 0.35
    x_syn = x_positions - width/2
    x_nonsyn = x_positions + width/2
    
    bars2_syn = ax2.bar(x_syn, syn_rates, width, label='Synonymous rate',
                       color='lightblue', alpha=0.5, edgecolor='blue', linewidth=0.5,)
    bars2_nonsyn = ax2.bar(x_nonsyn, nonsyn_rates, width, label='Non-synonymous rate',
                          color='lightcoral', alpha=0.5, edgecolor='red', linewidth=0.5,)

    # **NEW: Extract individual bar patches (needed for broken axes compatibility)**
    if isinstance(bars2_syn, BarContainer):
        bar_patches2_syn = bars2_syn.patches
    elif isinstance(bars2_syn, list) and len(bars2_syn) > 0 and isinstance(bars2_syn[0], BarContainer):
        bar_patches2_syn = [patch for container in bars2_syn for patch in container.patches]
    else:
        bar_patches2_syn = bars2_syn

    if isinstance(bars2_nonsyn, BarContainer):
        bar_patches2_nonsyn = bars2_nonsyn.patches
    elif isinstance(bars2_nonsyn, list) and len(bars2_nonsyn) > 0 and isinstance(bars2_nonsyn[0], BarContainer):
        bar_patches2_nonsyn = [patch for container in bars2_nonsyn for patch in container.patches]
    else:
        bar_patches2_nonsyn = bars2_nonsyn

    # **NEW: Add special styling for encompassing proteins and invalid sequences**
    for i, (valid, ptype) in enumerate(zip(is_valid, types)):
        if not valid:
            bar_patches2_syn[i].set_hatch('///')
            bar_patches2_nonsyn[i].set_hatch('///')
            bar_patches2_syn[i].set_alpha(0.5)
            bar_patches2_nonsyn[i].set_alpha(0.5)

        if ptype == 'EncompassingProtein':
            bar_patches2_syn[i].set_edgecolor('blue')
            bar_patches2_syn[i].set_linewidth(3)
            bar_patches2_syn[i].set_alpha(0.4)
            bar_patches2_nonsyn[i].set_edgecolor('red')
            bar_patches2_nonsyn[i].set_linewidth(3)
            bar_patches2_nonsyn[i].set_alpha(0.4)
    
    ax2.set_xticks(x_positions)
    ax2.set_xticklabels(labels, rotation=45, fontsize=10)
    ax2.set_ylabel('Mutation Rate (per site)', fontsize=12, fontweight='bold')
    #ax2.set_title('Synonymous vs Non-synonymous Rates', fontsize=14, fontweight='bold')
    if legend_outside:
        legend1 = ax2.legend(loc='upper left', bbox_to_anchor=(1, 1), fancybox=True, fontsize=10)
    else:
        legend1 = ax2.legend(loc='upper left', fancybox=True, fontsize=10, framealpha=0.3, edgecolor='grey', prop={'weight': 'bold'})
    ax2.add_artist(legend1)
    ax2.grid(False)
    
    # Add gene legend if multiple genes
    if len(results_dict) > 1:
        gene_legend_elements = [plt.Rectangle((0,0),1,1, facecolor=gene_colors[gene], alpha=0.7, label=gene) 
                               for gene in results_dict.keys()]
        ax1.legend(handles=gene_legend_elements + [plt.Line2D([0], [0], color='red', linestyle='--', label='Neutral (dN/dS=1)')], 
                   loc='upper right', fontsize=10, prop={'weight': 'bold'})
    
    # **NEW: Add type legend**
    if include_encompassing_proteins and any(t == 'EncompassingProtein' for t in types):
        from matplotlib.patches import Rectangle
        
        type_legend_elements = [
            Rectangle((0, 0), 1, 1, facecolor='lightblue', alpha=0.4, edgecolor='blue', linewidth=1, label='Micro-proteins'),
            Rectangle((0, 0), 1, 1, facecolor='lightblue', alpha=0.4, edgecolor='blue', linewidth=3, label='Encompassing proteins')
        ]
        if any(not valid for valid in is_valid):
            type_legend_elements.append(Rectangle((0, 0), 1, 1, facecolor='gray', alpha=0.5, hatch='///', label='Invalid sequences'))
        
        # Create the new legend
        if legend_outside:
            new_legend = ax2.legend(handles=type_legend_elements, 
                             bbox_to_anchor=(1, 1), fontsize=10)
        else:
            new_legend = ax2.legend(handles=type_legend_elements, 
                            loc='upper right', fontsize=10, framealpha=0.3, edgecolor='grey', prop={'weight': 'bold'})
        
        # Add it as an artist (this preserves any existing legend)
        ax2.add_artist(new_legend)
    
    plt.tight_layout()
    
    # **NEW: Updated summary statistics to include encompassing proteins**
    summary_text = f"Total items: {len(valid_data)}\n"
    summary_text += f"Genes analyzed: {', '.join(results_dict.keys())}\n"
    
    # Separate counts by type
    mp_count = sum(1 for t in types if t == 'MicroProtein')
    enc_count = sum(1 for t in types if t == 'EncompassingProtein')
    summary_text += f"Micro-proteins: {mp_count}, Encompassing proteins: {enc_count}\n"
    
    # Calculate selection statistics
    positive_selection = sum(1 for ratio in dnds_ratios if ratio > 1)
    negative_selection = sum(1 for ratio in dnds_ratios if ratio < 1)
    neutral_selection = sum(1 for ratio in dnds_ratios if ratio == 1)
    
    summary_text += f"Positive selection (dN/dS > 1): {positive_selection}\n"
    summary_text += f"Negative selection (dN/dS < 1): {negative_selection}\n"
    summary_text += f"Neutral selection (dN/dS = 1): {neutral_selection}\n"
    
    # Add significance statistics
    significant_count = sum(1 for is_sig in is_significant if is_sig)
    summary_text += f"Statistically significant: {significant_count}/{len(valid_data)}\n"
    
    # Add genetic code summary
    mito_count = sum(1 for gc in genetic_codes if gc == 'mitochondrial')
    nuclear_count = sum(1 for gc in genetic_codes if gc == 'nuclear')
    summary_text += f"Mitochondrial code: {mito_count}, Nuclear code: {nuclear_count}"
    sns.despine()

    if display_summary_text:
        fig.text(0.02, 0.02, summary_text, fontsize=9, verticalalignment='bottom', 
             bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8))
    
    if output_file:
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        logger.info(f"Combined plot saved to {output_file}")
    else:
        plt.show()
        
def create_summary_table(results_dict: Dict[str, Dict]) -> pd.DataFrame:
    """
    Create a summary table from multiple analyze_microprotein_dnds results.
    
    Parameters:
    -----------
    results_dict : Dict[str, Dict]
        Dictionary mapping gene names to their analysis results.
    
    Returns:
    --------
    pd.DataFrame
        Summary table with all micro-protein results.
    """
    all_results = []
    
    for gene_name, results in results_dict.items():
        if 'microprotein_results' not in results:
            continue
            
        for mp_name, mp_results in results['microprotein_results'].items():
            row = {
                'Gene': gene_name,
                'MicroProtein': mp_name,
                'AbsoluteStart': mp_results['coordinates']['absolute_start'],
                'AbsoluteEnd': mp_results['coordinates']['absolute_end'],
                'Strand': mp_results['coordinates']['strand'],
                'GeneticCode': mp_results['sequence_info']['genetic_code'],
                'Length': mp_results['sequence_info']['length'],
                'AminoAcidSequence': mp_results['sequence_info']['amino_acid_sequence'],
                'UniqueVariants': mp_results['sequence_info']['unique_variants'],
                'dN_dS_Ratio': mp_results['rates']['dnds_ratio'],
                'SynonymousRate': mp_results['rates']['synonymous_rate'],
                'NonsynonymousRate': mp_results['rates']['nonsynonymous_rate'],
                'Statistic': mp_results['statistical_tests']['statistic'],
                'P_Value': mp_results['statistical_tests']['p_value'],
                'Is_Significant': mp_results['statistical_tests']['is_significant'],
                'Test_Method': mp_results['statistical_tests']['test_method'],
                'IsValid': mp_results['validation']['is_valid'],
                'ValidationIssues': '; '.join(mp_results['validation']['issues']),
                'CoordinateAdjusted': mp_results['coordinates'].get('coordinate_adjusted', False)
            }
            
            # Add encompassing protein data if available
            if 'encompassing_protein_results' in results and mp_name in results['encompassing_protein_results']:
                enc_results = results['encompassing_protein_results'][mp_name]
                row.update({
                    'EncompassingProtein': enc_results['name'],
                    'EncompassingProtein_Length': enc_results['length'],
                    'EncompassingProtein_dN_dS_Ratio': enc_results['dnds_ratio'],
                    'EncompassingProtein_Statistic': enc_results.get('z_statistic', enc_results.get('statistic', 0.0)),
                    'EncompassingProtein_P_Value': enc_results['p_value'],
                    'EncompassingProtein_Is_Significant': enc_results['is_significant'],
                    'EncompassingProtein_Test_Method': enc_results.get('test_method', 'Unknown')
                })
            
            all_results.append(row)
    
    return pd.DataFrame(all_results)

if __name__ == "__main__":
    # Run example
    microproteins, genetic_codes = example_usage_with_new_features()