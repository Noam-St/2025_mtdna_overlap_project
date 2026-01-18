#!/usr/bin/env python3
"""
Adaptive cleaning function for dN/dS analysis.
Handles ambiguous bases and gaps based on their frequency.
"""

from typing import List, Tuple, Set
from dataclasses import dataclass
import sys

# IUPAC ambiguity codes
AMBIGUOUS_CODES = set('RYKMSWBDHVN')


@dataclass
class CleaningResult:
    """Results from the cleaning process"""
    cleaned_sequences: List[str]
    removed_sequence_indices: Set[int]
    removed_codon_positions: Set[int]  # 0-based codon indices
    original_seq_count: int
    original_codon_count: int
    final_seq_count: int
    final_codon_count: int
    iterations: int


def is_problematic(base: str) -> bool:
    """
    Check if a base is problematic (ambiguous or gap).
    """
    base_upper = base.upper()
    return base_upper in AMBIGUOUS_CODES or base_upper == '-' or base_upper == '.'


def clean_sequences_adaptive(
    sequences: List[str],
    threshold: float = 0.05,
    max_iterations: int = 10,
    verbose: bool = False
) -> CleaningResult:
    """
    Adaptively clean sequences based on frequency of problematic sites.
    
    Logic:
    - If a codon position has problems in >threshold fraction of sequences:
      → Remove that codon position from ALL sequences (systematic issue)
    - If a codon position has problems in ≤threshold fraction of sequences:
      → Remove the problematic sequences (low-quality sequences)
    
    Iterates until no more changes occur or max_iterations reached.
    
    Args:
        sequences: List of aligned sequences (all same length, multiple of 3)
        threshold: Frequency threshold (default 0.05 = 5%)
        max_iterations: Maximum cleaning iterations (default 10)
        verbose: Print progress information
    
    Returns:
        CleaningResult object with cleaned sequences and statistics
    
    Example:
        >>> seqs = ["ATGCCC---GGG", "ATGCCCAAAGGG", "ATGCCCAAAGGG", "ATGNNNNNNGG"]
        >>> result = clean_sequences_adaptive(seqs, threshold=0.25)
        >>> print(result.cleaned_sequences)
        ['ATGCCCGGG', 'ATGCCCGGG', 'ATGCCCGGG']  # Removed bad codon and bad sequence
    """
    
    # Validation
    if not sequences:
        raise ValueError("Empty sequence list provided")
    
    if not all(len(seq) == len(sequences[0]) for seq in sequences):
        raise ValueError("All sequences must have the same length")
    
    seq_length = len(sequences[0])
    if seq_length % 3 != 0:
        raise ValueError(f"Sequence length ({seq_length}) must be divisible by 3")
    
    original_seq_count = len(sequences)
    original_codon_count = seq_length // 3
    
    # Track what we've removed
    removed_sequence_indices = set()
    removed_codon_positions = set()
    
    # Working copies
    active_sequences = list(sequences)
    active_indices = set(range(len(sequences)))
    
    iteration = 0
    changes_made = True
    
    while changes_made and iteration < max_iterations:
        iteration += 1
        changes_made = False
        
        if verbose:
            print(f"\nIteration {iteration}:")
            print(f"  Active sequences: {len(active_sequences)}")
            print(f"  Sequence length: {len(active_sequences[0])} bp "
                  f"({len(active_sequences[0])//3} codons)")
        
        # Analyze each codon position
        codon_count = len(active_sequences[0]) // 3
        codons_to_remove = []
        sequences_to_remove = []
        
        for codon_idx in range(codon_count):
            pos_start = codon_idx * 3
            pos_end = pos_start + 3
            
            # Check which sequences have problems in this codon
            problematic_seq_indices = []
            
            for seq_idx, seq in enumerate(active_sequences):
                codon = seq[pos_start:pos_end]
                if any(is_problematic(base) for base in codon):
                    problematic_seq_indices.append(seq_idx)
            
            # Calculate frequency of problems
            if problematic_seq_indices:
                problem_frequency = len(problematic_seq_indices) / len(active_sequences)
                
                if problem_frequency > threshold:
                    # Common problem → remove codon from all sequences
                    codons_to_remove.append(codon_idx)
                    changes_made = True
                    
                    if verbose:
                        print(f"  Codon {codon_idx}: {problem_frequency:.1%} problematic "
                              f"(>{threshold:.1%}) → removing position")
                else:
                    # Rare problem → remove problematic sequences
                    sequences_to_remove.extend(problematic_seq_indices)
                    changes_made = True
                    
                    if verbose:
                        print(f"  Codon {codon_idx}: {problem_frequency:.1%} problematic "
                              f"(≤{threshold:.1%}) → removing {len(problematic_seq_indices)} sequence(s)")
        
        # Remove problematic sequences
        if sequences_to_remove:
            sequences_to_remove = sorted(set(sequences_to_remove), reverse=True)
            
            for seq_idx in sequences_to_remove:
                # Track original index
                original_idx = list(active_indices)[seq_idx]
                removed_sequence_indices.add(original_idx)
                
                # Remove from active set
                del active_sequences[seq_idx]
                active_indices = set(range(len(active_sequences)))
            
            if verbose:
                print(f"  Removed {len(sequences_to_remove)} sequence(s)")
        
        # Remove problematic codon positions
        if codons_to_remove:
            codons_to_remove = sorted(set(codons_to_remove), reverse=True)
            
            for codon_idx in codons_to_remove:
                removed_codon_positions.add(codon_idx)
            
            # Rebuild sequences without removed codons
            kept_codon_indices = [i for i in range(codon_count) 
                                 if i not in codons_to_remove]
            
            new_sequences = []
            for seq in active_sequences:
                new_seq = ''.join(
                    seq[idx*3:(idx*3)+3] for idx in kept_codon_indices
                )
                new_sequences.append(new_seq)
            
            active_sequences = new_sequences
            
            if verbose:
                print(f"  Removed {len(codons_to_remove)} codon position(s)")
        
        # Check if we still have data
        if not active_sequences or len(active_sequences[0]) == 0:
            raise ValueError("All sequences or all positions were removed! "
                           "Try a higher threshold or less stringent filtering.")
    
    if iteration >= max_iterations and changes_made:
        print(f"Warning: Reached maximum iterations ({max_iterations}). "
              f"Cleaning may be incomplete.", file=sys.stderr)
    
    # Final statistics
    final_seq_count = len(active_sequences)
    final_codon_count = len(active_sequences[0]) // 3 if active_sequences else 0
    
    result = CleaningResult(
        cleaned_sequences=active_sequences,
        removed_sequence_indices=removed_sequence_indices,
        removed_codon_positions=removed_codon_positions,
        original_seq_count=original_seq_count,
        original_codon_count=original_codon_count,
        final_seq_count=final_seq_count,
        final_codon_count=final_codon_count,
        iterations=iteration
    )
    
    if verbose:
        print_cleaning_summary(result)
    
    return result


def print_cleaning_summary(result: CleaningResult):
    """Print a summary of the cleaning results"""
    print(f"\n{'='*60}")
    print(f"ADAPTIVE CLEANING SUMMARY")
    print(f"{'='*60}")
    print(f"Iterations:                   {result.iterations:>10}")
    print(f"")
    print(f"Original sequences:           {result.original_seq_count:>10}")
    print(f"Final sequences:              {result.final_seq_count:>10}")
    print(f"Removed sequences:            {len(result.removed_sequence_indices):>10}")
    print(f"")
    print(f"Original codons:              {result.original_codon_count:>10}")
    print(f"Final codons:                 {result.final_codon_count:>10}")
    print(f"Removed codon positions:      {len(result.removed_codon_positions):>10}")
    print(f"")
    
    if result.final_seq_count > 0:
        seq_retention = (result.final_seq_count / result.original_seq_count) * 100
        codon_retention = (result.final_codon_count / result.original_codon_count) * 100
        print(f"Sequence retention:           {seq_retention:>9.1f}%")
        print(f"Codon retention:              {codon_retention:>9.1f}%")
        print(f"")
        print(f"Final alignment:              {result.final_seq_count} sequences × "
              f"{result.final_codon_count} codons")
        print(f"                              = {result.final_seq_count * result.final_codon_count} "
              f"total codons for analysis")
    
    print(f"{'='*60}")


# Example usage and testing
if __name__ == "__main__":
    import sys
    
    # Test case 1: Mix of common and rare problems
    print("Test 1: Mixed problems")
    print("-" * 60)
    test_sequences_1 = [
        "ATGCCCAAAGGG",  # Perfect
        "ATGCCCAAAGGG",  # Perfect
        "ATGCCCAAAGGG",  # Perfect
        "ATGCCCAAAGGG",  # Perfect
        "ATGNNNAAAGGG",  # Rare problem (1/6 = 16.7%)
        "ATGCCC---GGG",  # Common problem position (6/6 = 100%)
    ]
    
    result1 = clean_sequences_adaptive(test_sequences_1, threshold=0.20, verbose=True)
    print("\nCleaned sequences:")
    for i, seq in enumerate(result1.cleaned_sequences):
        print(f"  Seq {i}: {seq}")
    
    
    # Test case 2: Rare problems only
    print("\n\n\nTest 2: Only rare problems")
    print("-" * 60)
    test_sequences_2 = [
        "ATGCCCAAAGGGTTTTTT",  # Perfect
        "ATGCCCAAAGGGTTTTTT",  # Perfect
        "ATGCCCAAAGGGTTTTTT",  # Perfect
        "ATGCCCAAAGGGTTTTTT",  # Perfect
        "ATGNNNAAAGGGTTTTTT",  # Problem in codon 1 (1/10 = 10%)
        "ATGCCCNNNGGGTTTTTT",  # Problem in codon 2 (1/10 = 10%)
        "ATGCCCAAAGGGNNNNNN",  # Problem in codon 5 (1/10 = 10%)
        "ATGCCCAAAGGGTTTTTT",  # Perfect
        "ATGCCCAAAGGGTTTTTT",  # Perfect
        "ATGCCCAAAGGGTTTTTT",  # Perfect
    ]
    
    result2 = clean_sequences_adaptive(test_sequences_2, threshold=0.15, verbose=True)
    print("\nCleaned sequences:")
    for i, seq in enumerate(result2.cleaned_sequences):
        print(f"  Seq {i}: {seq}")
    
    
    # Test case 3: Common problems
    print("\n\n\nTest 3: Common problems")
    print("-" * 60)
    test_sequences_3 = [
        "ATGCCCAAAGGG",
        "ATGCCCAAAGGG",
        "ATG---AAAGGG",  # 50% have gap at codon 2
        "ATG---AAAGGG",
        "ATGCCCAAA---",  # 50% have gap at codon 4
        "ATGCCCAAA---",
    ]
    
    result3 = clean_sequences_adaptive(test_sequences_3, threshold=0.40, verbose=True)
    print("\nCleaned sequences:")
    for i, seq in enumerate(result3.cleaned_sequences):
        print(f"  Seq {i}: {seq}")