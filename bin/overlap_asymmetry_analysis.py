"""
Step 2: Overlapping Reading Frame Asymmetry Analysis

This script provides functions to analyze asymmetric evolution in overlapping 
reading frames, following Szklarczyk et al. 2007 methodology.

Key Concept:
- In overlapping regions, the SAME DNA encodes TWO different proteins
- A mutation that is synonymous in one frame may be nonsynonymous in the other
- Asymmetric evolution (different dN/dS ratios) suggests both frames are functional

Usage:
    results = analyze_overlap_asymmetry(
        df=your_dataframe,
        canonical_gene='ATP8',
        canonical_seq_col='ATP8',
        alternative_gene='ATP6', 
        alternative_seq_col='ATP8_ATP6',  # or use 'ATP6' if available
        overlap_start=161,
        overlap_end=207,
        canonical_genetic_code='mitochondrial',
        alternative_genetic_code='mitochondrial'
    )

Author: Analysis for mitochondrial microproteins
"""

import pandas as pd
import numpy as np
from typing import Dict, List, Tuple, Optional
from scipy import stats
import matplotlib.pyplot as plt
import seaborn as sns
from collections import Counter
import warnings
import clean_amb_and_gaps as cag
from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor, as_completed
from functools import partial
from statsmodels.stats.contingency_tables import mcnemar

# Import the existing dN/dS analysis functions
from enhanced_mtdna_dnds_analysis_biopython import (
    calculate_synonymous_nonsynonymous,
    analyze_protein_dnds,
    translate_sequence,
    MITOCHONDRIAL_GENETIC_CODE,
    NUCLEAR_GENETIC_CODE
)

from importlib import reload
import overlap_sub_model
reload(overlap_sub_model)
from overlap_sub_model import ARFomeNormalized

# Set style for plots
sns.set_style("whitegrid")
plt.rcParams['figure.figsize'] = (12, 8)

# Define nucleotide transition/transversion categories
TRANSITIONS = {('A', 'G'), ('G', 'A'), ('C', 'T'), ('T', 'C')}
TRANSVERSIONS = {('A', 'C'), ('C', 'A'), ('A', 'T'), ('T', 'A'), 
                 ('G', 'C'), ('C', 'G'), ('G', 'T'), ('T', 'G')}


def infer_frameshift(canonical_cds: str, alternative_cds: str, return_start_index: bool = False) -> str:
    """
    Infers the frameshift of an alternative CDS relative to a canonical CDS.
    
    The function attempts to locate the start of the alternative CDS within 
    the canonical CDS. The remainder of the start index divided by 3 
    determines the shift.
    
    Args:
        canonical_cds (str): The reference nucleotide sequence.
        alternative_cds (str): The variant nucleotide sequence.
        return_start_index (bool): Whether to return the start index of the alternative CDS in the canonical CDS.
        
    Returns:
        str: "+1", "+2", "0" (in-frame), or "No match found".
    """
    
    # 1. Normalize inputs (uppercase and remove whitespace)
    seq_can = canonical_cds.strip().upper()
    seq_alt = alternative_cds.strip().upper()
    
    # 2. Define a 'seed' length. 
    # We use the start of the alternative CDS to anchor the alignment.
    # We don't match the whole string to allow for downstream splicing differences.
    # 15 nucleotides (5 codons) is usually sufficient for uniqueness without being too strict.
    seed_length = 15
    
    if len(seq_alt) < seed_length:
        # If the alternative is tiny, use its whole length
        seed = seq_alt
    else:
        seed = seq_alt[:seed_length]
        
    # 3. Find the position of the alternative start within the canonical sequence
    start_index = seq_can.find(seed)
    result = "No match found"
    # 4. Infer logic
    frame_offset = start_index % 3
    
    if frame_offset == 1:
        result = 1
    elif frame_offset == 2:
        result = 2
    else:
        result = 0
    
    if return_start_index:
        return result, start_index
    else:
        return result

def classify_substitution(old_base: str, new_base: str) -> str:
    """
    Classify a nucleotide substitution as transition or transversion.
    
    Parameters:
    -----------
    old_base : str
        Original nucleotide
    new_base : str
        Substituted nucleotide
        
    Returns:
    --------
    str : 'transition', 'transversion', or 'none'
    """
    if old_base == new_base:
        return 'none'
    
    pair = (old_base.upper(), new_base.upper())
    if pair in TRANSITIONS:
        return 'transition'
    elif pair in TRANSVERSIONS:
        return 'transversion'
    else:
        return 'unknown'


def analyze_overlap_asymmetry(df: pd.DataFrame,
                              canonical_gene: str,
                              canonical_seq_col: str,
                              alternative_gene: str,
                              alternative_seq_col: str,
                              overlap_start: Optional[int] = None,
                              overlap_end: Optional[int] = None,
                              canonical_genetic_code: str = 'mitochondrial',
                              alternative_genetic_code: str = 'mitochondrial',
                              dnds_method: str = 'NG86',
                              use_full_sequences: bool = False,
                              frame_shift: str = 'infer',
                              include_dual_coding_test: bool = True,
                              alternative_strand: str = '+',
                              verbose: bool = True,
                              threshold = 0.05,
                              ref_id: Optional[str] = None) -> Dict:
    """
    Analyze asymmetric evolution in overlapping reading frames using ARFomeNormalized model.
    
    This function:
    1. Prepares sequences (using full canonical sequence for context).
    2. Uses ARFomeNormalized to calculate substitution rates and betaSTOP.
    3. Derives dN/dS metrics for compatibility with existing plots.
    4. Tests for significant asymmetry.
    
    Parameters:
    -----------
    df : pd.DataFrame
        DataFrame with sequence data
    canonical_gene : str
        Name of canonical/main gene
    canonical_seq_col : str
        Column name containing canonical gene sequences
    alternative_gene : str
        Name of alternative gene/microprotein
    alternative_seq_col : str
        Column name containing alternative gene sequences (used for initial check/inference)
    overlap_start : int, optional
        Start position of overlap in canonical gene (1-based)
    overlap_end : int, optional
        End position of overlap in canonical gene (1-based)
    canonical_genetic_code : str
        Genetic code for canonical gene
    alternative_genetic_code : str
        Genetic code for alternative gene
    dnds_method : str
        Method for dN/dS calculation (kept for API compatibility, but analysis uses new model)
    use_full_sequences : bool
        If True, use full sequences instead of overlap region
    frame_shift : str or int (default=infer)
        Reading frame shift between canonical and alternative (1 or 2)
    include_dual_coding_test : bool
        If True, perform betaSTOP analysis (always done in new model)
    alternative_strand : str
        Strand of the alternative gene relative to canonical ('+' or '-'). Default is '+'.
    verbose : bool
        Print progress messages
    threshold : float
        Threshold to remove the sequence from the analysis (passed to cag)
    

    Returns:
    --------
    Dict
        Comprehensive results including dN/dS, asymmetry metrics, and dual-coding evidence.
    """
    if verbose:
        print("=" * 80)
        print(f"ASYMMETRY ANALYSIS (ARFomeNormalized): {canonical_gene} vs {alternative_gene}")
        print(f"Alternative Strand: {alternative_strand}")
        print("=" * 80)
    
    # Determine frame shift if not provided

    if frame_shift not in [1, 2]:
        if verbose:
            print("Inferring frame shift between canonical and alternative sequences...")
        sample_can_seq = df[canonical_seq_col].dropna().iloc[0]
        sample_alt_seq = df[alternative_seq_col].dropna().iloc[0]
        
        # If opposite strand, we need to handle that in inference
        # infer_frameshift expects sequences that match.
        # If alt is opposite, we reverse complement it to match canonical strand for alignment
        seq_to_infer = sample_alt_seq
        if alternative_strand == '-':
            from Bio.Seq import Seq
            seq_to_infer = str(Seq(sample_alt_seq).reverse_complement())
        
        inferred_shift, start_index = infer_frameshift(sample_can_seq, seq_to_infer, return_start_index=True)
        if inferred_shift == "No match found":
            warnings.warn("Could not infer frame shift; defaulting to +1")
            frame_shift = 1
        elif frame_shift == 'infer':
            frame_shift = inferred_shift
        if verbose:
            print(f"Inferred frame shift: +{frame_shift}")
            print(f"Inferred start index: {start_index}")
    
    # Prepare sequences
    # We use the full canonical sequence to provide context for the new model
    canonical_seqs_full = df[canonical_seq_col].dropna().copy()
    
    # Clean up the sequences using cag
    cleanup_canonical = cag.clean_sequences_adaptive(canonical_seqs_full.to_list(), threshold = threshold)
    canonical_seqs_full = pd.Series(cleanup_canonical.cleaned_sequences)
    
    # Determine ARF start index and extract alternative sequences
    if use_full_sequences:
        if verbose:
            print("\nUsing full sequences (no overlap region specified)")
        arf_start_index = 0
        # If using full sequences, we assume canonical and alternative are same length/aligned
        alternative_seqs = canonical_seqs_full.copy()
        overlap_len = len(canonical_seqs_full.iloc[0])

    elif overlap_start is None and overlap_end is None:
        # If no overlap region specified, infer it from canonical sequences
        if verbose:
            print("\nInferring overlap region from canonical sequences")
            arf_start_index = start_index
            overlap_len = len(sample_alt_seq)
            overlap_end = arf_start_index + overlap_len
            overlap_start = arf_start_index + 1

    if overlap_start is not None and overlap_end is not None and not use_full_sequences:
        if verbose:
            print(f"\nExtracting overlap region: positions {overlap_start}-{overlap_end}")
        
        # Convert 1-based to 0-based index
        arf_start_index = overlap_start - 1
        
        # Extract alternative sequences from canonical using the overlap coordinates
        # This ensures the alternative sequence maps exactly to the canonical one
        # and provides the correct context for the model
        def extract_alt(seq):
            if pd.isna(seq) or len(seq) < overlap_end:
                return ""
            return seq[arf_start_index:overlap_end]
            
        alternative_seqs = canonical_seqs_full.apply(extract_alt)
        
        # If opposite strand, the extracted sequence from canonical needs to be reverse complemented
        # to represent the alternative gene sequence
        if alternative_strand == '-':
            from Bio.Seq import Seq
            alternative_seqs = alternative_seqs.apply(lambda s: str(Seq(s).reverse_complement()) if s else s)
            
        overlap_len = overlap_end - overlap_start + 1

    # Filter valid sequences
    valid_idx = (canonical_seqs_full.str.len() > 0) & (alternative_seqs.str.len() > 0)
    canonical_seqs_full = canonical_seqs_full[valid_idx].tolist()
    alternative_seqs = alternative_seqs[valid_idx].tolist()
    
    if len(canonical_seqs_full) == 0:
        warnings.warn("No valid sequences found for analysis")
        return {'error': 'No valid sequences'}

    # Map genetic code to ID
    # 'mitochondrial' -> 2 (Vertebrate Mitochondrial)
    # 'standard'/'nuclear' -> 1 (Standard)
    genetic_code_id = 2
    if canonical_genetic_code in ['standard', 'nuclear'] or alternative_genetic_code in ['standard', 'nuclear']:
        genetic_code_id = 1
    
    if verbose:
        print(f"Using Genetic Code ID: {genetic_code_id}")
        print(f"ARF Start Index: {arf_start_index}")
        print(f"Frame Shift: +{frame_shift}")

    # Report the reference sequence if it is set
    ref_seq = None
    if ref_id is not None:
        ref_seq = df[df['ID'] == ref_id][canonical_seq_col].values[0]
        if pd.isna(ref_seq) or ref_seq == '':
            ref_seq = None
        else:
            if verbose:
                print(f"Using reference sequence ID: {ref_id}")
    if ref_seq is None:
        ref_seq = canonical_seqs_full[0]
        if verbose:
            print("Using first canonical sequence as reference")

    # Run ARFomeNormalized Analysis
    arf_model = ARFomeNormalized(genetic_code_id=genetic_code_id)
    model_results = arf_model.run_analysis(
        canonical_seqs=canonical_seqs_full,
        alternative_seqs=alternative_seqs,
        arf_start_index=arf_start_index,
        frame_shift=frame_shift,
        alternative_strand=alternative_strand,
        verbose=verbose,
        ref_seq = ref_seq
    )

    # Extract results
    rates = model_results['raw_rates']
    observed = model_results['observed']
    potential = model_results['potential']
    beta_stop = model_results['beta_stop']
    n_sequences = len(canonical_seqs_full)  # Track number of sequences for frequency calculations
    
    # Calculate dN/dS metrics using Biopython's cal_dn_ds with NG86 method
    # Following the approach from enhanced_mtdna_dnds_analysis_biopython

    def safe_div(n, d):
        return n / d if d > 0 else 0.0

    def calculate_frame_dnds(sequences, genetic_code, method='NG86'):
        """
        Calculate dN/dS for a set of sequences using Biopython's cal_dn_ds.

        Parameters:
        -----------
        sequences : list
            List of DNA sequences
        genetic_code : str
            'mitochondrial' or 'nuclear'
        method : str
            Biopython method: 'NG86', 'LWL85', 'ML', 'YN00'

        Returns:
        --------
        tuple: (dN, dS, omega, obs_syn, obs_nonsyn)
        """
        from collections import Counter

        # Use most common sequence as reference
        variant_counts = Counter(sequences)
        reference_seq = variant_counts.most_common(1)[0][0]
        total_sequences = len(sequences)

        synonymous_diffs = []
        nonsynonymous_diffs = []
        synonymous_sites_list = []
        nonsynonymous_sites_list = []

        for variant_seq, count in variant_counts.items():
            if variant_seq == reference_seq:
                continue  # Skip reference sequence

            frequency = count / total_sequences

            # Calculate differences using Biopython
            diff_results = calculate_synonymous_nonsynonymous(
                reference_seq, variant_seq,
                strand='+',  # Strand is already handled in sequence preparation
                genetic_code=genetic_code,
                method=method
            )

            # Weight by frequency
            synonymous_diffs.append(diff_results['synonymous_differences'] * frequency)
            nonsynonymous_diffs.append(diff_results['nonsynonymous_differences'] * frequency)
            synonymous_sites_list.append(diff_results['synonymous_sites'] * frequency)
            nonsynonymous_sites_list.append(diff_results['nonsynonymous_sites'] * frequency)

        # Calculate summary statistics
        total_syn_diffs = sum(synonymous_diffs) if synonymous_diffs else 0
        total_nonsyn_diffs = sum(nonsynonymous_diffs) if nonsynonymous_diffs else 0
        total_syn_sites = sum(synonymous_sites_list) if synonymous_sites_list else 1
        total_nonsyn_sites = sum(nonsynonymous_sites_list) if nonsynonymous_sites_list else 1

        # Calculate rates
        dS = total_syn_diffs / total_syn_sites if total_syn_sites > 0 else 0
        dN = total_nonsyn_diffs / total_nonsyn_sites if total_nonsyn_sites > 0 else 0
        omega = safe_div(dN, dS)

        return dN, dS, omega, total_syn_diffs, total_nonsyn_diffs

    # Canonical
    obs_syn_can = observed.get('SS', 0) + observed.get('SN', 0)
    pot_syn_can = potential.get('SS', 0) + potential.get('SN', 0)
    obs_nonsyn_can = observed.get('NS', 0) + observed.get('NN', 0)
    pot_nonsyn_can = potential.get('NS', 0) + potential.get('NN', 0)
    
    dS_can = safe_div(obs_syn_can, pot_syn_can)
    dN_can = safe_div(obs_nonsyn_can, pot_nonsyn_can)
    omega_can = safe_div(dN_can, dS_can)
    
    # Alternative
    obs_syn_alt = observed.get('SS', 0) + observed.get('NS', 0)
    pot_syn_alt = potential.get('SS', 0) + potential.get('NS', 0)
    obs_nonsyn_alt = observed.get('SN', 0) + observed.get('NN', 0)
    pot_nonsyn_alt = potential.get('SN', 0) + potential.get('NN', 0)
    
    dS_alt = safe_div(obs_syn_alt, pot_syn_alt)
    dN_alt = safe_div(obs_nonsyn_alt, pot_nonsyn_alt)
    omega_alt = safe_div(dN_alt, dS_alt)
    
    # Asymmetry Ratio
    asymmetry_ratio = safe_div(dN_alt, dN_can)
    omega_ratio = safe_div(omega_alt, omega_can)
    
    # McNemar's Test for Paired Categorical Data
    # Tests for asymmetry in mutation effects between overlapping frames
    # Table structure: [[SS, SN], [NS, NN]]
    # where rows = canonical frame, columns = alternative frame
    # SS = Synonymous in both, SN = Syn in canonical/Nonsyn in alt
    # NS = Nonsyn in canonical/Syn in alt, NN = Nonsyn in both

    SS = observed.get('SS', 0)  # Synonymous in both frames
    SN = observed.get('SN', 0)  # Syn in canonical, Nonsyn in alternative
    NS = observed.get('NS', 0)  # Nonsyn in canonical, Syn in alternative
    NN = observed.get('NN', 0)  # Nonsynonymous in both frames

    # Report the McNemar input if verbose
    if verbose:
        print("\nMcNemar's Test Contingency Table (Paired Data):")
        print("                Alternative Frame")
        print("                Syn      Nonsyn")
        print(f"Canonical Syn   {SS:<8} {SN:<8}")
        print(f"         Nonsyn {NS:<8} {NN:<8}")
        print(f"\nOff-diagonal counts (test focus): SN={SN}, NS={NS}")

    mcnemar_table = [
        [SS, SN],
        [NS, NN]
    ]

    # McNemar's test with continuity correction (using statsmodels)
    result = mcnemar(mcnemar_table, exact=False, correction=True)

    fisher_results = {
        'p_value': result.pvalue,
        'significant': result.pvalue < 0.05,
        'statistic': result.statistic
    }
    

# -------------------------------------------------------------------------
    # Dual Coding Evidence: Beta STOP Analysis
    # -------------------------------------------------------------------------
    # Rationale:
    # beta_stop measures the "acceptance rate" of Stop codons in the alternative frame.
    # beta_stop = (Rate of STOP mutations) / (Rate of Neutral Baseline mutations)
    #
    #   - If beta_stop ~ 1.0: Stops accumulate at the background mutation rate (Neutral).
    #   - If beta_stop << 1.0: Stops are being actively purged (Purifying Selection).
    #
    # Statistical Test (Poisson Depletion):
    # We test if the observed number of STOP mutations is significantly LOWER
    # than expected by chance, given the sequence context and mutation rate.
    #
    #   H0 (Null): Stops occur randomly at the baseline neutral rate.
    #   p-value:   Probability of observing k or fewer stops under H0.
    #
    # Interpretation:
    #   - p < 0.05: Significant depletion of Stop codons. Strong evidence FOR dual coding.
    #   - p > 0.05: Observed stops are consistent with random noise. No evidence.
    # -------------------------------------------------------------------------
    
    n_stop_contexts = int(observed.get('STOP_F1', 0))
    baseline_rate = rates[model_results['baseline_used']]
    
    # Calculate Expected Stops under Null Hypothesis (Neutral Evolution)
    # Expected = (Total Sites where Stops COULD occur) * (Background Mutation Rate)
    expected_stops = potential.get('STOP_F1', 0) * baseline_rate
    
    # Calculate one-sided p-value for depletion (Observation <= Expected)
    if expected_stops > 0:
        p_val_beta = stats.poisson.cdf(n_stop_contexts, expected_stops)
    else:
        p_val_beta = 1.0
        
    dual_coding_results = {
        'summary': {
            'betaSTOP': beta_stop,
            'p_value': p_val_beta, # This is p-value for depletion of stops
            'Ti/Tv': model_results['titv_used'],
            'stop_contexts': n_stop_contexts,
            'total_contexts': sum(potential.values()),
            'baseline_used': model_results['baseline_used']
        },
        'dual_coding_evidence': 'HIGH CONFIDENCE' if beta_stop < 0.5 else 'LOW CONFIDENCE', # Simple heuristic
        'evidence_reasons': [f"betaSTOP = {beta_stop:.4f}"]
    }

    # Compile results
    results = {
        'region_name': f"{canonical_gene}_{alternative_gene}",
        'canonical_gene': canonical_gene,
        'alternative_gene': alternative_gene,
        'overlap_coords': {
            'start': overlap_start,
            'end': overlap_end,
            'length': overlap_len,
            'frame_shift': frame_shift
        },
        'canonical_frame': {
            'dN': dN_can,
            'dS': dS_can,
            'omega': omega_can,
            'N_subs': obs_nonsyn_can,
            'S_subs': obs_syn_can,
            'genetic_code': canonical_genetic_code
        },
        'alternative_frame': {
            'dN': dN_alt,
            'dS': dS_alt,
            'omega': omega_alt,
            'N_subs': obs_nonsyn_alt,
            'S_subs': obs_syn_alt,
            'genetic_code': alternative_genetic_code
        },
        'asymmetry': {
            'ratio_dN_alt_vs_can': asymmetry_ratio,
            'ratio_omega_alt_vs_can': omega_ratio,
            'fisher_test': fisher_results
        },
        'dual_coding_evidence': dual_coding_results,
        'interpretation': {
            'faster_evolving_frame': alternative_gene if asymmetry_ratio > 1 else canonical_gene,
            'asymmetry_detected': fisher_results['significant'],
            'evidence_for_dual_function': fisher_results['significant'] or (beta_stop < 0.5)
        },
        'arf_results': {
            'rates': rates,
            'potential': potential,
            'observed': observed,
            'n_sequences': n_sequences,
            'mutation_frequencies': model_results.get('mutation_frequencies', {})
        }  # Add raw results from ARFomeNormalized model
    }
    
    if verbose:
        print_asymmetry_summary(results)
        
    return results


def print_asymmetry_summary(results: Dict):
    """Print a formatted summary of asymmetry analysis results."""
    print("\n" + "=" * 80)
    print("RESULTS SUMMARY")
    print("=" * 80)
    
    can = results['canonical_frame']
    alt = results['alternative_frame']
    asym = results['asymmetry']
    
    print(f"\nRegion: {results['region_name']}")
    if results['overlap_coords']['length']:
        print(f"Overlap length: {results['overlap_coords']['length']} bp "
              f"({results['overlap_coords']['length']//3} codons)")
    
    print(f"\n{'Frame':<20} {'dN':>10} {'dS':>10} {'omega (dN/dS)':>12} {'N_subs':>8} {'S_subs':>8}")
    print("-" * 80)
    print(f"{results['canonical_gene']:<20} {can['dN']:>10.4f} {can['dS']:>10.4f} "
          f"{can['omega']:>12.4f} {can['N_subs']:>8.1f} {can['S_subs']:>8.1f}")
    print(f"{results['alternative_gene']:<20} {alt['dN']:>10.4f} {alt['dS']:>10.4f} "
          f"{alt['omega']:>12.4f} {alt['N_subs']:>8.1f} {alt['S_subs']:>8.1f}")
    
    print(f"\n{'Asymmetry Metrics':<40} {'Value':>15}")
    print("-" * 80) 
    print(f"{'dN ratio (alt/can):':<40} {asym['ratio_dN_alt_vs_can']:>15.2f}")
    print(f"{'omega ratio (alt/can):':<40} {asym['ratio_omega_alt_vs_can']:>15.2f}")
    print(f"{'Fisher exact test p-value:':<40} {asym['fisher_test']['p_value']:>15.4f}")
    print(f"{'Significant asymmetry (p<0.05):':<40} {asym['fisher_test']['significant']:>15}")
    
    # NEW: Print dual-coding evidence
    if 'dual_coding_evidence' in results and results['dual_coding_evidence']:
        dc = results['dual_coding_evidence']
        if 'error' not in dc:
            print(f"\n{'Dual-Coding Evidence (ARFomeNormalized)':<40}")
            print("-" * 80)
            if 'summary' in dc:
                print(f"{'betaSTOP (Obs/Baseline):':<40} {dc['summary']['betaSTOP']:>15.6f}")
                print(f"{'p-value (Depletion):':<40} {dc['summary']['p_value']:>15.4f}")
                print(f"{'Baseline Used:':<40} {dc['summary']['baseline_used']:>15}")
                print(f"{'Stop-introducing contexts:':<40} {dc['summary']['stop_contexts']:>15}")
                print(f"{'Total contexts analyzed:':<40} {dc['summary']['total_contexts']:>15}")
            
            if 'dual_coding_evidence' in dc:
                print(f"\n{'Evidence Level:':<40} {dc['dual_coding_evidence']}")
            
            if 'evidence_reasons' in dc and dc['evidence_reasons']:
                print(f"\nEvidence based on:")
                for reason in dc['evidence_reasons']:
                    print(f"  • {reason}")
    
    print(f"\n{'Interpretation':<40}")
    print("-" * 80)
    print(f"Faster evolving frame: {results['interpretation']['faster_evolving_frame']}")
    print(f"Evidence for dual functionality: "
          f"{'YES' if results['interpretation']['evidence_for_dual_function'] else 'NO'}")
    
    if asym['fisher_test']['significant']:
        print("\n⚠️  SIGNIFICANT ASYMMETRY DETECTED!")
        print("This suggests both reading frames are under selection,")
        print("providing evidence for dual coding functionality.")
    
    # NEW: Highlight betaSTOP results
    if 'dual_coding_evidence' in results and results['dual_coding_evidence']:
        dc = results['dual_coding_evidence']
        if 'error' not in dc and 'summary' in dc:
            beta = dc['summary']['betaSTOP']
            if beta < 0.5:
                print("\n✓ LOW betaSTOP DETECTED!")
                print(f"betaSTOP = {beta:.4f} (normalized rate of stop introduction)")
                print("This is STRONG EVIDENCE for functional alternative reading frame.")
                print("Stop-introducing substitutions are strongly selected against.")


def compare_multiple_overlaps(df: pd.DataFrame,
                              overlap_configs: List[Dict],
                              dnds_method: str = 'NG86',
                              include_dual_coding_test: bool = True) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    Analyze multiple overlapping regions and compile results.

    Parameters:
    -----------
    df : pd.DataFrame
        DataFrame with sequence data
    overlap_configs : List[Dict]
        List of configuration dictionaries for each overlap
        Each dict should contain: canonical_gene, canonical_seq_col,
        alternative_gene, alternative_seq_col, overlap_start, overlap_end
    dnds_method : str
        Method for dN/dS calculation
    include_dual_coding_test : bool
        If True, perform betaSTOP and codon substitution analysis

    Returns:
    --------
    Tuple[pd.DataFrame, pd.DataFrame]
        - Summary DataFrame with one row per region
        - Detailed mutations DataFrame with one row per individual mutation
    """
    all_results = []
    all_mutations = []  # Track individual mutations across all regions
    
    for config in overlap_configs:
        print(f"\n{'='*80}")
        print(f"Analyzing: {config['canonical_gene']} vs {config['alternative_gene']}")
        print(f"{'='*80}")
        
        # Add frame_shift if not specified
        if 'frame_shift' not in config:
            config['frame_shift'] = 'infer'
        
        results = analyze_overlap_asymmetry(
            df=df,
            dnds_method=dnds_method,
            include_dual_coding_test=include_dual_coding_test,
            **config
        )
        
        # Extract key metrics for table
        row = {
            'Region': results['region_name'],
            'Canonical_Frame': results['canonical_gene'],
            'Alternative_Frame': results['alternative_gene'],
            'Overlap_Length_bp': results['overlap_coords']['length'],
            'Frame_Shift': results['overlap_coords'].get('frame_shift', 1),
            'dN_canonical': results['canonical_frame']['dN'],
            'dS_canonical': results['canonical_frame']['dS'],
            'omega_canonical': results['canonical_frame']['omega'],
            'dN_alternative': results['alternative_frame']['dN'],
            'dS_alternative': results['alternative_frame']['dS'],
            'omega_alternative': results['alternative_frame']['omega'],
            'Asymmetry_Ratio_dN': results['asymmetry']['ratio_dN_alt_vs_can'],
            'Asymmetry_Ratio_omega': results['asymmetry']['ratio_omega_alt_vs_can'],
            'Fisher_P_value': results['asymmetry']['fisher_test']['p_value'],
            'Significant_Asymmetry': results['asymmetry']['fisher_test']['significant'],
            'N_subs_canonical': results['canonical_frame']['N_subs'],
            'S_subs_canonical': results['canonical_frame']['S_subs'],
            'N_subs_alternative': results['alternative_frame']['N_subs'],
            'S_subs_alternative': results['alternative_frame']['S_subs']
        }

        # Add all results from the ARFomeNormalized model
        if 'arf_results' in results:
            arf_results = results['arf_results']
            rates = arf_results.get('rates', {})
            observed = arf_results.get('observed', {})
            potential = arf_results.get('potential', {})
            n_sequences = arf_results.get('n_sequences', np.nan)
            mutation_frequencies = arf_results.get('mutation_frequencies', {})

            # Calculate mutation frequency statistics for each category
            def calc_freq_stats(freqs_list):
                if not freqs_list:
                    return {'n_mutations': 0, 'mean_freq': 0, 'max_freq': 0, 'total_affected': 0}
                frequencies = [m['frequency_pct'] for m in freqs_list]
                counts = [m['count'] for m in freqs_list]
                return {
                    'n_mutations': len(freqs_list),
                    'mean_freq': np.mean(frequencies),
                    'max_freq': np.max(frequencies),
                    'total_affected': sum(counts)  # Total individuals with any mutation in this category
                }

            ss_stats = calc_freq_stats(mutation_frequencies.get('SS', []))
            sn_stats = calc_freq_stats(mutation_frequencies.get('SN', []))
            ns_stats = calc_freq_stats(mutation_frequencies.get('NS', []))
            nn_stats = calc_freq_stats(mutation_frequencies.get('NN', []))
            stop_stats = calc_freq_stats(mutation_frequencies.get('STOP_F1', []))

            row.update({
                'N_Sequences': n_sequences,
                'Rate_SS': rates.get('SS', np.nan),
                'Rate_SN': rates.get('SN', np.nan),
                'Rate_NS': rates.get('NS', np.nan),
                'Rate_NN': rates.get('NN', np.nan),
                'Rate_STOP_F1': rates.get('STOP_F1', np.nan),
                'Observed_SS': observed.get('SS', np.nan),
                'Observed_SN': observed.get('SN', np.nan),
                'Observed_NS': observed.get('NS', np.nan),
                'Observed_NN': observed.get('NN', np.nan),
                'Observed_STOP_F1': observed.get('STOP_F1', np.nan),
                'Potential_SS': potential.get('SS', np.nan),
                'Potential_SN': potential.get('SN', np.nan),
                'Potential_NS': potential.get('NS', np.nan),
                'Potential_NN': potential.get('NN', np.nan),
                'Potential_STOP_F1': potential.get('STOP_F1', np.nan),
                # Mutation frequency statistics
                'SS_N_Mutations': ss_stats['n_mutations'],
                'SS_Mean_Freq_Pct': ss_stats['mean_freq'],
                'SS_Max_Freq_Pct': ss_stats['max_freq'],
                'SS_Total_Affected': ss_stats['total_affected'],
                'SN_N_Mutations': sn_stats['n_mutations'],
                'SN_Mean_Freq_Pct': sn_stats['mean_freq'],
                'SN_Max_Freq_Pct': sn_stats['max_freq'],
                'SN_Total_Affected': sn_stats['total_affected'],
                'NS_N_Mutations': ns_stats['n_mutations'],
                'NS_Mean_Freq_Pct': ns_stats['mean_freq'],
                'NS_Max_Freq_Pct': ns_stats['max_freq'],
                'NS_Total_Affected': ns_stats['total_affected'],
                'NN_N_Mutations': nn_stats['n_mutations'],
                'NN_Mean_Freq_Pct': nn_stats['mean_freq'],
                'NN_Max_Freq_Pct': nn_stats['max_freq'],
                'NN_Total_Affected': nn_stats['total_affected'],
                'STOP_N_Mutations': stop_stats['n_mutations'],
                'STOP_Mean_Freq_Pct': stop_stats['mean_freq'],
                'STOP_Max_Freq_Pct': stop_stats['max_freq'],
                'STOP_Total_Affected': stop_stats['total_affected']
            })

            # Store detailed mutation frequencies in row for later access
            row['mutation_frequency_details'] = mutation_frequencies

            # Collect individual mutations for the detailed DataFrame
            region_name = results['region_name']
            for mutation_type in ['SS', 'SN', 'NS', 'NN', 'STOP_F1']:
                if mutation_type in mutation_frequencies:
                    for mut in mutation_frequencies[mutation_type]:
                        all_mutations.append({
                            'Region': region_name,
                            'Canonical_Frame': results['canonical_gene'],
                            'Alternative_Frame': results['alternative_gene'],
                            'Mutation_Category': mutation_type,
                            'Position': mut['position'],
                            'Mutation': mut['mutation'],
                            'Count': mut['count'],
                            'Frequency_Pct': mut['frequency_pct'],
                            'N_Sequences': n_sequences
                        })

        # Add dual-coding evidence metrics
        if include_dual_coding_test and 'dual_coding_evidence' in results:
            dc = results['dual_coding_evidence']
            if 'error' not in dc:
                row.update({
                    'Beta_STOP': dc['summary']['betaSTOP'] if 'summary' in dc else np.nan,
                    'Beta_STOP_P_value': dc['summary']['p_value'] if 'summary' in dc else np.nan,
                    'Ti_Tv_Ratio': dc['summary']['Ti/Tv'] if 'summary' in dc else np.nan,
                    'Stop_Contexts': dc['summary']['stop_contexts'] if 'summary' in dc else 0,
                    'Total_Contexts': dc['summary']['total_contexts'] if 'summary' in dc else 0,
                    'Dual_Coding_Evidence': dc.get('dual_coding_evidence', 'UNKNOWN'),
                    'Evidence_Score': dc.get('evidence_score', 0)
                })
            else:
                row.update({
                    'Beta_STOP': np.nan,
                    'Beta_STOP_P_value': np.nan,
                    'Ti_Tv_Ratio': np.nan,
                    'Stop_Contexts': 0,
                    'Total_Contexts': 0,
                    'Dual_Coding_Evidence': 'ERROR',
                    'Evidence_Score': 0
                })
        
        all_results.append(row)
    
    results_df = pd.DataFrame(all_results)
    mutations_df = pd.DataFrame(all_mutations) if all_mutations else pd.DataFrame()

    # Sort mutations DataFrame by region and frequency
    if not mutations_df.empty:
        mutations_df = mutations_df.sort_values(
            ['Region', 'Mutation_Category', 'Frequency_Pct'],
            ascending=[True, True, False]
        )

    # Print summary
    print("\n" + "="*80)
    print("SUMMARY OF ALL ANALYSES")
    print("="*80)

    if include_dual_coding_test and 'Beta_STOP' in results_df.columns:
        # Count high-confidence dual-coding regions
        high_conf = (results_df['Dual_Coding_Evidence'] == 'HIGH CONFIDENCE').sum()
        mod_conf = (results_df['Dual_Coding_Evidence'] == 'MODERATE CONFIDENCE').sum()

        print(f"\nTotal regions analyzed: {len(results_df)}")
        print(f"High confidence dual-coding: {high_conf}")
        print(f"Moderate confidence dual-coding: {mod_conf}")
        print(f"Regions with betaSTOP = 0 (p>0.05): {(results_df['Beta_STOP_P_value'] > 0.05).sum()}")
        print(f"Regions with significant asymmetry: {results_df['Significant_Asymmetry'].sum()}")

    # Print mutation frequency summary if available
    if 'SS_N_Mutations' in results_df.columns:
        print("\n" + "="*80)
        print("MUTATION FREQUENCY SUMMARY")
        print("="*80)
        for _, row in results_df.iterrows():
            print(f"\n{row['Region']}:")
            print(f"  Sequences analyzed: {row['N_Sequences']}")
            print(f"  SS: {row['SS_N_Mutations']} mutations (mean freq: {row['SS_Mean_Freq_Pct']:.2f}%, max: {row['SS_Max_Freq_Pct']:.2f}%)")
            print(f"  SN: {row['SN_N_Mutations']} mutations (mean freq: {row['SN_Mean_Freq_Pct']:.2f}%, max: {row['SN_Max_Freq_Pct']:.2f}%)")
            print(f"  NS: {row['NS_N_Mutations']} mutations (mean freq: {row['NS_Mean_Freq_Pct']:.2f}%, max: {row['NS_Max_Freq_Pct']:.2f}%)")
            print(f"  NN: {row['NN_N_Mutations']} mutations (mean freq: {row['NN_Mean_Freq_Pct']:.2f}%, max: {row['NN_Max_Freq_Pct']:.2f}%)")
            print(f"  STOP: {row['STOP_N_Mutations']} mutations (mean freq: {row['STOP_Mean_Freq_Pct']:.2f}%, max: {row['STOP_Max_Freq_Pct']:.2f}%)")

        if not mutations_df.empty:
            print(f"\nDetailed mutations DataFrame contains {len(mutations_df)} individual mutations")
            print(f"Use: summary_df, mutations_df = compare_multiple_overlaps(...)")
            print(f"     mutations_df has columns: {', '.join(mutations_df.columns.tolist())}")

    return results_df, mutations_df


def export_mutation_frequencies(results_df: pd.DataFrame, output_file: Optional[str] = None) -> Dict[str, pd.DataFrame]:
    """
    Export detailed mutation frequencies for each region to separate DataFrames.

    Parameters:
    -----------
    results_df : pd.DataFrame
        Results DataFrame from compare_multiple_overlaps containing mutation_frequency_details
    output_file : str, optional
        If provided, saves each region's mutation frequencies to a CSV file

    Returns:
    --------
    Dict[str, pd.DataFrame]
        Dictionary mapping region names to DataFrames containing detailed mutation frequencies
    """
    mutation_dfs = {}

    for _, row in results_df.iterrows():
        region = row['Region']
        if 'mutation_frequency_details' not in row:
            continue

        mutation_details = row['mutation_frequency_details']
        all_mutations = []

        for mutation_type in ['SS', 'SN', 'NS', 'NN', 'STOP_F1']:
            if mutation_type in mutation_details:
                for mut in mutation_details[mutation_type]:
                    all_mutations.append({
                        'Region': region,
                        'Mutation_Type': mutation_type,
                        'Position': mut['position'],
                        'Mutation': mut['mutation'],
                        'Count': mut['count'],
                        'Frequency_Pct': mut['frequency_pct']
                    })

        if all_mutations:
            df = pd.DataFrame(all_mutations)
            df = df.sort_values(['Mutation_Type', 'Frequency_Pct'], ascending=[True, False])
            mutation_dfs[region] = df

            if output_file:
                # Save to CSV with region name
                safe_region = region.replace('/', '_').replace('\\', '_')
                filename = f"{output_file}_{safe_region}_mutations.csv"
                df.to_csv(filename, index=False)
                print(f"Saved mutation frequencies for {region} to {filename}")

    return mutation_dfs


def calculate_nucleotide_diversity(bases_list: List[str]) -> float:
    """
    Calculate nucleotide diversity (pi) for a position.

    Pi = (n/(n-1)) * sum(p_i * p_j) for all pairs of different nucleotides
    where p_i is the frequency of nucleotide i and n is the sample size.

    The n/(n-1) factor is a bias correction for finite sample sizes.

    Parameters:
    -----------
    bases_list : List[str]
        List of bases at a specific position

    Returns:
    --------
    float : Nucleotide diversity (0 to 1)
    """
    # Filter valid bases
    valid_bases = [b for b in bases_list if b in 'ATGC']

    if len(valid_bases) < 2:
        return 0.0

    # Count frequencies
    n = len(valid_bases)
    base_counts = Counter(valid_bases)

    # Calculate diversity (gene diversity / heterozygosity)
    diversity = 0.0
    for base_i, count_i in base_counts.items():
        for base_j, count_j in base_counts.items():
            if base_i != base_j:
                p_i = count_i / n
                p_j = count_j / n
                diversity += p_i * p_j

    # Apply bias correction for finite sample size
    # This converts gene diversity to nucleotide diversity (pi)
    diversity *= n / (n - 1)

    return diversity

def calculate_nucleotide_diversity_unweighted(bases_list: List[str]) -> float:
    """
    Calculate whether a position is segregating (polymorphic) - unweighted.

    Returns 1.0 if the position has more than one allele, 0.0 otherwise.
    When averaged across positions, this gives the proportion of segregating sites.

    Parameters:
    -----------
    bases_list : List[str]
        List of bases at a specific position

    Returns:
    --------
    float : 1.0 if segregating, 0.0 if monomorphic
    """
    # Filter valid bases
    valid_bases = [b for b in bases_list if b in 'ATGC']

    if len(valid_bases) < 2:
        return 0.0

    # Check if position is segregating (has more than one unique allele)
    unique_bases = set(valid_bases)

    return 1.0 if len(unique_bases) > 1 else 0.0

def calculate_codon_position_titv(sequences: List[str], position: int = 0) -> Dict:
    """
    Calculate Ti/Tv ratio and nucleotide diversity for a specific codon position.

    Parameters:
    -----------
    sequences : List[str]
        List of aligned DNA sequences (should be in-frame, length divisible by 3)
    position : int
        Codon position (0=first, 1=second, 2=third)

    Returns:
    --------
    Dict containing:
        - ti_count: Number of transitions observed
        - tv_count: Number of transversions observed
        - ti_tv_ratio: Transition/transversion ratio (with +0.5 pseudocount) - known as Haldane-Anscombe correction
        - ti_tv_ratio_raw: Raw ratio without pseudocount
        - total_substitutions: Total number of substitutions analyzed
        - nucleotide_diversity: Average nucleotide diversity across positions
    """
    if not sequences or len(sequences) < 2:
        return {
            'ti_count': 0,
            'tv_count': 0,
            'ti_tv_ratio': 1.0,  # (0+1)/(0+1)
            'ti_tv_ratio_raw': np.nan,
            'total_substitutions': 0,
            'nucleotide_diversity': 0.0
        }

    # Clean sequences and ensure they're in frame
    clean_seqs = []
    for seq in sequences:
        seq = seq.upper().strip()
        # Trim to multiple of 3
        seq = seq[:len(seq) - (len(seq) % 3)] # ACGAC - > ACG
        if len(seq) >= 3:
            clean_seqs.append(seq)

    if len(clean_seqs) < 2:
        return {
            'ti_count': 0,
            'tv_count': 0,
            'ti_tv_ratio': 1.0,
            'ti_tv_ratio_raw': np.nan,
            'total_substitutions': 0,
            'nucleotide_diversity': 0.0
        }

    # Extract bases at the specified codon position
    position_bases = []
    for seq in clean_seqs:
        bases = [seq[i] for i in range(position, len(seq), 3) if i < len(seq)]
        position_bases.append(bases)

    # Calculate nucleotide diversity for each position
    n_positions = len(position_bases[0])
    diversities = []

    for pos_idx in range(n_positions):
        bases_at_pos = [seq_bases[pos_idx] for seq_bases in position_bases
                        if pos_idx < len(seq_bases)]
        diversity = calculate_nucleotide_diversity(bases_at_pos)
        diversities.append(diversity)

    avg_diversity = np.mean(diversities) if diversities else 0.0

    # Count transitions and transversions
    ti_count = 0
    tv_count = 0

    # Use the first sequence as reference
    ref_bases = position_bases[0]

    # Collect unique substitutions: (position_index, alt_base)

    for alt_bases in position_bases[1:]:
        min_len = min(len(ref_bases), len(alt_bases))
        for i in range(min_len):
            ref_base = ref_bases[i]
            alt_base = alt_bases[i]

            # Skip if same or if either is ambiguous
            if ref_base == alt_base or ref_base not in 'ATGC' or alt_base not in 'ATGC':
                continue

            sub_type = classify_substitution(ref_base, alt_base)
            if sub_type == 'transition':
                ti_count += 1
            elif sub_type == 'transversion':
                tv_count += 1

    total_subs = ti_count + tv_count

    # Raw ratio (can be inf or nan)
    ti_tv_ratio_raw = ti_count / tv_count if tv_count > 0 else (np.inf if ti_count > 0 else np.nan)

    # Ratio with pseudocount (add 0.5 to both to avoid division by zero)
    ti_tv_ratio = (ti_count + 0.5) / (tv_count + 0.5)

    return {
        'ti_count': ti_count,
        'tv_count': tv_count,
        'ti_tv_ratio': ti_tv_ratio,
        'ti_tv_ratio_raw': ti_tv_ratio_raw,
        'total_substitutions': total_subs,
        'nucleotide_diversity': avg_diversity
    }


def _process_random_region(region_seqs: List[str]) -> Dict:
    """
    Helper function to process a single random region (for parallel processing).

    Parameters:
    -----------
    region_seqs : List[str]
        List of sequences for a random region

    Returns:
    --------
    Dict with Ti/Tv results for all three positions, or None if invalid
    """
    # Ensure in-frame
    region_seqs = [seq[:len(seq) - (len(seq) % 3)] for seq in region_seqs if len(seq) >= 3]

    if len(region_seqs) < 2:
        return None

    # Calculate Ti/Tv for each position
    first_result = calculate_codon_position_titv(region_seqs, position=0)
    second_result = calculate_codon_position_titv(region_seqs, position=1)
    third_result = calculate_codon_position_titv(region_seqs, position=2)

    # Combined first + second
    ti_12 = first_result['ti_count'] + second_result['ti_count']
    tv_12 = first_result['tv_count'] + second_result['tv_count']
    titv_12 = (ti_12 + 0.5) / (tv_12 + 0.5)  # With pseudocount
    titv_12_raw = ti_12 / tv_12 if tv_12 > 0 else (np.inf if ti_12 > 0 else np.nan)

    # Average diversity across first and second positions
    diversity_12 = (first_result['nucleotide_diversity'] + second_result['nucleotide_diversity']) / 2

    return {
        'first': first_result,
        'second': second_result,
        'third': third_result,
        'titv_12': titv_12,
        'titv_12_raw': titv_12_raw,
        'diversity_12': diversity_12
    }


def extract_random_regions(sequences: List[str],
                           region_length: int,
                           excluded_ranges: List[Tuple[int, int]],
                           n_iterations: int = 1000,
                           seed: Optional[int] = None) -> List[List[str]]:
    """
    Extract random regions from sequences, avoiding specified excluded ranges.

    Parameters:
    -----------
    sequences : List[str]
        List of full-length sequences (e.g., rRNA genes)
    region_length : int
        Length of random regions to extract (should be divisible by 3)
    excluded_ranges : List[Tuple[int, int]]
        List of (start, end) tuples for regions to exclude (0-based, end-exclusive)
    n_iterations : int
        Number of random regions to extract
    seed : Optional[int]
        Random seed for reproducibility

    Returns:
    --------
    List[List[str]]
        List of random region sequence sets, each containing sequences from all samples
    """
    if seed is not None:
        np.random.seed(seed)

    # Ensure region_length is divisible by 3
    region_length = region_length - (region_length % 3)

    # Find valid positions (positions where we can start a region)
    # Get minimum sequence length
    min_len = min(len(seq) for seq in sequences)

    # Create array of valid start positions
    valid_positions = []
    for start_pos in range(0, min_len - region_length + 1):
        end_pos = start_pos + region_length

        # Check if this range overlaps with any excluded range
        is_valid = True
        for excl_start, excl_end in excluded_ranges:
            # Check for overlap
            if not (end_pos <= excl_start or start_pos >= excl_end):
                is_valid = False
                break

        if is_valid:
            valid_positions.append(start_pos)

    if len(valid_positions) == 0:
        warnings.warn("No valid positions found for random region extraction")
        return []

    # Sample random positions
    n_samples = min(n_iterations, len(valid_positions))
    sampled_positions = np.random.choice(valid_positions, size=n_samples, replace=True)

    # Extract regions
    random_regions = []
    for start_pos in sampled_positions:
        region_seqs = [seq[start_pos:start_pos + region_length] for seq in sequences]
        random_regions.append(region_seqs)

    return random_regions


def analyze_rrna_mdp_asymmetry(df: pd.DataFrame,
                               rrna_seq_col: str,
                               mdp_configs: List[Dict[str, any]],
                               n_iterations: int = 1000,
                               seed: Optional[int] = 42,
                               threshold: float = 0.05,
                               verbose: bool = True) -> Dict:
    """
    Analyze Ti/Tv asymmetry between MDPs and random regions within rRNA genes.

    This function tests whether MDPs embedded in rRNA genes show protein-like
    substitution patterns (first + second codon position Ti/Tv < third position Ti/Tv)
    compared to randomly selected rRNA regions that don't overlap with MDPs.

    IMPORTANT:
    - rRNA sequences are NOT cleaned (they are non-coding)
    - MDP sequences ARE cleaned using clean_sequences_adaptive (they are protein-coding)
    - Random regions from rRNA are NOT cleaned (for fair comparison with non-coding baseline)

    Parameters:
    -----------
    df : pd.DataFrame
        DataFrame with sequence data
    rrna_seq_col : str
        Column name containing rRNA gene sequences (will NOT be cleaned)
    mdp_configs : List[Dict]
        List of MDP configuration dictionaries, each containing:
        - 'name': Name of the MDP (e.g., 'SHLP3')
        - 'seq_col': Column name with MDP sequences (optional)
        - 'start': Start position in rRNA (0-based)
        - 'end': End position in rRNA (0-based, exclusive)
    n_iterations : int
        Number of random regions to sample for comparison (default=1000)
    seed : Optional[int]
        Random seed for reproducibility (default=42)
    threshold : float
        Threshold for clean_sequences_adaptive for MDP sequences only (default=0.05)
    verbose : bool
        Print progress messages (default=True)

    Returns:
    --------
    Dict containing:
        - mdp_results: Dict with Ti/Tv ratios for each MDP at each codon position
        - random_results: Dict with Ti/Tv distributions for random regions
        - comparison: Statistical comparison of patterns
        - parameters: Analysis parameters used

    Example:
    --------
    >>> mdp_configs = [
    ...     {'name': 'SHLP3', 'seq_col': 'SHLP3_seq', 'start': 1701, 'end': 1819},
    ...     {'name': 'SHLP2', 'seq_col': 'SHLP2_seq', 'start': 2086, 'end': 2168}
    ... ]
    >>> results = analyze_rrna_mdp_asymmetry(df, 'RNR2', mdp_configs)
    """
    if verbose:
        print("=" * 80)
        print("rRNA-MDP Ti/Tv ASYMMETRY ANALYSIS")
        print("=" * 80)
        print(f"\nAnalyzing {len(mdp_configs)} MDPs within {rrna_seq_col}")

    # Extract rRNA sequences (no cleaning - rRNA is not protein-coding)
    rrna_seqs = df[rrna_seq_col].dropna().tolist()

    if len(rrna_seqs) < 2:
        return {'error': 'Insufficient rRNA sequences for analysis'}

    if verbose:
        print(f"\nUsing {len(rrna_seqs)} rRNA sequences")
        print("Note: rRNA sequences are NOT cleaned (not protein-coding)")
        print("      Only MDP sequences (protein-coding) will be cleaned")

    # Analyze each MDP
    mdp_results = {}

    for config in mdp_configs:
        mdp_name = config['name']
        mdp_seq_col = config.get('seq_col')
        start = config.get('start')
        end = config.get('end')

        if verbose:
            print(f"\n{'='*60}")
            print(f"Analyzing MDP: {mdp_name}")
            print(f"Seq column: {mdp_seq_col if mdp_seq_col else 'N/A'}")
            print(f"{'='*60}")

        # Extract MDP sequences from rRNA sequences
        if mdp_seq_col in df.columns:
            # Use pre-extracted MDP sequences from dataframe
            mdp_seqs_raw = df[mdp_seq_col].dropna().tolist()
            if verbose:
                print(f"Using pre-extracted sequences from column '{mdp_seq_col}'")
        elif start is not None and end is not None:
            mdp_seqs_raw = [seq[start:end] for seq in rrna_seqs if len(seq) >= end]
            if verbose:
                print(f"Extracted from rRNA positions {start}-{end}")
        else:
            warnings.warn(f"Cannot extract sequences for {mdp_name}: missing coordinates or sequence column")
            continue

        if len(mdp_seqs_raw) < 2:
            if verbose:
                print(f"  Skipping {mdp_name}: insufficient sequences (n={len(mdp_seqs_raw)})")
            continue

        # Clean MDP sequences (they are protein-coding)
        if verbose:
            print(f"  Cleaning {len(mdp_seqs_raw)} MDP sequences (threshold={threshold})...")

        cleaned_mdp = cag.clean_sequences_adaptive(mdp_seqs_raw, threshold=threshold)
        mdp_seqs = cleaned_mdp.cleaned_sequences

        if len(mdp_seqs) < 2:
            if verbose:
                print(f"  Skipping {mdp_name}: insufficient sequences after cleaning (n={len(mdp_seqs)})")
            continue

        # Ensure sequences are in frame (divisible by 3)
        mdp_seqs = [seq[:len(seq) - (len(seq) % 3)] for seq in mdp_seqs if len(seq) >= 3]

        if len(mdp_seqs) < 2:
            if verbose:
                print(f"  Skipping {mdp_name}: insufficient sequences after frame adjustment")
            continue

        if verbose:
            from Bio.Seq import Seq

            print(f"  Analyzing {len(mdp_seqs)} sequences, length={len(mdp_seqs[0])} bp")

            # Report first and last 3 codons for verification
            seq_length = len(mdp_seqs[0])
            n_codons = seq_length // 3

            if n_codons >= 6:
                # Get first 3 codons from first sequence
                first_3_codons = [mdp_seqs[0][i*3:(i+1)*3] for i in range(3)]
                # Get last 3 codons from first sequence
                last_3_codons = [mdp_seqs[0][i*3:(i+1)*3] for i in range(n_codons-3, n_codons)]

                # Translate codons to amino acids using BioPython
                first_3_aa = [str(Seq(codon).translate()) for codon in first_3_codons]
                last_3_aa = [str(Seq(codon).translate()) for codon in last_3_codons]

                print(f"  First 3 codons: {' '.join(first_3_codons)} -> {' '.join(first_3_aa)}")
                print(f"  Last 3 codons:  {' '.join(last_3_codons)} -> {' '.join(last_3_aa)}")
            else:
                # If fewer than 6 codons, just show all
                all_codons = [mdp_seqs[0][i*3:(i+1)*3] for i in range(n_codons)]
                all_aa = [str(Seq(codon).translate()) for codon in all_codons]
                print(f"  All {n_codons} codons: {' '.join(all_codons)} -> {' '.join(all_aa)}")

        # Calculate Ti/Tv for each codon position
        mdp_titv = {}
        for pos in [0, 1, 2]:
            pos_name = ['first', 'second', 'third'][pos]
            result = calculate_codon_position_titv(mdp_seqs, position=pos)
            mdp_titv[pos_name] = result

            if verbose:
                if result['total_substitutions'] > 0:
                    print(f"  {pos_name:6s} codon position: Ti/Tv = {result['ti_tv_ratio']:6.3f} "
                          f"(Ti={result['ti_count']:3d}, Tv={result['tv_count']:3d}, "
                          f"Total={result['total_substitutions']:4d})")
                else:
                    print(f"  {pos_name:6s} codon position: No substitutions found")

        # Calculate combined first + second position Ti/Tv
        ti_12 = mdp_titv['first']['ti_count'] + mdp_titv['second']['ti_count']
        tv_12 = mdp_titv['first']['tv_count'] + mdp_titv['second']['tv_count']
        titv_12 = (ti_12 + 0.5) / (tv_12 + 0.5)  # With pseudocount
        titv_12_raw = ti_12 / tv_12 if tv_12 > 0 else (np.inf if ti_12 > 0 else np.nan)

        # Average diversity across first and second positions
        diversity_12 = (mdp_titv['first']['nucleotide_diversity'] +
                       mdp_titv['second']['nucleotide_diversity']) / 2

        if verbose:
            print(f"  Combined (1st+2nd): Ti/Tv = {titv_12:6.3f} (Ti={ti_12:3d}, Tv={tv_12:3d})")
            print(f"  Nucleotide diversity (1st+2nd): {diversity_12:6.4f}")
            print(f"  Nucleotide diversity (3rd): {mdp_titv['third']['nucleotide_diversity']:6.4f}")

        mdp_results[mdp_name] = {
            'first_pos': mdp_titv['first'],
            'second_pos': mdp_titv['second'],
            'third_pos': mdp_titv['third'],
            'combined_12_titv': titv_12,
            'combined_12_titv_raw': titv_12_raw,
            'combined_12_diversity': diversity_12,
            'length': len(mdp_seqs[0]),
            'n_sequences': len(mdp_seqs),
            'start': start,
            'end': end
        }

    # Calculate global Ti/Tv ratios across all MDPs
    if len(mdp_results) > 0:
        global_mdp_ti_first = sum(res['first_pos']['ti_count'] for res in mdp_results.values())
        global_mdp_tv_first = sum(res['first_pos']['tv_count'] for res in mdp_results.values())
        global_mdp_ti_second = sum(res['second_pos']['ti_count'] for res in mdp_results.values())
        global_mdp_tv_second = sum(res['second_pos']['tv_count'] for res in mdp_results.values())
        global_mdp_ti_third = sum(res['third_pos']['ti_count'] for res in mdp_results.values())
        global_mdp_tv_third = sum(res['third_pos']['tv_count'] for res in mdp_results.values())

        # Combined first + second
        global_mdp_ti_12 = global_mdp_ti_first + global_mdp_ti_second
        global_mdp_tv_12 = global_mdp_tv_first + global_mdp_tv_second

        # Calculate global Ti/Tv ratios with pseudocount
        global_mdp_titv_first = (global_mdp_ti_first + 0.5) / (global_mdp_tv_first + 0.5)
        global_mdp_titv_second = (global_mdp_ti_second + 0.5) / (global_mdp_tv_second + 0.5)
        global_mdp_titv_third = (global_mdp_ti_third + 0.5) / (global_mdp_tv_third + 0.5)
        global_mdp_titv_12 = (global_mdp_ti_12 + 0.5) / (global_mdp_tv_12 + 0.5)

        # Calculate RoA for diversities (average)
        global_mdp_diversity_first = np.mean([res['first_pos']['nucleotide_diversity'] for res in mdp_results.values()])
        global_mdp_diversity_second = np.mean([res['second_pos']['nucleotide_diversity'] for res in mdp_results.values()])
        global_mdp_diversity_third = np.mean([res['third_pos']['nucleotide_diversity'] for res in mdp_results.values()])
        global_mdp_diversity_12 = np.mean([res['combined_12_diversity'] for res in mdp_results.values()])

        # Check if global pattern is protein-like
        global_mdp_protein_pattern = global_mdp_titv_12 < global_mdp_titv_third

        if verbose:
            print(f"\n{'='*80}")
            print("GLOBAL MDP STATISTICS (ALL MDPs COMBINED)")
            print(f"{'='*80}")
            print(f"  First position:  Ti/Tv = {global_mdp_titv_first:6.3f} (Ti={global_mdp_ti_first:4d}, Tv={global_mdp_tv_first:4d})")
            print(f"  Second position: Ti/Tv = {global_mdp_titv_second:6.3f} (Ti={global_mdp_ti_second:4d}, Tv={global_mdp_tv_second:4d})")
            print(f"  Third position:  Ti/Tv = {global_mdp_titv_third:6.3f} (Ti={global_mdp_ti_third:4d}, Tv={global_mdp_tv_third:4d})")
            print(f"  Combined (1st+2nd): Ti/Tv = {global_mdp_titv_12:6.3f} (Ti={global_mdp_ti_12:4d}, Tv={global_mdp_tv_12:4d})")
            print(f"  RoA Diversity (1st+2nd): {global_mdp_diversity_12:6.4f}")
            print(f"  RoA Diversity (3rd): {global_mdp_diversity_third:6.4f}")
            if global_mdp_protein_pattern:
                print(f"  ✓ Global pattern shows protein-coding signature (Ti/Tv 1st+2nd < 3rd)")
            else:
                print(f"  ✗ Global pattern does NOT show protein-coding signature")

        global_mdp_stats = {
            'first_pos_titv': global_mdp_titv_first,
            'second_pos_titv': global_mdp_titv_second,
            'third_pos_titv': global_mdp_titv_third,
            'combined_12_titv': global_mdp_titv_12,
            'first_pos_ti': global_mdp_ti_first,
            'first_pos_tv': global_mdp_tv_first,
            'second_pos_ti': global_mdp_ti_second,
            'second_pos_tv': global_mdp_tv_second,
            'third_pos_ti': global_mdp_ti_third,
            'third_pos_tv': global_mdp_tv_third,
            'combined_12_ti': global_mdp_ti_12,
            'combined_12_tv': global_mdp_tv_12,
            'roa_diversity_first': global_mdp_diversity_first,
            'roa_diversity_second': global_mdp_diversity_second,
            'roa_diversity_third': global_mdp_diversity_third,
            'roa_diversity_12': global_mdp_diversity_12,
            'shows_protein_pattern': global_mdp_protein_pattern
        }
    else:
        global_mdp_stats = None

    # Prepare excluded ranges (all MDP positions)
    excluded_ranges = [(config['start'], config['end'])
                       for config in mdp_configs
                       if 'start' in config and 'end' in config and
                       config['start'] is not None and config['end'] is not None]

    if verbose:
        print(f"\n{'='*80}")
        print("RANDOM REGION ANALYSIS")
        print(f"{'='*80}")
        print(f"Extracting {n_iterations} random regions for comparison...")
        print(f"Excluded ranges (MDP positions): {excluded_ranges}")

    # Analyze random regions
    random_results = {}

    # For each MDP, sample random regions of the same length
    for mdp_name, mdp_data in mdp_results.items():
        region_length = mdp_data['length']

        if verbose:
            print(f"\n{'='*60}")
            print(f"Random regions for {mdp_name} (length={region_length} bp)")
            print(f"{'='*60}")

        # Extract random regions from cleaned rRNA sequences
        random_regions = extract_random_regions(
            rrna_seqs,
            region_length,
            excluded_ranges,
            n_iterations,
            seed
        )

        if len(random_regions) == 0:
            if verbose:
                print(f"  No valid random regions found for {mdp_name}")
            continue

        # Process random regions in parallel using ThreadPoolExecutor
        if verbose:
            print(f"  Processing {len(random_regions)} random regions in parallel...")

        random_titv_first = []
        random_titv_second = []
        random_titv_third = []
        random_titv_12 = []
        random_diversity_12 = []
        random_diversity_3 = []

        # Use ThreadPoolExecutor for I/O-bound operations (faster for this use case)
        with ThreadPoolExecutor(max_workers=8) as executor:
            # Submit all tasks
            future_to_idx = {executor.submit(_process_random_region, region_seqs): idx
                            for idx, region_seqs in enumerate(random_regions)}

            # Collect results as they complete
            for future in as_completed(future_to_idx):
                result = future.result()

                if result is None:
                    continue

                # Store results with pseudocount ratios
                if result['first']['total_substitutions'] > 0:
                    random_titv_first.append(result['first']['ti_tv_ratio'])
                if result['second']['total_substitutions'] > 0:
                    random_titv_second.append(result['second']['ti_tv_ratio'])
                if result['third']['total_substitutions'] > 0:
                    random_titv_third.append(result['third']['ti_tv_ratio'])

                random_titv_12.append(result['titv_12'])
                random_diversity_12.append(result['diversity_12'])
                random_diversity_3.append(result['third']['nucleotide_diversity'])

        # Filter out infinite values (pseudocount should prevent this, but just in case)
        random_titv_first = [x for x in random_titv_first if np.isfinite(x)]
        random_titv_second = [x for x in random_titv_second if np.isfinite(x)]
        random_titv_third = [x for x in random_titv_third if np.isfinite(x)]
        random_titv_12 = [x for x in random_titv_12 if np.isfinite(x)]
        random_diversity_12 = [x for x in random_diversity_12 if np.isfinite(x)]
        random_diversity_3 = [x for x in random_diversity_3 if np.isfinite(x)]

        if verbose:
            print(f"  Valid random regions: {len(random_titv_12)}/{len(random_regions)}")
            if len(random_titv_first) > 0:
                print(f"  First position Ti/Tv:  {np.mean(random_titv_first):6.3f} ± {np.std(random_titv_first):5.3f}")
            if len(random_titv_second) > 0:
                print(f"  Second position Ti/Tv: {np.mean(random_titv_second):6.3f} ± {np.std(random_titv_second):5.3f}")
            if len(random_titv_third) > 0:
                print(f"  Third position Ti/Tv:  {np.mean(random_titv_third):6.3f} ± {np.std(random_titv_third):5.3f}")
            if len(random_titv_12) > 0:
                print(f"  Combined (1st+2nd):    {np.mean(random_titv_12):6.3f} ± {np.std(random_titv_12):5.3f}")
            if len(random_diversity_12) > 0:
                print(f"  Diversity (1st+2nd):   {np.mean(random_diversity_12):6.4f} ± {np.std(random_diversity_12):5.4f}")
            if len(random_diversity_3) > 0:
                print(f"  Diversity (3rd):       {np.mean(random_diversity_3):6.4f} ± {np.std(random_diversity_3):5.4f}")

        random_results[mdp_name] = {
            'first_pos_titv': random_titv_first,
            'second_pos_titv': random_titv_second,
            'third_pos_titv': random_titv_third,
            'combined_12_titv': random_titv_12,
            'combined_12_diversity': random_diversity_12,
            'third_diversity': random_diversity_3,
            'n_valid_regions': len(random_titv_12)
        }

    # Calculate global Ti/Tv ratios across all random regions
    # Need to re-process to get Ti/Tv counts (not just ratios)
    if len(random_results) > 0:
        if verbose:
            print(f"\n{'='*80}")
            print("GLOBAL RANDOM REGION STATISTICS (ALL RANDOM REGIONS COMBINED)")
            print(f"{'='*80}")
            print("  Aggregating Ti/Tv counts from all random regions...")

        # Collect all Ti/Tv counts from all random regions across all MDPs
        all_random_ti_first = []
        all_random_tv_first = []
        all_random_ti_second = []
        all_random_tv_second = []
        all_random_ti_third = []
        all_random_tv_third = []
        all_random_diversity_first = []
        all_random_diversity_second = []
        all_random_diversity_third = []
        all_random_diversity_12 = []

        # Re-process random regions to get Ti/Tv counts
        for mdp_name, mdp_data in mdp_results.items():
            if mdp_name not in random_results:
                continue

            region_length = mdp_data['length']

            # Re-extract random regions (using same parameters)
            random_regions = extract_random_regions(
                rrna_seqs,
                region_length,
                excluded_ranges,
                n_iterations,
                seed
            )

            # Process each random region to get Ti/Tv counts
            with ThreadPoolExecutor(max_workers=8) as executor:
                future_to_idx = {executor.submit(_process_random_region, region_seqs): idx
                                for idx, region_seqs in enumerate(random_regions)}

                for future in as_completed(future_to_idx):
                    result = future.result()
                    if result is None:
                        continue

                    all_random_ti_first.append(result['first']['ti_count'])
                    all_random_tv_first.append(result['first']['tv_count'])
                    all_random_ti_second.append(result['second']['ti_count'])
                    all_random_tv_second.append(result['second']['tv_count'])
                    all_random_ti_third.append(result['third']['ti_count'])
                    all_random_tv_third.append(result['third']['tv_count'])
                    all_random_diversity_first.append(result['first']['nucleotide_diversity'])
                    all_random_diversity_second.append(result['second']['nucleotide_diversity'])
                    all_random_diversity_third.append(result['third']['nucleotide_diversity'])
                    all_random_diversity_12.append(result['diversity_12'])

        # Sum all Ti and Tv counts
        global_random_ti_first = sum(all_random_ti_first)
        global_random_tv_first = sum(all_random_tv_first)
        global_random_ti_second = sum(all_random_ti_second)
        global_random_tv_second = sum(all_random_tv_second)
        global_random_ti_third = sum(all_random_ti_third)
        global_random_tv_third = sum(all_random_tv_third)

        # Combined first + second
        global_random_ti_12 = global_random_ti_first + global_random_ti_second
        global_random_tv_12 = global_random_tv_first + global_random_tv_second

        # Calculate global Ti/Tv ratios with pseudocount
        global_random_titv_first = (global_random_ti_first + 0.5) / (global_random_tv_first + 0.5)
        global_random_titv_second = (global_random_ti_second + 0.5) / (global_random_tv_second + 0.5)
        global_random_titv_third = (global_random_ti_third + 0.5) / (global_random_tv_third + 0.5)
        global_random_titv_12 = (global_random_ti_12 + 0.5) / (global_random_tv_12 + 0.5)

        # Calculate RoA for diversities (average)
        global_random_diversity_first = np.mean(all_random_diversity_first)
        global_random_diversity_second = np.mean(all_random_diversity_second)
        global_random_diversity_third = np.mean(all_random_diversity_third)
        global_random_diversity_12 = np.mean(all_random_diversity_12)

        # Check if global pattern is protein-like
        global_random_protein_pattern = global_random_titv_12 < global_random_titv_third

        if verbose:
            print(f"  Processed {len(all_random_ti_first)} total random regions")
            print(f"  First position:  Ti/Tv = {global_random_titv_first:6.3f} (Ti={global_random_ti_first:4d}, Tv={global_random_tv_first:4d})")
            print(f"  Second position: Ti/Tv = {global_random_titv_second:6.3f} (Ti={global_random_ti_second:4d}, Tv={global_random_tv_second:4d})")
            print(f"  Third position:  Ti/Tv = {global_random_titv_third:6.3f} (Ti={global_random_ti_third:4d}, Tv={global_random_tv_third:4d})")
            print(f"  Combined (1st+2nd): Ti/Tv = {global_random_titv_12:6.3f} (Ti={global_random_ti_12:4d}, Tv={global_random_tv_12:4d})")
            print(f"  RoA Diversity (1st+2nd): {global_random_diversity_12:6.4f}")
            print(f"  RoA Diversity (3rd): {global_random_diversity_third:6.4f}")
            if global_random_protein_pattern:
                print(f"  ✓ Global pattern shows protein-coding signature (Ti/Tv 1st+2nd < 3rd)")
            else:
                print(f"  ✗ Global pattern does NOT show protein-coding signature")

        global_random_stats = {
            'first_pos_titv': global_random_titv_first,
            'second_pos_titv': global_random_titv_second,
            'third_pos_titv': global_random_titv_third,
            'combined_12_titv': global_random_titv_12,
            'first_pos_ti': global_random_ti_first,
            'first_pos_tv': global_random_tv_first,
            'second_pos_ti': global_random_ti_second,
            'second_pos_tv': global_random_tv_second,
            'third_pos_ti': global_random_ti_third,
            'third_pos_tv': global_random_tv_third,
            'combined_12_ti': global_random_ti_12,
            'combined_12_tv': global_random_tv_12,
            'roa_diversity_first': global_random_diversity_first,
            'roa_diversity_second': global_random_diversity_second,
            'roa_diversity_third': global_random_diversity_third,
            'roa_diversity_12': global_random_diversity_12,
            'shows_protein_pattern': global_random_protein_pattern,
            'n_regions': len(all_random_ti_first)
        }
    else:
        global_random_stats = None

    # Statistical comparison
    if verbose:
        print(f"\n{'='*80}")
        print("STATISTICAL COMPARISON")
        print(f"{'='*80}")

    comparison_results = {}

    for mdp_name in mdp_results.keys():
        if mdp_name not in random_results:
            continue

        mdp = mdp_results[mdp_name]
        random = random_results[mdp_name]

        if verbose:
            print(f"\n{mdp_name}:")

        # Test if MDP shows protein-like pattern: Ti/Tv(1st+2nd) < Ti/Tv(3rd)
        mdp_titv_12 = mdp['combined_12_titv']
        mdp_titv_3 = mdp['third_pos']['ti_tv_ratio']
        mdp_protein_pattern = mdp_titv_12 < mdp_titv_3 if np.isfinite(mdp_titv_12) and np.isfinite(mdp_titv_3) else False

        if verbose:
            print(f"  MDP Ti/Tv(1st+2nd) = {mdp_titv_12:.3f}, Ti/Tv(3rd) = {mdp_titv_3:.3f}")
            print(f"  MDP shows protein pattern: {mdp_protein_pattern}")

        # For random regions, calculate proportion showing protein pattern
        random_protein_pattern_count = 0
        for i in range(len(random['combined_12_titv'])):
            if (i < len(random['third_pos_titv']) and
                random['combined_12_titv'][i] < random['third_pos_titv'][i]):
                random_protein_pattern_count += 1

        prop_random_protein = (random_protein_pattern_count / len(random['combined_12_titv'])
                               if len(random['combined_12_titv']) > 0 else 0)

        # Statistical test: is the Ti/Tv difference (3rd - 1st+2nd) in MDP significantly
        # larger than in random regions?
        # H0: MDP difference = random difference
        # H1: MDP difference > random difference (protein-like pattern is stronger in MDP)
        if (len(random['combined_12_titv']) > 0 and
            len(random['third_pos_titv']) > 0 and
            np.isfinite(mdp_titv_12) and
            np.isfinite(mdp_titv_3)):

            # Calculate Ti/Tv difference for MDP
            mdp_titv_diff = mdp_titv_3 - mdp_titv_12

            # Calculate Ti/Tv differences for each random region
            random_titv_diffs = []
            for i in range(min(len(random['combined_12_titv']), len(random['third_pos_titv']))):
                if np.isfinite(random['combined_12_titv'][i]) and np.isfinite(random['third_pos_titv'][i]):
                    diff = random['third_pos_titv'][i] - random['combined_12_titv'][i]
                    random_titv_diffs.append(diff)

            if len(random_titv_diffs) > 0:
                random_mean_diff = np.mean(random_titv_diffs)
                random_median_diff = np.median(random_titv_diffs)

                # Empirical p-value for one-sided test (H1: MDP diff > random diff)
                # p-value = proportion of random differences >= MDP difference
                n_equal_or_greater = sum(1 for rd in random_titv_diffs if rd >= mdp_titv_diff)
                p_value_diff = n_equal_or_greater / len(random_titv_diffs)

                # Also calculate summary stats for the individual Ti/Tv values (for reporting)
                random_mean_12 = np.mean(random['combined_12_titv'])
                random_median_12 = np.median(random['combined_12_titv'])
                random_mean_3 = np.mean(random['third_pos_titv'])
            else:
                random_mean_diff = np.nan
                random_median_diff = np.nan
                p_value_diff = np.nan
                random_mean_12 = np.nan
                random_median_12 = np.nan
                random_mean_3 = np.nan
        else:
            mdp_titv_diff = np.nan
            random_mean_diff = np.nan
            random_median_diff = np.nan
            p_value_diff = np.nan
            random_mean_12 = np.nan
            random_median_12 = np.nan
            random_mean_3 = np.nan

        # Calculate diversity metrics
        mdp_diversity_12 = mdp['combined_12_diversity']
        mdp_diversity_3 = mdp['third_pos']['nucleotide_diversity']
        random_mean_diversity_12 = np.mean(random['combined_12_diversity']) if random['combined_12_diversity'] else np.nan
        random_mean_diversity_3 = np.mean(random['third_diversity']) if random['third_diversity'] else np.nan

        if verbose:
            print(f"  Random regions showing protein pattern: {random_protein_pattern_count}/{len(random['combined_12_titv'])} ({prop_random_protein:.1%})")
            print(f"  Random Ti/Tv(1st+2nd): mean={random_mean_12:.3f}, median={random_median_12:.3f}")
            if not np.isnan(random_mean_3):
                print(f"  Random Ti/Tv(3rd): mean={random_mean_3:.3f}")
            print(f"  MDP Ti/Tv difference (3rd - 1st+2nd): {mdp_titv_diff:.3f}")
            if not np.isnan(random_mean_diff):
                print(f"  Random Ti/Tv difference: mean={random_mean_diff:.3f}, median={random_median_diff:.3f}")
            print(f"  MDP diversity (1st+2nd): {mdp_diversity_12:.4f}, Random mean: {random_mean_diversity_12:.4f}")
            print(f"  MDP diversity (3rd): {mdp_diversity_3:.4f}, Random mean: {random_mean_diversity_3:.4f}")
            if not np.isnan(p_value_diff):
                print(f"  p-value (one-sided): {p_value_diff:.4f}")
                if p_value_diff < 0.05:
                    print(f"  ✓ MDP Ti/Tv difference significantly greater than random (p < 0.05)")

        comparison_results[mdp_name] = {
            'mdp_shows_protein_pattern': mdp_protein_pattern,
            'mdp_titv_12': mdp_titv_12,
            'mdp_titv_3': mdp_titv_3,
            'mdp_titv_diff': mdp_titv_diff,
            'mdp_diversity_12': mdp_diversity_12,
            'mdp_diversity_3': mdp_diversity_3,
            'random_protein_pattern_proportion': prop_random_protein,
            'random_mean_titv_12': random_mean_12,
            'random_median_titv_12': random_median_12,
            'random_mean_titv_3': random_mean_3,
            'random_mean_titv_diff': random_mean_diff,
            'random_median_titv_diff': random_median_diff,
            'random_mean_diversity_12': random_mean_diversity_12,
            'random_mean_diversity_3': random_mean_diversity_3,
            'p_value': p_value_diff,
            'significant': p_value_diff < 0.05 if not np.isnan(p_value_diff) else False
        }

    if verbose:
        print(f"\n{'='*80}")
        print("SUMMARY")
        print(f"{'='*80}")
        n_protein_pattern = sum(1 for comp in comparison_results.values() if comp['mdp_shows_protein_pattern'])
        n_significant = sum(1 for comp in comparison_results.values() if comp['significant'])

        print(f"\nMDPs showing protein-like pattern: {n_protein_pattern}/{len(comparison_results)}")
        print(f"MDPs significantly different from random: {n_significant}/{len(comparison_results)}")

        for mdp_name, comp in comparison_results.items():
            status = ""
            if comp['mdp_shows_protein_pattern']:
                status += "✓ Protein pattern"
            if comp['significant']:
                status += " | ✓ Significant (p={:.4f})".format(comp['p_value'])

            print(f"\n{mdp_name}: {status if status else 'No protein pattern'}")
            print(f"  Ti/Tv ratio: {comp['mdp_titv_12']:.3f} (1st+2nd) vs {comp['mdp_titv_3']:.3f} (3rd)")
            print(f"  Random protein pattern: {comp['random_protein_pattern_proportion']:.1%}")

    return {
        'mdp_results': mdp_results,
        'random_results': random_results,
        'comparison': comparison_results,
        'global_mdp_stats': global_mdp_stats,
        'global_random_stats': global_random_stats,
        'parameters': {
            'n_iterations': n_iterations,
            'seed': seed,
            'rrna_gene': rrna_seq_col,
            'threshold': threshold,
            'n_mdps_analyzed': len(mdp_results),
            'excluded_ranges': excluded_ranges
        }
    }


def analyze_canonical_gene_titv(df: pd.DataFrame,
                                gene_configs: List[Dict[str, any]],
                                threshold: float = 0.05,
                                verbose: bool = True) -> Dict:
    """
    Analyze Ti/Tv ratios and nucleotide diversity for canonical protein coding genes.

    This function calculates Ti/Tv ratios for each codon position (first, second, third)
    and nucleotide diversity for canonical genes like ND5, ND1, COX1, etc.
    This serves as a baseline to compare against overlapping reading frames or MDPs.

    Parameters:
    -----------
    df : pd.DataFrame
        DataFrame with sequence data
    gene_configs : List[Dict]
        List of gene configuration dictionaries, each containing:
        - 'name': Name of the gene (e.g., 'ND5', 'ND1', 'COX1')
        - 'seq_col': Column name with gene sequences
    threshold : float
        Threshold for clean_sequences_adaptive (default=0.05)
    verbose : bool
        Print progress messages (default=True)

    Returns:
    --------
    Dict containing:
        - gene_results: Dict with Ti/Tv ratios and diversity for each gene at each codon position
        - parameters: Analysis parameters used

    Example:
    --------
    >>> gene_configs = [
    ...     {'name': 'ND5', 'seq_col': 'ND5'},
    ...     {'name': 'ND1', 'seq_col': 'ND1'},
    ...     {'name': 'COX1', 'seq_col': 'COX1'}
    ... ]
    >>> results = analyze_canonical_gene_titv(df, gene_configs)
    """
    if verbose:
        print("=" * 80)
        print("CANONICAL GENE Ti/Tv AND DIVERSITY ANALYSIS")
        print("=" * 80)
        print(f"\nAnalyzing {len(gene_configs)} canonical protein coding genes")

    gene_results = {}

    for config in gene_configs:
        gene_name = config['name']
        seq_col = config['seq_col']

        if verbose:
            print(f"\n{'='*60}")
            print(f"Analyzing Gene: {gene_name}")
            print(f"Sequence column: {seq_col}")
            print(f"{'='*60}")

        # Check if column exists
        if seq_col not in df.columns:
            warnings.warn(f"Column '{seq_col}' not found in DataFrame. Skipping {gene_name}")
            continue

        # Extract gene sequences
        gene_seqs_raw = df[seq_col].dropna().tolist()

        if len(gene_seqs_raw) < 2:
            if verbose:
                print(f"  Skipping {gene_name}: insufficient sequences (n={len(gene_seqs_raw)})")
            continue

        # Clean gene sequences (they are protein-coding)
        if verbose:
            print(f"  Cleaning {len(gene_seqs_raw)} sequences (threshold={threshold})...")

        cleaned_gene = cag.clean_sequences_adaptive(gene_seqs_raw, threshold=threshold)
        gene_seqs = cleaned_gene.cleaned_sequences

        if len(gene_seqs) < 2:
            if verbose:
                print(f"  Skipping {gene_name}: insufficient sequences after cleaning (n={len(gene_seqs)})")
            continue

        # Ensure sequences are in frame (divisible by 3)
        gene_seqs = [seq[:len(seq) - (len(seq) % 3)] for seq in gene_seqs if len(seq) >= 3]

        if len(gene_seqs) < 2:
            if verbose:
                print(f"  Skipping {gene_name}: insufficient sequences after frame adjustment")
            continue

        if verbose:
            from Bio.Seq import Seq

            print(f"  Analyzing {len(gene_seqs)} sequences, length={len(gene_seqs[0])} bp")

            # Report first and last 3 codons for verification
            seq_length = len(gene_seqs[0])
            n_codons = seq_length // 3

            if n_codons >= 6:
                # Get first 3 codons from first sequence
                first_3_codons = [gene_seqs[0][i*3:(i+1)*3] for i in range(3)]
                # Get last 3 codons from first sequence
                last_3_codons = [gene_seqs[0][i*3:(i+1)*3] for i in range(n_codons-3, n_codons)]

                # Translate codons to amino acids using BioPython using the standard mitochondrial code
                first_3_aa = [str(Seq(codon).translate(table="Vertebrate Mitochondrial")) for codon in first_3_codons]
                last_3_aa = [str(Seq(codon).translate(table="Vertebrate Mitochondrial")) for codon in last_3_codons]

                print(f"  First 3 codons: {' '.join(first_3_codons)} -> {' '.join(first_3_aa)}")
                print(f"  Last 3 codons:  {' '.join(last_3_codons)} -> {' '.join(last_3_aa)}")
            else:
                # If fewer than 6 codons, just show all
                all_codons = [gene_seqs[0][i*3:(i+1)*3] for i in range(n_codons)]
                all_aa = [str(Seq(codon).translate()) for codon in all_codons]
                print(f"  All {n_codons} codons: {' '.join(all_codons)} -> {' '.join(all_aa)}")

        # Calculate Ti/Tv for each codon position
        gene_titv = {}
        for pos in [0, 1, 2]:
            pos_name = ['first', 'second', 'third'][pos]
            result = calculate_codon_position_titv(gene_seqs, position=pos)
            gene_titv[pos_name] = result

            if verbose:
                if result['total_substitutions'] > 0:
                    print(f"  {pos_name:6s} codon position: Ti/Tv = {result['ti_tv_ratio']:6.3f} "
                          f"(Ti={result['ti_count']:3d}, Tv={result['tv_count']:3d}, "
                          f"Total={result['total_substitutions']:4d}, "
                          f"Diversity={result['nucleotide_diversity']:6.4f})")
                else:
                    print(f"  {pos_name:6s} codon position: No substitutions found")

        # Calculate combined first + second position Ti/Tv
        ti_12 = gene_titv['first']['ti_count'] + gene_titv['second']['ti_count']
        tv_12 = gene_titv['first']['tv_count'] + gene_titv['second']['tv_count']
        titv_12 = (ti_12 + .5) / (tv_12 + .5)  # With 0.5 pseudocount
        titv_12_raw = ti_12 / tv_12 if tv_12 > 0 else (np.inf if ti_12 > 0 else np.nan)

        # Average diversity across first and second positions
        diversity_12 = (gene_titv['first']['nucleotide_diversity'] +
                       gene_titv['second']['nucleotide_diversity']) / 2

        if verbose:
            print(f"  Combined (1st+2nd): Ti/Tv = {titv_12:6.3f} (Ti={ti_12:3d}, Tv={tv_12:3d})")
            print(f"  Nucleotide diversity (1st+2nd): {diversity_12:6.4f}")
            print(f"  Nucleotide diversity (3rd): {gene_titv['third']['nucleotide_diversity']:6.4f}")

        # Check for protein-like pattern (Ti/Tv ratio: 1st+2nd < 3rd)
        protein_pattern = titv_12 < gene_titv['third']['ti_tv_ratio'] if (
            np.isfinite(titv_12) and np.isfinite(gene_titv['third']['ti_tv_ratio'])
        ) else False

        if verbose:
            if protein_pattern:
                print(f"  ✓ Shows expected protein-coding pattern (Ti/Tv 1st+2nd < 3rd)")
            else:
                print(f"  ⚠ Does NOT show expected protein-coding pattern")

        gene_results[gene_name] = {
            'first_pos': gene_titv['first'],
            'second_pos': gene_titv['second'],
            'third_pos': gene_titv['third'],
            'combined_12_titv': titv_12,
            'combined_12_titv_raw': titv_12_raw,
            'combined_12_diversity': diversity_12,
            'length': len(gene_seqs[0]),
            'n_sequences': len(gene_seqs),
            'n_codons': len(gene_seqs[0]) // 3,
            'shows_protein_pattern': protein_pattern,
            'titv_diff': gene_titv['third']['ti_tv_ratio'] - titv_12 if (
                np.isfinite(titv_12) and np.isfinite(gene_titv['third']['ti_tv_ratio'])
            ) else np.nan
        }

    # Calculate global Ti/Tv ratios across all genes
    if len(gene_results) > 0:
        global_gene_ti_first = sum(res['first_pos']['ti_count'] for res in gene_results.values())
        global_gene_tv_first = sum(res['first_pos']['tv_count'] for res in gene_results.values())
        global_gene_ti_second = sum(res['second_pos']['ti_count'] for res in gene_results.values())
        global_gene_tv_second = sum(res['second_pos']['tv_count'] for res in gene_results.values())
        global_gene_ti_third = sum(res['third_pos']['ti_count'] for res in gene_results.values())
        global_gene_tv_third = sum(res['third_pos']['tv_count'] for res in gene_results.values())

        # Combined first + second
        global_gene_ti_12 = global_gene_ti_first + global_gene_ti_second
        global_gene_tv_12 = global_gene_tv_first + global_gene_tv_second

        # Calculate global Ti/Tv ratios with pseudocount
        global_gene_titv_first = (global_gene_ti_first + 0.5) / (global_gene_tv_first + 0.5)
        global_gene_titv_second = (global_gene_ti_second + 0.5) / (global_gene_tv_second + 0.5)
        global_gene_titv_third = (global_gene_ti_third + 0.5) / (global_gene_tv_third + 0.5)
        global_gene_titv_12 = (global_gene_ti_12 + 0.5) / (global_gene_tv_12 + 0.5)

        # Calculate RoA for diversities (average)
        global_gene_diversity_first = np.mean([res['first_pos']['nucleotide_diversity'] for res in gene_results.values()])
        global_gene_diversity_second = np.mean([res['second_pos']['nucleotide_diversity'] for res in gene_results.values()])
        global_gene_diversity_third = np.mean([res['third_pos']['nucleotide_diversity'] for res in gene_results.values()])
        global_gene_diversity_12 = np.mean([res['combined_12_diversity'] for res in gene_results.values()])

        # Check if global pattern is protein-like
        global_gene_protein_pattern = global_gene_titv_12 < global_gene_titv_third

        global_gene_stats = {
            'first_pos_titv': global_gene_titv_first,
            'second_pos_titv': global_gene_titv_second,
            'third_pos_titv': global_gene_titv_third,
            'combined_12_titv': global_gene_titv_12,
            'first_pos_ti': global_gene_ti_first,
            'first_pos_tv': global_gene_tv_first,
            'second_pos_ti': global_gene_ti_second,
            'second_pos_tv': global_gene_tv_second,
            'third_pos_ti': global_gene_ti_third,
            'third_pos_tv': global_gene_tv_third,
            'combined_12_ti': global_gene_ti_12,
            'combined_12_tv': global_gene_tv_12,
            'roa_diversity_first': global_gene_diversity_first,
            'roa_diversity_second': global_gene_diversity_second,
            'roa_diversity_third': global_gene_diversity_third,
            'roa_diversity_12': global_gene_diversity_12,
            'shows_protein_pattern': global_gene_protein_pattern
        }
    else:
        global_gene_stats = None

    # Summary statistics
    if verbose:
        print(f"\n{'='*80}")
        print("SUMMARY")
        print(f"{'='*80}")

        if len(gene_results) > 0:
            n_protein_pattern = sum(1 for res in gene_results.values() if res['shows_protein_pattern'])
            print(f"\nGenes analyzed: {len(gene_results)}")
            print(f"Genes showing protein-coding pattern: {n_protein_pattern}/{len(gene_results)}")

            # Summary table
            print(f"\n{'Gene':<12} {'Ti/Tv(1+2)':>12} {'Ti/Tv(3)':>12} {'Div(1+2)':>12} {'Div(3)':>12} {'Pattern':>10}")
            print("-" * 80)
            for gene_name, res in gene_results.items():
                pattern = "✓" if res['shows_protein_pattern'] else "✗"
                print(f"{gene_name:<12} {res['combined_12_titv']:>12.3f} "
                      f"{res['third_pos']['ti_tv_ratio']:>12.3f} "
                      f"{res['combined_12_diversity']:>12.4f} "
                      f"{res['third_pos']['nucleotide_diversity']:>12.4f} "
                      f"{pattern:>10}")

            # Average values
            avg_titv_12 = np.mean([res['combined_12_titv'] for res in gene_results.values()
                                   if np.isfinite(res['combined_12_titv'])])
            avg_titv_3 = np.mean([res['third_pos']['ti_tv_ratio'] for res in gene_results.values()
                                  if np.isfinite(res['third_pos']['ti_tv_ratio'])])
            avg_div_12 = np.mean([res['combined_12_diversity'] for res in gene_results.values()
                                  if np.isfinite(res['combined_12_diversity'])])
            avg_div_3 = np.mean([res['third_pos']['nucleotide_diversity'] for res in gene_results.values()
                                 if np.isfinite(res['third_pos']['nucleotide_diversity'])])

            print("-" * 80)
            print(f"{'AVERAGE':<12} {avg_titv_12:>12.3f} {avg_titv_3:>12.3f} "
                  f"{avg_div_12:>12.4f} {avg_div_3:>12.4f}")

            # Print global statistics
            print(f"\n{'='*80}")
            print("GLOBAL STATISTICS (ALL GENES COMBINED)")
            print(f"{'='*80}")
            print(f"  First position:  Ti/Tv = {global_gene_stats['first_pos_titv']:6.3f} (Ti={global_gene_stats['first_pos_ti']:4d}, Tv={global_gene_stats['first_pos_tv']:4d})")
            print(f"  Second position: Ti/Tv = {global_gene_stats['second_pos_titv']:6.3f} (Ti={global_gene_stats['second_pos_ti']:4d}, Tv={global_gene_stats['second_pos_tv']:4d})")
            print(f"  Third position:  Ti/Tv = {global_gene_stats['third_pos_titv']:6.3f} (Ti={global_gene_stats['third_pos_ti']:4d}, Tv={global_gene_stats['third_pos_tv']:4d})")
            print(f"  Combined (1st+2nd): Ti/Tv = {global_gene_stats['combined_12_titv']:6.3f} (Ti={global_gene_stats['combined_12_ti']:4d}, Tv={global_gene_stats['combined_12_tv']:4d})")
            print(f"  RoA Diversity (1st+2nd): {global_gene_stats['roa_diversity_12']:6.4f}")
            print(f"  RoA Diversity (3rd): {global_gene_stats['roa_diversity_third']:6.4f}")
            if global_gene_stats['shows_protein_pattern']:
                print(f"  ✓ Global pattern shows protein-coding signature (Ti/Tv 1st+2nd < 3rd)")
            else:
                print(f"  ✗ Global pattern does NOT show protein-coding signature")
        else:
            print("\nNo genes successfully analyzed.")

    return {
        'gene_results': gene_results,
        'global_gene_stats': global_gene_stats,
        'parameters': {
            'threshold': threshold,
            'n_genes_analyzed': len(gene_results),
            'gene_names': list(gene_results.keys())
        }
    }


def test_codon_position_diversity_difference(results: Dict,
                                             test_type: str = 'wilcoxon',
                                             verbose: bool = True) -> Dict:
    """
    Test if nucleotide diversity differs significantly between 1st+2nd vs 3rd codon positions.

    Uses paired tests since each gene has both diversity values.

    Parameters:
    -----------
    results : Dict
        Output from analyze_canonical_gene_titv()
    test_type : str
        'wilcoxon' (default, non-parametric) or 'ttest' (parametric)
    verbose : bool
        Print results

    Returns:
    --------
    Dict with test statistics and p-value containing:
        - test_name: Name of the statistical test used
        - test_type: Type of test ('wilcoxon' or 'ttest')
        - n_genes: Number of genes included in the test
        - statistic: Test statistic value
        - p_value: P-value from the test
        - significant: Boolean, True if p < 0.05
        - diversity_12_mean: Mean diversity for 1st+2nd positions
        - diversity_12_median: Median diversity for 1st+2nd positions
        - diversity_3_mean: Mean diversity for 3rd position
        - diversity_3_median: Median diversity for 3rd position
        - mean_difference: Mean of (3rd - [1st+2nd])
        - median_difference: Median of (3rd - [1st+2nd])
        - ci_95_lower: Lower bound of 95% CI for mean difference
        - ci_95_upper: Upper bound of 95% CI for mean difference
        - cohens_d: Effect size (Cohen's d)
        - effect_size_interpretation: Interpretation of effect size
        - n_3rd_greater: Number of genes where 3rd > 1st+2nd
        - n_3rd_lesser: Number of genes where 3rd < 1st+2nd
        - gene_names: List of gene names included
        - diversity_12: List of 1st+2nd diversity values
        - diversity_3: List of 3rd position diversity values
        - differences: List of differences (3rd - [1st+2nd])

    Example:
    --------
    >>> results = analyze_canonical_gene_titv(df, gene_configs)
    >>> test_results = test_codon_position_diversity_difference(results,
    ...                                                         test_type='wilcoxon',
    ...                                                         verbose=True)
    >>> if test_results['significant']:
    ...     print(f"Significant difference found (p={test_results['p_value']:.4f})")
    """
    from scipy import stats
    import numpy as np

    if 'gene_results' not in results:
        raise ValueError("Expected results from analyze_canonical_gene_titv()")

    gene_results = results['gene_results']

    # Extract paired diversity values
    diversity_12 = []
    diversity_3 = []
    gene_names = []

    for gene_name, gene_data in gene_results.items():
        div_12 = gene_data['combined_12_diversity']
        div_3 = gene_data['third_pos']['nucleotide_diversity']

        # Only include if both values are finite
        if np.isfinite(div_12) and np.isfinite(div_3):
            diversity_12.append(div_12)
            diversity_3.append(div_3)
            gene_names.append(gene_name)

    diversity_12 = np.array(diversity_12)
    diversity_3 = np.array(diversity_3)
    differences = diversity_3 - diversity_12  # 3rd - (1st+2nd)

    n_genes = len(diversity_12)

    if n_genes < 3:
        return {
            'test': None,
            'error': f'Insufficient genes for testing (n={n_genes})'
        }

    # Perform statistical test
    if test_type == 'wilcoxon':
        statistic, p_value = stats.wilcoxon(diversity_12, diversity_3,
                                            alternative='two-sided')
        test_name = 'Wilcoxon signed-rank test'
    elif test_type == 'ttest':
        statistic, p_value = stats.ttest_rel(diversity_12, diversity_3)
        test_name = 'Paired t-test'
    else:
        raise ValueError(f"Unknown test_type: {test_type}")

    # Calculate effect size (Cohen's d for paired data)
    mean_diff = np.mean(differences)
    std_diff = np.std(differences, ddof=1)
    cohens_d = mean_diff / std_diff if std_diff > 0 else np.nan

    # Calculate confidence interval for mean difference (95%)
    se_diff = std_diff / np.sqrt(n_genes)
    ci_95 = stats.t.interval(0.95, n_genes - 1, loc=mean_diff, scale=se_diff)

    # Descriptive statistics
    median_12 = np.median(diversity_12)
    median_3 = np.median(diversity_3)
    mean_12 = np.mean(diversity_12)
    mean_3 = np.mean(diversity_3)

    # Determine effect size interpretation
    if abs(cohens_d) >= 0.8:
        effect_size_interpretation = "large"
    elif abs(cohens_d) >= 0.5:
        effect_size_interpretation = "medium"
    elif abs(cohens_d) >= 0.2:
        effect_size_interpretation = "small"
    else:
        effect_size_interpretation = "negligible"

    if verbose:
        print("=" * 80)
        print("CODON POSITION DIVERSITY COMPARISON")
        print("=" * 80)
        print(f"\nSample size: {n_genes} genes")
        print(f"\n1st+2nd Position Diversity:")
        print(f"  Mean:   {mean_12:.6f}")
        print(f"  Median: {median_12:.6f}")
        print(f"  Range:  {np.min(diversity_12):.6f} - {np.max(diversity_12):.6f}")

        print(f"\n3rd Position Diversity:")
        print(f"  Mean:   {mean_3:.6f}")
        print(f"  Median: {median_3:.6f}")
        print(f"  Range:  {np.min(diversity_3):.6f} - {np.max(diversity_3):.6f}")

        print(f"\nDifference (3rd - [1st+2nd]):")
        print(f"  Mean difference: {mean_diff:+.6f}")
        print(f"  95% CI: [{ci_95[0]:+.6f}, {ci_95[1]:+.6f}]")
        print(f"  Median difference: {np.median(differences):+.6f}")

        print(f"\n{test_name}:")
        print(f"  Test statistic: {statistic:.4f}")
        print(f"  P-value: {p_value:.6f}")

        if p_value < 0.001:
            sig_str = "*** (p < 0.001)"
        elif p_value < 0.01:
            sig_str = "** (p < 0.01)"
        elif p_value < 0.05:
            sig_str = "* (p < 0.05)"
        else:
            sig_str = "n.s."
        print(f"  Significance: {sig_str}")

        print(f"\nEffect size (Cohen's d): {cohens_d:.3f}")
        print(f"  Interpretation: {effect_size_interpretation}")

        # Per-gene breakdown
        print(f"\n{'Gene':<12} {'Div(1+2)':>12} {'Div(3)':>12} {'Diff':>12} {'Direction':>12}")
        print("-" * 80)
        for i, gene in enumerate(gene_names):
            diff = differences[i]
            direction = "3rd > 1+2" if diff > 0 else ("3rd < 1+2" if diff < 0 else "Equal")
            print(f"{gene:<12} {diversity_12[i]:>12.6f} {diversity_3[i]:>12.6f} "
                  f"{diff:>+12.6f} {direction:>12}")

        print(f"\nGenes with 3rd > 1st+2nd: {np.sum(differences > 0)}/{n_genes}")
        print(f"Genes with 3rd < 1st+2nd: {np.sum(differences < 0)}/{n_genes}")

    return {
        'test_name': test_name,
        'test_type': test_type,
        'n_genes': n_genes,
        'statistic': statistic,
        'p_value': p_value,
        'significant': p_value < 0.05,
        'diversity_12_mean': mean_12,
        'diversity_12_median': median_12,
        'diversity_3_mean': mean_3,
        'diversity_3_median': median_3,
        'mean_difference': mean_diff,
        'median_difference': np.median(differences),
        'ci_95_lower': ci_95[0],
        'ci_95_upper': ci_95[1],
        'cohens_d': cohens_d,
        'effect_size_interpretation': effect_size_interpretation,
        'n_3rd_greater': int(np.sum(differences > 0)),
        'n_3rd_lesser': int(np.sum(differences < 0)),
        'gene_names': gene_names,
        'diversity_12': diversity_12.tolist(),
        'diversity_3': diversity_3.tolist(),
        'differences': differences.tolist()
    }


def test_mdp_diversity_permutation(results: Dict,
                                   n_permutations: int = 10000,
                                   seed: int = None,
                                   verbose: bool = True) -> Dict:
    """
    Perform one-sided permutation test on MDP nucleotide diversity.

    Tests whether the difference in nucleotide diversity (3rd - [1st+2nd]) for MDPs
    is significantly LARGER than for random rRNA regions using a permutation test.

    H0: The diversity difference for MDPs equals that of random regions
    H1: The diversity difference for MDPs is GREATER than random regions
        (indicating stronger codon position bias, expected for protein-coding)

    Supports both single and combined dictionaries:
    - Single: output from analyze_rrna_mdp_asymmetry()
    - Combined: {'rnr1_results': results1, 'rnr2_results': results2, ...}

    Parameters:
    -----------
    results : Dict
        Output from analyze_rrna_mdp_asymmetry() or combined dictionary with multiple result sets
    n_permutations : int
        Number of permutations (default=10000)
    seed : int, optional
        Random seed for reproducibility
    verbose : bool
        Print detailed results

    Returns:
    --------
    Dict containing for each MDP:
        - mdp_diversity_diff: Observed difference (3rd - [1st+2nd]) for MDP
        - random_diversity_diffs: List of differences for random regions
        - random_mean_diff: Mean difference for random regions
        - random_median_diff: Median difference for random regions
        - p_value_empirical: Empirical p-value (proportion random >= observed)
        - p_value_permutation: Permutation test p-value
        - significant_empirical: Boolean (p_empirical < 0.05)
        - significant_permutation: Boolean (p_permutation < 0.05)
        - n_permutations: Number of permutations used

    For combined dictionaries, returns a nested structure matching the input.

    Example:
    --------
    >>> # Single result
    >>> results = analyze_rrna_mdp_asymmetry(df, 'RNR2', mdp_configs)
    >>> test_results = test_mdp_diversity_permutation(results, n_permutations=10000)
    >>>
    >>> # Combined results
    >>> combined = {'rnr1_results': results1, 'rnr2_results': results2}
    >>> test_results = test_mdp_diversity_permutation(combined, n_permutations=10000)
    """
    from scipy import stats
    import numpy as np

    # Detect if this is a combined dictionary
    is_combined = not any(key in results for key in ['mdp_results', 'random_results'])

    if is_combined:
        # Process each sub-dictionary separately and return combined results
        if verbose:
            print("=" * 80)
            print("MDP DIVERSITY PERMUTATION TEST (COMBINED)")
            print("=" * 80)

        combined_test_results = {}
        for source_key, source_results in results.items():
            if verbose:
                print(f"\n{'='*80}")
                print(f"Processing: {source_key}")
                print(f"{'='*80}")

            # Recursively call this function on each sub-result
            sub_test_results = test_mdp_diversity_permutation(
                source_results,
                n_permutations=n_permutations,
                seed=seed,
                verbose=verbose
            )
            combined_test_results[source_key] = sub_test_results

        return combined_test_results

    # Single dictionary - original processing
    if 'mdp_results' not in results or 'random_results' not in results:
        raise ValueError("Expected results from analyze_rrna_mdp_asymmetry()")

    mdp_results = results['mdp_results']
    random_results = results['random_results']

    if seed is not None:
        np.random.seed(seed)

    if verbose:
        print("=" * 80)
        print("MDP DIVERSITY PERMUTATION TEST")
        print("=" * 80)
        print(f"\nTesting H1: MDP diversity difference (3rd - [1st+2nd]) > Random")
        print(f"Number of permutations: {n_permutations}")

    test_results = {}

    for mdp_name in mdp_results.keys():
        if mdp_name not in random_results:
            if verbose:
                print(f"\n{mdp_name}: No random results available, skipping")
            continue

        mdp = mdp_results[mdp_name]
        random = random_results[mdp_name]

        # Get MDP diversity values
        mdp_div_12 = mdp['combined_12_diversity']
        mdp_div_3 = mdp['third_pos']['nucleotide_diversity']
        mdp_diversity_diff = mdp_div_3 - mdp_div_12

        # Get random region diversity values
        random_div_12_list = random['combined_12_diversity']
        random_div_3_list = random['third_diversity']

        # Calculate differences for each random region
        random_diversity_diffs = []
        for i in range(min(len(random_div_12_list), len(random_div_3_list))):
            if np.isfinite(random_div_12_list[i]) and np.isfinite(random_div_3_list[i]):
                diff = random_div_3_list[i] - random_div_12_list[i]
                random_diversity_diffs.append(diff)

        if len(random_diversity_diffs) == 0:
            if verbose:
                print(f"\n{mdp_name}: No valid random differences, skipping")
            continue

        random_diversity_diffs = np.array(random_diversity_diffs)
        random_mean_diff = np.mean(random_diversity_diffs)
        random_median_diff = np.median(random_diversity_diffs)

        # Empirical p-value (simple comparison to random distribution)
        n_equal_or_greater = np.sum(random_diversity_diffs >= mdp_diversity_diff)
        p_value_empirical = n_equal_or_greater / len(random_diversity_diffs)

        # Permutation test
        # Pool all diversity values (MDP + all random regions)
        # For 1st+2nd position
        all_div_12 = np.array([mdp_div_12] + list(random_div_12_list))
        # For 3rd position
        all_div_3 = np.array([mdp_div_3] + list(random_div_3_list))

        # Ensure same length (in case random lists have different lengths)
        min_len = min(len(all_div_12), len(all_div_3))
        all_div_12 = all_div_12[:min_len]
        all_div_3 = all_div_3[:min_len]

        # Calculate observed test statistic
        # Test statistic: difference between MDP's (3rd - [1st+2nd]) and mean random (3rd - [1st+2nd])
        observed_stat = mdp_diversity_diff - random_mean_diff

        # Permutation distribution
        perm_stats = []
        for _ in range(n_permutations):
            # Randomly select one observation to be "MDP" (rest are "random")
            perm_idx = np.random.randint(0, min_len)

            # Calculate "MDP" difference for this permutation
            perm_mdp_diff = all_div_3[perm_idx] - all_div_12[perm_idx]

            # Calculate mean "random" difference for this permutation
            perm_random_mask = np.ones(min_len, dtype=bool)
            perm_random_mask[perm_idx] = False
            perm_random_diffs = all_div_3[perm_random_mask] - all_div_12[perm_random_mask]
            perm_random_mean = np.mean(perm_random_diffs)

            # Test statistic for this permutation
            perm_stat = perm_mdp_diff - perm_random_mean
            perm_stats.append(perm_stat)

        perm_stats = np.array(perm_stats)

        # One-sided p-value: proportion of permutations with stat >= observed
        p_value_permutation = np.sum(perm_stats >= observed_stat) / n_permutations

        if verbose:
            print(f"\n{'-'*60}")
            print(f"MDP: {mdp_name}")
            print(f"{'-'*60}")
            print(f"  MDP diversity (1st+2nd): {mdp_div_12:.6f}")
            print(f"  MDP diversity (3rd):     {mdp_div_3:.6f}")
            print(f"  MDP difference:          {mdp_diversity_diff:+.6f}")
            print(f"\n  Random regions (n={len(random_diversity_diffs)}):")
            print(f"    Mean difference:   {random_mean_diff:+.6f}")
            print(f"    Median difference: {random_median_diff:+.6f}")
            print(f"    Std difference:    {np.std(random_diversity_diffs):.6f}")
            print(f"\n  Empirical p-value:   {p_value_empirical:.4f} {'*' if p_value_empirical < 0.05 else ''}")
            print(f"  Permutation p-value: {p_value_permutation:.4f} {'*' if p_value_permutation < 0.05 else ''}")

            if p_value_permutation < 0.05:
                print(f"  ✓ MDP shows significantly stronger codon position bias (p < 0.05)")
            else:
                print(f"  ✗ No significant difference from random regions")

        test_results[mdp_name] = {
            'mdp_diversity_12': mdp_div_12,
            'mdp_diversity_3': mdp_div_3,
            'mdp_diversity_diff': mdp_diversity_diff,
            'random_diversity_diffs': random_diversity_diffs.tolist(),
            'random_mean_diff': random_mean_diff,
            'random_median_diff': random_median_diff,
            'random_std_diff': np.std(random_diversity_diffs),
            'p_value_empirical': p_value_empirical,
            'p_value_permutation': p_value_permutation,
            'significant_empirical': p_value_empirical < 0.05,
            'significant_permutation': p_value_permutation < 0.05,
            'n_permutations': n_permutations,
            'n_random_regions': len(random_diversity_diffs),
            'observed_test_statistic': observed_stat
        }

    if verbose:
        print(f"\n{'='*80}")
        print("SUMMARY")
        print(f"{'='*80}")
        n_sig_empirical = sum(1 for res in test_results.values() if res['significant_empirical'])
        n_sig_permutation = sum(1 for res in test_results.values() if res['significant_permutation'])
        print(f"MDPs with significant diversity bias (empirical test):   {n_sig_empirical}/{len(test_results)}")
        print(f"MDPs with significant diversity bias (permutation test): {n_sig_permutation}/{len(test_results)}")

    return test_results


if __name__ == "__main__":
    print("Overlapping Reading Frame Asymmetry Analysis Module")
    print("=" * 80)
    print("\nThis module provides functions for analyzing asymmetric evolution")
    print("in overlapping reading frames, following Szklarczyk et al. 2007.")
    print("\nMain function: analyze_overlap_asymmetry()")
    print("\nSee docstrings for usage examples.")