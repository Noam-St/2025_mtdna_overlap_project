import numpy as np
from collections import Counter, defaultdict
from Bio.Data import CodonTable
from Bio.Seq import Seq

def inspect_alignment_quality(sequences, n_rows=10, n_cols=50):
    """
    Visual check of the first few rows and columns of the alignment.
    """
    print(f"--- Inspecting Alignment (First {n_rows} seqs, first {n_cols} bases) ---")
    for i, seq in enumerate(sequences[:n_rows]):
        print(f"Seq {i}: {seq[:n_cols]}")

    # Check for length consistency
    lengths = set(len(s) for s in sequences)
    print(f"\nSequence Lengths present in data: {lengths}")
    if len(lengths) > 1:
        print("CRITICAL WARNING: Sequences have different lengths! They are not aligned.")
    else:
        print("Lengths are consistent.")

class ARFomeNormalized:
    """
    Analyzes selection pressure in overlapping reading frames by calculating
    normalized substitution rates (similar to dN/dS).
    
    Refactored to support:
    1. Dynamic Ti/Tv ratio estimation from input data.
    2. Weighted potential sites (denominator) based on Ti/Tv.
    3. Counting segregating sites (removing haplogroup frequency bias).
    4. Opposite strand overlaps.
    """
    
    def __init__(self, genetic_code_id=2):
        # Load mitochondrial code (id=2 for vertebrate mito)
        self.table = CodonTable.unambiguous_dna_by_id[genetic_code_id]
        self.forward_table = self.table.forward_table
        self.bases = ['A', 'C', 'G', 'T']
        
        # Define Transitions for fast lookup
        self.transitions = {('A', 'G'), ('G', 'A'), ('C', 'T'), ('T', 'C')}

    def _classify_mutation(self, codon_f0, new_codon_f0, 
                          flanking_u, flanking_d, frame_shift, alternative_strand='+', verbose = False):
        """
        Classifies a mutation into SS, SN, NS, NN, or STOP based on its
        effect on both the Canonical (F0) and Alternative (F1) frames.
        
        Handles opposite strand overlaps by reverse complementing the context.
        """
        # Make sure that the codons only contain valid bases
        assert all(base in self.bases for base in codon_f0), f"Invalid base in codon_f0: {codon_f0}"
        assert all(base in self.bases for base in new_codon_f0), f"Invalid base in new_codon_f0: {new_codon_f0}"

        # 1. Analyze Frame 0 (Canonical) - Always Forward
        aa_f0_old = self.forward_table.get(codon_f0, '*')
        aa_f0_new = self.forward_table.get(new_codon_f0, '*')
        
        if aa_f0_new == '*': return 'STOP_F0' 
        is_syn_f0 = (aa_f0_old == aa_f0_new)

        if verbose:
            print(f"\nFrame 0 (Canonical): {codon_f0} -> {new_codon_f0}")
            print(f"Frame 0 (Canonical): {aa_f0_old} -> {aa_f0_new}")
            print(f"Frame 0 (Canonical): {is_syn_f0}")

        # 2. Analyze Frame 1 (Alternative)
        full_seq_old = flanking_u + codon_f0 + flanking_d
        full_seq_new = flanking_u + new_codon_f0 + flanking_d

    
        # --- SAME STRAND LOGIC (Original) ---
        # Identify mutation index inside the F0 codon
        mutation_idx_in_codon = -1
        for i in range(3):
            if codon_f0[i] != new_codon_f0[i]:
                mutation_idx_in_codon = i
                break
        
        # Mutation not found in codon
        if mutation_idx_in_codon == -1: return 'SS' 
        
        mutation_abs_idx = 2 + mutation_idx_in_codon
        
        # Align Frame 1 based on shift
        if frame_shift == 0:
            f1_start = mutation_abs_idx
        elif frame_shift == 1:
            f1_start = (mutation_abs_idx // 3) * 3  
        elif frame_shift == 2:
            f1_start = ((mutation_abs_idx - 1) // 3) * 3 + 1

        # Reverse Complement the codon if opposite strand
        if alternative_strand == '-':
            codon_f1_old = Seq(full_seq_old[f1_start : f1_start+3]).reverse_complement()
            codon_f1_new = Seq(full_seq_new[f1_start : f1_start+3]).reverse_complement() 
        else:
            codon_f1_old = Seq(full_seq_old[f1_start : f1_start+3])
            codon_f1_new = Seq(full_seq_new[f1_start : f1_start+3])
        
        assert len(codon_f1_old) == 3, f"Invalid codon length in codon_f1_old: {len(codon_f1_old)}"
        assert len(codon_f1_new) == 3, f"Invalid codon length in codon_f1_new: {len(codon_f1_new)}"
        
        if len(codon_f1_old) < 3: return 'UNKNOWN'
        # Translate F1
        aa_f1_old = self.forward_table.get(codon_f1_old, '*')
        aa_f1_new = self.forward_table.get(codon_f1_new, '*')

        # Translate F1
        aa_f1_old = self.forward_table.get(codon_f1_old, '*')
        aa_f1_new = self.forward_table.get(codon_f1_new, '*')

        # --- CORRECTED LOGIC ---
        if aa_f1_new == '*': 
            if is_syn_f0:
                # Case 1: Synonymous in Canonical, Stop in Alt.
                # This is a "Pure Stop" - the exact signal we want to test.
                return 'STOP_F1'  
            else:
                # Case 2: Nonsynonymous in Canonical, Stop in Alt.
                # The canonical gene filters this out. It acts like a double-hit.
                # We classify as NN (Nonsynonymous in both) to remove it from beta_stop
                # while keeping the mutation counts logically consistent.
                return 'NN'       
        
        is_syn_f1 = (aa_f1_old == aa_f1_new)

        if verbose:
            print(f"\nFrame 1 (Alternative): {codon_f1_old} -> {codon_f1_new}")
            print(f"Frame 1 (Alternative): {aa_f1_old} -> {aa_f1_new}")
            print(f"Frame 1 (Alternative): {is_syn_f1}")

        # 3. Determine Category
        if is_syn_f0 and is_syn_f1: return 'SS'    
        if is_syn_f0 and not is_syn_f1: return 'SN' 
        if not is_syn_f0 and is_syn_f1: return 'NS' 
        if not is_syn_f0 and not is_syn_f1: return 'NN' 
        return 'UNKNOWN'

    def calculate_potential_sites(self, consensus_str, window_start, window_end, 
                                  frame_shift, titv_ratio, alternative_strand='+', verbose = False):
        """
        Calculates the denominator for rates.
        NOW WEIGHTED: Transitions add 'titv_ratio', Transversions add 1.0.
        """
        potential_counts = defaultdict(float)
        start_codon_idx = (window_start // 3) * 3
        
        for i in range(start_codon_idx, window_end, 3):
            if i < window_start or i+3 > window_end: continue
            
            codon = consensus_str[i:i+3]
            u = consensus_str[max(0, i-2):i]
            d = consensus_str[i+3:i+5]
            
            if len(u) < 2 or len(d) < 2: continue
            
            # Simulate every possible mutation
            for pos in range(3):
                original_base = codon[pos]
                for mut_base in self.bases:
                    if mut_base == original_base: continue
                    
                    # Create mutated codon
                    mut_codon_list = list(codon)
                    mut_codon_list[pos] = mut_base
                    mut_codon = "".join(mut_codon_list)
                    
                    m_type = self._classify_mutation(codon, mut_codon, u, d, 
                                                   frame_shift, alternative_strand, verbose=False if verbose <= 1 else verbose-1)
                    
                    # --- WEIGHTING LOGIC ---
                    # Check if this theoretical mutation is a Transition or Transversion
                    # Since transitions are much more likely to occur than transversions, the potential for transitions is weighted down by a factor of 'titv_ratio' (potentials are denominators when calculating rates)
                    is_transition = (original_base, mut_base) in self.transitions
                    weight = titv_ratio if is_transition else 1.0
                    
                    potential_counts[m_type] += weight
                    
        return potential_counts

    def get_observed_counts(self, canonical_list, alternative_list, 
                           arf_start_index, frame_shift, count_unique_sites=True,
                           alternative_strand='+', verbose = False, ref_seq = None):
        """
        Counts actual mutations.
        Returns: counts, consensus_sequence, and raw Ti/Tv counts for estimation.
        """
        mat_can = np.array([list(s.upper()) for s in canonical_list])
        _, len_can = mat_can.shape
        
        # Consensus generation
        if ref_seq is not None:
            cons_can_str = ref_seq.upper()
        else:
            cons_can = [Counter(mat_can[:,i]).most_common(1)[0][0] for i in range(len_can)]
            cons_can_str = "".join(cons_can)
        
        window_start = arf_start_index
        # Use actual alt length for window end
        window_end = arf_start_index + len(alternative_list[0]) 
        
        observed_counts = defaultdict(int)
        mutation_details = defaultdict(list)  # Track individual mutations with their counts
        n_transitions = 0
        n_transversions = 0

        # Valid DNA bases for Ti/Tv Calculation
        valid_dna = {'A', 'C', 'G', 'T'}

        start_codon_idx = (window_start // 3) * 3

        if verbose:
            print(f"\nConsensus Sequence: {cons_can_str}")
            print(f"Alternative Sequence: {cons_can_str[start_codon_idx:window_end]}")
            print(f"Window Start: {window_start}")
            print(f"Window End: {window_end}")
            print(f"Start Codon Index: {start_codon_idx}")
        
        for i in range(start_codon_idx, window_end, 3):
            if i < window_start or i+3 > window_end: 
                if verbose:
                    print(f"Skipping index {i} (out of bounds)")
                continue

            consensus_codon = cons_can_str[i:i+3]
            u = cons_can_str[max(0, i-2):i] 
            d = cons_can_str[i+3:i+5]       
            
            # Skip edge cases with insufficient context
            if len(u) < 2 or len(d) < 2: 
                if verbose:
                    print(f"Skipping index {i} ({consensus_codon}={variant_codon}) (insufficient context)")
                continue 
            
            variants_at_codon = mat_can[:, i:i+3]
            unique_rows, counts = np.unique(variants_at_codon, axis=0, return_counts=True)
            if verbose:
                pass
                #print(f"Unique Rows: {unique_rows}")
                #print(f"Counts: {counts}")
            
            for row, count in zip(unique_rows, counts): # Each unique codon variant
                variant_codon = "".join(row)
                if variant_codon == consensus_codon: 
                    if verbose:
                        print(f"Skipping index {i} ({consensus_codon}={variant_codon}) (same as consensus)")
                    continue
                
                bases_old = []
                bases_new = []
                # Check single nucleotide mutations - This needs to work for cases where there is more than one mutation in a codon.
                diffs = sum(1 for a, b in zip(consensus_codon, variant_codon) if a != b)
                if diffs >= 1:
                    # Identify specific base change
                    for k in range(3):
                        if consensus_codon[k] != variant_codon[k]:
                            bases_old.append(consensus_codon[k])
                            bases_new.append(variant_codon[k])
                    
                    # --- FIXED Ti/Tv LOGIC ---
                    # Only count strictly valid DNA bases. Ignore N, -, R, Y, etc.
                    if all(base_old in valid_dna and base_new in valid_dna for base_old, base_new in zip(bases_old, bases_new)):
                        for base_old, base_new in zip(bases_old, bases_new):
                            if (base_old, base_new) in self.transitions:
                                n_transitions += 1
                            else:
                                n_transversions += 1
                    # -------------------------
                    # This method classifies the observed mutation in both the canonical and the alternative frames.
                    m_type = self._classify_mutation(consensus_codon, variant_codon,
                                                   u, d, frame_shift, alternative_strand, verbose=False if verbose <= 1 else verbose-1)

                    if verbose > 1:
                        print(f"\nIndex {i}: {consensus_codon} -> {variant_codon} | Type: {m_type} | Count: {count}")
                    increment = 1 if count_unique_sites else count
                    observed_counts[m_type] += increment

                    # Track individual mutation details (always use actual count for frequency calculations)
                    mutation_details[m_type].append({
                        'position': i,
                        'consensus_codon': consensus_codon,
                        'variant_codon': variant_codon,
                        'mutation': f"{consensus_codon}->{variant_codon}",
                        'count': count  # Number of individuals with this specific mutation
                    })

        return observed_counts, cons_can_str, n_transitions, n_transversions, mutation_details

    def run_analysis(self, canonical_seqs, alternative_seqs, 
                    arf_start_index, frame_shift, default_titv=10.0,
                    alternative_strand='+', verbose = False, ref_seq=None):
        """
        Main execution. Automatically calculates Ti/Tv from data to weight the denominator.
        """
        print(f"Analyzing Nested ARF (Start: {arf_start_index}, Shift: {frame_shift}, Strand: {alternative_strand})...")
        
        # 1. Get Observed Numerators & Ti/Tv stats
        observed, consensus_str, n_ti, n_tv, mutation_details = self.get_observed_counts(
            canonical_seqs, alternative_seqs, arf_start_index, frame_shift,
            count_unique_sites=True, # Correct for population data
            alternative_strand=alternative_strand,
            verbose = verbose, ref_seq=ref_seq
        )
        
        # 2. Calculate Dynamic Ti/Tv Ratio
        if n_tv > 0:
            estimated_titv = n_ti / n_tv
            print(f"Dynamic Ti/Tv estimated from data: {estimated_titv:.2f} (Ti={n_ti}, Tv={n_tv})")
        else:
            estimated_titv = default_titv
            print(f"Warning: No transversions observed. Using default Ti/Tv: {estimated_titv}")

        # 3. Get Potential Denominators (Weighted by Ti/Tv)
        window_end = arf_start_index + len(alternative_seqs[0])
        potential = self.calculate_potential_sites(
            consensus_str, arf_start_index, window_end, frame_shift, 
            titv_ratio=estimated_titv, alternative_strand=alternative_strand
        )
        
        if verbose:
            print(f"Potential sites: {potential}")
            print(f"Observed sites: {observed}")
        
        # 4. Calculate Rates (Observed / Potential)
        rates = {}
        for k in ['SS', 'SN', 'NS', 'NN', 'STOP_F1']:
            if potential[k] > 0:
                rates[k] = observed[k] / potential[k]
            else:
                rates[k] = 0.0
        
        # 5. Beta Stop Calculation
        # Use SS as baseline, but switch to SN if SS is too low
        # Description - Beta Stop is calculated as the ratio of the observed stop mutations to the potential SS or SN mutations.
        # If the potential SS is too low, then the SN is used as the baseline.
        baseline_rate = rates['SS']
        baseline_type = 'SS'
        
        if potential['SS'] < 10 or observed['SS'] < 2:
            if rates['SN'] > 0:
                print("NOTICE: SS counts low. Using SN (Synonymous Canonical) as neutral baseline.")
                baseline_rate = rates['SN']
                baseline_type = 'SN'
        
        
        beta_stop = (rates['STOP_F1'] / baseline_rate) if baseline_rate > 0 else 0.0

        # Calculate individual mutation frequencies (percentage of population with each mutation)
        n_sequences = len(canonical_seqs)
        mutation_frequencies = {}
        for m_type, mutations in mutation_details.items():
            mutation_frequencies[m_type] = [
                {
                    'mutation': mut['mutation'],
                    'position': mut['position'],
                    'count': mut['count'],
                    'frequency_pct': (mut['count'] / n_sequences) * 100
                }
                for mut in mutations
            ]

        return {
            'beta_stop': beta_stop,
            'titv_used': estimated_titv,
            'baseline_used': baseline_type,
            'raw_rates': rates,
            'observed': dict(observed),
            'potential': dict(potential),
            'mutation_frequencies': mutation_frequencies,
            'n_sequences': n_sequences
        }