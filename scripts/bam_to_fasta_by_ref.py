#!/usr/bin/env python3
"""
Convert BAM alignments to aligned FASTA format preserving reference coordinates.
Each read becomes a sequence with the same length as the reference.
Sequences with insertions are skipped by default to maintain coordinate integrity.
"""

import sys
import argparse
import pysam
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


def parse_arguments():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(
        description='Convert BAM alignments to FASTA format with reference coordinate preservation.',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  %(prog)s -b aligned.bam -r reference.fasta -o output.fasta
  %(prog)s -b aligned.bam -r reference.fasta -o output.fasta -v --include-reference
  %(prog)s -b COX1_aligned.bam -r COX1.fasta -o COX1_sequences.fasta

Notes:
  - All output sequences have the same length as the reference
  - Position i in output corresponds to position i in reference
  - Deletions are represented as gaps ('-')
  - Sequences with insertions are SKIPPED by default
  - Perfect for maintaining gene coordinates for dN/dS analysis
        """
    )
    
    parser.add_argument(
        '-b', '--bam',
        required=True,
        metavar='FILE',
        help='Input BAM file'
    )
    
    parser.add_argument(
        '-r', '--reference',
        required=True,
        metavar='FILE',
        help='Reference FASTA file used for alignment'
    )
    
    parser.add_argument(
        '-o', '--output',
        required=True,
        metavar='FILE',
        help='Output aligned FASTA file'
    )
    
    parser.add_argument(
        '--include-reference',
        action='store_true',
        help='Include reference sequence in output'
    )
    
    parser.add_argument(
        '--allow-insertions',
        action='store_true',
        help='Allow sequences with insertions (insertions will be removed, not recommended)'
    )
    
    parser.add_argument(
        '-v', '--verbose',
        action='store_true',
        help='Print verbose output'
    )
    
    return parser.parse_args()


def has_insertions(read):
    """
    Check if read has any insertions relative to reference.
    Returns True if insertions are present.
    """
    if read.cigartuples is None:
        return False
    
    for op, length in read.cigartuples:
        if op == 1:  # I - insertion
            return True
    return False


def create_aligned_sequence(read, ref_length):
    """
    Create aligned sequence maintaining reference coordinates.
    Returns sequence of length ref_length with gaps for deletions.
    Insertions are discarded to maintain coordinate correspondence.
    """
    # Initialize with gaps
    aligned_seq = ['-'] * ref_length
    
    ref_pos = read.reference_start
    query_pos = 0
    query_seq = read.query_sequence
    
    if query_seq is None:
        return None
    
    for op, length in read.cigartuples:
        if op == 0 or op == 7 or op == 8:  # M, =, X (match/mismatch)
            for i in range(length):
                if ref_pos + i < ref_length:
                    aligned_seq[ref_pos + i] = query_seq[query_pos + i]
            ref_pos += length
            query_pos += length
            
        elif op == 1:  # I - insertion (skip, no reference position)
            query_pos += length
            
        elif op == 2:  # D - deletion (gap in query, advance reference)
            ref_pos += length
            
        elif op == 4:  # S - soft clipping (skip)
            query_pos += length
            
        elif op == 5:  # H - hard clipping (already removed from sequence)
            pass
            
        elif op == 3:  # N - skipped region (e.g., spliced alignment)
            ref_pos += length
    
    return ''.join(aligned_seq)


def bam_to_aligned_fasta(bam_file, ref_fasta, output_fasta, 
                         include_reference=False, allow_insertions=False, 
                         verbose=False):
    """
    Convert BAM alignments to FASTA format preserving reference coordinates.
    """
    # Load reference sequences
    if verbose:
        print(f"Loading reference from {ref_fasta}...", file=sys.stderr)
    
    try:
        ref_dict = SeqIO.to_dict(SeqIO.parse(ref_fasta, "fasta"))
    except Exception as e:
        print(f"Error: Could not read reference file: {e}", file=sys.stderr)
        sys.exit(1)
    
    if len(ref_dict) == 0:
        print(f"Error: No sequences found in reference file", file=sys.stderr)
        sys.exit(1)
    
    if verbose:
        print(f"Loaded {len(ref_dict)} reference sequence(s)", file=sys.stderr)
        for ref_name, ref_rec in ref_dict.items():
            print(f"  {ref_name}: {len(ref_rec.seq)} bp", file=sys.stderr)
    
    # Get reference info (assuming single reference per gene)
    if len(ref_dict) > 1:
        print(f"Warning: Multiple references found in FASTA. Using first one.", 
              file=sys.stderr)
    
    ref_name = list(ref_dict.keys())[0]
    ref_record = ref_dict[ref_name]
    ref_length = len(ref_record.seq)
    
    # Open BAM file
    try:
        bam = pysam.AlignmentFile(bam_file, "rb")
    except Exception as e:
        print(f"Error: Could not open BAM file: {e}", file=sys.stderr)
        sys.exit(1)
    
    output_records = []
    
    # Add reference if requested
    if include_reference:
        output_records.append(SeqRecord(
            ref_record.seq,
            id=ref_name,
            description="reference"
        ))
        if verbose:
            print(f"Including reference sequence in output", file=sys.stderr)
    
    if verbose:
        print(f"\nProcessing alignments...", file=sys.stderr)
        if not allow_insertions:
            print(f"  Sequences with insertions will be SKIPPED", file=sys.stderr)
        else:
            print(f"  Sequences with insertions will be included (insertions removed)", 
                  file=sys.stderr)
    
    # Tracking statistics
    processed = 0
    skipped_unmapped = 0
    skipped_insertions = 0
    skipped_no_sequence = 0
    skipped_ref_mismatch = 0
    
    for read in bam.fetch(until_eof=True):
        # Skip unmapped reads
        if read.is_unmapped or read.cigartuples is None:
            skipped_unmapped += 1
            continue
        
        # Check reference name matches
        read_ref_name = read.reference_name
        if read_ref_name not in ref_dict:
            if verbose and skipped_ref_mismatch < 5:
                print(f"Warning: Reference '{read_ref_name}' in BAM not found in FASTA. "
                      f"Skipping read {read.query_name}", file=sys.stderr)
            skipped_ref_mismatch += 1
            continue
        
        # Skip sequences with insertions (unless allowed)
        if not allow_insertions and has_insertions(read):
            if verbose and skipped_insertions < 5:  # Only show first 5
                print(f"  Skipping {read.query_name}: contains insertions", file=sys.stderr)
            skipped_insertions += 1
            continue
        
        # Create aligned sequence
        aligned_seq = create_aligned_sequence(read, ref_length)
        
        if aligned_seq is None:
            if verbose:
                print(f"Warning: No sequence for read {read.query_name}. Skipping.", 
                      file=sys.stderr)
            skipped_no_sequence += 1
            continue
        
        # Create record
        aligned_record = SeqRecord(
            Seq(aligned_seq),
            id=read.query_name,
            description=f"aligned_to_{read_ref_name} span={read.reference_start+1}-{read.reference_end}"
        )
        output_records.append(aligned_record)
        processed += 1
        
        if verbose and processed % 100 == 0:
            print(f"  Processed {processed} reads...", file=sys.stderr)
    
    bam.close()
    
    # Write output
    if verbose:
        print(f"\nWriting {len(output_records)} sequences to {output_fasta}...", 
              file=sys.stderr)
    
    if len(output_records) == 0:
        print(f"Warning: No sequences to write! Check your BAM file and filters.", 
              file=sys.stderr)
    
    try:
        SeqIO.write(output_records, output_fasta, "fasta")
    except Exception as e:
        print(f"Error: Could not write output file: {e}", file=sys.stderr)
        sys.exit(1)
    
    # Summary statistics
    total_reads = (processed + skipped_unmapped + skipped_insertions + 
                   skipped_no_sequence + skipped_ref_mismatch)
    
    print(f"\n{'='*60}", file=sys.stderr)
    print(f"SUMMARY", file=sys.stderr)
    print(f"{'='*60}", file=sys.stderr)
    print(f"Total reads in BAM:           {total_reads:>10}", file=sys.stderr)
    print(f"Successfully processed:       {processed:>10}", file=sys.stderr)
    print(f"", file=sys.stderr)
    
    if total_reads - processed > 0:
        print(f"Skipped (total):              {total_reads - processed:>10}", file=sys.stderr)
        if skipped_unmapped > 0:
            print(f"  - Unmapped/invalid:         {skipped_unmapped:>10}", file=sys.stderr)
        if skipped_insertions > 0:
            print(f"  - Contains insertions:      {skipped_insertions:>10}", file=sys.stderr)
        if skipped_no_sequence > 0:
            print(f"  - No sequence data:         {skipped_no_sequence:>10}", file=sys.stderr)
        if skipped_ref_mismatch > 0:
            print(f"  - Reference mismatch:       {skipped_ref_mismatch:>10}", file=sys.stderr)
    
    print(f"{'='*60}", file=sys.stderr)
    
    if skipped_insertions > 0:
        pct_skipped = (skipped_insertions / total_reads) * 100 if total_reads > 0 else 0
        print(f"\nNote: {skipped_insertions} sequences ({pct_skipped:.1f}%) were skipped due to insertions.", 
              file=sys.stderr)
        print(f"This ensures all output sequences maintain exact reference coordinates.", 
              file=sys.stderr)
    
    # Verify all sequences have same length
    if len(output_records) > 0:
        lengths = set(len(rec.seq) for rec in output_records)
        if len(lengths) == 1:
            print(f"\n✓ All sequences have identical length: {lengths.pop()} bp", 
                  file=sys.stderr)
        else:
            print(f"\n✗ WARNING: Sequences have inconsistent lengths: {lengths}", 
                  file=sys.stderr)


def main():
    """Main entry point."""
    args = parse_arguments()
    
    bam_to_aligned_fasta(
        bam_file=args.bam,
        ref_fasta=args.reference,
        output_fasta=args.output,
        include_reference=args.include_reference,
        allow_insertions=args.allow_insertions,
        verbose=args.verbose
    )


if __name__ == '__main__':
    main()