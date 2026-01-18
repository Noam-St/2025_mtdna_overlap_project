#!/usr/bin/env python3
"""
FASTA File Splitter

This script splits a FASTA file into multiple chunks, each containing
a specified number of sequences.

Author: Assistant
"""

import argparse
import os
import sys
from pathlib import Path


def parse_fasta(fasta_file):
    """
    Parse a FASTA file and yield (header, sequence) tuples.
    
    Args:
        fasta_file (str): Path to the FASTA file
        
    Yields:
        tuple: (header, sequence) pairs
    """
    with open(fasta_file, 'r') as f:
        header = None
        sequence = []
        
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                # If we have a previous sequence, yield it
                if header is not None:
                    yield (header, ''.join(sequence))
                # Start new sequence
                header = line[1:]  # Remove '>' character
                sequence = []
            else:
                # Add to current sequence
                sequence.append(line)
        
        # Yield the last sequence
        if header is not None:
            yield (header, ''.join(sequence))


def count_sequences(fasta_file):
    """
    Count the total number of sequences in a FASTA file.
    
    Args:
        fasta_file (str): Path to the FASTA file
        
    Returns:
        int: Number of sequences
    """
    count = 0
    with open(fasta_file, 'r') as f:
        for line in f:
            if line.startswith('>'):
                count += 1
    return count


def write_fasta_chunk(sequences, output_file):
    """
    Write sequences to a FASTA file.
    
    Args:
        sequences (list): List of (header, sequence) tuples
        output_file (str): Output file path
    """
    with open(output_file, 'w') as f:
        for header, sequence in sequences:
            f.write(f">{header}\n")
            # Write sequence in lines of 80 characters (standard FASTA format)
            for i in range(0, len(sequence), 80):
                f.write(sequence[i:i+80] + '\n')


def split_fasta(input_file, sequences_per_chunk, output_prefix, output_dir):
    """
    Split a FASTA file into chunks.
    
    Args:
        input_file (str): Input FASTA file path
        sequences_per_chunk (int): Number of sequences per chunk
        output_prefix (str): Prefix for output files
        output_dir (str): Output directory
    """
    # Create output directory if it doesn't exist
    os.makedirs(output_dir, exist_ok=True)
    
    # Count total sequences
    total_sequences = count_sequences(input_file)
    total_chunks = (total_sequences + sequences_per_chunk - 1) // sequences_per_chunk
    
    print(f"Input file: {input_file}")
    print(f"Total sequences: {total_sequences}")
    print(f"Sequences per chunk: {sequences_per_chunk}")
    print(f"Total chunks to create: {total_chunks}")
    print(f"Output directory: {output_dir}")
    print("-" * 50)
    
    # Parse and split the FASTA file
    chunk_num = 1
    current_chunk = []
    sequences_processed = 0
    
    for header, sequence in parse_fasta(input_file):
        current_chunk.append((header, sequence))
        sequences_processed += 1
        
        # If chunk is full or we've processed all sequences
        if len(current_chunk) == sequences_per_chunk or sequences_processed == total_sequences:
            # Create output filename
            output_file = os.path.join(output_dir, f"{output_prefix}_chunk_{chunk_num:03d}.fasta")
            
            # Write chunk to file
            write_fasta_chunk(current_chunk, output_file)
            
            print(f"Created chunk {chunk_num}: {output_file} ({len(current_chunk)} sequences)")
            
            # Reset for next chunk
            current_chunk = []
            chunk_num += 1
    
    print("-" * 50)
    print(f"Successfully split {total_sequences} sequences into {chunk_num - 1} chunks")


def main():
    """Main function to parse arguments and execute the splitting."""
    parser = argparse.ArgumentParser(
        description="Split a FASTA file into multiple chunks with specified number of sequences per chunk",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  %(prog)s -i sequences.fasta -n 100
  %(prog)s -i input.fasta -n 50 -o output_dir -p split_seqs
  %(prog)s --input large_file.fasta --sequences-per-chunk 500 --output-dir chunks/ --prefix genome
        """
    )
    
    parser.add_argument('-i', '--input', 
                       required=True,
                       help='Input FASTA file path')
    
    parser.add_argument('-n', '--sequences-per-chunk',
                       type=int,
                       required=True,
                       help='Number of sequences per chunk')
    
    parser.add_argument('-o', '--output-dir',
                       default='fasta_chunks',
                       help='Output directory for chunk files (default: fasta_chunks)')
    
    parser.add_argument('-p', '--prefix',
                       default='chunk',
                       help='Prefix for output files (default: chunk)')
    
    parser.add_argument('--version',
                       action='version',
                       version='%(prog)s 1.0')
    
    args = parser.parse_args()
    
    # Validate input file
    if not os.path.isfile(args.input):
        print(f"Error: Input file '{args.input}' does not exist or is not a file", file=sys.stderr)
        sys.exit(1)
    
    # Validate sequences per chunk
    if args.sequences_per_chunk <= 0:
        print("Error: Number of sequences per chunk must be positive", file=sys.stderr)
        sys.exit(1)
    
    # Check if input file is empty
    if os.path.getsize(args.input) == 0:
        print("Error: Input file is empty", file=sys.stderr)
        sys.exit(1)
    
    try:
        split_fasta(args.input, args.sequences_per_chunk, args.prefix, args.output_dir)
    except Exception as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
