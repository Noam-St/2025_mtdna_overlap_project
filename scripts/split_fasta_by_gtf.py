#!/usr/bin/env python3
"""
Split a FASTA file into separate files based on GTF annotations.
Each gene annotation produces one FASTA file containing that region.
"""

import sys
import argparse
import re
from pathlib import Path
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


def parse_arguments():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(
        description='Split FASTA file into separate files based on GTF gene annotations.',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  %(prog)s -g annotations.gtf -f reference.fasta -o genes_output/
  %(prog)s -g annotations.gtf -f reference.fasta -o genes_output/ -t CDS
  %(prog)s -g annotations.gtf -f reference.fasta -o genes_output/ -v --prefix mt_

Notes:
  - GTF coordinates are 1-based, inclusive on both ends
  - Handles strand information (reverse complements minus strand features)
  - Output files are named: {prefix}{gene_name}.fasta
  - If multiple features have the same gene name, only the first is used
        """
    )
    
    parser.add_argument(
        '-g', '--gtf',
        required=True,
        metavar='FILE',
        help='Input GTF/GFF annotation file'
    )
    
    parser.add_argument(
        '-f', '--fasta',
        required=True,
        metavar='FILE',
        help='Input FASTA file (reference sequence)'
    )
    
    parser.add_argument(
        '-o', '--output-dir',
        required=True,
        metavar='DIR',
        help='Output directory for gene FASTA files'
    )
    
    parser.add_argument(
        '-t', '--feature-type',
        default='gene',
        metavar='TYPE',
        help='Feature type to extract (e.g., gene, CDS, exon). Default: gene'
    )
    
    parser.add_argument(
        '--prefix',
        default='',
        metavar='STR',
        help='Prefix for output filenames. Default: none'
    )
    
    parser.add_argument(
        '--no-reverse-complement',
        action='store_true',
        help='Do not reverse complement minus strand features'
    )
    
    parser.add_argument(
        '-v', '--verbose',
        action='store_true',
        help='Print verbose output'
    )
    
    return parser.parse_args()


def parse_gtf_attributes(attr_string):
    """
    Parse GTF attributes field to extract gene_id and gene_name.
    Handles both GTF and GFF3 formats.
    """
    attributes = {}
    
    # Try GTF format first (semicolon separated, quoted values)
    # Example: gene_id "COX1"; gene_name "COX1";
    gtf_pattern = r'(\w+)\s+"([^"]+)"'
    for match in re.finditer(gtf_pattern, attr_string):
        key, value = match.groups()
        attributes[key] = value
    
    # If no matches, try GFF3 format (semicolon separated, key=value)
    # Example: ID=gene:COX1;Name=COX1
    if not attributes:
        gff_pattern = r'(\w+)=([^;]+)'
        for match in re.finditer(gff_pattern, attr_string):
            key, value = match.groups()
            attributes[key] = value
    
    return attributes


def parse_gtf(gtf_file, feature_type='gene', verbose=False):
    """
    Parse GTF file and extract gene annotations.
    Returns list of tuples: (seqname, start, end, strand, gene_name)
    """
    annotations = []
    seen_genes = set()
    skipped = 0
    
    if verbose:
        print(f"Parsing GTF file: {gtf_file}", file=sys.stderr)
        print(f"Looking for feature type: {feature_type}", file=sys.stderr)
    
    try:
        with open(gtf_file, 'r') as f:
            for line_num, line in enumerate(f, 1):
                # Skip comments and empty lines
                if line.startswith('#') or not line.strip():
                    continue
                
                fields = line.strip().split('\t')
                
                if len(fields) < 9:
                    if verbose:
                        print(f"Warning: Line {line_num} has fewer than 9 fields. Skipping.", 
                              file=sys.stderr)
                    skipped += 1
                    continue
                
                seqname = fields[0]
                feature = fields[2]
                start = int(fields[3])  # GTF is 1-based
                end = int(fields[4])    # GTF is 1-based, inclusive
                strand = fields[6]
                attributes_str = fields[8]
                
                # Filter by feature type
                if feature != feature_type:
                    continue
                
                # Parse attributes
                attributes = parse_gtf_attributes(attributes_str)
                
                # Try to get gene name (try multiple possible fields)
                gene_name = None
                for key in ['gene_name', 'gene_id', 'Name', 'ID', 'locus_tag']:
                    if key in attributes:
                        gene_name = attributes[key]
                        # Clean up gene name (remove prefixes like "gene:")
                        gene_name = gene_name.split(':')[-1]
                        break
                
                if not gene_name:
                    if verbose:
                        print(f"Warning: Line {line_num} has no gene name. Skipping.", 
                              file=sys.stderr)
                    skipped += 1
                    continue
                
                # Skip duplicates (keep first occurrence)
                if gene_name in seen_genes:
                    if verbose:
                        print(f"Warning: Duplicate gene '{gene_name}' found. Using first occurrence.", 
                              file=sys.stderr)
                    skipped += 1
                    continue
                
                seen_genes.add(gene_name)
                annotations.append((seqname, start, end, strand, gene_name))
                
                if verbose and len(annotations) % 10 == 0:
                    print(f"  Parsed {len(annotations)} annotations...", file=sys.stderr)
    
    except FileNotFoundError:
        print(f"Error: GTF file '{gtf_file}' not found.", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(f"Error parsing GTF file: {e}", file=sys.stderr)
        sys.exit(1)
    
    if verbose:
        print(f"Found {len(annotations)} unique {feature_type} features", file=sys.stderr)
        if skipped > 0:
            print(f"Skipped {skipped} entries", file=sys.stderr)
    
    return annotations


def extract_and_write_genes(fasta_file, annotations, output_dir, prefix='',
                            reverse_complement=True, verbose=False):
    """
    Extract gene sequences from FASTA and write to separate files.
    """
    # Create output directory
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)
    
    if verbose:
        print(f"\nLoading reference sequences from {fasta_file}...", file=sys.stderr)
    
    # Load reference sequences
    try:
        seq_dict = SeqIO.to_dict(SeqIO.parse(fasta_file, "fasta"))
    except FileNotFoundError:
        print(f"Error: FASTA file '{fasta_file}' not found.", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(f"Error reading FASTA file: {e}", file=sys.stderr)
        sys.exit(1)
    
    if verbose:
        print(f"Loaded {len(seq_dict)} sequence(s)", file=sys.stderr)
        for seq_name in seq_dict:
            print(f"  {seq_name}: {len(seq_dict[seq_name].seq)} bp", file=sys.stderr)
    
    # Process each annotation
    if verbose:
        print(f"\nExtracting and writing gene sequences...", file=sys.stderr)
    
    extracted = 0
    skipped = 0
    
    for seqname, start, end, strand, gene_name in annotations:
        # Check if sequence exists in FASTA
        if seqname not in seq_dict:
            print(f"Warning: Sequence '{seqname}' not found in FASTA file. Skipping {gene_name}.", 
                  file=sys.stderr)
            skipped += 1
            continue
        
        # Get the full sequence
        full_seq = seq_dict[seqname].seq
        
        # Extract subsequence (GTF is 1-based, Python is 0-based)
        # GTF end is inclusive, so we use start-1:end
        gene_seq = full_seq[start-1:end]
        
        # Handle strand
        if strand == '-' and reverse_complement:
            gene_seq = gene_seq.reverse_complement()
            strand_info = "minus_strand_revcomp"
        elif strand == '-':
            strand_info = "minus_strand"
        else:
            strand_info = "plus_strand"
        
        # Create SeqRecord
        gene_record = SeqRecord(
            gene_seq,
            id=gene_name,
            description=f"{seqname}:{start}-{end}({strand}) length={len(gene_seq)}bp {strand_info}"
        )
        
        # Write to file
        output_file = output_path / f"{prefix}{gene_name}.fasta"
        try:
            SeqIO.write(gene_record, output_file, "fasta")
            extracted += 1
            
            if verbose:
                print(f"  Wrote {gene_name}: {len(gene_seq)} bp -> {output_file}", 
                      file=sys.stderr)
        except Exception as e:
            print(f"Error writing {gene_name}: {e}", file=sys.stderr)
            skipped += 1
    
    # Summary
    print(f"\n{'='*60}", file=sys.stderr)
    print(f"SUMMARY", file=sys.stderr)
    print(f"{'='*60}", file=sys.stderr)
    print(f"Total annotations:       {len(annotations):>10}", file=sys.stderr)
    print(f"Successfully extracted:  {extracted:>10}", file=sys.stderr)
    print(f"Skipped:                 {skipped:>10}", file=sys.stderr)
    print(f"Output directory:        {output_dir}", file=sys.stderr)
    print(f"{'='*60}", file=sys.stderr)


def main():
    """Main entry point."""
    args = parse_arguments()
    
    # Parse GTF
    annotations = parse_gtf(
        args.gtf,
        feature_type=args.feature_type,
        verbose=args.verbose
    )
    
    if not annotations:
        print(f"Error: No annotations of type '{args.feature_type}' found in GTF file.", 
              file=sys.stderr)
        sys.exit(1)
    
    # Extract and write genes
    extract_and_write_genes(
        args.fasta,
        annotations,
        args.output_dir,
        prefix=args.prefix,
        reverse_complement=not args.no_reverse_complement,
        verbose=args.verbose
    )


if __name__ == '__main__':
    main()