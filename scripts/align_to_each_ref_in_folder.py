#!/usr/bin/env python3
"""
Align input sequences to multiple gene references using minimap2.
Creates one alignment file per gene reference.
"""

import sys
import argparse
import subprocess
from pathlib import Path
import shutil


def parse_arguments():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(
        description='Align sequences to multiple gene references using minimap2.',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  %(prog)s -r gene_refs/ -i sequences.fasta -o alignments/
  %(prog)s -r gene_refs/ -i sequences.fasta -o alignments/ -x asm20 -t 20
  %(prog)s -r gene_refs/ -i seq1.fasta seq2.fasta -o alignments/ --sam
  %(prog)s -r gene_refs/ -i sequences.fasta -o alignments/ --minimap2-args "--secondary=no --eqx"

Notes:
  - One alignment file is created per gene reference
  - Output files are named: {gene_name}_aligned.bam (or .sam)
  - BAM files are automatically sorted and indexed
  - Requires minimap2 and samtools in PATH
        """
    )
    
    parser.add_argument(
        '-r', '--reference-dir',
        required=True,
        metavar='DIR',
        help='Directory containing gene reference FASTA files'
    )
    
    parser.add_argument(
        '-i', '--input',
        required=True,
        nargs='+',
        metavar='FILE',
        help='Input FASTA/FASTQ file(s) to align (can specify multiple files)'
    )
    
    parser.add_argument(
        '-o', '--output-dir',
        required=True,
        metavar='DIR',
        help='Output directory for alignment files'
    )
    
    parser.add_argument(
        '-x', '--preset',
        default='asm5',
        metavar='STR',
        help='Minimap2 preset (e.g., asm5, asm10, asm20, map-ont, map-pb). Default: asm5'
    )
    
    parser.add_argument(
        '-t', '--threads',
        type=int,
        default=4,
        metavar='INT',
        help='Number of threads for minimap2. Default: 4'
    )
    
    parser.add_argument(
        '--minimap2-args',
        default='--secondary=no --eqx',
        metavar='STR',
        help='Additional arguments to pass to minimap2. Default: "--secondary=no --eqx"'
    )
    
    parser.add_argument(
        '--sam',
        action='store_true',
        help='Output SAM format instead of BAM (not recommended for large files)'
    )
    
    parser.add_argument(
        '--no-index',
        action='store_true',
        help='Do not index BAM files (saves time if index not needed)'
    )
    
    parser.add_argument(
        '-v', '--verbose',
        action='store_true',
        help='Print verbose output'
    )
    
    return parser.parse_args()


def check_dependencies(verbose=False):
    """Check if required tools are available."""
    required = ['minimap2', 'samtools']
    missing = []
    
    for tool in required:
        if shutil.which(tool) is None:
            missing.append(tool)
        elif verbose:
            print(f"Found {tool}: {shutil.which(tool)}", file=sys.stderr)
    
    if missing:
        print(f"Error: Required tool(s) not found in PATH: {', '.join(missing)}", 
              file=sys.stderr)
        print("Please install missing tools and ensure they are in your PATH.", 
              file=sys.stderr)
        sys.exit(1)


def find_reference_files(ref_dir, verbose=False):
    """Find all FASTA files in reference directory."""
    ref_path = Path(ref_dir)
    
    if not ref_path.exists():
        print(f"Error: Reference directory '{ref_dir}' does not exist.", file=sys.stderr)
        sys.exit(1)
    
    if not ref_path.is_dir():
        print(f"Error: '{ref_dir}' is not a directory.", file=sys.stderr)
        sys.exit(1)
    
    # Find FASTA files (common extensions)
    fasta_files = []
    for ext in ['*.fasta', '*.fa', '*.fna', '*.fasta.gz', '*.fa.gz']:
        fasta_files.extend(ref_path.glob(ext))
    
    if not fasta_files:
        print(f"Error: No FASTA files found in '{ref_dir}'", file=sys.stderr)
        sys.exit(1)
    
    if verbose:
        print(f"Found {len(fasta_files)} reference file(s):", file=sys.stderr)
        for f in sorted(fasta_files):
            print(f"  {f.name}", file=sys.stderr)
    
    return sorted(fasta_files)


def check_input_files(input_files, verbose=False):
    """Verify input files exist."""
    for input_file in input_files:
        if not Path(input_file).exists():
            print(f"Error: Input file '{input_file}' does not exist.", file=sys.stderr)
            sys.exit(1)
    
    if verbose:
        print(f"\nInput file(s):", file=sys.stderr)
        for f in input_files:
            print(f"  {f}", file=sys.stderr)


def run_minimap2(ref_file, input_files, output_file, preset, threads, 
                 extra_args, verbose=False):
    """
    Run minimap2 alignment.
    Returns True on success, False on failure.
    """
    # Build minimap2 command
    cmd = [
        'minimap2',
        '-ax', preset,
        '-t', str(threads),
        '-a'  # Output SAM format
    ]
    
    # Add extra arguments if provided
    if extra_args:
        cmd.extend(extra_args.split())
    
    # Add reference and input files
    cmd.append(str(ref_file))
    cmd.extend(input_files)
    
    if verbose:
        print(f"Running: {' '.join(cmd)}", file=sys.stderr)
    
    try:
        with open(output_file, 'w') as outf:
            result = subprocess.run(
                cmd,
                stdout=outf,
                stderr=subprocess.PIPE,
                text=True,
                check=True
            )
        
        if verbose and result.stderr:
            print(result.stderr, file=sys.stderr)
        
        return True
    
    except subprocess.CalledProcessError as e:
        print(f"Error: minimap2 failed for {ref_file.name}", file=sys.stderr)
        print(e.stderr, file=sys.stderr)
        return False
    except Exception as e:
        print(f"Error running minimap2: {e}", file=sys.stderr)
        return False


def sam_to_bam(sam_file, bam_file, threads, verbose=False):
    """Convert SAM to sorted BAM using samtools."""
    cmd = [
        'samtools', 'sort',
        '-@', str(threads),
        '-o', str(bam_file),
        str(sam_file)
    ]
    
    if verbose:
        print(f"Converting to BAM: {' '.join(cmd)}", file=sys.stderr)
    
    try:
        result = subprocess.run(
            cmd,
            stderr=subprocess.PIPE,
            text=True,
            check=True
        )
        
        if verbose and result.stderr:
            print(result.stderr, file=sys.stderr)
        
        # Remove SAM file after successful conversion
        Path(sam_file).unlink()
        return True
    
    except subprocess.CalledProcessError as e:
        print(f"Error: samtools sort failed", file=sys.stderr)
        print(e.stderr, file=sys.stderr)
        return False
    except Exception as e:
        print(f"Error converting to BAM: {e}", file=sys.stderr)
        return False


def index_bam(bam_file, verbose=False):
    """Index BAM file using samtools."""
    cmd = ['samtools', 'index', str(bam_file)]
    
    if verbose:
        print(f"Indexing: {' '.join(cmd)}", file=sys.stderr)
    
    try:
        result = subprocess.run(
            cmd,
            stderr=subprocess.PIPE,
            text=True,
            check=True
        )
        
        if verbose and result.stderr:
            print(result.stderr, file=sys.stderr)
        
        return True
    
    except subprocess.CalledProcessError as e:
        print(f"Error: samtools index failed", file=sys.stderr)
        print(e.stderr, file=sys.stderr)
        return False
    except Exception as e:
        print(f"Error indexing BAM: {e}", file=sys.stderr)
        return False


def align_to_references(ref_files, input_files, output_dir, preset, threads,
                       extra_args, output_sam=False, no_index=False, verbose=False):
    """Align input sequences to each reference file."""
    # Create output directory
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)
    
    successful = 0
    failed = 0
    
    print(f"\n{'='*60}", file=sys.stderr)
    print(f"Starting alignments with minimap2", file=sys.stderr)
    print(f"{'='*60}", file=sys.stderr)
    print(f"Preset: {preset}", file=sys.stderr)
    print(f"Threads: {threads}", file=sys.stderr)
    print(f"Extra args: {extra_args}", file=sys.stderr)
    print(f"Output format: {'SAM' if output_sam else 'BAM'}", file=sys.stderr)
    print(f"{'='*60}\n", file=sys.stderr)
    
    for i, ref_file in enumerate(ref_files, 1):
        gene_name = ref_file.stem  # Filename without extension
        
        print(f"[{i}/{len(ref_files)}] Processing {gene_name}...", file=sys.stderr)
        
        if output_sam:
            # Output directly to SAM
            output_file = output_path / f"{gene_name}_aligned.sam"
            
            if verbose:
                print(f"  Output: {output_file}", file=sys.stderr)
            
            success = run_minimap2(
                ref_file, input_files, output_file,
                preset, threads, extra_args, verbose
            )
            
            if success:
                successful += 1
                print(f"  ✓ Complete: {output_file}", file=sys.stderr)
            else:
                failed += 1
                print(f"  ✗ Failed", file=sys.stderr)
        
        else:
            # Output to SAM, then convert to BAM
            temp_sam = output_path / f"{gene_name}_aligned.temp.sam"
            output_bam = output_path / f"{gene_name}_aligned.bam"
            
            if verbose:
                print(f"  Output: {output_bam}", file=sys.stderr)
            
            # Run minimap2
            success = run_minimap2(
                ref_file, input_files, temp_sam,
                preset, threads, extra_args, verbose
            )
            
            if not success:
                failed += 1
                print(f"  ✗ Alignment failed", file=sys.stderr)
                continue
            
            # Convert to BAM
            success = sam_to_bam(temp_sam, output_bam, threads, verbose)
            
            if not success:
                failed += 1
                print(f"  ✗ BAM conversion failed", file=sys.stderr)
                continue
            
            # Index BAM (unless disabled)
            if not no_index:
                success = index_bam(output_bam, verbose)
                if not success:
                    print(f"  ⚠ Warning: BAM indexing failed (alignment file still usable)", 
                          file=sys.stderr)
            
            successful += 1
            print(f"  ✓ Complete: {output_bam}", file=sys.stderr)
        
        print()  # Empty line between genes
    
    # Summary
    print(f"{'='*60}", file=sys.stderr)
    print(f"ALIGNMENT SUMMARY", file=sys.stderr)
    print(f"{'='*60}", file=sys.stderr)
    print(f"Total references:     {len(ref_files):>10}", file=sys.stderr)
    print(f"Successful:           {successful:>10}", file=sys.stderr)
    print(f"Failed:               {failed:>10}", file=sys.stderr)
    print(f"Output directory:     {output_dir}", file=sys.stderr)
    print(f"{'='*60}", file=sys.stderr)
    
    if failed > 0:
        sys.exit(1)


def main():
    """Main entry point."""
    args = parse_arguments()
    
    # Check dependencies
    check_dependencies(args.verbose)
    
    # Find reference files
    ref_files = find_reference_files(args.reference_dir, args.verbose)
    
    # Check input files
    check_input_files(args.input, args.verbose)
    
    # Run alignments
    align_to_references(
        ref_files=ref_files,
        input_files=args.input,
        output_dir=args.output_dir,
        preset=args.preset,
        threads=args.threads,
        extra_args=args.minimap2_args,
        output_sam=args.sam,
        no_index=args.no_index,
        verbose=args.verbose
    )


if __name__ == '__main__':
    main()