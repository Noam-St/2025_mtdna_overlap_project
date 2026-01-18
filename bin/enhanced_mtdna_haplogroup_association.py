#!/usr/bin/env python3
"""
Enhanced Mitochondrial Haplogroup-Mutation Association Analysis

This script performs comprehensive association analysis between high frequency 
mutations in mitochondrial genes and haplogroups, with statistical testing,
absolute position annotation, and haplogroup marker detection.

Author: Assistant for PhD Genetics Student
Focus: Mito-nuclear co-regulation studies

ENHANCEMENTS:
- Absolute position annotation using GenBank coordinates
- Haplogroup marker detection and validation
- Sequence alignment for gene finding
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from collections import Counter, defaultdict
from scipy import stats
from scipy.stats import chi2_contingency, fisher_exact
import warnings
from typing import Dict, List, Tuple, Optional, Any, Union
import logging
from statsmodels.stats.multitest import multipletests
import re

# NEW IMPORTS FOR ENHANCED FUNCTIONALITY
from Bio import SeqIO, Entrez
from Bio.Align import PairwiseAligner
import urllib.request
import tempfile
import os

# Set up logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

# Default path for haplogroup markers file
DEFAULT_HAPLOGROUP_MARKERS_PATH = os.path.join(
    os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
    'data', 'haplogroups', 'haplogroup_markers.xlsx'
)

def load_haplogroup_markers(filepath: Optional[str] = None) -> Dict[str, List[int]]:
    """
    Load haplogroup markers from an Excel file.

    Parameters:
    -----------
    filepath : str, optional
        Path to the haplogroup markers Excel file. If not provided,
        uses the default path: data/haplogroups/haplogroup_markers.xlsx

    Returns:
    --------
    Dict[str, List[int]]: Dictionary mapping haplogroup names to lists of marker positions
    """
    if filepath is None:
        filepath = DEFAULT_HAPLOGROUP_MARKERS_PATH

    if not os.path.exists(filepath):
        logger.warning(f"Haplogroup markers file not found: {filepath}. Using empty dictionary.")
        return {}

    try:
        df = pd.read_excel(filepath)

        # Validate required columns
        required_cols = ['Top Level Haplogroup', 'HG Markers']
        missing_cols = [col for col in required_cols if col not in df.columns]
        if missing_cols:
            logger.warning(f"Missing columns in haplogroup markers file: {missing_cols}. Using empty dictionary.")
            return {}

        # Filter rows with valid haplogroup data
        valid_df = df[df['Top Level Haplogroup'].notna()].copy()

        markers_dict = {}
        for _, row in valid_df.iterrows():
            haplogroup = str(row['Top Level Haplogroup']).strip()
            markers_str = row['HG Markers']

            if pd.isna(markers_str) or not isinstance(markers_str, str):
                markers_dict[haplogroup] = []
                continue

            # Parse markers: extract position numbers from strings like "1048T", "3516a", "249d", "(310C)"
            positions = []
            for marker in markers_str.split(','):
                marker = marker.strip().strip('()')
                # Extract leading digits as the position
                match = re.match(r'^(\d+)', marker)
                if match:
                    positions.append(int(match.group(1)))

            markers_dict[haplogroup] = positions

        logger.info(f"Loaded haplogroup markers for {len(markers_dict)} haplogroups from {filepath}")
        return markers_dict

    except Exception as e:
        logger.warning(f"Error loading haplogroup markers from {filepath}: {e}. Using empty dictionary.")
        return {}

# Load haplogroup markers from file (can be reloaded with custom path using load_haplogroup_markers())
KNOWN_HAPLOGROUP_MARKERS = load_haplogroup_markers()

def filter_sequence_to_acgt(sequence: str) -> str:
    """
    Filter a sequence to contain only A, C, G, T bases.

    Parameters:
    -----------
    sequence : str
        Input DNA sequence

    Returns:
    --------
    str: Filtered sequence containing only A, C, G, T (uppercase)
    """
    if pd.isna(sequence) or not isinstance(sequence, str):
        return ''
    # Convert to uppercase and keep only A, C, G, T
    return re.sub(r'[^ACGT]', '', sequence.upper())


def filter_sequences_to_same_length(sequences: pd.Series,
                                    method: str = 'mode',
                                    min_length: int = 20) -> Tuple[pd.Series, int]:
    """
    Filter sequences to ensure they are all the same length.

    Parameters:
    -----------
    sequences : pd.Series
        Series of DNA sequences
    method : str
        Method to determine target length:
        - 'mode': Use the most common sequence length
        - 'max': Keep only sequences of the maximum length
        - 'min': Keep only sequences of the minimum length (above min_length)
    min_length : int
        Minimum acceptable sequence length

    Returns:
    --------
    Tuple of (filtered_sequences, target_length)
    """
    # First filter out empty/NaN sequences
    valid_mask = sequences.notna() & (sequences.str.len() >= min_length)
    valid_sequences = sequences[valid_mask]

    if len(valid_sequences) == 0:
        logger.warning("No valid sequences found after length filtering")
        return pd.Series(dtype=str), 0

    # Calculate sequence lengths
    lengths = valid_sequences.str.len()

    if method == 'mode':
        # Use most common length
        length_counts = lengths.value_counts()
        target_length = length_counts.index[0]
    elif method == 'max':
        target_length = lengths.max()
    elif method == 'min':
        target_length = lengths.min()
    else:
        raise ValueError(f"Unknown method: {method}")

    # Filter to target length
    length_mask = lengths == target_length
    filtered_sequences = valid_sequences[length_mask]

    n_filtered = len(valid_sequences) - len(filtered_sequences)
    if n_filtered > 0:
        logger.info(f"Filtered out {n_filtered} sequences with length != {target_length}")

    logger.info(f"Keeping {len(filtered_sequences)} sequences of length {target_length}")

    return filtered_sequences, target_length


class MtDNAHaplogroupAssociation:
    """
    Enhanced analysis of mtDNA mutations and haplogroup associations with
    absolute position annotation and haplogroup marker detection.

    Supports two modes:
    - 'full_mtdna': Analysis using columns ending in '_seq' for gene sequences
                    and optional 'sequence' column for full mtDNA
    - 'gene': Analysis using specific gene columns (like 'RNR1', 'CYTB') with
              gene-relative positions that are converted to absolute mtDNA positions
    """
    
    def __init__(self, hs_pop_seq: pd.DataFrame, haplogroup_file: str,
                 frequency_threshold: float = 0.01, id_column: str = 'SampleID',
                 haplogroup_column: str = 'Haplogroup',
                 reference_id: str = None, marker_threshold: float = 0.8,
                 accession_id: str = "NC_012920.1",
                 mode: str = 'full_mtdna',
                 sequence_column: Union[str, Dict[str, int]] = None,
                 filter_sequences: bool = True,
                 length_filter_method: str = 'mode',
                 min_sequence_length: int = 100):
        """
        Initialize the enhanced association analysis.

        Parameters:
        -----------
        hs_pop_seq : pd.DataFrame
            DataFrame with sequences. Format depends on mode:
            - full_mtdna mode: columns ending in '_seq' for gene sequences
            - gene mode: gene-specific columns (e.g., 'RNR1', 'CYTB')
        haplogroup_file : str
            Path to CSV file containing haplogroup assignments
        frequency_threshold : float
            Minimum frequency threshold for mutations (default: 0.01 = 1%)
        id_column : str
            Column name for sample IDs in haplogroup file
        haplogroup_column : str
            Column name for haplogroups in haplogroup file
        reference_id : str, optional
            Reference sequence ID or GenBank accession
        marker_threshold : float
            Threshold for haplogroup marker detection (default: 0.8 = 80%)
        accession_id : str
            Default GenBank accession ID for reference genome
        mode : str
            Analysis mode:
            - 'full_mtdna': Look for columns ending in '_seq' (default)
            - 'gene': Use specific gene columns with mtDNA start positions
        sequence_column : str or Dict[str, int], optional
            For gene mode:
            - str: Single gene column name (e.g., 'RNR1')
            - Dict: Mapping of gene column names to mtDNA start positions
                   e.g., {'RNR1': 1671, 'CYTB': 14747}
            For full_mtdna mode: ignored (uses '_seq' columns)
        filter_sequences : bool
            Whether to filter sequences to valid ACGT bases and same length
        length_filter_method : str
            Method for length filtering ('mode', 'max', 'min')
        min_sequence_length : int
            Minimum acceptable sequence length
        """
        self.hs_pop_seq = hs_pop_seq.copy()
        self.frequency_threshold = frequency_threshold
        self.marker_threshold = marker_threshold
        self.id_column = id_column
        self.haplogroup_column = haplogroup_column
        self.reference_id = reference_id
        self.accession_id = accession_id
        self.mode = mode
        self.sequence_column = sequence_column
        self.filter_sequences = filter_sequences
        self.length_filter_method = length_filter_method
        self.min_sequence_length = min_sequence_length

        # Gene start positions for absolute position calculation
        self.gene_start_positions = {}

        # NEW: Gene coordinates storage
        self.gene_coordinates = {}
        self.reference_sequence = None

        # Load and process haplogroup data
        self.haplogroup_df = pd.read_csv(haplogroup_file)

        # Debug information about data structures
        logger.info(f"Sequence dataframe shape: {self.hs_pop_seq.shape}")
        logger.info(f"Sequence dataframe index type: {type(self.hs_pop_seq.index[0]) if len(self.hs_pop_seq) > 0 else 'empty'}")
        logger.info(f"Haplogroup dataframe shape: {self.haplogroup_df.shape}")
        logger.info(f"Haplogroup ID column '{id_column}' type: {self.haplogroup_df[id_column].dtype}")
        logger.info(f"Analysis mode: {mode}")

        # Check if sequence dataframe has SampleID column or uses index
        self._setup_sample_ids()
        self._process_haplogroups()

        # Setup gene columns based on mode
        self._setup_gene_columns()

        # Filter sequences if requested
        if filter_sequences:
            self._filter_all_sequences()

        # NEW: Initialize reference sequence and gene coordinates (only for full_mtdna mode)
        if mode == 'full_mtdna':
            self._initialize_reference_and_coordinates()

        # Storage for results
        self.mutation_results = {}
        self.association_results = {}
        self.summary_stats = {}
        self.haplogroup_markers = {}  # NEW: Store detected markers

        logger.info(f"Initialized enhanced analysis with {len(self.gene_columns)} genes and "
                   f"{len(self.haplogroup_df)} samples with haplogroup information")
        logger.info(f"Gene columns found: {self.gene_names}")

    def _setup_gene_columns(self):
        """Setup gene columns based on mode."""
        if self.mode == 'gene':
            # Gene mode: use specified sequence columns
            if self.sequence_column is None:
                raise ValueError("sequence_column must be provided in gene mode")

            if isinstance(self.sequence_column, dict):
                # Multi-gene mode: dict of {gene_col: start_position}
                self.gene_columns = list(self.sequence_column.keys())
                self.gene_start_positions = self.sequence_column.copy()
                logger.info(f"Gene mode with {len(self.gene_columns)} gene columns")
            else:
                # Single gene mode
                self.gene_columns = [self.sequence_column]
                self.gene_start_positions = {}  # Will need to be set manually or via GenBank
                logger.info(f"Gene mode with single column: {self.sequence_column}")

            self.gene_names = self.gene_columns  # In gene mode, column names are gene names

            # Validate that columns exist
            for col in self.gene_columns:
                if col not in self.hs_pop_seq.columns:
                    # Try column with suffix '_seq'
                    if col + '_seq' in self.hs_pop_seq.columns:
                        self.gene_columns[self.gene_columns.index(col)] = col + '_seq'
                    else:
                        raise ValueError(f"Gene column '{col}' not found in DataFrame. "
                                   f"Available columns: {list(self.hs_pop_seq.columns)}")

        else:
            # Full mtDNA mode: look for columns ending in '_seq'
            self.gene_columns = [col for col in self.hs_pop_seq.columns if col.endswith('_seq')]
            self.gene_names = [col.replace('_seq', '') for col in self.gene_columns]
            logger.info(f"Full mtDNA mode with {len(self.gene_columns)} _seq columns")

    def _filter_all_sequences(self):
        """Filter all gene sequence columns to valid ACGT and same length."""
        logger.info("Filtering sequences to valid ACGT bases and uniform length...")

        for gene_col in self.gene_columns:
            if gene_col not in self.hs_pop_seq.columns:
                continue

            # Get original count
            original_count = self.hs_pop_seq[gene_col].notna().sum()

            # Filter to ACGT only
            self.hs_pop_seq[gene_col] = self.hs_pop_seq[gene_col].apply(filter_sequence_to_acgt)

            # Replace empty strings with NaN
            self.hs_pop_seq.loc[self.hs_pop_seq[gene_col] == '', gene_col] = np.nan

            # Get valid sequences for length filtering
            valid_seqs = self.hs_pop_seq[gene_col].dropna()
            if len(valid_seqs) == 0:
                logger.warning(f"No valid sequences remaining in {gene_col} after ACGT filtering")
                continue

            # Filter to same length
            filtered_seqs, target_length = filter_sequences_to_same_length(
                valid_seqs,
                method=self.length_filter_method,
                min_length=self.min_sequence_length
            )

            # Update dataframe - set non-matching length sequences to NaN
            if target_length > 0:
                length_mask = self.hs_pop_seq[gene_col].str.len() != target_length
                self.hs_pop_seq.loc[length_mask, gene_col] = np.nan

            final_count = self.hs_pop_seq[gene_col].notna().sum()
            logger.info(f"{gene_col}: {original_count} -> {final_count} sequences "
                       f"(target length: {target_length})")
        
    def _initialize_reference_and_coordinates(self):
        """NEW: Initialize reference sequence and get gene coordinates."""
        try:
            # Get reference sequence
            if self.reference_id:
                if self.reference_id in self.hs_pop_seq.index:
                    # Reference ID points to a row in our data
                    self.reference_sequence = self.hs_pop_seq.loc[self.reference_id, 'sequence']
                    logger.info(f"Using reference sequence from sample {self.reference_id}")
                else:
                    # Assume it's a GenBank accession
                    self.accession_id = self.reference_id
                    if not self.reference_id.endswith('.1'):
                        self.accession_id += '.1'
                    logger.info(f"Using GenBank accession {self.accession_id}")
            else:
                # Use most common sequence
                self.reference_sequence = self.hs_pop_seq['sequence'].value_counts().index[0]
                logger.info("Using most common sequence as reference")
            
            # Get gene coordinates from GenBank
            self.gene_coordinates = self._get_gene_coordinates_from_genbank()
            logger.info(f"Retrieved coordinates for {len(self.gene_coordinates)} genes")
            
        except Exception as e:
            logger.warning(f"Could not initialize reference/coordinates: {e}")
            self.gene_coordinates = {}
    
    def _get_gene_coordinates_from_genbank(self):
        """NEW: Fetch gene coordinates from GenBank file."""
        gene_coords = {}
        
        try:
            # Set email for Entrez (required)
            Entrez.email = "your_email@example.com"  # Replace with actual email
            
            # Fetch the GenBank record
            handle = Entrez.efetch(db="nucleotide", id=self.accession_id, rettype="gb", retmode="text")
            
            # Create temporary file to store GenBank data
            with tempfile.NamedTemporaryFile(mode='w', suffix='.gb', delete=False) as temp_file:
                temp_file.write(handle.read())
                temp_file_path = temp_file.name
            handle.close()
            
            # Parse the GenBank file
            with open(temp_file_path, 'r') as file:
                for record in SeqIO.parse(file, "genbank"):
                    if not self.reference_sequence:
                        self.reference_sequence = str(record.seq)
                    
                    for feature in record.features:
                        if feature.type == "gene":
                            gene_name = None
                            if "gene" in feature.qualifiers:
                                gene_name = feature.qualifiers["gene"][0]
                            elif "locus_tag" in feature.qualifiers:
                                gene_name = feature.qualifiers["locus_tag"][0]
                            
                            if gene_name:
                                # Convert to 0-based coordinates
                                start = int(feature.location.start)
                                end = int(feature.location.end)
                                strand = '+' if feature.location.strand == 1 else '-'
                                
                                # Map common gene name variations
                                gene_name = self._standardize_gene_name(gene_name)
                                gene_coords[gene_name] = (start, end, strand)
            
            # Clean up temporary file
            os.unlink(temp_file_path)
            
            logger.info(f"Successfully retrieved coordinates for genes: {list(gene_coords.keys())}")
            
        except Exception as e:
            logger.warning(f"Could not fetch GenBank data: {e}")
            logger.info("Will use sequence alignment as fallback")
        
        return gene_coords
    
    def _standardize_gene_name(self, gene_name):
        """NEW: Standardize gene names to match our data."""
        # Common gene name mappings
        name_mappings = {
            'CYTB': 'cytb',
            'COX1': 'cox1', 'COI': 'cox1',
            'COX2': 'cox2', 'COII': 'cox2',
            'COX3': 'cox3', 'COIII': 'cox3',
            'ND1': 'nad1', 'NADH1': 'nad1',
            'ND2': 'nad2', 'NADH2': 'nad2',
            'ND3': 'nad3', 'NADH3': 'nad3',
            'ND4': 'nad4', 'NADH4': 'nad4',
            'ND4L': 'nad4L', 'NADH4L': 'nad4L',
            'ND5': 'nad5', 'NADH5': 'nad5',
            'ND6': 'nad6', 'NADH6': 'nad6',
            'ATP6': 'atp6',
            'ATP8': 'atp8',
            '12S': 'rnr1', 'RNR1': 'rnr1',
            '16S': 'rnr2', 'RNR2': 'rnr2'
        }
        
        return name_mappings.get(gene_name.upper(), gene_name.lower())
    
    def _find_gene_by_alignment(self, gene_sequence, gene_name):
        """NEW: Find gene coordinates using sequence alignment."""
        if not self.reference_sequence:
            logger.warning("No reference sequence available for alignment")
            return None
        
        try:
            aligner = PairwiseAligner()
            aligner.match_score = 2
            aligner.mismatch_score = -1
            aligner.open_gap_score = -2
            aligner.extend_gap_score = -0.5
            
            alignments = aligner.align(self.reference_sequence, gene_sequence)
            
            if alignments:
                best_alignment = alignments[0]
                start = best_alignment.aligned[0][0][0]  # Start position in reference
                end = best_alignment.aligned[0][0][1]    # End position in reference
                
                logger.info(f"Found {gene_name} by alignment at positions {start}-{end}")
                return (start, end, '+')
            
        except Exception as e:
            logger.warning(f"Alignment failed for {gene_name}: {e}")
        
        return None
    
    def debug_data_structure(self):
        """
        Debug function to examine data structure and ID matching.
        """
        print("="*60)
        print("DATA STRUCTURE DEBUG INFORMATION")
        print("="*60)
        
        # Sequence dataframe info
        print(f"Sequence dataframe shape: {self.hs_pop_seq.shape}")
        print(f"Sequence dataframe columns: {list(self.hs_pop_seq.columns)}")
        print(f"Gene columns found: {self.gene_columns}")
        
        if len(self.hs_pop_seq) > 0:
            print(f"First few SampleIDs: {self.hs_pop_seq['SampleID'].head().tolist()}")
            print(f"SampleID data type: {self.hs_pop_seq['SampleID'].dtype}")
        
        # Haplogroup dataframe info  
        print(f"\nHaplogroup dataframe shape: {self.haplogroup_df.shape}")
        print(f"Haplogroup columns: {list(self.haplogroup_df.columns)}")
        
        if len(self.haplogroup_df) > 0:
            print(f"First few haplogroup IDs: {self.haplogroup_df[self.id_column].head().tolist()}")
            print(f"Haplogroup ID data type: {self.haplogroup_df[self.id_column].dtype}")
        
        # Check ID overlap
        if len(self.hs_pop_seq) > 0 and len(self.haplogroup_df) > 0:
            seq_ids = set(self.hs_pop_seq['SampleID'].astype(str))
            hg_ids = set(self.haplogroup_df[self.id_column].astype(str))
            
            overlap = seq_ids.intersection(hg_ids)
            print(f"\nID matching:")
            print(f"Sequence IDs: {len(seq_ids)}")
            print(f"Haplogroup IDs: {len(hg_ids)}")
            print(f"Overlapping IDs: {len(overlap)}")
            
            if len(overlap) > 0:
                print(f"Sample overlapping IDs: {list(overlap)[:10]}")
            else:
                print("WARNING: No overlapping IDs found!")
                print(f"Sample seq IDs: {list(seq_ids)[:5]}")
                print(f"Sample hg IDs: {list(hg_ids)[:5]}")
        
        # Haplogroup distribution
        if 'MajorHaplogroup' in self.haplogroup_df.columns:
            hg_counts = self.haplogroup_df['MajorHaplogroup'].value_counts()
            print(f"\nMajor haplogroup distribution:")
            print(hg_counts)
        
        # NEW: Gene coordinates info
        print(f"\nGene coordinates retrieved: {len(self.gene_coordinates)}")
        if self.gene_coordinates:
            for gene, coords in list(self.gene_coordinates.items())[:5]:
                print(f"  {gene}: {coords[0]}-{coords[1]} ({coords[2]})")
        
        print("="*60)
        
    def _setup_sample_ids(self):
        """Setup sample ID handling for proper merging."""
        # Check what we have for sample identification
        has_sampleid_col = 'SampleID' in self.hs_pop_seq.columns
        
        if not has_sampleid_col:
            # Check if index looks like sample IDs
            if len(self.hs_pop_seq) > 0:
                index_sample = str(self.hs_pop_seq.index[0])
                # If index looks numeric, create SampleID column from index
                self.hs_pop_seq['SampleID'] = self.hs_pop_seq.index.astype(str)
                logger.info("Created SampleID column from dataframe index")
            else:
                logger.warning("Empty sequence dataframe")
        else:
            # Ensure SampleID is string type for consistent merging
            self.hs_pop_seq['SampleID'] = self.hs_pop_seq['SampleID'].astype(str)
            logger.info("Using existing SampleID column")
        
        # Show sample of sample IDs for debugging
        if len(self.hs_pop_seq) > 0:
            sample_ids = self.hs_pop_seq['SampleID'].head(5).tolist()
            logger.info(f"Sample sequence IDs: {sample_ids}")
        
        # Also show sample of haplogroup IDs
        if len(self.haplogroup_df) > 0:
            hg_sample_ids = self.haplogroup_df[self.id_column].head(5).tolist()
            logger.info(f"Sample haplogroup IDs: {hg_sample_ids}")
            
            # Check for overlap
            seq_ids_set = set(self.hs_pop_seq['SampleID'].astype(str))
            hg_ids_set = set(self.haplogroup_df[self.id_column].astype(str))
            overlap = len(seq_ids_set.intersection(hg_ids_set))
            logger.info(f"ID overlap between datasets: {overlap} samples")
        
    def _process_haplogroups(self):
        """Process haplogroup data to major lineage level."""
        # Create major haplogroup categories based on the summary image
        def get_major_haplogroup(hg):
            if pd.isna(hg):
                return 'Unknown'
            hg = str(hg).strip()
            
            # African lineages (L)
            if hg.startswith('L'):
                if hg.startswith('L0'):
                    return 'L0'
                elif hg.startswith('L1'):
                    return 'L1'
                elif hg.startswith('L2'):
                    return 'L2'
                elif hg.startswith('L3'):
                    return 'L3'
                elif hg.startswith('L4'):
                    return 'L4'
                elif hg.startswith('L5'):
                    return 'L5'
                elif hg.startswith('L6'):
                    return 'L6'
                else:
                    return 'L_other'
            
            # Asian lineages (M)
            elif hg.startswith(('M', 'C', 'D', 'E', 'G', 'Q', 'Z')):
                if hg.startswith('M') and not hg.startswith(('M7', 'M8', 'M9')):
                    return 'M'
                elif hg.startswith('C'):
                    return 'C'
                elif hg.startswith('D'):
                    return 'D'
                elif hg.startswith('E'):
                    return 'E'
                elif hg.startswith('G'):
                    return 'G'
                elif hg.startswith('Q'):
                    return 'Q'
                elif hg.startswith('Z'):
                    return 'Z'
                else:
                    return 'M_other'
            
            # Eurasian lineages (N)
            elif hg.startswith(('N', 'A', 'B', 'F', 'H', 'I', 'J', 'K', 'P', 'R', 'S', 'T', 'U', 'V', 'W', 'X', 'Y', 'HV')):
                if hg.startswith('H') and not hg.startswith('HV'):
                    return 'H'
                elif hg.startswith('U'):
                    return 'U'
                elif hg.startswith('B'):
                    return 'B'
                elif hg.startswith('J'):
                    return 'J'
                elif hg.startswith('T'):
                    return 'T'
                elif hg.startswith('K'):
                    return 'K'
                elif hg.startswith('F'):
                    return 'F'
                elif hg.startswith('A'):
                    return 'A'
                elif hg.startswith('R'):
                    return 'R'
                elif hg.startswith('HV'):
                    return 'HV'
                elif hg.startswith('N') and len(hg) > 1:
                    return 'N'
                elif hg.startswith('I'):
                    return 'I'
                elif hg.startswith('V'):
                    return 'V'
                elif hg.startswith('X'):
                    return 'X'
                elif hg.startswith('W'):
                    return 'W'
                elif hg.startswith('P'):
                    return 'P'
                elif hg.startswith('Y'):
                    return 'Y'
                elif hg.startswith('S'):
                    return 'S'
                elif hg == 'O':
                    return 'O'
                else:
                    return 'N_other'
            else:
                return 'Unknown'
        
        self.haplogroup_df['MajorHaplogroup'] = self.haplogroup_df[self.haplogroup_column].apply(get_major_haplogroup)
        
        # Get haplogroup counts for filtering
        hg_counts = self.haplogroup_df['MajorHaplogroup'].value_counts()
        logger.info(f"Haplogroup distribution:\n{hg_counts}")
        
        # Filter out very rare haplogroups (less than 10 samples)
        min_samples = 10
        rare_hgs = hg_counts[hg_counts < min_samples].index.tolist()
        if rare_hgs:
            logger.info(f"Excluding rare haplogroups with <{min_samples} samples: {rare_hgs}")
            self.haplogroup_df = self.haplogroup_df[~self.haplogroup_df['MajorHaplogroup'].isin(rare_hgs)]
    
    def identify_mutations(self, gene_name: str, reference_id: Optional[str] = None) -> Dict[str, Any]:
        """
        ENHANCED: Identify high frequency mutations with absolute position annotation.

        Parameters:
        -----------
        gene_name : str
            Name of the gene to analyze
        reference_id : str, optional
            Reference sequence ID. If None, uses most common sequence.

        Returns:
        --------
        Dict containing mutation information with absolute positions
        """
        # Determine sequence column based on mode
        if self.mode == 'gene':
            # In gene mode, column name is the gene name directly
            seq_column = gene_name
        else:
            # In full_mtdna mode, column name is gene_name + '_seq'
            seq_column = f"{gene_name}_seq"

        if seq_column not in self.hs_pop_seq.columns:
            raise ValueError(f"Column {seq_column} not found in sequence data")

        logger.info(f"Analyzing mutations in {gene_name} (column: {seq_column})")

        # Get sequences and remove NaN values
        sequences = self.hs_pop_seq[seq_column].dropna()
        n_total = len(sequences)

        if n_total == 0:
            return {'error': f'No sequences found for {gene_name}'}

        # Determine reference sequence
        if reference_id is None:
            seq_counts = sequences.value_counts()
            reference_seq = seq_counts.index[0]
        else:
            if reference_id not in sequences.index:
                raise ValueError(f"Reference ID {reference_id} not found in sequence data")
            reference_seq = sequences.loc[reference_id]

        # Get gene start position for absolute positioning
        if self.mode == 'gene':
            # In gene mode, use gene_start_positions dict
            gene_start_abs = self.gene_start_positions.get(gene_name, 0)
            if gene_start_abs == 0:
                logger.warning(f"No start position defined for {gene_name}, using relative positions")
            gene_coords = (gene_start_abs, gene_start_abs + len(reference_seq), '+')
        else:
            # In full_mtdna mode, get coordinates from GenBank or alignment
            gene_coords = self._get_gene_coordinates(gene_name, reference_seq)
            gene_start_abs = gene_coords[0] if gene_coords else 0

        # Find mutations
        mutations = []
        mutation_carriers = defaultdict(list)
        relative_positions = {}  # Store relative positions for gene mode

        for idx, seq in sequences.items():
            if seq != reference_seq and len(seq) == len(reference_seq):
                # Find positions where sequences differ
                for pos, (ref_base, seq_base) in enumerate(zip(reference_seq, seq)):
                    if ref_base != seq_base and seq_base in 'ACGT':  # Only count valid base changes
                        # Calculate absolute position
                        abs_position = gene_start_abs + pos + 1  # 1-based
                        mutation = f"{ref_base}{abs_position}{seq_base}"  # Use absolute position
                        relative_position = pos + 1  # 1-based relative position

                        mutations.append(mutation)
                        mutation_carriers[mutation].append(idx)
                        relative_positions[mutation] = relative_position

        # Calculate mutation frequencies
        mutation_counts = Counter(mutations)
        high_freq_mutations = {}

        for mutation, count in mutation_counts.items():
            frequency = count / n_total
            if frequency >= self.frequency_threshold:
                high_freq_mutations[mutation] = {
                    'count': count,
                    'frequency': frequency,
                    'carriers': mutation_carriers[mutation],
                    'absolute_position': int(''.join(filter(str.isdigit, mutation))),
                    'relative_position': relative_positions.get(mutation),  # Gene-relative position
                    'gene_coordinates': gene_coords,
                    'gene_column': seq_column  # Track which column was used
                }

        logger.info(f"Found {len(high_freq_mutations)} high frequency mutations in {gene_name}")

        return {
            'gene_name': gene_name,
            'gene_column': seq_column,
            'n_sequences': n_total,
            'reference_length': len(reference_seq),
            'total_mutations': len(mutation_counts),
            'high_freq_mutations': high_freq_mutations,
            'frequency_threshold': self.frequency_threshold,
            'gene_coordinates': gene_coords,
            'gene_start_position': gene_start_abs,
            'mode': self.mode
        }
    
    def _get_gene_coordinates(self, gene_name, gene_sequence):
        """NEW: Get gene coordinates, using alignment as fallback."""
        # First try to get from GenBank data
        if gene_name in self.gene_coordinates:
            return self.gene_coordinates[gene_name]
        
        # Try common name variations
        name_variants = [gene_name.upper(), gene_name.lower(), 
                        gene_name.replace('nad', 'ND'), gene_name.replace('ND', 'nad')]
        
        for variant in name_variants:
            if variant in self.gene_coordinates:
                logger.info(f"Found coordinates for {gene_name} using variant {variant}")
                return self.gene_coordinates[variant]
        
        # Fallback to sequence alignment
        logger.info(f"Using sequence alignment to find coordinates for {gene_name}")
        coords = self._find_gene_by_alignment(gene_sequence, gene_name)
        
        if coords:
            self.gene_coordinates[gene_name] = coords
            return coords
        
        logger.warning(f"Could not determine coordinates for {gene_name}")
        return (0, len(gene_sequence), '+')  # Default to relative positioning
    
    def test_mutation_haplogroup_association(self, gene_name: str, mutation: str, 
                                           carriers: List, test_method: str = 'auto') -> Dict[str, Any]:
        """
        Test association between a specific mutation and haplogroups.
        
        Parameters:
        -----------
        gene_name : str
            Name of the gene
        mutation : str
            Mutation identifier (e.g., 'A8993G' with absolute position)
        carriers : List
            List of sample indices carrying the mutation
        test_method : str
            Statistical test method ('chi2', 'fisher', 'auto')
            
        Returns:
        --------
        Dict containing test results
        """
        # Get sample IDs for carriers
        # Note: carriers contains DataFrame index values (from sequences.items() iteration),
        # not integer positions, so we use .loc[] for label-based lookup
        carrier_ids = []
        for idx in carriers:
            try:
                if idx in self.hs_pop_seq.index:
                    carrier_id = str(self.hs_pop_seq.loc[idx, 'SampleID'])
                    carrier_ids.append(carrier_id)
                elif 'SampleID' in self.hs_pop_seq.columns and str(idx) in self.hs_pop_seq['SampleID'].astype(str).values:
                    # idx might already be a SampleID
                    carrier_ids.append(str(idx))
            except (IndexError, KeyError) as e:
                logger.warning(f"Could not get SampleID for index {idx}: {e}")
                continue
        
        if not carrier_ids:
            return {'error': 'No valid carrier IDs found'}
        
        # Create mutation status for all samples
        seq_df = self.hs_pop_seq.copy()
        seq_df['has_mutation'] = seq_df['SampleID'].isin(carrier_ids)
        
        # Prepare haplogroup data for merging with consistent data types
        haplogroup_merge_df = self.haplogroup_df.copy()
        haplogroup_merge_df[self.id_column] = haplogroup_merge_df[self.id_column].astype(str)
        
        # Merge with haplogroup data
        try:
            merged_df = seq_df.merge(
                haplogroup_merge_df[[self.id_column, 'MajorHaplogroup']], 
                left_on='SampleID',
                right_on=self.id_column, 
                how='inner'
            )
        except Exception as e:
            logger.error(f"Merge failed for {gene_name}:{mutation}: {e}")
            return {'error': f'Failed to merge data: {str(e)}'}
        
        if merged_df.empty:
            return {'error': 'No matching samples between sequence and haplogroup data'}
        
        # Create contingency table
        try:
            contingency = pd.crosstab(merged_df['MajorHaplogroup'], merged_df['has_mutation'])
        except Exception as e:
            return {'error': f'Failed to create contingency table: {str(e)}'}
        
        # Ensure we have both mutation carriers and non-carriers
        if contingency.shape[1] < 2:
            return {'error': 'No variation in mutation status'}
        
        # Statistical testing
        result = {
            'gene': gene_name,
            'mutation': mutation,
            'contingency_table': contingency,
            'n_carriers': contingency[True].sum() if True in contingency.columns else 0,
            'n_total': contingency.sum().sum()
        }
        
        # Choose test method
        if test_method == 'auto':
            # Use chi-square if all expected frequencies > 5, otherwise Fisher's exact
            try:
                expected_freq = stats.contingency.expected_freq(contingency)
                min_expected = expected_freq.min().min()
                test_method = 'chi2' if min_expected >= 5 else 'fisher'
            except:
                test_method = 'fisher'  # Default to Fisher's if calculation fails
        
        try:
            if test_method == 'chi2':
                chi2, p_value, dof, expected = chi2_contingency(contingency)
                result.update({
                    'test_method': 'chi2',
                    'statistic': chi2,
                    'p_value': p_value,
                    'degrees_of_freedom': dof
                })
            elif test_method == 'fisher':
                # For Fisher's exact test with >2x2 tables, we need to use a different approach
                if contingency.shape == (2, 2):
                    odds_ratio, p_value = fisher_exact(contingency)
                    result.update({
                        'test_method': 'fisher_exact',
                        'odds_ratio': odds_ratio,
                        'p_value': p_value
                    })
                else:
                    # Use Freeman-Halton extension for larger tables (approximate)
                    chi2, p_value, dof, expected = chi2_contingency(contingency)
                    result.update({
                        'test_method': 'chi2_approx',
                        'statistic': chi2,
                        'p_value': p_value,
                        'degrees_of_freedom': dof,
                        'note': 'Fisher exact approximated with chi-square for >2x2 table'
                    })
            
            # Calculate effect sizes
            if contingency.shape[1] == 2:  # Binary mutation status
                mutation_freq_by_hg = {}
                for hg in contingency.index:
                    total_hg = contingency.loc[hg].sum()
                    carriers_hg = contingency.loc[hg, True] if True in contingency.columns else 0
                    mutation_freq_by_hg[hg] = carriers_hg / total_hg if total_hg > 0 else 0
                
                result['mutation_frequencies'] = mutation_freq_by_hg
                
        except Exception as e:
            result['error'] = f"Statistical test failed: {str(e)}"
        
        return result
    
    def detect_haplogroup_markers(self) -> Dict[str, Any]:
        """
        NEW: Detect haplogroup-specific markers and validate against known markers.
        
        Returns:
        --------
        Dict containing detected markers and validation results
        """
        logger.info(f"Detecting haplogroup markers (threshold: {self.marker_threshold*100}%)")
        
        if not self.association_results:
            logger.warning("No association results available. Run analysis first.")
            return {}
        
        detected_markers = {}
        validation_results = {}
        
        # For each haplogroup, find mutations with high frequency within that group
        for hg in self.haplogroup_df['MajorHaplogroup'].unique():
            if hg == 'Unknown':
                continue
                
            hg_markers = []
            
            # Check all mutations across all genes
            for gene, gene_data in self.association_results.items():
                if 'associations' not in gene_data:
                    continue
                    
                for mutation, assoc in gene_data['associations'].items():
                    if 'mutation_frequencies' in assoc:
                        hg_freq = assoc['mutation_frequencies'].get(hg, 0)
                        
                        # Check if this mutation is a potential marker for this haplogroup
                        if hg_freq >= self.marker_threshold:
                            # Also check that it's not too common in other haplogroups
                            other_freqs = [freq for other_hg, freq in assoc['mutation_frequencies'].items() 
                                         if other_hg != hg]
                            max_other_freq = max(other_freqs) if other_freqs else 0
                            
                            # Require significant difference
                            if hg_freq - max_other_freq >= 0.3:  # At least 30% difference
                                abs_pos = None
                                try:
                                    abs_pos = int(''.join(filter(str.isdigit, mutation)))
                                except:
                                    pass
                                
                                hg_markers.append({
                                    'mutation': mutation,
                                    'gene': gene,
                                    'frequency_in_hg': hg_freq,
                                    'max_frequency_other_hg': max_other_freq,
                                    'absolute_position': abs_pos,
                                    'p_value': assoc.get('p_value', 1.0)
                                })
            
            if hg_markers:
                detected_markers[hg] = sorted(hg_markers, key=lambda x: x['frequency_in_hg'], reverse=True)
                
                # Validate against known markers
                validation_results[hg] = self._validate_markers(hg, hg_markers)
        
        self.haplogroup_markers = {
            'detected': detected_markers,
            'validation': validation_results
        }
        
        logger.info(f"Detected markers for {len(detected_markers)} haplogroups")
        
        return self.haplogroup_markers
    
    def _validate_markers(self, haplogroup, detected_markers):
        """NEW: Validate detected markers against known haplogroup markers."""
        validation = {
            'known_markers_found': [],
            'novel_markers': [],
            'validation_score': 0.0
        }
        
        if haplogroup not in KNOWN_HAPLOGROUP_MARKERS:
            validation['note'] = f"No known markers in database for {haplogroup}"
            return validation
        
        known_positions = set(KNOWN_HAPLOGROUP_MARKERS[haplogroup])
        detected_positions = set()
        
        for marker in detected_markers:
            if marker['absolute_position']:
                detected_positions.add(marker['absolute_position'])
                
                # Check if this matches a known marker
                if marker['absolute_position'] in known_positions:
                    validation['known_markers_found'].append({
                        'position': marker['absolute_position'],
                        'mutation': marker['mutation'],
                        'frequency': marker['frequency_in_hg']
                    })
                else:
                    validation['novel_markers'].append({
                        'position': marker['absolute_position'],
                        'mutation': marker['mutation'],
                        'frequency': marker['frequency_in_hg']
                    })
        
        # Calculate validation score
        overlap = len(known_positions.intersection(detected_positions))
        total_known = len(known_positions)
        validation['validation_score'] = overlap / total_known if total_known > 0 else 0.0
        
        return validation
    
    def run_full_analysis(self, genes: Optional[List[str]] = None, 
                         correct_multiple_testing: bool = True,
                         correction_method: str = 'fdr_bh',
                         test_method: str = 'auto',
                         detect_markers: bool = True) -> Dict[str, Any]:
        """
        ENHANCED: Run complete association analysis with marker detection.
        
        Parameters:
        -----------
        genes : List[str], optional
            List of genes to analyze. If None, analyzes all available genes.
        correct_multiple_testing : bool
            Whether to apply multiple testing correction
        correction_method : str
            Method for multiple testing correction ('fdr_bh', 'bonferroni', etc.)
        test_method : str
            Statistical test method ('auto', 'chi2', 'fisher')
        detect_markers : bool
            Whether to run haplogroup marker detection
            
        Returns:
        --------
        Dict containing complete analysis results with markers
        """
        if genes is None:
            genes = self.gene_names
        
        logger.info(f"Running enhanced analysis for {len(genes)} genes")
        
        all_results = {}
        all_p_values = []
        mutation_info = []
        
        for gene in genes:
            if gene not in self.gene_names:
                logger.warning(f"Gene {gene} not found in data")
                continue
            
            # Identify mutations (now with absolute positions)
            mutation_data = self.identify_mutations(gene, reference_id=self.reference_id)
            if 'error' in mutation_data:
                logger.warning(f"Error in {gene}: {mutation_data['error']}")
                continue
            
            gene_results = {
                'mutation_data': mutation_data,
                'associations': {}
            }
            
            # Test associations for each high frequency mutation
            for mutation, mut_info in mutation_data['high_freq_mutations'].items():
                assoc_result = self.test_mutation_haplogroup_association(
                    gene, mutation, mut_info['carriers'], test_method = test_method
                )
                
                if 'error' not in assoc_result:
                    gene_results['associations'][mutation] = assoc_result
                    all_p_values.append(assoc_result.get('p_value', 1.0))
                    mutation_info.append({
                        'gene': gene,
                        'mutation': mutation,
                        'p_value': assoc_result.get('p_value', 1.0)
                    })
            
            all_results[gene] = gene_results
        
        # Multiple testing correction
        if correct_multiple_testing and all_p_values:
            rejected, p_corrected, alpha_sidak, alpha_bonf = multipletests(
                all_p_values, method=correction_method
            )
            
            # Update results with corrected p-values
            idx = 0
            for gene in all_results:
                for mutation in all_results[gene]['associations']:
                    if idx < len(p_corrected):
                        all_results[gene]['associations'][mutation]['p_corrected'] = p_corrected[idx]
                        all_results[gene]['associations'][mutation]['significant_corrected'] = rejected[idx]
                        idx += 1
        
        # Generate summary statistics
        summary = self._generate_summary(all_results, correction_method if correct_multiple_testing else None)
        
        self.mutation_results = all_results
        self.association_results = all_results
        self.summary_stats = summary
        
        # NEW: Detect haplogroup markers
        if detect_markers:
            markers = self.detect_haplogroup_markers()
            summary['haplogroup_markers'] = markers
        
        return {
            'results': all_results,
            'summary': summary,
            'correction_method': correction_method if correct_multiple_testing else None,
            'haplogroup_markers': self.haplogroup_markers if detect_markers else {}  # NEW
        }
    
    def _generate_summary(self, results: Dict, correction_method: Optional[str] = None) -> Dict[str, Any]:
        """Generate summary statistics from analysis results."""
        summary = {
            'total_genes': len(results),
            'total_mutations': 0,
            'significant_associations': 0,
            'genes_with_associations': 0,
            'by_gene': {},
            'most_significant': [],
            'haplogroup_effects': defaultdict(list)
        }
        
        all_significant = []
        
        for gene, gene_data in results.items():
            if 'associations' not in gene_data:
                continue
                
            gene_summary = {
                'n_mutations': len(gene_data['mutation_data']['high_freq_mutations']),
                'n_associations_tested': len(gene_data['associations']),
                'significant_uncorrected': 0,
                'significant_corrected': 0
            }
            
            summary['total_mutations'] += gene_summary['n_mutations']
            
            for mutation, assoc in gene_data['associations'].items():
                p_val = assoc.get('p_value', 1.0)
                
                if p_val < 0.05:
                    gene_summary['significant_uncorrected'] += 1
                    
                if correction_method and assoc.get('significant_corrected', False):
                    gene_summary['significant_corrected'] += 1
                    
                all_significant.append({
                    'gene': gene,
                    'mutation': mutation,
                    'p_value': p_val,
                    'p_corrected': assoc.get('p_corrected', 1.0)
                })
            
            if gene_summary['significant_uncorrected'] > 0:
                summary['genes_with_associations'] += 1
            
            summary['by_gene'][gene] = gene_summary
        
        # Sort by significance
        all_significant.sort(key=lambda x: x['p_corrected'] if correction_method else x['p_value'])
        summary['most_significant'] = all_significant[:10]  # Top 10
        
        summary['significant_associations'] = len([x for x in all_significant if x['p_value'] < 0.05])
        
        return summary
    
    def plot_haplogroup_markers(self, figsize: Tuple[int, int] = (16, 10),
                               save_path: Optional[str] = None) -> plt.Figure:
        """
        NEW: Plot detected haplogroup markers with validation information.
        
        Parameters:
        -----------
        figsize : Tuple[int, int]
            Figure size
        save_path : str, optional
            Path to save the figure
            
        Returns:
        --------
        matplotlib Figure object
        """
        if not self.haplogroup_markers:
            logger.warning("No haplogroup markers detected. Run detect_haplogroup_markers() first.")
            return None
        
        detected = self.haplogroup_markers.get('detected', {})
        validation = self.haplogroup_markers.get('validation', {})
        
        if not detected:
            logger.warning("No markers to plot")
            return None
        
        # Prepare data for plotting
        plot_data = []
        
        for hg, markers in detected.items():
            for marker in markers:
                val_info = validation.get(hg, {})
                is_known = any(km['position'] == marker['absolute_position'] 
                             for km in val_info.get('known_markers_found', []))
                
                plot_data.append({
                    'Haplogroup': hg,
                    'Position': marker['absolute_position'] or 0,
                    'Frequency': marker['frequency_in_hg'],
                    'Gene': marker['gene'],
                    'Mutation': marker['mutation'],
                    'Known_Marker': is_known,
                    'Validation_Score': val_info.get('validation_score', 0)
                })
        
        if not plot_data:
            return None
        
        plot_df = pd.DataFrame(plot_data)
        
        # Create figure with subplots
        fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=figsize)
        fig.suptitle('Haplogroup Marker Analysis', fontsize=16, fontweight='bold')
        
        # Plot 1: Marker frequencies by haplogroup
        marker_counts = plot_df.groupby('Haplogroup').size()
        bars1 = ax1.bar(marker_counts.index, marker_counts.values)
        ax1.set_title('Number of Detected Markers by Haplogroup')
        ax1.set_xlabel('Haplogroup')
        ax1.set_ylabel('Number of Markers')
        ax1.tick_params(axis='x', rotation=45)
        
        # Color bars by validation score
        for i, (hg, bar) in enumerate(zip(marker_counts.index, bars1)):
            val_score = plot_df[plot_df['Haplogroup'] == hg]['Validation_Score'].iloc[0]
            color = plt.cm.RdYlGn(val_score)  # Red for low, green for high validation
            bar.set_color(color)
        
        # Plot 2: Marker positions along mtDNA
        for hg in plot_df['Haplogroup'].unique():
            hg_data = plot_df[plot_df['Haplogroup'] == hg]
            known_markers = hg_data[hg_data['Known_Marker']]
            novel_markers = hg_data[~hg_data['Known_Marker']]
            
            if not known_markers.empty:
                ax2.scatter(known_markers['Position'], [hg] * len(known_markers), 
                           c='green', s=60, marker='o', label='Known' if hg == plot_df['Haplogroup'].iloc[0] else "", alpha=0.7)
            if not novel_markers.empty:
                ax2.scatter(novel_markers['Position'], [hg] * len(novel_markers), 
                           c='red', s=60, marker='^', label='Novel' if hg == plot_df['Haplogroup'].iloc[0] else "", alpha=0.7)
        
        ax2.set_title('Marker Positions in mtDNA Genome')
        ax2.set_xlabel('mtDNA Position')
        ax2.set_ylabel('Haplogroup')
        ax2.legend()
        
        # Plot 3: Frequency distribution
        ax3.hist(plot_df['Frequency'], bins=20, alpha=0.7, edgecolor='black')
        ax3.axvline(x=self.marker_threshold, color='red', linestyle='--', 
                   label=f'Threshold ({self.marker_threshold*100}%)')
        ax3.set_title('Distribution of Marker Frequencies')
        ax3.set_xlabel('Frequency within Haplogroup')
        ax3.set_ylabel('Number of Markers')
        ax3.legend()
        
        # Plot 4: Validation scores
        val_scores = plot_df.groupby('Haplogroup')['Validation_Score'].first()
        bars4 = ax4.bar(val_scores.index, val_scores.values)
        ax4.set_title('Validation Scores by Haplogroup')
        ax4.set_xlabel('Haplogroup')
        ax4.set_ylabel('Validation Score (Known Markers Found)')
        ax4.tick_params(axis='x', rotation=45)
        ax4.set_ylim(0, 1)
        
        # Color bars by score
        for bar, score in zip(bars4, val_scores.values):
            bar.set_color(plt.cm.RdYlGn(score))
        
        plt.tight_layout()
        
        if save_path:
            plt.savefig(save_path, dpi=300, bbox_inches='tight')
            logger.info(f"Haplogroup markers plot saved to {save_path}")
        
        return fig
    
    def plot_association_heatmap(self, max_genes: int = 20, max_mutations_per_gene: int = 10,
                                figsize: Tuple[int, int] = (15, 10), 
                                save_path: Optional[str] = None,
                                p_value_threshold: float = 0.05,
                                use_corrected_p: bool = True,
                                highlight_markers: bool = True) -> plt.Figure:
        """
        ENHANCED: Create a heatmap showing mutation-haplogroup associations with marker highlighting.
        
        Parameters:
        -----------
        max_genes : int
            Maximum number of genes to include
        max_mutations_per_gene : int
            Maximum mutations per gene to show
        figsize : Tuple[int, int]
            Figure size
        save_path : str, optional
            Path to save the figure
        p_value_threshold : float
            P-value threshold for highlighting significant associations
        use_corrected_p : bool
            Whether to use corrected p-values if available
        highlight_markers : bool
            Whether to highlight detected haplogroup markers
            
        Returns:
        --------
        matplotlib Figure object
        """
        if not self.association_results:
            raise ValueError("No association results available. Run analysis first.")
        
        # Prepare data for heatmap
        heatmap_data = []
        
        # NEW: Get marker positions for highlighting
        marker_positions = set()
        if highlight_markers and self.haplogroup_markers:
            for hg, markers in self.haplogroup_markers.get('detected', {}).items():
                for marker in markers:
                    if marker['absolute_position']:
                        marker_positions.add(marker['absolute_position'])
        
        gene_count = 0
        for gene, gene_data in self.association_results.items():
            if gene_count >= max_genes:
                break
                
            if 'associations' not in gene_data:
                continue
            
            # Sort mutations by significance
            mutations = []
            for mutation, assoc in gene_data['associations'].items():
                p_key = 'p_corrected' if (use_corrected_p and 'p_corrected' in assoc) else 'p_value'
                p_val = assoc.get(p_key, 1.0)
                mutations.append((mutation, p_val, assoc))
            
            mutations.sort(key=lambda x: x[1])  # Sort by p-value
            
            mut_count = 0
            for mutation, p_val, assoc in mutations[:max_mutations_per_gene]:
                if 'mutation_frequencies' in assoc:
                    # NEW: Check if this is a detected marker
                    abs_pos = None
                    try:
                        abs_pos = int(''.join(filter(str.isdigit, mutation)))
                    except:
                        pass
                    is_marker = abs_pos in marker_positions if abs_pos else False
                    
                    for hg, freq in assoc['mutation_frequencies'].items():
                        heatmap_data.append({
                            'Gene': gene,
                            'Mutation': f"{gene}:{mutation}",
                            'Haplogroup': hg,
                            'Frequency': freq,
                            'P_value': p_val,
                            'Significant': p_val < p_value_threshold,
                            'Is_Marker': is_marker,  # NEW
                            'Absolute_Position': abs_pos,  # NEW
                            '-log10_P': -np.log10(max(p_val, 1e-300))
                        })
                mut_count += 1
            
            if mut_count > 0:
                gene_count += 1
        
        if not heatmap_data:
            raise ValueError("No data available for heatmap")
        
        heatmap_df = pd.DataFrame(heatmap_data)
        
        # Create pivot table for heatmap
        pivot_freq = heatmap_df.pivot_table(
            index='Mutation', columns='Haplogroup', 
            values='Frequency', fill_value=0
        )
        
        pivot_sig = heatmap_df.pivot_table(
            index='Mutation', columns='Haplogroup', 
            values='Significant', fill_value=False
        )
        
        # NEW: Create marker indicator
        pivot_marker = heatmap_df.pivot_table(
            index='Mutation', columns='Haplogroup', 
            values='Is_Marker', fill_value=False
        )
        
        # Create figure
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=figsize, gridspec_kw={'width_ratios': [4, 1]})
        
        # Main heatmap
        sns.heatmap(pivot_freq, 
                   cmap='Reds', 
                   annot=False,
                   cbar_kws={'label': 'Mutation Frequency'},
                   ax=ax1)
        
        # Add significance markers and haplogroup markers
        for i, mutation in enumerate(pivot_freq.index):
            for j, hg in enumerate(pivot_freq.columns):
                # Significance marker (blue asterisk)
                if pivot_sig.loc[mutation, hg]:
                    ax1.text(j + 0.5, i + 0.3, '*', 
                            ha='center', va='center', 
                            color='blue', fontsize=12, fontweight='bold')
                
                # NEW: Haplogroup marker (gold circle)
                if highlight_markers and pivot_marker.loc[mutation, hg]:
                    ax1.text(j + 0.5, i + 0.7, '●', 
                            ha='center', va='center', 
                            color='gold', fontsize=10, fontweight='bold')
        
        title_text = 'Mutation Frequencies by Haplogroup\n(* = significant, ● = haplogroup marker)'
        ax1.set_title(title_text, fontsize=14, fontweight='bold')
        ax1.set_xlabel('Haplogroup', fontweight='bold')
        ax1.set_ylabel('Mutation (Gene:Position)', fontweight='bold')
        
        # Rotate labels for better readability
        ax1.tick_params(axis='x', rotation=45)
        ax1.tick_params(axis='y', rotation=0)
        
        # Summary statistics subplot
        ax2.axis('off')
        
        # Calculate summary stats for display
        n_significant = heatmap_df['Significant'].sum()
        n_markers = heatmap_df['Is_Marker'].sum()  # NEW
        n_total = len(heatmap_df)
        
        summary_text = f"""
ASSOCIATION SUMMARY

Total Tests: {n_total}
Significant: {n_significant}
% Significant: {100*n_significant/n_total:.1f}%

Haplogroup Markers: {n_markers}
% Markers: {100*n_markers/n_total:.1f}%

Genes Analyzed: {len(pivot_freq.index)}
Haplogroups: {len(pivot_freq.columns)}

P-value threshold: {p_value_threshold}
Using {'corrected' if use_corrected_p else 'uncorrected'} p-values
Marker threshold: {self.marker_threshold*100:.0f}%

Most significant associations:
"""
        
        # Add top significant associations
        top_sig = heatmap_df[heatmap_df['Significant']].nsmallest(5, 'P_value')
        for _, row in top_sig.iterrows():
            marker_indicator = " (M)" if row['Is_Marker'] else ""
            summary_text += f"\n{row['Mutation']} in {row['Haplogroup']}{marker_indicator}"
            summary_text += f"\n  p = {row['P_value']:.2e}"
        
        ax2.text(0.05, 0.95, summary_text, transform=ax2.transAxes, 
                fontsize=10, verticalalignment='top', fontfamily='monospace',
                bbox=dict(boxstyle='round,pad=0.5', facecolor='lightgray', alpha=0.8))
        
        plt.tight_layout()
        
        if save_path:
            plt.savefig(save_path, dpi=300, bbox_inches='tight')
            logger.info(f"Enhanced heatmap saved to {save_path}")
        
        return fig
    
    def plot_manhattan_style(self, figsize: Tuple[int, int] = (16, 8),
                            save_path: Optional[str] = None,
                            use_corrected_p: bool = True,
                            color_by_gene: bool = True,
                            highlight_markers: bool = True) -> plt.Figure:
        """
        ENHANCED: Create a Manhattan-style plot with marker highlighting.
        
        Parameters:
        -----------
        figsize : Tuple[int, int]
            Figure size
        save_path : str, optional
            Path to save the figure
        use_corrected_p : bool
            Whether to use corrected p-values if available
        color_by_gene : bool
            Whether to color points by gene
        highlight_markers : bool
            Whether to highlight haplogroup markers
            
        Returns:
        --------
        matplotlib Figure object
        """
        if not self.association_results:
            raise ValueError("No association results available. Run analysis first.")
        
        # NEW: Get marker positions for highlighting
        marker_positions = set()
        if highlight_markers and self.haplogroup_markers:
            for hg, markers in self.haplogroup_markers.get('detected', {}).items():
                for marker in markers:
                    if marker['absolute_position']:
                        marker_positions.add(marker['absolute_position'])
        
        # Prepare data
        plot_data = []
        x_pos = 0
        gene_positions = {}
        
        for gene, gene_data in self.association_results.items():
            if 'associations' not in gene_data:
                continue
                
            gene_start = x_pos
            
            for mutation, assoc in gene_data['associations'].items():
                p_key = 'p_corrected' if (use_corrected_p and 'p_corrected' in assoc) else 'p_value'
                p_val = assoc.get(p_key, 1.0)
                
                # NEW: Check if this is a marker
                abs_pos = None
                try:
                    abs_pos = int(''.join(filter(str.isdigit, mutation)))
                except:
                    pass
                is_marker = abs_pos in marker_positions if abs_pos else False
                
                plot_data.append({
                    'Gene': gene,
                    'Mutation': mutation,
                    'Position': x_pos,
                    'P_value': p_val,
                    'Absolute_Position': abs_pos,  # NEW
                    'Is_Marker': is_marker,  # NEW
                    '-log10_P': -np.log10(max(p_val, 1e-300)),
                    'Significant': p_val < 0.05
                })
                x_pos += 1
            
            gene_positions[gene] = (gene_start + x_pos - 1) / 2  # Midpoint
        
        if not plot_data:
            raise ValueError("No data available for Manhattan plot")
        
        plot_df = pd.DataFrame(plot_data)
        
        # Create figure
        fig, ax = plt.subplots(figsize=figsize)
        
        if color_by_gene:
            # Color by gene, with special highlighting for markers
            genes = plot_df['Gene'].unique()
            colors = plt.cm.Set3(np.linspace(0, 1, len(genes)))
            gene_colors = dict(zip(genes, colors))
            
            for gene in genes:
                gene_data = plot_df[plot_df['Gene'] == gene]
                
                # Plot non-markers
                non_markers = gene_data[~gene_data['Is_Marker']]
                if not non_markers.empty:
                    ax.scatter(non_markers['Position'], non_markers['-log10_P'], 
                              c=[gene_colors[gene]], label=gene, alpha=0.7, s=50)
                
                # NEW: Plot markers with special highlighting
                if highlight_markers:
                    markers = gene_data[gene_data['Is_Marker']]
                    if not markers.empty:
                        ax.scatter(markers['Position'], markers['-log10_P'], 
                                  c=[gene_colors[gene]], s=100, marker='D', 
                                  edgecolors='gold', linewidths=2, alpha=0.9,
                                  label=f'{gene} (markers)' if not non_markers.empty else gene)
        else:
            # Color by significance and markers
            nonsig_nonmarker = plot_df[(~plot_df['Significant']) & (~plot_df['Is_Marker'])]
            sig_nonmarker = plot_df[plot_df['Significant'] & (~plot_df['Is_Marker'])]
            nonsig_marker = plot_df[(~plot_df['Significant']) & plot_df['Is_Marker']]
            sig_marker = plot_df[plot_df['Significant'] & plot_df['Is_Marker']]
            
            if not nonsig_nonmarker.empty:
                ax.scatter(nonsig_nonmarker['Position'], nonsig_nonmarker['-log10_P'], 
                          c='gray', alpha=0.6, s=30, label='Non-significant')
            if not sig_nonmarker.empty:
                ax.scatter(sig_nonmarker['Position'], sig_nonmarker['-log10_P'], 
                          c='red', alpha=0.8, s=50, label='Significant')
            if not nonsig_marker.empty:
                ax.scatter(nonsig_marker['Position'], nonsig_marker['-log10_P'], 
                          c='blue', alpha=0.8, s=80, marker='D', label='Marker (non-sig)')
            if not sig_marker.empty:
                ax.scatter(sig_marker['Position'], sig_marker['-log10_P'], 
                          c='gold', alpha=0.9, s=100, marker='D', 
                          edgecolors='black', linewidths=1, label='Marker (significant)')
        
        # Add significance threshold line
        ax.axhline(y=-np.log10(0.05), color='red', linestyle='--', alpha=0.7, 
                  label='p = 0.05')
        
        # Add gene labels
        for gene, pos in gene_positions.items():
            ax.text(pos, ax.get_ylim()[1] * 0.02, gene, 
                   ha='center', va='bottom', rotation=45, fontsize=8)
        
        ax.set_xlabel('Mutations (grouped by gene)', fontweight='bold')
        ax.set_ylabel('-log₁₀(p-value)', fontweight='bold')
        title = f"Mutation-Haplogroup Association Significance"
        if use_corrected_p:
            title += " (Multiple Testing Corrected)"
        if highlight_markers:
            title += "\n(◆ = Haplogroup Markers)"
        ax.set_title(title, fontweight='bold', fontsize=14)
        
        # Add grid
        ax.grid(True, alpha=0.3)
        
        # Legend
        if len(genes) <= 10 or not color_by_gene:
            ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
        
        plt.tight_layout()
        
        if save_path:
            plt.savefig(save_path, dpi=300, bbox_inches='tight')
            logger.info(f"Enhanced Manhattan plot saved to {save_path}")
        
        return fig
    
    def export_results(self, output_file: str, format: str = 'excel'):
        """
        ENHANCED: Export analysis results with marker information to file.
        
        Parameters:
        -----------
        output_file : str
            Output file path
        format : str
            Output format ('excel', 'csv')
        """
        if not self.association_results:
            raise ValueError("No results to export. Run analysis first.")
        
        # Prepare export data
        export_data = []
        
        for gene, gene_data in self.association_results.items():
            if 'associations' not in gene_data:
                continue
                
            for mutation, assoc in gene_data['associations'].items():
                # Basic info
                abs_pos = None
                try:
                    abs_pos = int(''.join(filter(str.isdigit, mutation)))
                except:
                    pass
                
                row = {
                    'Gene': gene,
                    'Mutation': mutation,
                    'Absolute_Position': abs_pos,  # NEW
                    'N_Carriers': assoc.get('n_carriers', 0),
                    'N_Total': assoc.get('n_total', 0),
                    'Test_Method': assoc.get('test_method', 'unknown'),
                    'P_Value': assoc.get('p_value', 1.0),
                    'P_Corrected': assoc.get('p_corrected', None),
                    'Significant_Uncorrected': assoc.get('p_value', 1.0) < 0.05,
                    'Significant_Corrected': assoc.get('significant_corrected', False)
                }
                
                # Add haplogroup-specific frequencies
                if 'mutation_frequencies' in assoc:
                    for hg, freq in assoc['mutation_frequencies'].items():
                        row[f'Freq_{hg}'] = freq
                
                export_data.append(row)
        
        export_df = pd.DataFrame(export_data)
        
        if format.lower() == 'excel':
            with pd.ExcelWriter(output_file, engine='openpyxl') as writer:
                export_df.to_excel(writer, sheet_name='Associations', index=False)
                
                # Add summary sheet
                summary_df = pd.DataFrame([self.summary_stats])
                summary_df.to_excel(writer, sheet_name='Summary', index=False)
                
                # NEW: Add haplogroup markers sheet
                if self.haplogroup_markers:
                    marker_data = []
                    for hg, markers in self.haplogroup_markers.get('detected', {}).items():
                        for marker in markers:
                            marker_row = {
                                'Haplogroup': hg,
                                'Mutation': marker['mutation'],
                                'Gene': marker['gene'],
                                'Frequency_in_Haplogroup': marker['frequency_in_hg'],
                                'Max_Frequency_Other_Haplogroup': marker['max_frequency_other_hg'],
                                'Absolute_Position': marker['absolute_position'],
                                'P_Value': marker['p_value']
                            }
                            
                            # Add validation info
                            val_info = self.haplogroup_markers.get('validation', {}).get(hg, {})
                            marker_row['Known_Marker'] = any(
                                km['position'] == marker['absolute_position'] 
                                for km in val_info.get('known_markers_found', [])
                            )
                            marker_row['Validation_Score'] = val_info.get('validation_score', 0)
                            
                            marker_data.append(marker_row)
                    
                    if marker_data:
                        markers_df = pd.DataFrame(marker_data)
                        markers_df.to_excel(writer, sheet_name='Haplogroup_Markers', index=False)
                
        elif format.lower() == 'csv':
            export_df.to_csv(output_file, index=False)
        
        logger.info(f"Enhanced results exported to {output_file}")


def run_haplogroup_association_analysis(hs_pop_seq: pd.DataFrame,
                                       haplogroup_file: str,
                                       genes: Optional[List[str]] = None,
                                       frequency_threshold: float = 0.01,
                                       marker_threshold: float = 0.8,
                                       output_dir: str = './haplogroup_analysis',
                                       generate_plots: bool = True,
                                       test_method: str = 'auto',
                                       reference_id: str = None,
                                       accession_id: str = "NC_012920.1",
                                       mode: str = 'full_mtdna',
                                       sequence_column: Union[str, Dict[str, int]] = None,
                                       filter_sequences: bool = True,
                                       length_filter_method: str = 'mode',
                                       min_sequence_length: int = 20) -> MtDNAHaplogroupAssociation:
    """
    ENHANCED: Convenience function to run complete haplogroup association analysis with markers.

    Parameters:
    -----------
    hs_pop_seq : pd.DataFrame
        DataFrame with sequence data. Format depends on mode:
        - full_mtdna mode: columns ending in '_seq' for gene sequences
        - gene mode: gene-specific columns (e.g., 'RNR1', 'CYTB')
    haplogroup_file : str
        Path to haplogroup CSV file
    genes : List[str], optional
        Specific genes to analyze. If None, analyzes all available genes.
    frequency_threshold : float
        Mutation frequency threshold
    marker_threshold : float
        Haplogroup marker detection threshold (default: 0.8 = 80%)
    output_dir : str
        Output directory for results and plots
    generate_plots : bool
        Whether to generate visualization plots
    test_method : str
        Statistical test method ('auto', 'chi2', 'fisher')
    reference_id : str, optional
        Reference sequence ID or GenBank accession
    accession_id : str
        Default GenBank accession ID
    mode : str
        Analysis mode:
        - 'full_mtdna': Look for columns ending in '_seq' (default)
        - 'gene': Use specific gene columns with mtDNA start positions
    sequence_column : str or Dict[str, int], optional
        For gene mode:
        - str: Single gene column name (e.g., 'RNR1')
        - Dict: Mapping of gene column names to mtDNA start positions
               e.g., {'RNR1': 1671, 'CYTB': 14747, 'ND4L-ND4': 10470}
        For full_mtdna mode: ignored (uses '_seq' columns)
    filter_sequences : bool
        Whether to filter sequences to valid ACGT bases and same length (default: True)
    length_filter_method : str
        Method for length filtering ('mode', 'max', 'min') (default: 'mode')
    min_sequence_length : int
        Minimum acceptable sequence length (default: 100)

    Returns:
    --------
    MtDNAHaplogroupAssociation object with completed enhanced analysis
    """
    import os

    # Create output directory
    os.makedirs(output_dir, exist_ok=True)

    # Initialize enhanced analysis
    analyzer = MtDNAHaplogroupAssociation(
        hs_pop_seq=hs_pop_seq,
        haplogroup_file=haplogroup_file,
        frequency_threshold=frequency_threshold,
        marker_threshold=marker_threshold,
        reference_id=reference_id,
        accession_id=accession_id,
        mode=mode,
        sequence_column=sequence_column,
        filter_sequences=filter_sequences,
        length_filter_method=length_filter_method,
        min_sequence_length=min_sequence_length
    )
    
    # Run enhanced analysis
    logger.info("Running complete enhanced haplogroup association analysis...")
    results = analyzer.run_full_analysis(
        genes=genes,
        correct_multiple_testing=True,
        correction_method='fdr_bh',
        test_method=test_method,
        detect_markers=True  # NEW
    )
    
    # Export enhanced results
    analyzer.export_results(
        os.path.join(output_dir, 'haplogroup_associations_enhanced.xlsx'),
        format='excel'
    )
    
    # Generate enhanced plots if requested
    if generate_plots:
        logger.info("Generating enhanced visualization plots...")
        
        try:
            # Enhanced heatmap with marker highlighting
            analyzer.plot_association_heatmap(
                save_path=os.path.join(output_dir, 'association_heatmap_enhanced.png'),
                figsize=(16, 12),
                highlight_markers=True
            )
        except Exception as e:
            logger.warning(f"Could not generate enhanced heatmap: {e}")
        
        try:
            # Enhanced Manhattan plot with marker highlighting
            analyzer.plot_manhattan_style(
                save_path=os.path.join(output_dir, 'manhattan_plot_enhanced.png'),
                figsize=(18, 8),
                highlight_markers=True
            )
        except Exception as e:
            logger.warning(f"Could not generate enhanced Manhattan plot: {e}")
        
        try:
            # NEW: Haplogroup markers plot
            analyzer.plot_haplogroup_markers(
                save_path=os.path.join(output_dir, 'haplogroup_markers.png'),
                figsize=(16, 10)
            )
        except Exception as e:
            logger.warning(f"Could not generate haplogroup markers plot: {e}")
    
    # Print enhanced summary
    print("\n" + "="*70)
    print("ENHANCED HAPLOGROUP ASSOCIATION ANALYSIS SUMMARY")
    print("="*70)
    
    summary = results['summary']
    print(f"Genes analyzed: {summary['total_genes']}")
    print(f"High frequency mutations found: {summary['total_mutations']}")
    print(f"Genes with significant associations: {summary['genes_with_associations']}")
    print(f"Total significant associations: {summary['significant_associations']}")
    
    # NEW: Haplogroup markers summary
    if 'haplogroup_markers' in summary and summary['haplogroup_markers']:
        markers = summary['haplogroup_markers']
        detected = markers.get('detected', {})
        validation = markers.get('validation', {})
        
        print(f"\nHAPLOGROUP MARKERS DETECTED:")
        print(f"Haplogroups with markers: {len(detected)}")
        
        total_markers = sum(len(m) for m in detected.values())
        print(f"Total markers detected: {total_markers}")
        
        validated_hgs = sum(1 for v in validation.values() 
                           if v.get('validation_score', 0) > 0)
        print(f"Haplogroups with validated markers: {validated_hgs}")
        
        # Show top markers for each haplogroup
        for hg, hg_markers in list(detected.items())[:5]:
            if hg_markers:
                top_marker = hg_markers[0]  # Highest frequency
                val_score = validation.get(hg, {}).get('validation_score', 0)
                print(f"  {hg}: {top_marker['mutation']} (freq: {top_marker['frequency_in_hg']:.1%}, "
                      f"validation: {val_score:.1%})")
    
    if summary['most_significant']:
        print(f"\nMOST SIGNIFICANT ASSOCIATIONS:")
        for i, assoc in enumerate(summary['most_significant'][:5], 1):
            print(f"{i}. {assoc['gene']}:{assoc['mutation']} (p = {assoc['p_value']:.2e})")
    
    print(f"\nResults saved to: {output_dir}")
    print("="*70)
    
    return analyzer


def test_data_compatibility(hs_pop_seq: pd.DataFrame, haplogroup_file: str,
                           id_column: str = 'SampleID',
                           mode: str = 'full_mtdna',
                           sequence_column: Union[str, Dict[str, int]] = None) -> bool:
    """
    Test data compatibility before running full analysis.

    Parameters:
    -----------
    hs_pop_seq : pd.DataFrame
        Sequence dataframe
    haplogroup_file : str
        Path to haplogroup CSV file
    id_column : str
        Column name for sample IDs in haplogroup file
    mode : str
        Analysis mode ('full_mtdna' or 'gene')
    sequence_column : str or Dict[str, int], optional
        For gene mode: gene column(s) to check

    Returns:
    --------
    bool: True if data is compatible, False otherwise
    """
    print(f"Testing data compatibility (mode: {mode})...")

    try:
        # Load haplogroup data
        hg_df = pd.read_csv(haplogroup_file)
        print(f"✓ Loaded haplogroup file with {len(hg_df)} samples")

        # Check sequence data
        print(f"✓ Sequence data has {len(hs_pop_seq)} samples")

        if mode == 'gene':
            # Gene mode: check for specified gene columns
            if sequence_column is None:
                print("❌ ERROR: sequence_column must be provided in gene mode!")
                return False

            if isinstance(sequence_column, dict):
                gene_cols = list(sequence_column.keys())
            else:
                gene_cols = [sequence_column]

            missing_cols = [col for col in gene_cols if col not in hs_pop_seq.columns]
            if missing_cols:
                print(f"❌ ERROR: Gene columns not found: {missing_cols}")
                print(f"   Available columns: {list(hs_pop_seq.columns)}")
                return False

            print(f"✓ Found {len(gene_cols)} gene columns: {gene_cols}")

            # Check sequence validity
            for col in gene_cols:
                valid_seqs = hs_pop_seq[col].dropna()
                if len(valid_seqs) == 0:
                    print(f"⚠️  WARNING: No valid sequences in column '{col}'")
                else:
                    lengths = valid_seqs.str.len()
                    print(f"✓ Column '{col}': {len(valid_seqs)} sequences, "
                          f"lengths {lengths.min()}-{lengths.max()}")
        else:
            # Full mtDNA mode: check for _seq columns
            gene_cols = [col for col in hs_pop_seq.columns if col.endswith('_seq')]
            if not gene_cols:
                print("⚠️  WARNING: No '_seq' columns found in full_mtdna mode")
            else:
                print(f"✓ Found {len(gene_cols)} gene columns: "
                      f"{[col.replace('_seq', '') for col in gene_cols]}")

            # Check for full sequence column (optional in full_mtdna mode)
            if 'sequence' in hs_pop_seq.columns:
                print("✓ Found 'sequence' column for absolute positioning")
            else:
                print("⚠️  WARNING: 'sequence' column not found - will use GenBank for coordinates")

        # Check SampleID situation
        seq_df = hs_pop_seq.copy()
        if 'SampleID' not in seq_df.columns:
            seq_df['SampleID'] = seq_df.index.astype(str)
            print("✓ Created SampleID column from index")
        else:
            seq_df['SampleID'] = seq_df['SampleID'].astype(str)
            print("✓ Using existing SampleID column")

        # Check ID overlap
        seq_ids = set(seq_df['SampleID'].astype(str))
        hg_ids = set(hg_df[id_column].astype(str))
        overlap = len(seq_ids.intersection(hg_ids))

        print(f"✓ ID overlap: {overlap} samples match between datasets")

        if overlap == 0:
            print("❌ ERROR: No matching sample IDs found!")
            print(f"Sample sequence IDs: {list(seq_ids)[:5]}")
            print(f"Sample haplogroup IDs: {list(hg_ids)[:5]}")
            return False
        elif overlap < 100:
            print(f"⚠️  WARNING: Only {overlap} samples overlap - this may limit analysis power")

        # Test a simple merge
        test_merge = seq_df[['SampleID']].merge(
            hg_df[[id_column]],
            left_on='SampleID',
            right_on=id_column,
            how='inner'
        )
        print(f"✓ Test merge successful with {len(test_merge)} samples")

        print("✅ Data compatibility test PASSED")
        print("  Ready for analysis with haplogroup marker detection!")
        return True

    except Exception as e:
        print(f"❌ Data compatibility test FAILED: {e}")
        return False


# Example usage with enhanced features
if __name__ == "__main__":
    print("Enhanced Mitochondrial Haplogroup Association Analysis module loaded!")
    print("\n" + "="*70)
    print("FEATURES:")
    print("="*70)
    print("- Absolute position annotation using GenBank coordinates")
    print("- Haplogroup marker detection and validation")
    print("- Enhanced visualizations with marker highlighting")
    print("- Sequence alignment fallback for gene finding")
    print("- ACGT sequence filtering (removes invalid characters)")
    print("- Uniform sequence length filtering")
    print("- Support for gene mode with multiple gene columns")
    print("\nUse test_data_compatibility() first, then run_haplogroup_association_analysis().")

    print("\n" + "="*70)
    print("MODE 1: Full mtDNA mode (default)")
    print("="*70)
    print("""
# For dataframes with '_seq' suffix columns (e.g., 'cytb_seq', 'cox1_seq')

analyzer = run_haplogroup_association_analysis(
    hs_pop_seq=hs_pop_seq,
    haplogroup_file='all_haplogroups_cleaned.csv',
    frequency_threshold=0.01,
    marker_threshold=0.8,
    output_dir='./haplogroup_analysis',
    mode='full_mtdna',  # default
    filter_sequences=True,  # Filter to valid ACGT
    length_filter_method='mode'  # Keep sequences of most common length
)
""")

    print("\n" + "="*70)
    print("MODE 2: Gene mode - Single gene column")
    print("="*70)
    print("""
# For dataframes with direct gene columns (e.g., 'RNR1', 'CYTB')

analyzer = run_haplogroup_association_analysis(
    hs_pop_seq=hs_pop_seq,
    haplogroup_file='all_haplogroups_cleaned.csv',
    frequency_threshold=0.01,
    marker_threshold=0.8,
    output_dir='./haplogroup_analysis',
    mode='gene',
    sequence_column='RNR1',  # Single gene column
    filter_sequences=True
)
""")

    print("\n" + "="*70)
    print("MODE 3: Gene mode - Multiple gene columns")
    print("="*70)
    print("""
# Analyze multiple gene columns with their mtDNA start positions
# This allows absolute position calculation for each gene

gene_columns = {
    'RNR1': 1671,       # 12S rRNA starts at mtDNA position 1671
    'RNR2': 3230,       # 16S rRNA starts at mtDNA position 3230
    'ND4L-ND4': 10470,  # ND4L-ND4 region starts at position 10470
    'CYTB': 14747,      # Cytochrome B starts at position 14747
}

analyzer = run_haplogroup_association_analysis(
    hs_pop_seq=hs_pop_seq,
    haplogroup_file='all_haplogroups_cleaned.csv',
    frequency_threshold=0.01,
    marker_threshold=0.8,
    output_dir='./haplogroup_analysis',
    mode='gene',
    sequence_column=gene_columns,  # Dict of {column: start_position}
    filter_sequences=True,
    length_filter_method='mode',
    min_sequence_length=100
)
""")

    print("\n" + "="*70)
    print("TESTING DATA COMPATIBILITY")
    print("="*70)
    print("""
# Test compatibility before running analysis

# For full_mtdna mode:
compatible = test_data_compatibility(
    hs_pop_seq=hs_pop_seq,
    haplogroup_file='all_haplogroups_cleaned.csv',
    mode='full_mtdna'
)

# For gene mode with multiple columns:
gene_columns = {'RNR1': 1671, 'CYTB': 14747}
compatible = test_data_compatibility(
    hs_pop_seq=hs_pop_seq,
    haplogroup_file='all_haplogroups_cleaned.csv',
    mode='gene',
    sequence_column=gene_columns
)
""")

    print("\n" + "="*70)
    print("OUTPUT FILES:")
    print("="*70)
    print("""
- haplogroup_associations_enhanced.xlsx: Complete results with:
  - Associations sheet: All mutation-haplogroup associations
  - Summary sheet: Analysis summary statistics
  - Haplogroup_Markers sheet: Detected haplogroup markers
- association_heatmap_enhanced.png: Heatmap of mutation frequencies
- manhattan_plot_enhanced.png: Manhattan-style significance plot
- haplogroup_markers.png: Haplogroup marker visualization
""")

    print("="*70)
    print("For direct class usage:")
    print("="*70)
    print("""
# Direct class instantiation for more control:

analyzer = MtDNAHaplogroupAssociation(
    hs_pop_seq=hs_pop_seq,
    haplogroup_file='all_haplogroups_cleaned.csv',
    frequency_threshold=0.01,
    marker_threshold=0.8,
    mode='gene',
    sequence_column={'RNR1': 1671, 'CYTB': 14747},
    filter_sequences=True
)

# Debug data structure
analyzer.debug_data_structure()

# Run analysis
results = analyzer.run_full_analysis(
    genes=None,  # None = all genes
    correct_multiple_testing=True,
    detect_markers=True
)

# Access results
print(analyzer.gene_names)  # List of gene names
print(analyzer.haplogroup_markers)  # Detected markers
""")

    # Example usage (uncomment and modify as needed):

    # # Test enhanced data compatibility first
    # compatible = test_data_compatibility(
    #     hs_pop_seq=hs_pop_seq,
    #     haplogroup_file='all_haplogroups_cleaned.csv',
    #     mode='gene',
    #     sequence_column={'RNR1': 1671, 'CYTB': 14747}
    # )
    #
    # if compatible:
    #     # Run enhanced analysis
    #     analyzer = run_haplogroup_association_analysis(
    #         hs_pop_seq=hs_pop_seq,
    #         haplogroup_file='all_haplogroups_cleaned.csv',
    #         frequency_threshold=0.01,
    #         marker_threshold=0.8,
    #         output_dir='./enhanced_haplogroup_analysis',
    #         mode='gene',
    #         sequence_column={'RNR1': 1671, 'CYTB': 14747},
    #         filter_sequences=True
    #     )
    #
    #     # Access detected markers
    #     if analyzer.haplogroup_markers:
    #         print("\nDetected haplogroup markers:")
    #         for hg, markers in analyzer.haplogroup_markers['detected'].items():
    #             print(f"{hg}: {len(markers)} markers")
    # else:
    #     print("Please check your data structure before running enhanced analysis")
