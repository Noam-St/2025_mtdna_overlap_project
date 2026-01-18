# 2025_mtdna_overlap_project

Mitochondrial DNA overlapping sequences analysis — code and materials for the manuscript "The functional importance of mitochondrial sequences with dual functions: overlapping protein-protein and protein-RNA gene sequences as test cases."

## Repository structure

### Analysis Modules (`bin/`)

#### Python Modules

| File | Description |
|------|-------------|
| `consts.py` | Project constants including mtDNA gene coordinates, codon tables, and reference data |
| `utils.py` | Utility functions for sequence handling, data processing, and common operations |
| `clean_amb_and_gaps.py` | Cleans ambiguous bases (N, R, Y, etc.) and gaps from FASTA sequences |
| `enhanced_mtdna_dnds_analysis_biopython.py` | Calculates dN/dS ratios for mtDNA genes using BioPython; includes mitochondrial genetic code definitions and codon-level analysis |
| `enhanced_mtdna_haplogroup_association.py` | Statistical analysis of haplogroup associations with genetic variants; performs chi-square tests with multiple testing correction |
| `enhanced_sequence_logo_creator.py` | Creates sequence logos for visualization of conservation patterns and mutation frequencies |
| `haplogroup_mutation_analysis_plot.py` | Plotting functions for haplogroup-specific mutation analysis and variant visualization |
| `overlap_asymmetry_analysis.py` | Core analysis module for asymmetric evolution in overlapping reading frames; implements methodology from Szklarczyk et al. (2007) for analyzing selection pressures |
| `overlap_asymmetry_plots.py` | Visualization functions for asymmetry analysis including scatter plots, heatmaps, and comparative displays |
| `overlap_sub_model.py` | Contains `ARFomeNormalized` class for calculating normalized substitution rates in overlapping reading frames |
| `regional_conservation_diversity.py` | Comprehensive module for per-position conservation analysis including Shannon entropy, nucleotide diversity (π), KL divergence from background frequencies, and autocorrelation for codon periodicity detection |
| `haplogroup_analyzer.py` | Haplogroup frequency table generator; processes haplogrep3 output and creates frequency tables in text/HTML formats |

#### Jupyter Notebooks

| File | Description |
|------|-------------|
| `01_preapre_overlap_data.ipynb` | Prepares overlap data by loading gene overlap positions and human sequence data, identifying overlapping regions for each gene pair |
| `02_run_dnds_analysis.ipynb` | Runner notebook for dN/dS ratio calculations across mtDNA genes |
| `03_analyze_per_position_conservation.ipynb` | Per-position conservation and diversity analysis across mtDNA sequences |
| `04_haplogroup_stats.ipynb` | Descriptive statistics and summaries of haplogroup distributions |
| `05_haplogroup_association_analysis.ipynb` | Analysis of haplogroup associations with genetic variants and statistical testing |
| `06_run_codon_dependency_test.ipynb` | Tests for codon usage dependencies and biases |

#### MitoRiboSeq Analysis (`mitoriboseq_analysis/bin/`)

| File | Description |
|------|-------------|
| `01_analyze_mitoriboseq.ipynb` | Detailed MitoRiboSeq analysis with visualization for ribosome profiling |
| `mito_riboseq_analysis.py` | Core functions for MitoRiboSeq data processing and analysis |
| `mito_riboseq_analysis_extension.py` | Extended analysis functions for MitoRiboSeq data |

### Pipeline Scripts (`scripts/`)

Shell scripts and Python CLI tools for data generation and preprocessing.

#### Shell Scripts

| File | Description |
|------|-------------|
| `run_minimap2_per_gene.csh` | Main pipeline script for generating per-gene alignments. Orchestrates the full alignment workflow using the Python scripts below |
| `run_haplogrep_simple.csh` | Runs haplogrep3 for haplogroup classification on FASTA chunks |

#### Python CLI Scripts

| File | Description |
|------|-------------|
| `split_fasta_by_gtf.py` | Splits full mtDNA sequences into individual gene regions using GTF/GFF3 annotations |
| `align_to_each_ref_in_folder.py` | Minimap2 wrapper for generating individual BAM files for each gene reference |
| `bam_to_fasta_by_ref.py` | Converts BAM alignments to aligned FASTA files for each genomic region |
| `split_fasta.py` | Splits large FASTA files into smaller chunks for parallel processing (used for haplogrep input) |

### Other Directories

- `data/`: Input data (not version controlled)
- `results/`: Output files (not version controlled)
- `docs/`: Documentation and supplementary information
- `figures/`: Generated figures and visualizations (not version controlled)

### Base Data Files

The following data files are required to run the analysis notebooks. These are included in the distribution package.

#### Main Data (`data/`)

| File | Description |
|------|-------------|
| `gene_overlap_positions.csv` | Coordinates of overlapping gene regions in human mtDNA, including Gene1, Gene2, Overlap_Start, and Overlap_End positions |
| `hs_row.parquet` | Human reference sequence data in Parquet format containing mtDNA sequence information |

**Data provenance:**

- **`gene_overlap_positions.csv`**: Adapted from the supplementary information in [Yen, K. et al. (2024) "Mitochondrial-derived microproteins: from discovery to function." *Trends in Genetics*](https://www.cell.com/trends/genetics/fulltext/S0168-9525(24)00292-0).

- **`hs_row.parquet`**: Generated from NCBI's human mtDNA RefSeq GenBank file ([NC_012920.1](https://www.ncbi.nlm.nih.gov/nuccore/NC_012920.1?report=genbank&to=16569)).

#### Per-Gene Alignments (`data/per_gene_alignment_fasta/`)

FASTA files containing multiple sequence alignments for each mtDNA gene across human samples:

- `ATP6_all.fasta`, `ATP8_all.fasta`, `ATP8_ATP6_all.fasta`
- `COX1_all.fasta`, `COX2_all.fasta`, `COX3_all.fasta`
- `CYTB_all.fasta`
- `ND1_all.fasta`, `ND2_all.fasta`, `ND3_all.fasta`, `ND4_all.fasta`, `ND4L_all.fasta`, `ND4L_ND4_all.fasta`, `ND5_all.fasta`, `ND6_all.fasta`
- `RNR1_all.fasta`, `RNR2_all.fasta`
- `TRNS2_TRNL2_ND5_all.fasta`

**Data provenance:**

These alignments were generated by running `scripts/run_minimap2_per_gene.csh`, which executes the following pipeline:

1. **Source data**: Human mtDNA sequences downloaded from NCBI (also available at [Mitomap Mitobank](https://www.mitomap.org/foswiki/bin/view/MITOMAP/Mitobank))
2. **Gene annotations**: `Homo_sapiens.gff3` generated from [NC_012920.1](https://www.ncbi.nlm.nih.gov/nuccore/NC_012920.1?report=genbank&to=16569)
3. **Split by gene**: `split_fasta_by_gtf.py` splits full mtDNA sequences into individual gene regions
4. **Alignment**: `align_to_each_ref_in_folder.py` (minimap2 wrapper) generates BAM files for each gene
5. **FASTA conversion**: `bam_to_fasta_by_ref.py` converts BAM alignments to aligned FASTA files per region

#### Haplogroup Data (`data/haplogroups/`)

| File | Description |
|------|-------------|
| `all_haplogroups.csv` | Haplogroup assignments for all samples |
| `haplogroup_markers.xlsx` | Reference table of haplogroup-defining markers and their positions. Downloaded from [Mitomap](https://www.mitomap.org/foswiki/bin/view/MITOMAP/GBFreqInfo) |

**Data provenance:**

- **`all_haplogroups.csv`**: Generated by running [haplogrep3](https://github.com/genepi/haplogrep3) (v3.2.1) on FASTA chunks created using `split_fasta.py`. The pipeline is orchestrated by `scripts/run_haplogrep_simple.csh`.

#### MitoRiboSeq Data (`mitoriboseq_analysis/data/hek_ini/`)

| File | Description |
|------|-------------|
| `mito_cumsum_table.csv` | Cumulative sum table of ribosome profiling data from HEK cells with initiation site mapping |

**Data provenance:**

Generated by running the [Ingolia Lab RiboSeq Snakemake pipeline](https://github.com/ingolia-lab/RiboSeq) ([Li, S.H.-J. et al.](https://github.com/ingolia-lab/RiboSeq)) on mitochondrial ribosome profiling data from [Wakigawa, T. et al.](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE237154). Analyzed samples:
- SAMN36415884
- SAMN36415885
- SAMN36415886

(GEO accession: [GSE237154](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE237154))

## Requirements

### Using Conda (recommended)

Create the conda environment from `environment.yml`:

```bash
conda env create -f environment.yml
conda activate mtdna_analysis
```

### Using pip

Install dependencies from `requirements.txt`:

```bash
pip install -r requirements.txt
```

### Key Dependencies

- **BioPython**: Sequence manipulation and codon analysis
- **NumPy/Pandas**: Data manipulation and analysis
- **SciPy**: Statistical functions
- **Matplotlib/Seaborn**: Visualization
- **Statsmodels**: Statistical modeling and multiple testing correction
- **statannotations**: Statistical annotations on plots
- **adjustText**: Automatic label adjustment in plots

## Notes

- This repository contains analysis code, intermediate results and utilities used for the manuscript. Large input data files are not tracked here — see `data/` for expected inputs.
- The analysis uses the vertebrate mitochondrial genetic code (NCBI translation table 2).
- The methodology for asymmetric evolution analysis in overlapping reading frames is based on Szklarczyk et al. (2007), *Nucleic Acids Research*, 35(10): 3384–3391.