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

#### Jupyter Notebooks

| File | Description |
|------|-------------|
| `analyze_per_position_conservation.ipynb` | Per-position conservation and diversity analysis across mtDNA sequences |
| `haplogroup_association_analysis.ipynb` | Analysis of haplogroup associations with genetic variants and statistical testing |
| `haplogroup_stats.ipynb` | Descriptive statistics and summaries of haplogroup distributions |
| `mitoriboseq_analysis.ipynb` | MitoRiboSeq data analysis for ribosome profiling |
| `run_codon_dependency_test.ipynb` | Tests for codon usage dependencies and biases |
| `run_dnds_analysis.ipynb` | Runner notebook for dN/dS ratio calculations |
| `mitoriboseq_analysis/analyze_mitoriboseq.ipynb` | Detailed MitoRiboSeq analysis with visualization |

### Pipeline Scripts (`scripts/`)

- `01_minimap2_alignment.sh`: Sequence alignment
- `02_variant_calling.sh`: Variant calling with lofreq
- `03_vcf_to_fasta.sh`: Convert VCF to FASTA
- `04_haplogrep.sh`: Haplogroup classification
- `05_compare_pipelines.py`: Compare mapping pipelines
- `06_extract_gene_regions.py`: Extract gene regions
- `07_hyphy_analysis.sh`: HyPhy selection analysis
- `08_parse_hyphy_results.py`: Parse HyPhy results
- `09_dnds_analysis.py`: Original dN/dS calculations

### Other Directories

- `data/`: Input data (not version controlled)
  - Sequences downloaded from NCBI

- `results/`: Output files (not version controlled)

- `docs/`: Documentation and supplementary information

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
- Update the project title in this file if the manuscript title changes.
