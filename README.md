# Mitochondrial DNA Overlapping Sequences Analysis

Analysis code for "The functional importance of mitochondrial sequences with dual functions: overlapping protein-protein and protein-RNA gene sequences as test cases."

## Repository Structure

- `scripts/`: Analysis scripts
  - `01_minimap2_alignment.sh`: Sequence alignment
  - `02_variant_calling.sh`: Variant calling with lofreq
  - `03_vcf_to_fasta.sh`: Convert VCF to FASTA
  - `04_haplogrep.sh`: Haplogroup classification
  - `05_compare_pipelines.py`: Compare mapping pipelines
  - `06_extract_gene_regions.py`: Extract gene regions
  - `07_hyphy_analysis.sh`: HyPhy selection analysis
  - `08_parse_hyphy_results.py`: Parse HyPhy results
  - `09_dnds_analysis.py`: Original dN/dS calculations
  - (additional scripts...)

- `data/`: Input data (not version controlled)
  - Sequences downloaded from NCBI

- `results/`: Output files (not version controlled)

- `docs/`: Documentation and supplementary information

## Requirements
```bash
conda env create -f environment.yml
