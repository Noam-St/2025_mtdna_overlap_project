# 2025_mtdna_overlap_project

Mitochondrial DNA overlapping sequences analysis — code and materials for the manuscript "The functional importance of mitochondrial sequences with dual functions: overlapping protein-protein and protein-RNA gene sequences as test cases."

## Repository structure

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

Create the conda environment from `environment.yml`:

```bash
conda env create -f environment.yml
```

## Notes

- This repository contains analysis code, intermediate results and utilities used for the manuscript. Large input data files are not tracked here — see `data/` for expected inputs.
- Update the project title in this file if the manuscript title changes.

If you'd like, I can also:

- Add a short usage example for the main scripts.
- Create a minimal `CONTRIBUTING.md` or `CITATION.cff`.

