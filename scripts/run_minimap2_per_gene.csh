#!/bin/env conda run -n noam_env_min
#$ -cwd
#$ -o "/gpfs0/biores/users/mishmarlab/Noam/multi_function_review/scripts/stdout/"
#$ -e "/gpfs0/biores/users/mishmarlab/Noam/multi_function_review/scripts/stderr/"
#$ -N minimap2_per_gene
#$ -q bioinfo.q
#$ -S /bin/bash
#$ -V
#$ -pe shared 20


python useful_tools/split_fasta_by_gtf.py -g multi_function_review/data/Homo_sapiens.gff3 -f sample_downloader/src/mtDNA_ref/Homo_sapiens/Homo_sapiens.fasta -o multi_function_review/data/rcs_fasta_per_gene -t gene -v

python useful_tools/align_to_each_ref_in_folder.py -r "/gpfs0/biores/users/mishmarlab/Noam/multi_function_review/data/rcs_fasta_per_gene/just_cytb/" -i "/gpfs0/biores/users/mishmarlab/Noam/multi_function_review/data/all_hs_mtdna.fasta" -o /gpfs0/biores/users/mishmarlab/Noam/multi_function_review/data/all_hs_aligned_per_gene/ -x asm20 -t 20 -v

python useful_tools/bam_to_fasta_by_ref.py -b /gpfs0/biores/users/mishmarlab/Noam/multi_function_review/data/all_hs_aligned_per_gene/ND4L_ND4_aligned.bam -r ../Noam/multi_function_review/data/rcs_fasta_per_gene/ND4L_ND4.fasta -o multi_function_review/data/all_hs_aligned_per_gene_fasta/ND4L_ND4_all.fasta -v

python useful_tools/bam_to_fasta_by_ref.py -b /gpfs0/biores/users/mishmarlab/Noam/multi_function_review/data/all_hs_aligned_per_gene/ATP8_ATP6_aligned.bam -r ../Noam/multi_function_review/data/rcs_fasta_per_gene/ATP8_ATP6.fasta -o multi_function_review/data/all_hs_aligned_per_gene_fasta/ATP8_ATP6_all.fasta -v

python useful_tools/bam_to_fasta_by_ref.py -b /gpfs0/biores/users/mishmarlab/Noam/multi_function_review/data/all_hs_aligned_per_gene/ATP8_ATP6_aligned.bam -r ../Noam/multi_function_review/data/rcs_fasta_per_gene/ATP8_ATP6.fasta -o multi_function_review/data/all_hs_aligned_per_gene_fasta/ATP8_ATP6_all.fasta -v

python useful_tools/bam_to_fasta_by_ref.py -b /gpfs0/biores/users/mishmarlab/Noam/multi_function_review/data/all_hs_aligned_per_gene/TRNS2_TRNL2_ND5_aligned.bam -r ../Noam/multi_function_review/data/rcs_fasta_per_gene/TRNS2_TRNL2_ND5.fasta -o multi_function_review/data/all_hs_aligned_per_gene_fasta/TRNS2_TRNL2_ND5_all.fasta -v
python useful_tools/bam_to_fasta_by_ref.py -b /gpfs0/biores/users/mishmarlab/Noam/multi_function_review/data/all_hs_aligned_per_gene/RNR1_aligned.bam -r ../Noam/multi_function_review/data/rcs_fasta_per_gene/RNR1.fasta -o multi_function_review/data/all_hs_aligned_per_gene_fasta/RNR1_all.fasta -v
python useful_tools/bam_to_fasta_by_ref.py -b /gpfs0/biores/users/mishmarlab/Noam/multi_function_review/data/all_hs_aligned_per_gene/RNR2_aligned.bam -r ../Noam/multi_function_review/data/rcs_fasta_per_gene/RNR2.fasta -o multi_function_review/data/all_hs_aligned_per_gene_fasta/RNR2_all.fasta -v
python useful_tools/bam_to_fasta_by_ref.py -b /gpfs0/biores/users/mishmarlab/Noam/multi_function_review/data/all_hs_aligned_per_gene/COX1_aligned.bam -r ../Noam/multi_function_review/data/rcs_fasta_per_gene/COX1.fasta -o multi_function_review/data/all_hs_aligned_per_gene_fasta/COX1_all.fasta -v
python useful_tools/bam_to_fasta_by_ref.py -b /gpfs0/biores/users/mishmarlab/Noam/multi_function_review/data/all_hs_aligned_per_gene/CYTB_aligned.bam -r ../Noam/multi_function_review/data/rcs_fasta_per_gene/CYTB.fasta -o multi_function_review/data/all_hs_aligned_per_gene_fasta/CYTB_all.fasta -v
python useful_tools/bam_to_fasta_by_ref.py -b /gpfs0/biores/users/mishmarlab/Noam/multi_function_review/data/all_hs_aligned_per_gene/ND4_aligned.bam -r ../Noam/multi_function_review/data/rcs_fasta_per_gene/ND4.fasta -o multi_function_review/data/all_hs_aligned_per_gene_fasta/ND4_all.fasta -v

python useful_tools/bam_to_fasta_by_ref.py -b /gpfs0/biores/users/mishmarlab/Noam/multi_function_review/data/all_hs_aligned_per_gene/ATP6_aligned.bam -r ../Noam/multi_function_review/data/rcs_fasta_per_gene/ATP6.fasta -o multi_function_review/data/all_hs_aligned_per_gene_fasta/ATP6_all.fasta -v

python useful_tools/bam_to_fasta_by_ref.py -b /gpfs0/biores/users/mishmarlab/Noam/multi_function_review/data/all_hs_aligned_per_gene/ATP8_aligned.bam -r ../Noam/multi_function_review/data/rcs_fasta_per_gene/ATP8.fasta -o multi_function_review/data/all_hs_aligned_per_gene_fasta/ATP8_all.fasta -v

python useful_tools/bam_to_fasta_by_ref.py -b /gpfs0/biores/users/mishmarlab/Noam/multi_function_review/data/all_hs_aligned_per_gene/ATP8_aligned.bam -r ../Noam/multi_function_review/data/rcs_fasta_per_gene/ATP8.fasta -o multi_function_review/data/all_hs_aligned_per_gene_fasta/ATP8_all.fasta -v
python useful_tools/bam_to_fasta_by_ref.py -b /gpfs0/biores/users/mishmarlab/Noam/multi_function_review/data/all_hs_aligned_per_gene/ND1_aligned.bam -r ../Noam/multi_function_review/data/rcs_fasta_per_gene/ND1.fasta -o multi_function_review/data/all_hs_aligned_per_gene_fasta/ND1_all.fasta -v

python useful_tools/bam_to_fasta_by_ref.py -b /gpfs0/biores/users/mishmarlab/Noam/multi_function_review/data/all_hs_aligned_per_gene/ND2_aligned.bam -r ../Noam/multi_function_review/data/rcs_fasta_per_gene/ND2.fasta -o multi_function_review/data/all_hs_aligned_per_gene_fasta/ND2_all.fasta -v

python useful_tools/bam_to_fasta_by_ref.py -b /gpfs0/biores/users/mishmarlab/Noam/multi_function_review/data/all_hs_aligned_per_gene/ND3_aligned.bam -r ../Noam/multi_function_review/data/rcs_fasta_per_gene/ND3.fasta -o multi_function_review/data/all_hs_aligned_per_gene_fasta/ND3_all.fasta -v

python useful_tools/bam_to_fasta_by_ref.py -b /gpfs0/biores/users/mishmarlab/Noam/multi_function_review/data/all_hs_aligned_per_gene/ND5_aligned.bam -r ../Noam/multi_function_review/data/rcs_fasta_per_gene/ND5.fasta -o multi_function_review/data/all_hs_aligned_per_gene_fasta/ND5_all.fasta -v

python useful_tools/bam_to_fasta_by_ref.py -b /gpfs0/biores/users/mishmarlab/Noam/multi_function_review/data/all_hs_aligned_per_gene/ND6_aligned.bam -r ../Noam/multi_function_review/data/rcs_fasta_per_gene/ND6.fasta -o multi_function_review/data/all_hs_aligned_per_gene_fasta/ND6_all.fasta -v

python useful_tools/bam_to_fasta_by_ref.py -b /gpfs0/biores/users/mishmarlab/Noam/multi_function_review/data/all_hs_aligned_per_gene/COX2_aligned.bam -r ../Noam/multi_function_review/data/rcs_fasta_per_gene/COX2.fasta -o multi_function_review/data/all_hs_aligned_per_gene_fasta/COX2_all.fasta -v

python useful_tools/bam_to_fasta_by_ref.py -b /gpfs0/biores/users/mishmarlab/Noam/multi_function_review/data/all_hs_aligned_per_gene/COX3_aligned.bam -r ../Noam/multi_function_review/data/rcs_fasta_per_gene/COX3.fasta -o multi_function_review/data/all_hs_aligned_per_gene_fasta/COX3_all.fasta -v