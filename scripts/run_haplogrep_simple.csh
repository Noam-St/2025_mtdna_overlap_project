#!/bin/env conda run -n noam_env_min
#$ -cwd
#$ -o "/gpfs0/biores/users/mishmarlab/Noam/multi_function_review/scripts/stdout/"
#$ -e "/gpfs0/biores/users/mishmarlab/Noam/multi_function_review/scripts/stderr/"
#$ -N haplogrep
#$ -q bioinfo.q
#$ -S /bin/bash
#$ -V
#$ -pe shared 20

for file in haplogrep /gpfs0/biores/users/mishmarlab/Noam/multi_function_review/data/fasta_chunks/*fasta
do
filename=$(basename $file)
haplogrep classify --format=fasta --input=$file --output=/gpfs0/biores/users/mishmarlab/Noam/multi_function_review/data/all_hs_chunks/$filename.txt
done
