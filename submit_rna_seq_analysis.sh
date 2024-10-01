#!/bin/bash
#SBATCH -c 1
#SBATCH -p shared
#SBATCH --time=0-4:00:00
#SBATCH --mem=4GB
#SBATCH -o %j.out
#SBATCH -e %j.err
Rscript ./rna_seq_analysis.R $result_dir
