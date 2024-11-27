#!/bin/bash
#SBATCH -c 1
#SBATCH -p shared
#SBATCH --time=0-1:00:00
#SBATCH --mem=4GB
#SBATCH -o %j.out
#SBATCH -e %j.err
Rscript ./download_and_process_data.R 
sbatch submit_fast_rna_seq_analysis.sh
sbatch submit_classical_rna_seq_analysis.sh
Rscript ./analyze_results.R