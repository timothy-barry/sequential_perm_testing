#!/bin/bash
#SBATCH -c 1
#SBATCH -p short
#SBATCH --time=0-1:00:00
#SBATCH --mem=4GB
#SBATCH -o %j.out
#SBATCH -e %j.err

module load gcc/14.2.0
module load R/4.4.2
module load conda/miniforge3/24.11.3-0

Rscript ./download_and_process_data.R
sbatch --wait submit_fast_rna_seq_analysis.sh
N_GROUPS="${N_GROUPS:-50}"
sbatch --wait --array=1-"${N_GROUPS}" submit_classical_rna_seq_analysis.sh
Rscript ./analyze_results.R
