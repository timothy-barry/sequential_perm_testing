#!/bin/bash
#SBATCH -c 1
#SBATCH -p shared
#SBATCH --time=0-8:00:00
#SBATCH --mem=4GB
#SBATCH --array=1-200
#SBATCH -o %j.out
#SBATCH -e %j.err
Rscript ./run_classical_perm_test.R ${SLURM_ARRAY_TASK_ID}
