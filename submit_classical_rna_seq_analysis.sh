#!/bin/bash
#SBATCH -c 1
#SBATCH -p short
#SBATCH --time=0-7:00:00
#SBATCH --mem=4GB
#SBATCH -o %j.out
#SBATCH -e %j.err
Rscript ./run_classical_perm_test.R ${SLURM_ARRAY_TASK_ID}
