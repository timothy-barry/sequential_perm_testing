#!/bin/bash
#SBATCH -c 1
#SBATCH -p shared
#SBATCH --time=0-1:00:00
#SBATCH --mem=4GB
#SBATCH -o %j.out
#SBATCH -e %j.err
Rscript ./run_fast_wilcox_methods.R
