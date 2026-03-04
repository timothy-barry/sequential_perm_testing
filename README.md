# sequential_perm_testing

This code replicates the genomics data analysis reported in the paper "Multiple testing with sequential permutation e-values."
It compares four Mann-Whitney procedures on real RNA-seq data: asymptotic MW, classical finite-sample permutation MW, Besag-Clifford adaptive permutation MW, and anytime-valid adaptive permutation MW.

The script `run_all.sh` reproduces the entire analysis. Two steps of this script --- namely, `sbatch --wait submit_fast_rna_seq_analysis.sh` and `sbatch --wait submit_classical_rna_seq_analysis.sh` --- assume the user is on a SLURM cluster. The fast script runs asymptotic + anytime-valid MW, while the array script runs classical + Besag-Clifford MW in parallel.

To reproduce these scripts, first add the following function to your `~/.Rprofile` file:
```
.get_config_path <- function(dir_name) {
  cmd <- paste0("source ~/.research_config; echo $", dir_name)
  system(command = cmd, intern = TRUE)
}
```
Next, create a file `~/.research_config` and add the following line to this file:

```
LOCAL_SEQUENTIAL_TEST_DATA_DIR="~/my_data_directory"
```
Here, `~/my_data_directory` should be replaced with the location on your machine in which you would like to store the data and results.

The analysis scripts require `adaptiveperm`. The asymptotic MW benchmark is implemented directly in `run_fast_wilcox_methods.R` via `stats::wilcox.test`.
