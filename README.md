# sequential_perm_testing

This code replicates the genomics data analysis reported in the paper "Multiple testing with sequential permutation e-values."

The script `run_all.sh` reproduces the entire analysis. Two steps of this script --- namely, `sbatch submit_fast_rna_seq_analysis.sh` and `sbatch submit_classical_rna_seq_analysis.sh` --- assume the user is on a SLURM cluster.

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