library(adaptiveperm)
offsite_dir <- .get_config_path("LOCAL_SEQUENTIAL_TEST_DATA_DIR")
args <- commandArgs(trailingOnly = TRUE)
i <- as.integer(args[1])
n_groups <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_COUNT", unset = NA))

# load the data
dat <- readRDS(paste0(offsite_dir, "gtex/processed/adipose_processed.rds"))
x <- dat$labels_binary
Y_list <- dat$Y_list
B <- as.integer(5 * round(length(Y_list)/0.1))

if (is.na(n_groups)) {
  stop("`SLURM_ARRAY_TASK_COUNT` is missing. Submit this script as a SLURM array job.")
}
if (is.na(i) || i < 1L || i > n_groups) {
  stop("Invalid array task index.")
}

# extract the hypotheses to test
set.seed(1)
ids <- sample(rep(seq(1L, n_groups), length.out = length(Y_list)))
Y_list_sub <- Y_list[ids == i]

cat(sprintf("[Task %d/%d] Running classical finite-sample permutation MW test...\n", i, n_groups))
finite_sample_classical_time <- system.time({
  finite_sample_classical_res <- run_permutation_test(
    Y_list = Y_list_sub,
    x = x,
    side = "two_tailed",
    alpha = 0.1,
    test_statistic = "MW",
    method = "classical",
    B = B)
})
cat(sprintf("[Task %d/%d] Finished classical finite-sample permutation MW test.\n", i, n_groups))

cat(sprintf("[Task %d/%d] Running Besag-Clifford finite-sample permutation MW test...\n", i, n_groups))
finite_sample_bc_time <- system.time({
  finite_sample_bc_res <- run_permutation_test(
    Y_list = Y_list_sub,
    x = x,
    side = "two_tailed",
    h = 15L,
    alpha = 0.1,
    test_statistic = "MW",
    method = "besag_clifford",
    B = B)
})
cat(sprintf("[Task %d/%d] Finished Besag-Clifford finite-sample permutation MW test.\n", i, n_groups))

# save result
finite_sample_classical_res$gene_id <- names(Y_list_sub)
finite_sample_bc_res$gene_id <- names(Y_list_sub)
saveRDS(object = list(finite_sample_classical_res = finite_sample_classical_res,
                      finite_sample_classical_time = finite_sample_classical_time),
        file = paste0(paste0(offsite_dir, "results/classical_res_", i, ".rds")))
saveRDS(object = list(finite_sample_bc_res = finite_sample_bc_res,
                      finite_sample_bc_time = finite_sample_bc_time),
        file = paste0(paste0(offsite_dir, "results/bc_res_", i, ".rds")))
