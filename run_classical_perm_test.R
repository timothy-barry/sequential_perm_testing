library(robustDESeq)
offsite_dir <- .get_config_path("LOCAL_SEQUENTIAL_TEST_DATA_DIR")
args <- commandArgs(trailingOnly = TRUE)
i <- as.integer(args[1])

# load the data
dat <- readRDS(paste0(offsite_dir, "gtex/processed/adipose_processed.rds"))
x <- dat$labels_binary
Y_list <- dat$Y_list
#B <- as.integer(round(10 * length(Y_list)/0.1))
B <- as.integer(5 * round(length(Y_list)/0.1))

# extract the hypotheses to test
set.seed(1)
n_groups <- 200L
ids <- sample(rep(seq(1L, n_groups), length.out = length(Y_list)))
Y_list_sub <- Y_list[ids == i]

print("Running classical finite sample test")
finite_sample_classical_time <- system.time({finite_sample_classical_res <- run_mann_whitney_test_permutations(
  Y_list = Y_list_sub,
  x = x,
  Z = NULL,
  side = "two_tailed",
  alpha = 0.1,
  B = B,
  adaptive_permutation_test = FALSE,
  order_result_df = FALSE)})

# save result
finite_sample_classical_res$gene_id <- names(Y_list_sub)
saveRDS(object = list(finite_sample_classical_res = finite_sample_classical_res,
                      finite_sample_classical_time = finite_sample_classical_time),
        file = paste0(paste0(offsite_dir, "results/classical_res_", i, ".rds")))
