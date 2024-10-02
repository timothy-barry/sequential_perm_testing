library(robustDESeq)
offsite_dir <- .get_config_path("LOCAL_SEQUENTIAL_TEST_DATA_DIR")

# load the data
dat <- readRDS(paste0(offsite_dir, "gtex/processed/adipose_processed.rds"))
x <- dat$labels_binary
Y_list <- dat$Y_list

# run and time the asymptotic mw test
print("Running asymptotic test")
asymptotic_time <- system.time({asymptotic_res <- run_mann_whitney_test_asymptotic(Y_list = Y_list,
                                                                                   x = x,
                                                                                   Z = NULL,
                                                                                   side = "two_tailed",
                                                                                   implementation = "r",
                                                                                   alpha = 0.1,
                                                                                   order_result_df = FALSE)})
asymptotic_res$gene_id <- names(Y_list)

# run and time the adaptive mw test
print("Running adaptive finite sample test")
finite_sample_adaptive_time <- system.time({finite_sample_adaptive_res <- run_mann_whitney_test_permutations(Y_list = Y_list,
                                                                                                    x = x,
                                                                                                    Z = NULL,
                                                                                                    side = "two_tailed",
                                                                                                    h = 15,
                                                                                                    alpha = 0.1,
                                                                                                    adaptive_permutation_test = TRUE,
                                                                                                    order_result_df = FALSE)})
finite_sample_adaptive_res$gene_id <- names(Y_list)

# save the result
res_list <- list(asymptotic_time = asymptotic_time,
                 finite_sample_adaptive_time = finite_sample_adaptive_time,
                 asymptotic_res = asymptotic_res,
                 finite_sample_adaptive_res = finite_sample_adaptive_res)
saveRDS(object = res_list, file = paste0(offsite_dir, "results/fast_results.rds"))
