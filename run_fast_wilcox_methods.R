library(adaptiveperm)
offsite_dir <- .get_config_path("LOCAL_SEQUENTIAL_TEST_DATA_DIR")

run_mann_whitney_test_asymptotic <- function(Y_list, x, Z, side = "two_tailed", implementation = "r", alpha = 0.1) {
  if (!(implementation %in% c("r"))) stop("`implementation` not recognized.")
  alternative <- switch(side, two_tailed = "two.sided", right = "greater", left = "less")
  p_values <- sapply(X = Y_list, FUN = function(y) {
    s_trt <- y[x == 1L]
    s_cntrl <- y[x == 0L]
    fit <- stats::wilcox.test(x = s_trt, y = s_cntrl, alternative = alternative, correct = TRUE, exact = FALSE)
    fit$p.value
  })
  rejected <- stats::p.adjust(p_values, method = "BH") < alpha
  data.frame(p_value = p_values, rejected = rejected)
}

# load the data
dat <- readRDS(paste0(offsite_dir, "gtex/processed/adipose_processed.rds"))
x <- dat$labels_binary
Y_list <- dat$Y_list

# run and time the asymptotic mw test
cat("[1/2] Running asymptotic MW test...\n")
asymptotic_time <- system.time({asymptotic_res <- run_mann_whitney_test_asymptotic(Y_list = Y_list,
                                                                                   x = x,
                                                                                   Z = NULL,
                                                                                   side = "two_tailed",
                                                                                   alpha = 0.1,
                                                                                   implementation = "r")})
cat("[1/2] Finished asymptotic MW test.\n")
asymptotic_res$gene_id <- names(Y_list)

# run and time the adaptive mw test
cat("[2/2] Running anytime-valid finite-sample permutation MW test...\n")
finite_sample_anytime_time <- system.time({
  finite_sample_anytime_res <- run_permutation_test(
    Y_list = Y_list,
    x = x,
    side = "two_tailed",
    h = 15L,
    alpha = 0.1,
    test_statistic = "MW",
    method = "anytime")
})
cat("[2/2] Finished anytime-valid finite-sample permutation MW test.\n")
finite_sample_anytime_res$gene_id <- names(Y_list)

# save the result
res_list <- list(asymptotic_time = asymptotic_time,
                 finite_sample_anytime_time = finite_sample_anytime_time,
                 asymptotic_res = asymptotic_res,
                 finite_sample_anytime_res = finite_sample_anytime_res)
saveRDS(object = res_list, file = paste0(offsite_dir, "results/fast_results.rds"))
