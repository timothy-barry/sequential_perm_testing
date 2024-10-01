library(robustDESeq)
library(airway)
data("airway")
offsite_dir <- .get_config_path("LOCAL_SEQUENTIAL_TEST_DATA_DIR")

# obtain the expression matrix
expression_matrix <- assays(airway)$counts
# compute the lib sizes
lib_sizes <- colSums(expression_matrix)
# determine the genes with zero expression across samples; remove those
zero_gene <- apply(X = expression_matrix, MARGIN = 1, FUN = function(r) all(r == 0))
expression_matrix <- expression_matrix[!zero_gene,]
# normalize the counts by dividing by the library size
expression_matrix_norm <- apply(X = expression_matrix, MARGIN = 1, FUN = function(r) {
  (r/lib_sizes) * 10000
}) |> t()
# extract each of the genes
Y_list <- apply(X = expression_matrix_norm, MARGIN = 1, identity, simplify = FALSE)
Y_list <- Y_list[1:500]
# obtain the treatment idx
x <- as.integer(colData(airway)$dex == "trt")

asymptotic_time <- system.time({asymptotic_res <- run_mann_whitney_test_asymptotic(Y_list = Y_list,
                                                                                   x = x,
                                                                                   Z = NULL,
                                                                                   side = "two_tailed",
                                                                                   implementation = "custom",
                                                                                   alpha = 0.1)})
finite_sample_adaptive_time <- system.time({finite_sample_adaptive_res <- run_mann_whitney_test_permutations(Y_list = Y_list,
                                                                                                    x = x,
                                                                                                    Z = NULL,
                                                                                                    side = "two_tailed",
                                                                                                    h = 15,
                                                                                                    alpha = 0.1,
                                                                                                    adaptive_permutation_test = TRUE)})
finite_sample_classical_time <- system.time({finite_sample_classical_res <- run_mann_whitney_test_permutations(Y_list = Y_list,
                                                                                                     x = x,
                                                                                                     Z = NULL,
                                                                                                     side = "two_tailed",
                                                                                                     alpha = 0.1,
                                                                                                     adaptive_permutation_test = FALSE)})
res_list <- list(asymptotic_time = asymptotic_time,
                 finite_sample_adaptive_time = finite_sample_adaptive_time,
                 finite_sample_classical_time = finite_sample_classical_time,
                 asymptotic_res = asymptotic_res,
                 finite_sample_adaptive_res = finite_sample_adaptive_res,
                 finite_sample_classical_res = finite_sample_classical_res)
saveRDS(object = res_list, file = paste0(offsite_dir, "airway_result.rds"))
