library(lubridate)
offsite_dir <- .get_config_path("LOCAL_SEQUENTIAL_TEST_DATA_DIR")
result_dir <- paste0(offsite_dir, "results/")

# slow res
fs <- list.files(result_dir, pattern = "classical_")
# ensure all files present
f_names <- paste0("classical_res_", seq(1, 200), ".rds")
all(sort(fs) == sort(f_names))

# load the files and process the results
res_list <- lapply(paste0(result_dir, f_names), FUN = function(f) readRDS(f))
result_df <- lapply(res_list, FUN = function(r) r$finite_sample_classical_res) |>
  data.table::rbindlist()
result_df$rejected <- result_df$hyp_idx <- NULL
result_df$rejected <- p.adjust(p = result_df$p_value, method = "BH") < 0.1
total_s <- sapply(res_list, FUN = function(r) r$finite_sample_classical_time[["elapsed"]]) |> sum()
seconds_to_period(total_s)
classical_result_df <- result_df |> dplyr::arrange(gene_id)

# load the file for the fast methods
fast_res <- readRDS(paste0(offsite_dir, "results/fast_results.rds"))
asymptotic_result_df <- fast_res$asymptotic_res |> dplyr::arrange(gene_id)
asymptotic_total_s <- fast_res$asymptotic_time[["elapsed"]]
seconds_to_period(asymptotic_total_s)
adaptive_result_df <- fast_res$finite_sample_adaptive_res |> dplyr::arrange(gene_id)
adaptive_total_s <- fast_res$finite_sample_adaptive_time[["elapsed"]]
seconds_to_period(adaptive_total_s)

# compare the results
all(classical_result_df$gene_id == adaptive_result_df$gene_id) & 
  all(classical_result_df$gene_id == asymptotic_result_df$gene_id)
mean(classical_result_df$rejected)
mean(asymptotic_result_df$rejected)
mean(adaptive_result_df$rejected)

# what is the agreement? Among discoveries made by EITHER method A or method B, what percentage are common to both?
compute_discovery_overlap <- function(df_a, df_b) {
  tab <- table(df_a$rejected, df_b$rejected)
  n_union_discoveries <- sum(tab) - tab[1,1]
  n_intersection_discoveries <- tab[2,2]
  n_intersection_discoveries/n_union_discoveries
}
compute_discovery_overlap(classical_result_df, adaptive_result_df)
compute_discovery_overlap(classical_result_df, asymptotic_result_df)
