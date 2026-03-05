library(lubridate)

offsite_dir <- .get_config_path("LOCAL_SEQUENTIAL_TEST_DATA_DIR")
result_dir <- paste0(offsite_dir, "results/")
alpha <- 0.1

format_runtime <- function(total_s) {
  d <- floor(total_s/86400)
  h <- floor((total_s %% 86400)/3600)
  m <- floor((total_s %% 3600)/60)
  s <- total_s %% 60
  if (d > 0) {
    sprintf("%dd %02dh %02dm %05.2fs", d, h, m, s)
  } else if (h > 0) {
    sprintf("%dh %02dm %05.2fs", h, m, s)
  } else {
    sprintf("%dm %05.2fs", m, s)
  }
}

load_chunked_method <- function(result_dir, file_pattern, result_field, time_field, alpha) {
  f_names <- sort(list.files(result_dir, pattern = file_pattern))
  if (length(f_names) == 0L) stop(sprintf("No files found for pattern: %s", file_pattern))

  res_list <- lapply(paste0(result_dir, f_names), readRDS)
  df <- data.table::rbindlist(lapply(res_list, function(r) r[[result_field]]))
  df$rejected <- NULL
  if ("hyp_idx" %in% names(df)) df$hyp_idx <- NULL
  df$rejected <- stats::p.adjust(p = df$p_value, method = "BH") < alpha
  df <- dplyr::arrange(df, gene_id)

  total_s <- sum(sapply(res_list, function(r) r[[time_field]][["elapsed"]]))
  list(df = df, total_s = total_s, files = f_names)
}

classical <- load_chunked_method(
  result_dir = result_dir,
  file_pattern = "^classical_res_[0-9]+\\.rds$",
  result_field = "finite_sample_classical_res",
  time_field = "finite_sample_classical_time",
  alpha = alpha
)

bc <- load_chunked_method(
  result_dir = result_dir,
  file_pattern = "^bc_res_[0-9]+\\.rds$",
  result_field = "finite_sample_bc_res",
  time_field = "finite_sample_bc_time",
  alpha = alpha
)

if (length(classical$files) != length(bc$files)) stop("Classical and BC chunk counts differ.")

fast_res <- readRDS(paste0(result_dir, "fast_results.rds"))
asymptotic_df <- dplyr::arrange(fast_res$asymptotic_res, gene_id)
anytime_df <- dplyr::arrange(fast_res$finite_sample_anytime_res, gene_id)
asymptotic_total_s <- fast_res$asymptotic_time[["elapsed"]]
anytime_total_s <- fast_res$finite_sample_anytime_time[["elapsed"]]

method_dfs <- list(
  classical = classical$df,
  besag_clifford = bc$df,
  anytime = anytime_df,
  asymptotic = asymptotic_df
)

ref_gene_ids <- method_dfs$classical$gene_id
for (name in names(method_dfs)) {
  df <- method_dfs[[name]]
  if (anyDuplicated(df$gene_id) > 0L) stop(sprintf("Duplicate gene_id found in %s results.", name))
  if (!all(df$gene_id == ref_gene_ids)) stop(sprintf("Gene ordering mismatch for %s results.", name))
  if (!all(df$rejected %in% c(TRUE, FALSE))) stop(sprintf("Invalid rejected column in %s results.", name))
}

summary_table <- data.frame(
  method = c(
    "Asymptotic MW",
    "Anytime-valid permutation MW",
    "Classical permutation MW",
    "Besag-Clifford permutation MW"
  ),
  percent_significant = 100 * c(
    mean(asymptotic_df$rejected),
    mean(anytime_df$rejected),
    mean(classical$df$rejected),
    mean(bc$df$rejected)
  ),
  runtime_seconds = c(
    asymptotic_total_s,
    anytime_total_s,
    classical$total_s,
    bc$total_s
  )
)

summary_table$runtime <- vapply(summary_table$runtime_seconds, format_runtime, character(1))
summary_table$percent_significant <- round(summary_table$percent_significant, 4)
summary_table$runtime_seconds <- round(summary_table$runtime_seconds, 3)
summary_table <- summary_table |> dplyr::arrange(runtime_seconds)

print(summary_table)

out_rds <- paste0(result_dir, "method_summary_table.rds")
out_csv <- paste0(result_dir, "method_summary_table.csv")

tryCatch(saveRDS(summary_table, file = out_rds),
         error = function(e) warning(sprintf("Failed to write %s: %s", out_rds, e$message)))
tryCatch(utils::write.csv(summary_table, file = out_csv, row.names = FALSE),
         error = function(e) warning(sprintf("Failed to write %s: %s", out_csv, e$message)))
