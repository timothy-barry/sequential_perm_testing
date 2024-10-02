library(readr)
library(R.utils)
sequential_perm_dir <- .get_config_path("LOCAL_SEQUENTIAL_TEST_DATA_DIR")

# 1. set up offsite directory structure
dir.create(path = paste0(sequential_perm_dir, "gtex"), recursive = TRUE)
dir.create(path = paste0(sequential_perm_dir, "gtex/processed"))
dir.create(path = paste0(sequential_perm_dir, "results/"))

# 2. download raw data from Zenodo
download.file(url = "https://zenodo.org/records/8320659/files/permuted%20datasets_GTEx.zip?download=1",
              destfile = paste0(sequential_perm_dir, "gtex/gtex_raw.zip"))
unzip(zipfile = paste0(sequential_perm_dir, "gtex/gtex_raw.zip"), exdir = paste0(sequential_perm_dir, "gtex"))
file.rename(paste0(sequential_perm_dir, "gtex/permuted datasets_GTEx"),
            paste0(sequential_perm_dir, "gtex/raw"))

# 3. load, process, and save the data
labels_fp <- paste0(sequential_perm_dir, "gtex/raw/AdiposeSubcutaneous.AdiposeVisceralOmentum.conditions.tsv")
matrix_fp <- paste0(sequential_perm_dir, "gtex/raw/AdiposeSubcutaneous.AdiposeVisceralOmentum.gene_readCount.tsv")
labels <- read_tsv(labels_fp, col_types = "c", col_names = FALSE)
labels <- as.character(labels[1,])
expression_matrix <- read_tsv(matrix_fp, col_types = "c", col_names = TRUE)
gene_ids <- expression_matrix$Name_Description
expression_matrix <- expression_matrix[,-1] |> as.matrix()
rownames(expression_matrix) <- gene_ids
lib_sizes <- colSums(expression_matrix)
zero_gene <- apply(X = expression_matrix, MARGIN = 1, FUN = function(r) all(r == 0))
expression_matrix <- expression_matrix[!zero_gene,]
expression_matrix_norm <- apply(X = expression_matrix, MARGIN = 1, FUN = function(r) {
  r/lib_sizes * 1000000
}) |> t()
Y_list <- apply(X = expression_matrix_norm, MARGIN = 1, identity, simplify = FALSE)
tab <- table(labels)
trt_label <- names(which.min(tab))
labels_binary <- integer(length = length(labels))
labels_binary[labels == trt_label] <- 1L
saveRDS(object = list(labels_binary = labels_binary,
                      Y_list = Y_list),
        file = paste0(sequential_perm_dir, "gtex/processed/adipose_processed.rds"))
