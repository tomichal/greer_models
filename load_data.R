# load_data.R

load_data <- function() {
  setwd("~/dev/r/greer_models")

  ms_counts <- read_csv("./data/gene_count_lazar_1000.csv")
  counts <- ms_counts %>%
    as.data.frame() %>%
    column_to_rownames(var = "...1")

  meta_MS <- read_csv("./data/meta_data_lazar.csv")

  columns_to_keep <- meta_MS$study_id

  subset_data_ms_drop <- counts[, colnames(counts) %in% columns_to_keep]

  matrix_ms_predrop <- as.matrix(subset_data_ms_drop)

  # pc_genes <-read.delim("~/biomart_protein_coding_genes_withgeneinfo.tsv", stringsAsFactors = FALSE)

  list(
    counts = counts,
    meta_MS = meta_MS,
    matrix_ms_predrop = matrix_ms_predrop
    # pc_genes = pc_genes
  )
}

