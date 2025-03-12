# load_data.R

load_data <- function(metadata_filepath, counts_filepath, pc_genes_filepath) {
  ms_counts <- read_csv(counts_filepath)
  counts <- ms_counts %>%
    as.data.frame() %>%
    column_to_rownames(var = "...1")

  meta_MS <- read_csv(metadata_filepath)

  columns_to_keep <- meta_MS$study_id

  subset_data_ms_drop <- counts[, colnames(counts) %in% columns_to_keep]

  matrix_ms_predrop <- as.matrix(subset_data_ms_drop)

  pc_genes <- read_csv(pc_genes_filepath)

  list(
    # counts = counts,
    meta_MS = meta_MS,
    matrix_ms_predrop = matrix_ms_predrop,
    pc_genes = pc_genes
  )
}

