# load_data.R

load_data <- function() {
  # 1. Set the working directory
  setwd("~/dev/r/greer_models")

  # 2. Load the counts dataset
  ms_counts <- read_csv("./data/gene_count_lazar_1000.csv")
  counts <- ms_counts %>%
    as.data.frame() %>%
    column_to_rownames(var = "...1")

  # 3. Load the metadata
  meta_MS <- read_csv("./data/meta_data_lazar.csv")

  # 4. Filter counts columns to include only those in the metadata's study_id
  columns_to_keep <- meta_MS$study_id
  subset_data_ms_drop <- counts[, colnames(counts) %in% columns_to_keep]

  # 5. Convert the filtered data to a matrix
  matrix_ms_predrop <- as.matrix(subset_data_ms_drop)

  # 6. Return the results as a list for easy unpacking
  list(
    counts = counts,
    meta_MS = meta_MS,
    matrix_ms_predrop = matrix_ms_predrop
    # Uncomment the following line if pc_genes is necessary in the future:
    # pc_genes = read.delim("~/biomart_protein_coding_genes_withgeneinfo.tsv", stringsAsFactors = FALSE)
  )
}
