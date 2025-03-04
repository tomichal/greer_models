# prepare_data.R

prepare_data <- function(matrix_ms_predrop, meta_MS) {
  #dropping zero genes
  row_sums_combo <- rowSums(matrix_ms_predrop)
  zero_row_names_combo <- rownames(matrix_ms_predrop)[row_sums_combo == 0]
  print(matrix_ms_predrop)
  matrix_combo_ms_ns_drop <- matrix_ms_predrop[rowSums(matrix_ms_predrop) != 0,]

  #checking if there are zero sum rows
  row_sums_combo_drop <- rowSums(matrix_ms_predrop)
  zero_row_names_combo_after_drop <- rownames(matrix_ms_predrop)[rowSums(matrix_ms_predrop) == 0]
  print(zero_row_names_combo_after_drop)


  print(addmargins(table(meta_MS$pleocytosis)))
  meta_MS$pleocytosis[is.na(meta_MS$pleocytosis)] <- 3

  print(addmargins(table(meta_MS$enhancing)))
  meta_MS$enhancing[is.na(meta_MS$enhancing)] <- "Missing"

  print(addmargins(table(meta_MS$immunetherapyattimeoflp)))
  meta_MS$immunetherapyattimeoflp[is.na(meta_MS$immunetherapyattimeoflp)] <- 2

  # Return the processed data
  list(
    matrix_combo_ms_ns_drop = matrix_combo_ms_ns_drop,
    meta_MS = meta_MS
  )
}