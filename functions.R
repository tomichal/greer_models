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

run_model <- function(Y, xFormula, maxIterOptimize = 400, V = NULL) {

  run_zinbFit <- function(zeroInflation) {
    zinbFit(
      Y,
      K = 2,
      X = if (is.null(V)) {
        model.matrix(
          xFormula,
          data = colData(se_combo_drop_nieve)
        )
      } else {
        model.matrix(
          xFormula,
          data = colData(se_combo_drop_nieve),
          V = V
        )
      },
      BPPARAM = BiocParallel::SerialParam(),
      maxiter.optimize = maxIterOptimize,
      verbose = FALSE,
      zeroinflation = zeroInflation
    )
  }

  #ZINB
  zinb_nieve_cohort <- run_zinbFit(zeroInflation = TRUE)
  ZI_ll <- loglik(zinb_nieve_cohort, zinbSim(zinb_nieve_cohort)$counts)
  ZI_aic <- zinbAIC(zinb_nieve_cohort, t(assay(se_combo_drop_nieve)))
  ZI_df <- nParams(zinb_nieve_cohort)

  #neg binomial crude mdoel
  nb_nieve_cohort <- run_zinbFit(zeroInflation = FALSE)
  NB_ll <- loglik(nb_nieve_cohort, zinbSim(nb_nieve_cohort)$counts)
  NB_aic <- zinbAIC(nb_nieve_cohort, t(assay(se_combo_drop_nieve)))
  NB_df <- nParams(nb_nieve_cohort)

  list(
    ZI_ll = ZI_ll,
    ZI_aic = ZI_aic,
    ZI_df = ZI_df,
    NB_ll = NB_ll,
    NB_aic = NB_aic,
    NB_df = NB_df
  )
}

summarized_experiment <- function(matrix_combo_ms_ns_drop, meta_MS) {
  categorize_quartiles <- function(value, quartiles) {
    if (value <= quartiles[2]) {
      return("Q1")
    } else if (value <= quartiles[3]) {
      return("Q2")
    } else if (value <= quartiles[4]) {
      return("Q3")
    } else {
      return("Q4")
    }
  }

  se_combo_drop_nieve <- SummarizedExperiment(
    assays = list(counts = matrix_combo_ms_ns_drop), colData = data.frame(meta_MS)
  )

  #Add that quartile variable to the Summarized Experiment object
  se_combo_drop_nieve$Quartiletotal <- sapply(
    se_combo_drop_nieve$totaltranscripts,
    categorize_quartiles,

    #Make Quartiles of Number of Non Zero Genes and Total Transcripts and time from collection to sequencing just for ease of model convering
    quartiles = quantile(se_combo_drop_nieve$totaltranscripts, probs = c(0, 0.25, 0.5, 0.75, 1))
  )

  #Add that quartile variable to the Summarized Experiment object
  se_combo_drop_nieve$Quartilenonzero <- sapply(
    se_combo_drop_nieve$numnonzero,
    categorize_quartiles,

    #numnonzero quartiles
    quartiles = quantile(se_combo_drop_nieve$numnonzero, probs = c(0, 0.25, 0.5, 0.75, 1))
  )

  #Add that quartile variable to the Summarized Experiment object
  se_combo_drop_nieve$Quartilefreezer <- sapply(
    se_combo_drop_nieve$coll_to_seq,
    categorize_quartiles,

    #freezer time quartiles
    quartiles = quantile(se_combo_drop_nieve$coll_to_seq, probs = c(0, 0.25, 0.5, 0.75, 1))
  )

  return(se_combo_drop_nieve)
}


