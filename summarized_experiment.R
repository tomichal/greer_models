# summarized_experiment.R


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