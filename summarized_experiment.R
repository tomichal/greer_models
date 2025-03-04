# summarized_experiment.R

summarized_experiment <- function(meta_MS, matrix_combo_ms_ns_drop) {
  se_combo_drop_nieve <- SummarizedExperiment(
    assays = list(counts = matrix_combo_ms_ns_drop), colData = data.frame(meta_MS)
  )

  #Make Quartiles of Number of Non Zero Genes and Total Transcripts and time from collection to sequencing just for ease of model convering
  quartiles_total <- quantile(se_combo_drop_nieve$totaltranscripts, probs = c(0, 0.25, 0.5, 0.75, 1))

  categorize_quartile_total <- function(value, quartiles_total) {
    if (value <= quartiles_total[2]) {
      return("Q1")
    } else if (value <= quartiles_total[3]) {
      return("Q2")
    } else if (value <= quartiles_total[4]) {
      return("Q3")
    } else {
      return("Q4")
    }
  }

  #Add that quartile variable to the Summarized Experiment object
  se_combo_drop_nieve$Quartiletotal <- sapply(se_combo_drop_nieve$totaltranscripts, categorize_quartile_total, quartiles_total = quartiles_total)

  #numnonzero quartiles
  quartiles <- quantile(se_combo_drop_nieve$numnonzero, probs = c(0, 0.25, 0.5, 0.75, 1))

  categorize_quartile <- function(value, quartiles) {
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

  #Add that quartile variable to the Summarized Experiment object
  se_combo_drop_nieve$Quartilenonzero <- sapply(se_combo_drop_nieve$numnonzero, categorize_quartile, quartiles = quartiles)


  #freezer time quartiles
  quartiles <- quantile(se_combo_drop_nieve$coll_to_seq, probs = c(0, 0.25, 0.5, 0.75, 1))

  categorize_quartile <- function(value, quartiles) {
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

  #Add that quartile variable to the Summarized Experiment object
  se_combo_drop_nieve$Quartilefreezer <- sapply(se_combo_drop_nieve$coll_to_seq, categorize_quartile, quartiles = quartiles)

  return(se_combo_drop_nieve)
}