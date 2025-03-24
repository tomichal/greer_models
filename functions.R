#
# This is a file containing function definitions.
# It doesn't actually run any functions, it just defines them.
#

#Load in your data, first file path is meta, second is count, third is protein coding gene info
load_data <- function(metadata_filepath, counts_filepath) {
  meta_data <- read_csv(metadata_filepath) #this makes a metadata dataframe

  counts <- read_csv(counts_filepath) #this makes a count dataframe
  counts <- counts %>%
    as.data.frame() %>%
    column_to_rownames(var = "...1") #this makes a new count dataframe with gene names as row names

  columns_to_keep <- meta_data$study_id #this give you the sampleids in your meta data as a vector 

  subset_data_drop <- counts[, colnames(counts) %in% columns_to_keep] #this subsets the count data so that only sampleids in the meta data are kept in the count data 

  matrix_predrop <- as.matrix(subset_data_drop) #this makes count dataframe (with only meta-data sampleid) a count matric

  list(
    meta_data = meta_data,
    matrix_predrop = matrix_predrop
  )   #this is just lists what me made so we can access these things later when we run the function
}

prepare_data <- function(meta_data, matrix_predrop) {
  #dropping all zero genes 
  row_sums_pre <- rowSums(matrix_predrop) #makes a vectors of the rowSums (total gene counts)
  zero_row_names_combo <- rownames(matrix_predrop)[row_sums_pre == 0] #makes a vector of the gene names that have total zero counts
  matrix_after_drop <- matrix_predrop[rowSums(matrix_predrop) != 0,] #makes a count matrix after dropping all zero gene counts

  #checking if there are zero sum rows
  zero_row_names_after_drop <- rownames(matrix_after_drop)[rowSums(matrix_after_drop) == 0]
  print(zero_row_names_after_drop)


  #below will add a "missing value" to any variable you want to include in the models that has missing data
  #print(addmargins(table(meta_data$pleocytosis)))
  #meta_MS$pleocytosis[is.na(meta_data$pleocytosis)] <- 3

  #print(addmargins(table(meta_data$enhancing)))
  #meta_MS$enhancing[is.na(meta_data$enhancing)] <- "Missing"

  #print(addmargins(table(meta_data$immunetherapyattimeoflp)))
  #meta_MS$immunetherapyattimeoflp[is.na(meta_data$immunetherapyattimeoflp)] <- 2

  #Return the processed data
  list(
    matrix_after_drop = matrix_after_drop,
    meta_data = meta_data
  )
}

#Making a Summarized Experiment Object and adding quantile of continuous data to summarized experiment
summarized_experiment_with_quartiles <- function(metadata_file_path, count_matrix_file_path) {
  loaded_data <- load_data(metadata_file_path, count_matrix_file_path)
  prepared_data <- prepare_data(loaded_data$meta_data, loaded_data$matrix_predrop)

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

  se <- SummarizedExperiment(
    assays = list(counts = prepared_data$matrix_after_drop), colData = data.frame(prepared_data$meta_data)
  )

  for (col in  c(
    "totaltranscripts",
    "numnonzero",
    "coll_to_seq",
    "age_at_csf",
    "percentzero"
  )) {
    quartile_name <- paste0("quartile_", col)
    se[[quartile_name]] <- sapply(
      se[[col]],
      categorize_quartiles,

      quartiles = quantile(se[[col]], probs = c(0, 0.25, 0.5, 0.75, 1))
    )
  }

  return(se)
}
