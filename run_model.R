run_model <- function(yData, xDataFormula, maxIterOptimize = 25) {

  run_zinbFit <- function(zeroInflation) {
    zinbFit(
      yData,
      K = 2,
      X = model.matrix(
        xDataFormula,
        data = colData(se_combo_drop_nieve)
      ),
      BPPARAM = BiocParallel::SerialParam(),
      maxiter.optimize = maxIterOptimize,
      verbose = FALSE,
      zeroinflation = zeroInflation
    )
  }

  #ZINB
  zinb_nieve_cohort <- run_zinbFit(yData, xDataFormula, zeroInflation = TRUE)
  ZI_ll <- loglik(zinb_nieve_cohort, zinbSim(zinb_nieve_cohort)$counts)
  ZI_aic <- zinbAIC(zinb_nieve_cohort, t(assay(se_combo_drop_nieve)))
  ZI_df <- nParams(zinb_nieve_cohort)

  #neg binomial crude mdoel
  nb_nieve_cohort <- run_zinbFit(yData, xDataFormula, zeroInflation = FALSE)
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