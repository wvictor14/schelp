#' Summarize quality control data
#'
#' Calculates statistics across cells for # genes detected, # UMIs, and
#' mitochondrial gene fraction.
#'
#' Will also count the number of cells failing QC thresholds.
#'
#' Pass a grouped tibble to calculate statistics by group.
#'
#' @param df a tibble or data frame with the following columns:
#' n_gene, n_umi, qc_mtfraction, scrublet_call, qc_log10genesperumi
#'
#' @param ngene_th a numeric that is the lower bound of number of genes detected
#'  per cell
#' @param numi_th a numeric that is the lower bound of the number of UMIs
#' detected per cell
#' @param mtr_th ranges between 0-1. upper bound of the fraction of UMIs mapping
#'  to mitochondrial genes
#'
#' @return a tibble containing summarized qc info
#' @export
#'
summarize_qc <- function(df, ngene_th = 250, numi_th = 500, mtr_th = 0.2) {

  # stat functions
  stats_fn <- list(
    median = ~median(.x, na.rm = TRUE),
    min = ~min(.x, na.rm = TRUE),
    max = ~max(.x, na.rm = TRUE),
    mean = ~mean(.x, na.rm = TRUE)
  )


  summarize(
    df,

    # number of cells
    ncell = n(),

    # number of cells after qc
    ncell_after_qc_filter = sum(n_gene >= ngene_th &
                                  n_umi >= numi_th &
                                  qc_mtfraction <= mtr_th &
                                  scrublet_call == 'Singlet', na.rm =TRUE
    ),


    # number of cells filtered by thresholds
    ncell_lt_ngene_th = sum(n_gene < ngene_th, na.rm = TRUE),
    ncell_lt_numi_th = sum(n_umi < numi_th, na.rm = TRUE),
    ncell_gt_mitoratio_th = sum(qc_mtfraction > mtr_th, na.rm = TRUE),
    ncell_scrublet = sum(scrublet_call == 'Doublet', na.rm = TRUE),

    # some stats over qc measures
    across(c(n_umi, n_gene, qc_log10genesperumi, qc_mtfraction), stats_fn))
}
