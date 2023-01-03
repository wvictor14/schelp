#' calculate_qc_metrics
#'
#' @param scobj a seurat object
#'
#' @return a tibble
#' @export
#'
#' @examples
calculate_qc_metrics <- function(scobj) {
  metadata <- scobj@meta.data %>% tibble::as_tibble(rownames = "cellid")

  # extract QC metrics
  metadata <- metadata %>%

    # rename columns
    dplyr::rename(n_umi = nCount_RNA,
                  n_gene = nFeature_RNA) %>%

    # calculate log10genesperumi and mtRatio
    dplyr::mutate(
      qc_log10genesperumi = log10(n_gene) / log10(n_umi),
      qc_mtfraction = Seurat::PercentageFeatureSet(scobj, pattern = "^MT-")$nCount_RNA / 100,
      qc_rbfraction = Seurat::PercentageFeatureSet(scobj, pattern = "^RP[SL]")$nCount_RNA / 100)

  metadata
}
