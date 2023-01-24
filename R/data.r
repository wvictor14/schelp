#' genelists
#'
#' Sets of genes that are interesting
#'
#' Smillie 2019 has sets of genes that were able to separate epithelial,
#' stromal, and immune cells
#'
#' Ferraro 2014 has sets of genes that are differentially expressed in Tregs vs
#' Tconv cells.
#'
#' @format ## `genelists`
#' A list containing 5 sets of genes
#'
#' @source Smillie 2019 methods, Ferarro 2014 supplemental
"genelists"

#' Cell metadata
#'
#' An example subset of 1% of the original data 2623 cells are included with
#' the package. To run the package with the full data, overwrite metadata and
#' counts objects in the environment.
#'
#' @format ## `metadata`
#' tibble of cell metadata.
"metadata"

#' Gene expression matrix
#'
#' An example subset of  1% of the original data 2623 cells and 30 genes are
#' included with the package. To run the package with the full data, overwrite
#' metadata and counts objects in the environment.
#'
#' @format ## `counts`
#' Sparse dgc matrix of gene expression by cells
"counts"