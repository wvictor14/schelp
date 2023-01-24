#' extracts gene expression data into metadata
#'
#' @param metadata cell metadata
#' @param genes character vector of gene names
#' @param expr count matrix, cells in columns, genes in rows
#' @param cellid_col column in metadata that corresponds to rownames of expr
#'
#' @return tibble
#' @export
#'
extract_gene <- function(metadata, expr, genes, cellid_col = 'cellid') {

  if (sum(genes %in% rownames(expr)) == 0) stop('No genes were matched.')
  if (!all(genes %in% rownames(expr))) {
    warning(
      glue::glue(
        'These genes were not found in the data: {missing}',
        missing = paste0(
          setdiff(genes, rownames(expr)), collapse = ', '))
    )
  }

  # extract those genes that match
  genes <- intersect(genes, rownames(expr))
  extracted_data <- expr[genes,,drop = FALSE] %>%
    as.matrix() %>% t() %>% tibble::as_tibble(rownames = cellid_col)

  # if some cells in expr are not in metadata
  if (!all(colnames(expr) %in% metadata[[cellid_col]])) {
    cells_missing <- length(setdiff(colnames(expr), metadata[[cellid_col]]))
    warning(glue::glue("{cells_missing} cells in expr data are not found in metadata"))
  }

  # if some cells in metadata are not in counts
  if (!all(metadata[[cellid_col]] %in% colnames(expr))) {
    cells_missing <- length(setdiff(metadata[[cellid_col]], colnames(expr)))
    warning(glue::glue("{cells_missing} cells in metadata are not found in expr"))
  }

  metadata <- metadata %>%
    dplyr::left_join(extracted_data, by = cellid_col)

  return(metadata)
}