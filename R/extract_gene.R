#' extracts gene expression data into metadata
#'
#' First checks if gene names are in counts rownames
#' If genes already exist in metadata, message that they are being overwritten
#' rownames of counts will be joined on a cell identifier column in metadata
#' using left_join
#'
#' @param metadata cell metadata
#' @param cellid_col column in metadata containing cell identifiers, should
#' match rownames(counts). Default is "cellid"
#' @param genes character vector of gene names
#' @param expr count matrix, cells in columns, genes in rows
#'
#' @return tibble
#' @export
#'
extract_gene <- function(metadata, cellid_col = 'cellid', expr, genes) {
  # check if genes are already present in metadata
  if (any(genes %in% colnames(metadata))) {
    message('Genes exist in metadata columns')
    message(glue::glue(
      "Overwriting columns: {overlap}",
      overlap = glue::glue_collapse(intersect(genes, colnames(metadata)),
                                    ", ")
    ))

    # remove existing gene columns
    metadata <- metadata %>%
      dplyr::select(-tidyselect::any_of(genes))
  }

  # collapse gene list
  gene_pattern <- stringr::str_c("\\b(", stringr::str_c(genes, collapse = "|"), ")\\b")
  gene_matches <- stringr::str_subset(rownames(expr), gene_pattern)

  # print genes that did not match
  gene_no_match <- genes[!genes %in% gene_matches]
  if (length(gene_no_match) > 0) {
    print(stringr::str_c("These genes were not found: ",
                         stringr::str_c(gene_no_match, collapse = ", ")))
  }

  # if no matches found, stop
  stopifnot(length(gene_matches) > 0)

  extracted_data <- expr[gene_matches,,drop = FALSE] %>%
    as.matrix() %>% t() %>% tibble::as_tibble(rownames = 'cellid')

  # warn if not all cells in counts are found in metadata
  if (!all(extracted_data$cellid %in% metadata$cellid)) {
    cells_missing <- length(setdiff(extracted_data$cellid, metadata$cellid))
    warning(glue::glue("{cells_missing} cells in expr data are not found in metadata$cellid"))
  }

  # add to metadata
  metadata <- metadata %>%
    dplyr::left_join(extracted_data, by = 'cellid')

  return(metadata)
}
