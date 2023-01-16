#' extracts gene expression data into metadata
#'
#' @param metadata cell metadata
#' @param genes character vector of gene names
#' @param expr count matrix, cells in columns, genes in rows
#'
#' @return tibble
#' @export
#'
extract_gene <- function(metadata, expr, genes) {
  # check if genes are already present in metadata
  if (any(genes %in% colnames(metadata))) {
    message('Genes exist in metadata columns')
    message(glue::glue(
      "Overwriting columns: {overlap}",
      overlap = glue::glue_collapse(intersect(genes, colnames(metadata)),
                                    ", ")
    ))

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
                         stringr::str_c(gene_no_match, collapse = ", ")
    )
    )
  }

  stopifnot(length(gene_matches) > 0)

  extracted_data <- expr[gene_matches,,drop = FALSE] %>%
    as.matrix() %>% t() %>% tibble::as_tibble(rownames = 'cellid')

  if (!all(extracted_data$cellid %in% metadata$cellid)) {

    cells_missing <- length(setdiff(extracted_data$cellid, metadata$cellid))

    warning(glue::glue("{cells_missing} cells in expr data are not found in metadata$cellid"))
  }


  metadata <- metadata %>%
    dplyr::left_join(extracted_data, by = 'cellid')

  return(metadata)
}
