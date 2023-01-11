#' plots expression on umap coordinates
#'
#' Use extract_gene so that gene expression exists as columns of metadata.
#'
#' Color is scaled from the minimum to 99th percentile only if the difference in
#' 99th percentile and max is greater than 1 unit, and the 99th percentile isn't
#' equal to 0. This addresses the situation where a very small percentage (<1%)
#' cells have massively higher expression of a gene (common for small clusters),
#' which makes it hard to visualize the variability of the other 99% of the
#' data.
#'
#'
#' @param metadata cells metadata
#' @param genes character vector of gene names
#' @param zscore Whether expression should be zscore or not
#'
#' @return ggplot2
#' @export
#' @import ggplot2
plot_umap_expr <- function(metadata, genes,
                          UMAP1 = 'UMAP_1', UMAP2 = 'UMAP_2', zscore = TRUE,
                          pointsize = 1,
                          alpha = 1,
                          breaks = 4,
                          ncol = NULL,
                          nrow = NULL){

  label <- 'ln(TP10K)'

  metadata <- metadata %>%
    dplyr::select(cellid, .data[[UMAP1]], .data[[UMAP2]],
                  tidyselect::any_of(genes))

  # calculate z score
  if (zscore) {
    metadata <- metadata %>%
      tidyr::pivot_longer(cols = tidyselect::any_of(genes),
                          names_to = 'gene',
                          values_to = 'expression') %>%
      dplyr::group_by(gene)  %>%
      dplyr::mutate(expression_zscore = (expression - mean(expression)) / sd(expression))

    label <- glue::glue('Z score')
  }


  plots <- metadata  %>%

    # split by gene
    split(.$gene) %>%
    purrr::imap(
      .f = ~.x %>% {
        ggplot(.x, aes(x = .data[[UMAP1]], y = .data[[UMAP2]],
                       color = expression)) +
          scattermore::geom_scattermore(pointsize = pointsize,
                                        alpha = alpha) +
          theme(strip.placement = 'outside') +
          theme_bw() +
          theme(panel.border = element_blank(),
                panel.grid = element_blank(),
                axis.line = element_line(),
                aspect.ratio = 1) +
          scale_color_viridis_c(
            limits = c(
              # minimum of expression
              min(.x$expression),

              # if large (>1 unit) difference in max - 99th percentile
              # and 99p is not equal to minimum (happens if 99=0)
              # then set upper limit the 99th percentile
              # otherwise set to max (unchanged)
              #
              # This makes it easier to see when there is non-zero expression
              # in many (~99%) cells e.g. for expression in small clusters
              ifelse(max(.x$expression) - quantile(.x$expression, 0.99) > 1 &
                       quantile(.x$expression, 0.99) != min(.x$expression),
                     quantile(.x$expression, 0.99),
                     max(.x$expression)
              )
            ),
            breaks = scales::breaks_extended(n=breaks),

            oob = scales::squish,
            guide = guide_colorbar(barwidth = 0.5, barheight = 2.5,
                                   ticks.colour = 'black')
          ) +
          labs(color = label, title = .y,
               x = UMAP1, y = UMAP2)
      }
    )
  patchwork::wrap_plots(plots, ncol = ncol, nrow = nrow)
}

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

    metadata %>%
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
