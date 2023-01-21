#' plots a heatmap / dotplot
#'
#' @param metadata cells metadata
#' @param genes character vector of gene names
#' @param expr counts matrix
#' @param celltype_col column in metadata that contains celltype cluster labels
#' @param cellid_col column in metadata that contains cell IDs, matches rownames
#' of expr
#' @param facet_row variable to facet rows (cells) by
#' @param zscore zscores expression by row (gene)
#' @param fontsize_p a multiplier that controls fontsize.
#'
#' @return ggplot2
#' @export
plot_dotplot <- function(metadata, expr, genes, zscore = TRUE,
                         celltype_col,
                         cellid_col = 'cellid',
                         fontsize_p = 1,
                         facet_row = NULL) {
  # extract gene expr
  metadata <- extract_gene(
    metadata = metadata, expr = expr, cellid_col = cellid_col, genes = genes)

  # tidy
  metadata <- metadata %>%

    # pivot for plotting
    tidyr::pivot_longer(cols = tidyselect::any_of(genes),
                        names_to = 'gene',
                        values_to = 'expression')

  label <- 'Mean\nexpression\nlnTP10K'
  # z score across genes
  if (zscore) {
    metadata <- metadata %>%
      dplyr::group_by(gene)  %>%
      dplyr::mutate(expression =
                      (expression - mean(expression)) / sd(expression))
    label <- 'Mean\nexpression\nZ score'
  }

  metadata <- metadata  %>%
    dplyr::group_by({{facet_row}}, {{celltype_col}}, gene) %>%

    # summarize
    dplyr::summarize(p_expr = sum(expression > 0) / dplyr::n(),
                     mean_expr = mean(expression),
                     median_expr = median(expression))

  plot <- metadata %>%
    {
      ggplot(data= ., aes(y = {{celltype_col}},  x = gene,
                          color = mean_expr, size = p_expr))+
        geom_tile(fill = NA, # ensures all cell types are represented on y
                  aes(color = NULL, size = NULL) # unmaps color and size for tile
        ) +
        geom_point(data = dplyr::filter(., p_expr > 0.01) # removes zeroes
        ) +

        theme_bw(base_size = 10*fontsize_p) +
        theme(axis.ticks = element_blank(),
              axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
              axis.title.x = element_blank(),
              axis.title.y = element_blank(),
              aspect.ratio =  NULL,
              panel.grid.major = element_blank(),
              panel.spacing = unit(0,'lines'),
              panel.border = element_rect(color = 'grey', fill = NA)) +
        scale_color_viridis_c(option = 'mako', direction = -1,
                              name = label) +
        scale_size_binned(labels = scales::percent,
                          name = '% cells\nexpressing',
                          range = c(0.1, 3.5)
        ) +
        guides(size = guide_bins(show.limits = TRUE))
    }
  if (!rlang::quo_is_null(enquo(facet_row))) {
    plot <- plot +
      facet_grid(rows = vars({{facet_row}}), scale = 'free', space = 'free')
  }

  plot
}
