#' plots a heatmap / dotplot
#'
#' @param metadata cells metadata
#' @param genes character vector of gene names
#' @param expr counts matrix
#' @param celltype_col column in metadata that contains celltype cluster labels
#' @param cellid_col column in metadata that contains cell IDs, matches rownames
#' of expr. expects a character
#' @param facet_row variable to facet rows (cells) by
#' @param fontsize_p a multiplier that controls fontsize.
#' @param zscore zscores expression by row (gene)
#' @param zscore_max zscores over this will have same color. Only when zscore =
#' TRUE.
#' @param zscore_min zscores below this will be white.
#'
#'
#' @return ggplot2
#' @export
plot_dotplot <- function(metadata, expr, genes,
                         zscore = TRUE,
                         zscore_max = 3,
                         zscore_min = 0,
                         celltype_col,
                         cellid_col = 'cellid',
                         fontsize_p = 1,
                         facet_row = NULL) {
  # extract gene expr
  metadata <- extract_gene(
    metadata = metadata, expr = expr, cellid_col = cellid_col,
    genes = unlist(genes))

  # tidy
  metadata <- metadata %>%

    # pivot for plotting
    tidyr::pivot_longer(
      cols = tidyselect::any_of(unlist(genes, use.names = FALSE)),
      names_to = 'gene', values_to = 'expression')

  # z score across genes
  if (zscore) {
    metadata <- metadata %>%
      dplyr::group_by(gene)  %>%
      dplyr::mutate(expression = scale(expression))
    label <- 'Mean\nexpression\nZ score'
    limits <- c(zscore_min, zscore_max)
  } else {
    label <- 'Mean\nexpression\nlnTP10K'
    limits <- NULL
  }

  # calculate p_expr and mean expr
  metadata <- metadata  %>%
    dplyr::group_by({{facet_row}}, {{celltype_col}}, gene) %>%

    # summarize
    dplyr::summarize(p_expr = sum(expression > 0) / dplyr::n(),
                     mean_expr = mean(expression),
                     median_expr = median(expression))

  # add gene categories
  if (is.list(genes)) {
    genes <- genes %>% stack() %>%dplyr::rename(gene = values, gene_group = ind)

    metadata <- metadata %>%
      dplyr::left_join(genes, by = 'gene')
    facet_col <- vars(gene_group)
  } else {
    facet_col <- vars()
  }

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
              panel.border = element_rect(color = 'grey', fill = NA),
              strip.background = element_blank(),
              strip.placement = 'outside',
              strip.text = element_text(size = 10*fontsize_p)) +
        scale_color_viridis_c(option = 'mako', direction = -1,
                              name = label, limits = limits,
                              oob = scales::squish) +
        scale_size_binned(labels = scales::percent,
                          name = '% cells\nexpressing',
                          range = c(0.1, 3.5)
        ) +
        guides(size = guide_bins(show.limits = TRUE))
    }

  # add facet for cells
  if (!rlang::quo_is_null(enquo(facet_row))) {
    facet_row <- vars({{facet_row}})
  } else {
    facet_row <- vars()
  }
  plot <- plot +
    facet_grid(
      rows = facet_row, cols = facet_col,
      scales = 'free', space = 'free', switch = 'x')
  plot
}

#' @describeIn plot_dotplot puts sample groups on the x-axis and column facets
#' by gene
#' @param groups column of metadata to plot on x-axis (e.g. disease)
#' @export
plot_dotplot_groups <- function(metadata, expr, genes,
                         zscore = TRUE,
                         zscore_max = 3,
                         zscore_min = 0,
                         celltype_col,
                         cellid_col = 'cellid',
                         fontsize_p = 1,
                         facet_row = NULL,
                         groups) {
  # extract gene expr
  metadata <- extract_gene(
    metadata = metadata, expr = expr, cellid_col = cellid_col,
    genes = unlist(genes))

  # tidy
  metadata <- metadata %>%

    # pivot for plotting
    tidyr::pivot_longer(
      cols = tidyselect::any_of(unlist(genes, use.names = FALSE)),
      names_to = 'gene', values_to = 'expression')

  # z score across genes
  if (zscore) {
    metadata <- metadata %>%
      dplyr::group_by(gene)  %>%
      dplyr::mutate(expression = scale(expression))
    label <- 'Mean\nexpression\nZ score'
    limits <- c(zscore_min, zscore_max)
  } else {
    label <- 'Mean\nexpression\nlnTP10K'
    limits <- NULL
  }

  # calculate p_expr and mean expr
  metadata <- metadata  %>%
    dplyr::group_by({{facet_row}}, {{celltype_col}}, {{groups}}, gene) %>%

    # summarize
    dplyr::summarize(p_expr = sum(expression > 0) / dplyr::n(),
                     mean_expr = mean(expression),
                     median_expr = median(expression))

  # add gene categories
  if (is.list(genes)) {
    genes <- genes %>% stack() %>%dplyr::rename(gene = values, gene_group = ind)

    metadata <- metadata %>%
      dplyr::left_join(genes, by = 'gene')
    facet_col <- vars(gene_group, gene)
  } else {
    facet_col <- vars(gene)
  }

  plot <- metadata %>%
    {
      ggplot(data= ., aes(y = {{celltype_col}},  x = {{groups}},
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
              panel.border = element_rect(color = 'grey', fill = NA),
              strip.background = element_blank(),
              strip.placement = 'outside',
              strip.text = element_text(size = 10*fontsize_p)) +
        scale_color_viridis_c(option = 'mako', direction = -1,
                              name = label, limits = limits,
                              oob = scales::squish) +
        scale_size_binned(labels = scales::percent,
                          name = '% cells\nexpressing',
                          range = c(0.1, 3.5)
        ) +
        guides(size = guide_bins(show.limits = TRUE))
    }

  # add facet for cells
  if (!rlang::quo_is_null(enquo(facet_row))) {
    facet_row <- vars({{facet_row}})
  } else {
    facet_row <- vars()
  }
  plot <- plot +
    facet_grid(
      rows = facet_row, cols = facet_col,
      scales = 'free', space = 'free')
  plot
}
