#' plots expression by cell cluster
#'
#' Use extract_gene so that gene expression exists as columns of metadata.
#'
#' @param metadata cell metadata
#' @param genes character vector of gene names
#' @param expr count matrix, cells in columns, genes in rows
#'
#' @param cellid_col column in metadata containing cell identifiers, should
#' match rownames(counts). Default is "cellid"
#' @param celltype_col column in metadata referring to cell type identity
#'
#' @param type one of "mean_se", "bar_mean", "bar_mean_se", "bar_median",
#' "violin", "boxplot"
#' @param zscore Whether expression should be zscore or not
#' @param add_jitter adds jittered points default = FALSE,
#' @param jitterheight multiplier that controls jitter height. Default = 0.75
#' @param alpha controls the alpha of jittered points = 0.5
#' @param pointsize = 1,
#' @param expand_mult expand_mult and *_add increases limit of x axis to make
#' space for text labels. Default mult is 0
#' @param expand_add default is 0.45
#' @param fontsize_p Multiplier that controls relative fontsize. Default is 1
#' @param p_expr_highlight_thresh clusters over/under this proportion (0-1) are
#' highlighted
#' @param pexpr_label_size controls size of percent expression labels. Default
#' is 2.25
#'
#' @return ggplot2
#' @export
#' @import ggplot2
plot_genebycell <- function(metadata, cellid_col = 'cellid', expr, genes,
                            celltype_col,
                            type = 'boxplot',
                            zscore = FALSE,
                            add_jitter = FALSE,
                            jitterheight = 0.75,
                            alpha = 0.5,
                            pointsize = 1,
                            expand_mult = 0,
                            expand_add = 0.45,
                            fontsize_p = 1,
                            p_expr_highlight_thresh = 0.05,
                            pexpr_label_size = 2.25){
  label <- 'ln(TP10K)'

  # extract gene expr
  metadata <- extract_gene(
    metadata = metadata, expr = expr, cellid_col = cellid_col, genes = genes)

  # reshape to tidy format
  metadata <- metadata %>%
    dplyr::select(cellid, {{celltype_col}},
                  tidyselect::any_of(genes)) %>%
    tidyr::pivot_longer(cols = tidyselect::any_of(genes),
                        names_to = 'gene',
                        values_to = 'expression')


  # calculate % expression
  metadata <- metadata %>%
    dplyr::group_by( gene, {{celltype_col}}) %>%
    dplyr::summarize(p_expr = sum(expression > 0) / dplyr::n()) %>%
    dplyr::mutate(
      p_expr_highlight = ifelse(p_expr > p_expr_highlight_thresh, 'black', 'grey'),
      p_expr = p_expr %>% scales::percent(accuracy = 0.1) %>%
        stringr::str_pad(width = 5)) %>%
    dplyr::ungroup() %>%

    dplyr::arrange({{celltype_col}}) %>%
    dplyr::mutate(# refactor
      celltype_pexpr = paste0(
        {{celltype_col}} %>% stringr::str_pad(width = 2),
        glue::glue("({p_expr})") %>% stringr::str_pad(width = 8)
      ) %>%
        forcats::fct_inorder()
    ) %>%
    dplyr::left_join(metadata,
                     .,
                     by = c('gene', rlang::as_name(ensym(celltype_col))))

  # calculate z score
  if (zscore) {
    metadata <- metadata %>%
      dplyr::group_by(gene)  %>%
      dplyr::mutate(expression_zscore = (expression - mean(expression)) / sd(expression))

    label <- glue::glue('Z score {label}')
  }

  # start plot
  plot <- ggplot(metadata, aes(y = {{celltype_col}}, x = expression, color = p_expr_highlight)) +
    facet_grid(cols = vars(gene), scales = 'free') +
    geom_text(aes(label = p_expr),
              x = -0.25, size = pexpr_label_size, family = 'mono') +
    scale_color_identity()


  # select geom and summary function
  if (type == 'mean_se') {
    plot <- plot +
      stat_summary(fun.data = mean_se)
  }

  if (type == 'bar_mean') {
    plot <- plot +
      stat_summary(fun = mean, geom = 'bar')
    label <- glue::glue('mean {label}')
  }

  if (type == 'bar_mean_se') {
    plot <- plot +
      stat_summary(fun = mean, geom = 'bar') +
      stat_summary(fun.data = mean_se, geom = "errorbar")
    label <- glue::glue('mean {label}')
  }

  if (type == 'bar_median') {
    plot <- plot +
      stat_summary(fun = median, geom = 'bar')
    label <- glue::glue('median {label}')
  }

  if (type == 'violin') {
    plot <- plot +
      geom_violin()
  }
  if (type == 'boxplot') {
    plot <- plot +
      geom_boxplot(outlier.shape = ifelse(add_jitter, NA, 19))
  }

  if (add_jitter) {

    if (type == 'violin') {
      plot <- plot +
        ggforce::geom_sina(alpha = alpha, size = pointsize, adjust = 1*jitterheight)
    } else {
      plot <- plot +
        geom_jitter(alpha = alpha, size = pointsize, height = 0.2*jitterheight)
    }
  }

  plot <- plot +
    # add theme customization
    scale_x_continuous(expand = expansion(mult=expand_mult, add=expand_add)) +
    theme_bw(base_size = 10*fontsize_p) +
    theme(panel.grid.major.y = element_blank(),
          panel.grid.minor.y = element_blank(),
          panel.spacing = unit(0, 'lines'),
          axis.line = element_line(),
          axis.ticks = element_blank(),
          strip.background = element_blank(),
          strip.placement = 'outside') +
    labs(x = label)
  return(plot)
}
