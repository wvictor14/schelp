#' plots expression by cell cluster
#'
#' Use extract_gene so that gene expression exists as columns of metadata.
#'
#'
#'
#' @param metadata cells metadata
#' @param genes character vector of gene names
#' @param celltype column of metadata that denotes cell identity
#' @param type one of "mean_se", "bar_mean", "bar_mean_se", "bar_median",
#' "violin", "boxplot"
#' @param zscore Whether expression should be zscore or not
#' @param add_jitter adds jittered points default = FALSE,
#' @param jitterheight multiplier that controls jitter height. Default = 0.75
#' @param alpha controls the alpha of jittered points = 0.5
#' @param pointsize = 1,
#' @param ncol for facetwrap
#' @param nrow
#' @param fontsize_p Multiplier that controls relative fontsize. Default is 1
#'
#' @return ggplot2
#' @export
#' @import ggplot2
plot_bar_cellgene <- function(metadata, genes, celltype,
                              type = 'boxplot',
                              zscore = FALSE,
                              add_jitter = FALSE,
                              jitterheight = 0.75,
                              alpha = 0.5,
                              pointsize = 1,
                              ncol = NULL,
                              nrow = NULL,
                              fontsize_p = 1){


  label <- 'ln(TP10K)'

  metadata <- metadata %>%
    dplyr::select(cellid, {{celltype}},
                  tidyselect::any_of(genes)) %>%

    # reshape to tidy format
    tidyr::pivot_longer(cols = tidyselect::any_of(genes),
                        names_to = 'gene',
                        values_to = 'expression')


  # calculate % expression
  metadata <- metadata %>%
    dplyr::group_by( gene, {{celltype}}) %>%
    dplyr::summarize(p_expr = sum(expression > 0) / dplyr::n()) %>%
    ungroup() %>%

    arrange({{celltype}}) %>%
    mutate(p_expr = p_expr %>% scales::percent(accuracy = 0.1),
           # refactor
           celltype_pexpr = paste0(
             {{celltype}} %>% stringr::str_pad(width = 2),
             glue::glue("({p_expr})") %>% stringr::str_pad(width = 8)
             ) %>%
             forcats::fct_inorder()
             ) %>%
    dplyr::left_join(metadata,
                     .,
                     by = c('gene', rlang::as_name(ensym(celltype))))

  # calculate z score
  if (zscore) {
    metadata <- metadata %>%
      dplyr::group_by(gene)  %>%
      dplyr::mutate(expression_zscore = (expression - mean(expression)) / sd(expression))

    label <- glue::glue('Z score {label}')
  }

  # start plot
  plot <- ggplot(metadata, aes(y = celltype_pexpr, x = expression)) +
    facet_wrap(vars(gene), nrow = nrow, ncol = ncol, scale = 'free')

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
    theme(strip.placement = 'outside') +
    theme_bw(base_size = 10*fontsize_p) +
    theme(panel.grid.major.y = element_blank(),
          panel.grid.minor.y = element_blank(),
          axis.line = element_line(),
          strip.background = element_blank()) +
    labs(x = label)
  return(plot)
}
