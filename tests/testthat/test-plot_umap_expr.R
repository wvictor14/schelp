test_that("plot_umap_expr works with 1 gene", {
  counts <- scibd::counts
  metadata <- scibd::metadata

  plot_umap_expr(metadata = metadata, genes = 'FOXP3', expr = counts,
                 UMAP1 = 'UMAP1_all', UMAP2 = 'UMAP2_all')
  expect_true(TRUE) # if this runs then a plot was successfully generated
})

test_that("plot_umap_expr works with 2 genes", {
  counts <- scibd::counts
  metadata <- scibd::metadata
  genes <- c('FOXP3', 'NOX1')

  plot_umap_expr(metadata = metadata, genes = genes, expr = counts,
                 UMAP1 = 'UMAP1_all', UMAP2 = 'UMAP2_all')
  expect_true(TRUE)
})
