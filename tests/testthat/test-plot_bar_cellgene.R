test_that("plot_bar_cellgene works with 1 gene", {
  counts <- scibd::counts
  metadata <- scibd::metadata

  plot_bar_cellgene(metadata = metadata, genes = 'FOXP3', expr = counts,
                    celltype = ori.Cluster_group)
  plot_bar_cellgene(metadata = metadata, genes = 'NOX1', expr = counts,
                    celltype = ori.Cluster_group)
  plot_bar_cellgene(metadata = metadata, genes = 'MS4A10', expr = counts,
                    celltype = ori.Cluster_group)
  expect_true(TRUE)
})

test_that("plot_bar_cellgene works with 2 genes", {
  counts <- scibd::counts
  metadata <- scibd::metadata

  genes <- c('FOXP3', 'NOX1')
  plot_bar_cellgene(metadata = metadata, genes = genes, expr = counts,
                    celltype = ori.Cluster_group)

  genes <- c('MS4A10', 'NXPE1')
  plot_bar_cellgene(metadata = metadata, genes = genes, expr = counts,
                    celltype = ori.Cluster_group)

  expect_true(TRUE)
})
