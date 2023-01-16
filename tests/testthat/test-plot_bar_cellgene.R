test_that("plot_bar_cellgene works with 1 gene", {
  counts <- scibd::counts
  metadata <- scibd::metadata

  genes <- c('FOXP3')
  metadata$FOXP3 <- counts[genes,]
  plot_bar_cellgene(metadata = metadata, genes = genes, celltype = ori.Cluster_group)

  genes <- c('NOX1')
  metadata$NOX1 <- counts[genes,]
  plot_bar_cellgene(metadata = metadata, genes = genes, celltype = ori.Cluster_group)

  genes <- c('MS4A10')
  metadata$MS4A10 <- counts[genes,]
  plot_bar_cellgene(metadata = metadata, genes = genes, celltype = ori.Cluster_group)

  expect_true(TRUE)
})

test_that("plot_bar_cellgene works with 2 genes", {
  counts <- scibd::counts
  metadata <- scibd::metadata

  genes <- c('FOXP3', 'NOX1')
  metadata <- dplyr::bind_cols(metadata, t(as.matrix(counts[genes,])))
  plot_bar_cellgene(metadata = metadata, genes = genes, celltype = ori.Cluster_group)

  genes <- c('MS4A10', 'NXPE1')
  metadata <- dplyr::bind_cols(metadata, t(as.matrix(counts[genes,])))
  plot_bar_cellgene(metadata = metadata, genes = genes, celltype = ori.Cluster_group)
  expect_true(TRUE)
})
