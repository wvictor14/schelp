test_that("plot_genebycell works with 1 gene", {

  plot_genebycell(metadata = metadata, genes = 'FOXP3', expr = counts,
                    celltype_col = ori.Cluster_group)
  plot_genebycell(metadata = metadata, genes = 'NOX1', expr = counts,
                    celltype_col = ori.Cluster_group)
  plot_genebycell(metadata = metadata, genes = 'MS4A10', expr = counts,
                    celltype_col = ori.Cluster_group)
  expect_true(TRUE)
})

test_that("plot_genebycell works with 2 genes", {
  genes <- c('FOXP3', 'NOX1')
  plot_genebycell(metadata = metadata, genes = genes, expr = counts,
                    celltype_col = ori.Cluster_group)

  genes <- c('MS4A10', 'NXPE1')
  plot_genebycell(metadata = metadata, genes = genes, expr = counts,
                    celltype_col = ori.Cluster_group)

  expect_true(TRUE)
})
test_that("multiplication works", {
  expect_equal(2 * 2, 4)
})
