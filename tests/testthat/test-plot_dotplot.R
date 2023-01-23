test_that("plot_dotplot works with character vector gene input", {
  counts <- scibd::counts
  metadata <- scibd::metadata
  plot_dotplot(metadata, counts, 'FOXP3', celltype_col = ori.Cluster_group)
  plot_dotplot(metadata, counts, rownames(counts)[1:20],
               celltype_col = ori.Cluster_group)
  expect_true(TRUE)
})

test_that("plot_dotplot z score works", {
  counts <- scibd::counts
  metadata <- scibd::metadata
  plot_dotplot(metadata, counts, 'FOXP3', celltype_col = ori.Cluster_group)
  plot_dotplot(metadata, counts, rownames(counts)[10:20],
               zscore = FALSE,
               celltype_col = ori.Cluster_group)
  expect_true(TRUE)
})

test_that("plot_dotplot facet_row works", {
  counts <- scibd::counts
  metadata <- scibd::metadata
  plot_dotplot(metadata, counts, rownames(counts)[10:20],
               zscore = FALSE,
               celltype_col = ori.Cluster_group,
               facet_row = Compartment)
  expect_true(TRUE)
})

test_that("plot_dotplot works with gene list input", {
  counts <- scibd::counts
  metadata <- scibd::metadata
  genelist <- list(
    'Targets' = c('NOX1', 'NXPE1', 'MS4A10'),
    'Gene example' = c('FOXP3', 'FKBP6')
  )
  plot_dotplot(metadata, counts, genes = genelist,
               zscore = FALSE,
               celltype_col = ori.Cluster_group)
  expect_true(TRUE)
})

test_that("plot_dotplot works with gene list input + facet_row", {
  counts <- scibd::counts
  metadata <- scibd::metadata
  genelist <- list(
    'Targets' = c('NOX1', 'NXPE1', 'MS4A10'),
    'Gene example' = c('FOXP3', 'FKBP6')
  )
  plot_dotplot(metadata, counts, genes = genelist,
               zscore = FALSE,
               celltype_col = ori.Cluster_group,
               facet_row = Compartment)
  expect_true(TRUE)
})
