test_that("plot_dotplot works", {
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
