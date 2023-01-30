# plot_dotplot ----
test_that("plot_dotplot works with character vector gene input", {
  plot_dotplot(metadata, counts, 'FOXP3', celltype_col = ori.Cluster_group)
  plot_dotplot(metadata, counts, rownames(counts)[1:20],
               celltype_col = ori.Cluster_group)
  expect_true(TRUE)
})

test_that("plot_dotplot z score works", {
  plot_dotplot(metadata, counts, 'FOXP3', celltype_col = ori.Cluster_group)
  plot_dotplot(metadata, counts, rownames(counts)[10:20],
               zscore = FALSE,
               celltype_col = ori.Cluster_group)
  expect_true(TRUE)
})

test_that("plot_dotplot facet_row works", {
  plot_dotplot(metadata, counts, rownames(counts)[10:20],
               zscore = FALSE,
               celltype_col = ori.Cluster_group,
               facet_row = Compartment)
  expect_true(TRUE)
})

test_that("plot_dotplot works with gene list input", {
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

test_that("plot_dotplot max zscore", {
  genelist <- list(
    'Targets' = c('NOX1', 'NXPE1', 'MS4A10'),
    'Gene example' = c('FOXP3', 'FKBP6')
  )
  plot_dotplot(metadata, counts, genes = genelist,
               zscore = TRUE,
               zscore_max = 3,
               celltype_col = ori.Cluster_group,
               facet_row = Compartment)
  expect_true(TRUE)
})

# plot_dotplot_groups ----
test_that("plot_dotplot_groups works with named list input", {
  genelist <- list(
    'Targets' = c('NOX1', 'NXPE1', 'MS4A10'),
    'Gene example' = c('FOXP3', 'FKBP6')
  )
  plot_dotplot_groups(
    metadata, counts, genes = genelist,
    zscore = TRUE, zscore_max = 3, celltype_col = ori.Cluster_group,
    facet_row = Compartment, groups = Health_group_pretty)
  expect_true(TRUE)
})

test_that("plot_dotplot_groups works with multi gene character vector input", {
  genelist <- list(
    'Targets' = c('NOX1', 'NXPE1', 'MS4A10'),
    'Gene example' = c('FOXP3', 'FKBP6')
  )
  plot_dotplot_groups(
    metadata, counts, genes = unlist(genelist, use.names = FALSE),
    zscore = TRUE, zscore_max = 3, celltype_col = ori.Cluster_group,
    facet_row = Compartment, groups = Health_group_pretty)
  expect_true(TRUE)
})
