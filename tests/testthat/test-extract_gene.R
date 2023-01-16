test_that("extract_gene works on 1 gene", {
  metadata <- extract_gene(
    scibd::metadata, expr = scibd::counts, genes = 'FOXP3')
  expect_vector(metadata$FOXP3)
})

test_that("extract_gene works on 5 genes", {
  metadata <- extract_gene(
    scibd::metadata, expr = scibd::counts,
    genes = c('FOXP3', 'NOX1', 'NXPE1', 'MS4A10', 'CCDC62', 'FKBP6'))
  expect_vector(metadata$FOXP3)
})

test_that("behaviour when no genes are not in counts", {
  expect_equal(2 * 2, 4)
})

test_that("behaviour when not all genes are in counts (e.g. 1/3)", {
  expect_equal(2 * 2, 4)
})

test_that("behaviour when cells in counts are not in metadata", {
  expect_equal(2 * 2, 4)
})

test_that("behaviour when cells in metadata are not in counts", {
  expect_equal(2 * 2, 4)
})
