test_that("extract_gene works on 1 gene", {
  dim(counts)
  metadata <- extract_gene(
    metadata, expr = counts, genes = 'FOXP3')
  expect_vector(metadata$FOXP3)
})

test_that("extract_gene works on 5 genes", {
  metadata <- extract_gene(
    metadata, expr = counts,
    genes = c('FOXP3', 'NOX1', 'NXPE1', 'MS4A10', 'CCDC62', 'FKBP6'))
  expect_vector(metadata$FOXP3)
})

test_that("behaviour when no genes are not in counts", {
  expect_error(
    extract_gene(
      metadata, expr = counts, genes = 'NOTAGENE')
  )
})

test_that("behaviour when not all genes are in counts (e.g. 1/3)", {
  expect_warning(
    extract_gene(
      metadata, expr = counts, genes = c('FOXP3', 'NOTAGENE'))
  )
})

test_that("behaviour when cells in counts are not in metadata", {
  expect_warning(
    extract_gene(
      metadata, expr = counts[,1:6], genes = 'FOXP3')
  )
})

test_that("behaviour when cells in metadata are not in counts", {
  expect_warning(
    extract_gene(
      head(metadata), expr = counts, genes = 'FOXP3')
  )
})

