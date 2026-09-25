# The original code signalled several failures by printing a message and
# returning NULL, which surfaced far downstream as a confusing error. These
# tests pin the behaviour to real conditions.

test_that("normalizeMatrixN rejects an empty or non-tabular matrix", {
  session <- local_copyscat_session()
  expect_error(normalizeMatrixN(data.frame()), "empty")
  expect_error(normalizeMatrixN("not a matrix"), "data frame or matrix")
})

test_that("normalizeMatrixN errors when every bin is blacklisted", {
  refs <- make_test_references(binsPerChrom = 10L)
  session <- local_copyscat_session(refs)
  counts <- simulate_counts(refs, nCells = 20L, baseCount = 10)

  # A cutoff far above the data blacklists everything.
  expect_error(
    normalizeMatrixN(counts, blacklistCutoff = 1e9, blacklistProp = 0.1,
                     dividingFactor = 1),
    "blacklisted"
  )
})

test_that("normalizeMatrixN errors when no cell passes the empty-bin filter", {
  refs <- make_test_references(binsPerChrom = 10L)
  session <- local_copyscat_session(refs)
  counts <- simulate_counts(refs, nCells = 20L, baseCount = 5)

  expect_error(
    normalizeMatrixN(counts, blacklistCutoff = 1e6, blacklistProp = 0.99,
                     maxZero = 0, dividingFactor = 1),
    "maxZero|blacklisted"
  )
})

test_that("filterCells errors rather than returning an empty matrix", {
  refs <- make_test_references(binsPerChrom = 10L)
  session <- local_copyscat_session(refs)
  counts <- simulate_counts(refs, nCells = 30L, baseCount = 400)
  norm <- normalizeMatrixN(counts, blacklistCutoff = 10, dividingFactor = 1,
                           upperFilterQuantile = 0.99)
  collapsed <- collapseChrom3N(norm, minimumChromValue = 0, tssEnrich = 1,
                               minCPG = 100)

  # No cell can cover 10 000 segments.
  expect_error(
    filterCells(collapsed, minimumSegments = 10000, minDensity = 0),
    "No cell covered"
  )
})

test_that("collapseChrom3N errors when references and matrix disagree", {
  refs <- make_test_references(binsPerChrom = 10L)
  session <- local_copyscat_session(refs)
  counts <- simulate_counts(refs, nCells = 20L, baseCount = 400)
  norm <- normalizeMatrixN(counts, blacklistCutoff = 10, dividingFactor = 1,
                           upperFilterQuantile = 0.99)

  # Drop rows so the matrix no longer lines up with the loaded references.
  expect_error(collapseChrom3N(norm[1:5, ]), "matched by position")
})

test_that("identifyCNVClusters requires normalCells when centring on them", {
  refs <- make_test_references(binsPerChrom = 10L)
  session <- local_copyscat_session(refs)
  fake <- data.frame(chrom = c("chr1p", "chr1q"),
                     `CELL001-1` = c(1, 2), check.names = FALSE)

  expect_error(
    identifyCNVClusters(fake, list(), medianQuantileCutoff = -1,
                        normalCells = NULL),
    "normalCells"
  )
})

test_that("identifyCNVClusters rejects normalCells absent from the matrix", {
  refs <- make_test_references(binsPerChrom = 10L)
  session <- local_copyscat_session(refs)
  fake <- data.frame(chrom = c("chr1p", "chr1q"),
                     `CELL001-1` = c(1, 2), check.names = FALSE)

  expect_error(
    identifyCNVClusters(fake, list(), medianQuantileCutoff = -1,
                        normalCells = c("NOT_A_CELL-1")),
    "not columns of"
  )
})

test_that("annotateCNV4B rejects an empty or unmatched normal set", {
  session <- local_copyscat_session()
  cnvResults <- list(
    data.frame(Cells = c("CELL001-1", "CELL002-1"), chr1p = c(1, 2)),
    data.frame(Chrom = "chr1p", V1 = 0, V2 = 1)
  )
  expect_error(annotateCNV4B(cnvResults, character(0), saveOutput = FALSE),
               "empty")
  expect_error(annotateCNV4B(cnvResults, "NOBODY-1", saveOutput = FALSE),
               "None of")
})

test_that("smoothClusters validates its inputs", {
  session <- local_copyscat_session()
  expect_error(
    smoothClusters(data.frame(Barcode = "a"), inputCNVList = data.frame()),
    "at least 2 columns"
  )
  expect_error(
    smoothClusters(data.frame(Barcode = "a", clust = 1),
                   inputCNVList = data.frame(),
                   inputCNVClusterFile = file.path(tempdir(), "absent.csv")),
    "does not exist"
  )
})

test_that("getAlteredSegments rejects a malformed clusterResults argument", {
  session <- local_copyscat_session()
  expect_error(
    getAlteredSegments(data.frame(chrom = "chr1", pos = 1),
                       clusterResults = list(wrong = 1)),
    "identifyNonNeoplastic"
  )
})
