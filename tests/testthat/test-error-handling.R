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

test_that("normalizeMatrixN names empty barcodes instead of failing inside cpm()", {
  refs <- make_test_references(chroms = paste0("chr", 1:4), binsPerChrom = 20L)
  session <- local_copyscat_session(refs)

  counts <- simulate_counts(refs, nCells = 40L, baseCount = 400)
  counts["DEADCELL-1", ] <- 0

  # maxZero at or above the bin count disables the empty-bin filter, so the
  # barcode reaches edgeR::cpm(), whose own message ("library sizes should be
  # greater than zero") names neither the cell nor the parameter to change.
  err <- tryCatch(
    normalizeMatrixN(counts, maxZero = 10000, blacklistProp = 0.9,
                     blacklistCutoff = 125, dividingFactor = 1,
                     upperFilterQuantile = 0.9),
    error = conditionMessage
  )

  expect_match(err, "DEADCELL-1", fixed = TRUE)
  expect_match(err, "maxZero", fixed = TRUE)
  expect_match(err, "filter is disabled", fixed = TRUE)
  expect_false(grepl("library sizes should be greater than zero", err,
                     fixed = TRUE))
})

test_that("an active empty-bin filter drops zero-signal barcodes silently", {
  refs <- make_test_references(chroms = paste0("chr", 1:4), binsPerChrom = 20L)
  session <- local_copyscat_session(refs)

  counts <- simulate_counts(refs, nCells = 40L, baseCount = 400)
  counts["DEADCELL-1", ] <- 0

  # maxZero below the bin count keeps the filter working, so normalisation
  # completes and the empty barcode is simply absent from the result.
  norm <- normalizeMatrixN(counts, maxZero = 70, blacklistProp = 0.9,
                           blacklistCutoff = 125, dividingFactor = 1,
                           upperFilterQuantile = 0.9)
  expect_false("DEADCELL-1" %in% colnames(norm))
})

test_that("annotateCNV4B collapses the normal-cluster index to a vector", {
  # t() on the one-row normal_clusters frame yields an n-by-1 matrix. Older
  # dplyr coerced it inside mutate(); newer dplyr keeps it a matrix and
  # if_else() then fails with "`condition` must be a logical vector, not a
  # logical matrix". Checking the shape of the result keeps this test
  # independent of the installed dplyr version.
  session <- local_copyscat_session()

  # clusterCNV() returns Chrom plus V1..V6, which annotateCNV4B relies on when
  # it maps cluster indices onto factor levels 0..6.
  arms <- c("chr1p", "chr1q")
  cnvResults <- list(
    data.frame(
      Cells = c("N1-1", "N2-1", "T1-1", "T2-1"),
      chr1p = c(1, 1, 2, 2),
      chr1q = c(1, 1, 2, 2),
      stringsAsFactors = FALSE
    ),
    data.frame(
      Chrom = arms,
      V1 = c(0.10, 0.20), V2 = c(1.10, 1.20), V3 = c(0, 0),
      V4 = c(0, 0), V5 = c(0, 0), V6 = c(0, 0),
      stringsAsFactors = FALSE
    )
  )

  res <- annotateCNV4B(cnvResults, expectedNormals = c("N1-1", "N2-1"),
                       saveOutput = FALSE, filterResults = FALSE)

  expect_length(res, 3L)
  # zoffset must be a plain numeric vector, one value per arm, not a matrix.
  zoff <- res[[1]]$zoffset
  expect_null(dim(zoff))
  expect_type(zoff, "double")
})

test_that("annotateCNV4B rejects cluster results of the wrong shape", {
  session <- local_copyscat_session()

  # Three arms in the cluster means but only two columns of assignments.
  cnvResults <- list(
    data.frame(Cells = c("N1-1", "N2-1"), chr1p = c(1, 1), chr1q = c(1, 1),
               stringsAsFactors = FALSE),
    data.frame(Chrom = c("chr1p", "chr1q", "chr2p"),
               V1 = c(0.1, 0.2, 0.3), V2 = c(1.1, 1.2, 1.3),
               V3 = 0, V4 = 0, V5 = 0, V6 = 0,
               stringsAsFactors = FALSE)
  )

  expect_error(
    annotateCNV4B(cnvResults, expectedNormals = c("N1-1", "N2-1"),
                  saveOutput = FALSE),
    "one normal-cluster index per chromosome arm"
  )
})
