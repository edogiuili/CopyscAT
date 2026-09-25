# End-to-end tests: run the documented workflow on simulated data whose
# copy number truth is known, and assert the pipeline recovers it.

test_that("the normalise -> collapse -> filter pipeline recovers a known gain", {
  refs <- make_test_references(
    chroms = c("chr1", "chr2", "chr3", "chrX", "chrY"), binsPerChrom = 20L
  )
  session <- local_copyscat_session(refs)

  # Half the cells carry a 2x gain on chr2; everything else is diploid.
  counts <- simulate_counts(
    refs, nCells = 80L, baseCount = 400,
    gains = list(chr2 = 2), propAltered = 0.5
  )

  norm <- normalizeMatrixN(
    counts, logNorm = FALSE, maxZero = 2000, imputeZeros = FALSE,
    blacklistProp = 0.8, blacklistCutoff = 10, dividingFactor = 1,
    upperFilterQuantile = 0.99
  )

  expect_true(all(c("chrom", "pos", "blacklist") %in% colnames(norm)))
  expect_equal(nrow(norm), nrow(refs$bins))
  # `pos` must survive as a full coordinate, not a truncated fragment.
  expect_false(any(is.na(as.numeric(norm$pos))))

  collapsed <- collapseChrom3N(
    norm, summaryFunction = cutAverage, binExpand = 1,
    minimumChromValue = 0, logTrans = FALSE, tssEnrich = 1,
    logBase = 2, minCPG = 100, powVal = 0.73
  )

  # Centromeric bins are dropped, so each chromosome yields a p and a q arm.
  expect_true(all(c("chr1p", "chr1q", "chr2p", "chr2q") %in% collapsed$chrom))

  cellCols <- setdiff(colnames(collapsed), c("chrom", "medianNorm"))
  alteredCells <- grep("^CELL0(0[1-9]|[1-3][0-9]|40)-1$", cellCols, value = TRUE)
  alteredCells <- intersect(alteredCells, cellCols)

  chr2p <- unlist(collapsed[collapsed$chrom == "chr2p", alteredCells])
  chr1p <- unlist(collapsed[collapsed$chrom == "chr1p", alteredCells])

  # The gained arm must carry visibly more signal in the altered cells.
  expect_gt(mean(chr2p), 1.5 * mean(chr1p))
})

test_that("filterCells removes cells and reports honest counts", {
  refs <- make_test_references(binsPerChrom = 20L)
  session <- local_copyscat_session(refs)

  counts <- simulate_counts(refs, nCells = 60L, baseCount = 400)
  norm <- normalizeMatrixN(
    counts, blacklistCutoff = 10, dividingFactor = 1,
    upperFilterQuantile = 0.99
  )
  collapsed <- collapseChrom3N(
    norm, minimumChromValue = 0, tssEnrich = 1, minCPG = 100
  )

  filtered <- filterCells(collapsed, minimumSegments = 1, minDensity = 0,
                          signalSDcut = 3)

  expect_true(ncol(filtered) <= ncol(collapsed))
  expect_true("chrom" %in% colnames(filtered))
  # finalFilterCells must match what actually survived, not an earlier stage.
  expect_equal(scCNVCaller$finalFilterCells, ncol(filtered) - 1L)
})

test_that("computeCenters and scaleMatrix produce centred Z-scores", {
  refs <- make_test_references(binsPerChrom = 20L)
  session <- local_copyscat_session(refs)

  counts <- simulate_counts(refs, nCells = 60L, baseCount = 400)
  norm <- normalizeMatrixN(counts, blacklistCutoff = 10, dividingFactor = 1,
                           upperFilterQuantile = 0.99)
  collapsed <- collapseChrom3N(norm, minimumChromValue = 0, tssEnrich = 1,
                               minCPG = 100)
  filtered <- filterCells(collapsed, minimumSegments = 1, minDensity = 0,
                          signalSDcut = 3)

  centers <- computeCenters(filtered, summaryFunction = cutAverage)
  expect_length(centers, 2L)
  expect_named(centers[[1]], c("chrom", "Value"))
  expect_named(centers[[2]], c("chrom", "Value"))
  expect_equal(nrow(centers[[1]]), nrow(centers[[2]]))

  scaled <- scaleMatrix(filtered, centers, spread = FALSE)
  expect_true("chrom" %in% colnames(scaled))
  # A diploid simulation should sit near zero once scaled.
  numeric_vals <- unlist(scaled[, setdiff(colnames(scaled), "chrom")])
  expect_lt(abs(stats::median(numeric_vals, na.rm = TRUE)), 2)

  longForm <- scaleMatrix(filtered, centers, spread = TRUE)
  expect_true(all(c("Cell", "Density") %in% colnames(longForm)))
})

test_that("saveSummaryStats reflects the counts recorded during the run", {
  refs <- make_test_references(binsPerChrom = 20L)
  session <- local_copyscat_session(refs)

  counts <- simulate_counts(refs, nCells = 50L, baseCount = 400)
  norm <- normalizeMatrixN(counts, blacklistCutoff = 10, dividingFactor = 1,
                           upperFilterQuantile = 0.99)
  collapsed <- collapseChrom3N(norm, minimumChromValue = 0, tssEnrich = 1,
                               minCPG = 100)
  filterCells(collapsed, minimumSegments = 1, minDensity = 0, signalSDcut = 3)

  stats <- saveSummaryStats("_stats")
  expect_equal(stats$startingCellCount, 50)
  expect_gt(stats$cellsPassingFilter, 0)
  expect_gt(stats$meanReadsPerCell, 0)
  expect_true(file.exists(file.path(session$outputDir, "testsample_stats.tsv")))
})

test_that("graphCNVDistribution writes a PDF", {
  refs <- make_test_references(binsPerChrom = 20L)
  session <- local_copyscat_session(refs)

  counts <- simulate_counts(refs, nCells = 40L, baseCount = 400)
  norm <- normalizeMatrixN(counts, blacklistCutoff = 10, dividingFactor = 1,
                           upperFilterQuantile = 0.99)
  collapsed <- collapseChrom3N(norm, minimumChromValue = 0, tssEnrich = 1,
                               minCPG = 100)

  graphCNVDistribution(collapsed, outputSuffix = "_qc")
  expect_true(file.exists(
    file.path(session$outputDir, "testsample_qc_graph.pdf")
  ))
})
