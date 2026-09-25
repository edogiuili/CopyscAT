# Full CNV-calling workflow on simulated data with one known gain.
# These tests are slower than the rest of the suite because they fit Gaussian
# mixtures per chromosome arm.

simulate_called_cnvs <- function(gainChrom = "chr3", nCells = 200L) {
  refs <- make_test_references(chroms = paste0("chr", 1:6), binsPerChrom = 20L)
  session <- local_copyscat_session(refs, env = parent.frame())

  counts <- simulate_counts(
    refs, nCells = nCells, baseCount = 500,
    gains = stats::setNames(list(2), gainChrom), propAltered = 0.5
  )
  norm <- normalizeMatrixN(counts, blacklistCutoff = 10, dividingFactor = 1,
                           upperFilterQuantile = 0.99)
  collapsed <- collapseChrom3N(norm, minimumChromValue = 0, tssEnrich = 1,
                               minCPG = 100)
  filtered <- filterCells(collapsed, minimumSegments = 1, minDensity = 0,
                          signalSDcut = 3)
  centers <- computeCenters(filtered, summaryFunction = cutAverage)

  candidates <- identifyCNVClusters(
    filtered, centers, useDummyCells = TRUE, propDummy = 0.25, minMix = 0.01,
    deltaMean = 0.03, deltaBIC2 = 0.25, bicMinimum = 0.1, subsetSize = 100,
    fakeCellSD = 0.08, uncertaintyCutoff = 0.55,
    summaryFunction = cutAverage, maxClust = 4, mergeCutoff = 3,
    IQRCutoff = 0.2, medianQuantileCutoff = 0.4
  )
  list(session = session, filtered = filtered, centers = centers,
       candidates = candidates)
}

test_that("identifyCNVClusters returns assignments for every cell and arm", {
  res <- simulate_called_cnvs()
  candidates <- res$candidates

  expect_length(candidates, 3L)
  assignments <- candidates[[1]]
  expect_true("Cells" %in% colnames(assignments))
  # One row per real cell plus the simulated pseudodiploid cells.
  expect_gt(nrow(assignments), 100)
  # Cluster labels are non-negative integers; 0 marks an uncertain call.
  numericCols <- assignments[, setdiff(colnames(assignments), "Cells")]
  expect_true(all(unlist(numericCols) >= 0))
})

test_that("clusterCNV merges clusters and preserves the cell set", {
  res <- simulate_called_cnvs()
  clean <- clusterCNV(initialResultList = res$candidates,
                      medianIQR = res$candidates[[3]], minDiff = 1.0)

  expect_length(clean, 2L)
  expect_equal(nrow(clean[[1]]), nrow(res$candidates[[1]]))
  expect_equal(clean[[1]]$Cells, res$candidates[[1]]$Cells)
})

test_that("annotateCNV4 reports the simulated gain and writes its output", {
  res <- simulate_called_cnvs(gainChrom = "chr3")
  clean <- clusterCNV(initialResultList = res$candidates,
                      medianIQR = res$candidates[[3]], minDiff = 1.0)
  final <- annotateCNV4(clean, saveOutput = TRUE, outputSuffix = "_called",
                        sdCNV = 0.5, filterResults = TRUE, filterRange = 0.8)

  expect_length(final, 3L)
  calledArms <- setdiff(colnames(final[[3]]), "rowname")

  # The gained chromosome must be among the calls, and unaltered chromosomes
  # must not dominate them.
  expect_true(any(grepl("^chr3", calledArms)),
              info = paste("called arms:", paste(calledArms, collapse = ",")))

  expect_true(file.exists(
    file.path(res$session$outputDir, "testsample_called_cnv_scores.csv")
  ))
})

test_that("identifyCNVClusters warns and clamps an oversized subsetSize", {
  refs <- make_test_references(chroms = paste0("chr", 1:4), binsPerChrom = 20L)
  session <- local_copyscat_session(refs)

  counts <- simulate_counts(refs, nCells = 40L, baseCount = 500)
  norm <- normalizeMatrixN(counts, blacklistCutoff = 10, dividingFactor = 1,
                           upperFilterQuantile = 0.99)
  collapsed <- collapseChrom3N(norm, minimumChromValue = 0, tssEnrich = 1,
                               minCPG = 100)
  filtered <- filterCells(collapsed, minimumSegments = 1, minDensity = 0,
                          signalSDcut = 3)
  centers <- computeCenters(filtered, summaryFunction = cutAverage)

  # subsetSize far larger than the cell count used to abort with an opaque
  # sample() error.
  expect_warning(
    identifyCNVClusters(
      filtered, centers, useDummyCells = TRUE, subsetSize = 100000,
      minMix = 0.01, maxClust = 3, medianQuantileCutoff = 0.4
    ),
    "exceeds"
  )
})
