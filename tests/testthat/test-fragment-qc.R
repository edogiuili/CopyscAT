make_demo_matrix <- function(nCells = 20L, nBins = 40L, seed = 3) {
  set.seed(seed)
  m <- matrix(
    stats::rpois(nCells * nBins, lambda = 500),
    nrow = nCells,
    dimnames = list(
      paste0("CELL", sprintf("%03d", seq_len(nCells)), "-1"),
      paste0("chr1_", seq_len(nBins) * 1e6)
    )
  )
  as.data.frame(m)
}

test_that("signalPerCell summarises an in-memory matrix", {
  m <- make_demo_matrix()
  res <- signalPerCell(list(demo = m))

  expect_s3_class(res, "data.frame")
  expect_equal(nrow(res), nrow(m))
  expect_named(res, c("sample", "barcode", "signal_bp", "nonzero_bins"))
  expect_equal(res$signal_bp, unname(rowSums(m)))
  expect_true(all(res$sample == "demo"))
})

test_that("signalPerCell combines several samples and labels them", {
  m <- make_demo_matrix(nCells = 10L)
  res <- signalPerCell(list(first = m, second = m[1:4, ]))

  expect_equal(nrow(res), 14L)
  expect_equal(sort(unique(res$sample)), c("first", "second"))
  expect_equal(sum(res$sample == "second"), 4L)
})

test_that("signalPerCell reads matrices from disk", {
  m <- make_demo_matrix(nCells = 8L)
  path <- withr::local_tempfile(fileext = ".tsv")
  utils::write.table(cbind(Cell_id = rownames(m), m), path, sep = "\t",
                     quote = FALSE, row.names = FALSE)

  res <- signalPerCell(c(onDisk = path))
  expect_equal(nrow(res), 8L)
  expect_equal(res$signal_bp, unname(rowSums(m)))
})

test_that("signalPerCell rejects unnamed, empty or missing input", {
  m <- make_demo_matrix(nCells = 4L)
  expect_error(signalPerCell(list()), "empty")
  expect_error(signalPerCell(list(m)), "must be named")
  expect_error(signalPerCell(c(gone = "no_such_file_xyz.tsv")),
               "does not exist")
})

test_that("fragmentsPerCell counts barcodes and skips comment headers", {
  path <- withr::local_tempfile(fileext = ".tsv")
  writeLines("# cellranger-atac style header", path)
  frags <- data.frame(
    chrom = "chr1",
    start = seq(1, 900, by = 100),
    end = seq(101, 1000, by = 100),
    barcode = rep(c("AAA-1", "BBB-1", "CCC-1"), each = 3),
    n = 1
  )
  utils::write.table(frags, path, sep = "\t", quote = FALSE, row.names = FALSE,
                     col.names = FALSE, append = TRUE)

  res <- fragmentsPerCell(c(demo = path))
  expect_named(res, c("sample", "barcode", "n_fragments"))
  expect_equal(nrow(res), 3L)
  expect_true(all(res$n_fragments == 3))
  # the '#' line must not become a barcode
  expect_false(any(startsWith(res$barcode, "#")))
})

test_that("fragmentsPerCell applies minFragments", {
  path <- withr::local_tempfile(fileext = ".tsv")
  frags <- data.frame(
    chrom = "chr1", start = 1:10, end = 11:20,
    barcode = c(rep("DEEP-1", 8), rep("SHALLOW-1", 2)), n = 1
  )
  utils::write.table(frags, path, sep = "\t", quote = FALSE, row.names = FALSE,
                     col.names = FALSE)

  res <- fragmentsPerCell(c(demo = path), minFragments = 5)
  expect_equal(res$barcode, "DEEP-1")
})

test_that("fragmentsPerCell reads a CellRanger singlecell.csv", {
  path <- withr::local_tempfile(fileext = ".csv")
  utils::write.csv(
    data.frame(barcode = c("AAA-1", "BBB-1"), passed_filters = c(5000, 12000)),
    path, row.names = FALSE
  )

  res <- fragmentsPerCell(c(cr = path))
  expect_equal(nrow(res), 2L)
  expect_equal(res$n_fragments, c(5000, 12000))

  filtered <- fragmentsPerCell(c(cr = path), minFragments = 6000)
  expect_equal(filtered$barcode, "BBB-1")
})

test_that("fragmentsPerCell reports a singlecell.csv missing its columns", {
  path <- withr::local_tempfile(fileext = ".csv")
  utils::write.csv(data.frame(cell = "AAA-1", total = 10), path,
                   row.names = FALSE)
  expect_error(fragmentsPerCell(c(cr = path)), "passed_filters")
})

test_that("summariseCellDepth reports one row per sample", {
  m <- make_demo_matrix(nCells = 12L)
  res <- summariseCellDepth(signalPerCell(list(a = m, b = m[1:5, ])))

  expect_equal(nrow(res), 2L)
  expect_named(res, c("sample", "cells", "min", "q25", "median", "q75", "max"))
  expect_equal(res$cells[res$sample == "a"], 12L)
  expect_equal(res$cells[res$sample == "b"], 5L)
  expect_true(all(res$min <= res$median & res$median <= res$max))
})

test_that("the depth column is detected, or named explicitly", {
  m <- make_demo_matrix(nCells = 6L)
  sig <- signalPerCell(list(a = m))

  # signal_bp is found without being named
  expect_equal(summariseCellDepth(sig)$cells, 6L)
  # an explicit column also works
  expect_equal(summariseCellDepth(sig, "nonzero_bins")$cells, 6L)
  # and an unknown one is reported
  expect_error(summariseCellDepth(sig, "not_a_column"), "not a column")
})

test_that("plotting functions return ggplot objects", {
  m <- make_demo_matrix(nCells = 15L)
  sig <- signalPerCell(list(a = m, b = m[1:6, ]))

  expect_s3_class(plotCellDistribution(sig), "ggplot")
  expect_s3_class(plotCellDistribution(sig, cutoff = 1e4), "ggplot")
  expect_s3_class(plotCellKnee(sig), "ggplot")
  expect_s3_class(plotCellKnee(sig, cutoff = 1e4), "ggplot")
})

test_that("plotCellDistribution warns when the log scale drops cells", {
  m <- make_demo_matrix(nCells = 5L)
  sig <- signalPerCell(list(a = m))
  sig$signal_bp[1] <- 0

  expect_warning(plotCellDistribution(sig), "dropped for the log scale")
  expect_silent(plotCellDistribution(sig, logScale = FALSE))
})
