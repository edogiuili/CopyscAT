test_that("initialiseEnvironment populates the session environment", {
  refs <- make_test_references()
  session <- local_copyscat_session(refs)

  expect_equal(scCNVCaller$binSize, refs$binSize)
  expect_equal(scCNVCaller$cellSuffix, "-1")
  expect_equal(nrow(scCNVCaller$cpg_data), nrow(refs$bins))
  expect_equal(nrow(scCNVCaller$cytoband_data), nrow(refs$bins))
  expect_named(scCNVCaller$chrom_sizes, c("chrom", "length"))
  expect_named(scCNVCaller$cpg_data,
               c("chrom", "start", "end", "cpg_density"))
})

test_that("initialiseEnvironment rejects missing reference files", {
  refs <- make_test_references()
  expect_error(
    initialiseEnvironment(
      genomeFile = file.path(tempdir(), "definitely_absent.tsv"),
      cytobandFile = refs$cytobandFile, cpgFile = refs$cpgFile
    ),
    "does not exist"
  )
  reset_copyscat_session()
})

test_that("initialiseEnvironment rejects invalid trim quantiles", {
  refs <- make_test_references()
  expect_error(
    initialiseEnvironment(
      genomeFile = refs$genomeFile, cytobandFile = refs$cytobandFile,
      cpgFile = refs$cpgFile, lowerTrim = 0.9, upperTrim = 0.2
    ),
    "lowerTrim"
  )
  reset_copyscat_session()
})

test_that("initialiseEnvironment warns when references disagree on bin count", {
  refs <- make_test_references()
  shortCpg <- file.path(refs$dir, "short_cpg.tsv")
  utils::write.table(
    refs$bins[seq_len(5), c("chrom", "start", "end", "cpg")],
    shortCpg, sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE
  )
  expect_warning(
    initialiseEnvironment(
      genomeFile = refs$genomeFile, cytobandFile = refs$cytobandFile,
      cpgFile = shortCpg
    ),
    "same bins"
  )
  reset_copyscat_session()
})

test_that("setOutputFile appends a trailing slash and rejects missing dirs", {
  refs <- make_test_references()
  out <- tempfile("out_no_slash"); dir.create(out)
  initialiseEnvironment(refs$genomeFile, refs$cytobandFile, refs$cpgFile)
  setOutputFile(out, "sample")
  expect_match(scCNVCaller$locPrefix, "/$")
  expect_equal(scCNVCaller$outPrefix, "sample")

  expect_error(
    setOutputFile(file.path(tempdir(), "no_such_dir_xyz"), "s"),
    "does not exist"
  )
  reset_copyscat_session()
})

test_that("functions refuse to run before the session is initialised", {
  reset_copyscat_session()
  expect_error(saveSummaryStats(), "not initialised")
})

test_that("readInputTable round-trips a matrix and honours `sep`", {
  refs <- make_test_references()
  session <- local_copyscat_session(refs)
  counts <- simulate_counts(refs, nCells = 5L)

  tsv <- tempfile(fileext = ".tsv")
  utils::write.table(
    cbind(Cell_id = rownames(counts), counts),
    tsv, sep = "\t", quote = FALSE, row.names = FALSE
  )
  back <- readInputTable(tsv)
  expect_equal(nrow(back), nrow(counts))
  expect_equal(ncol(back), ncol(counts))
  expect_equal(rownames(back), rownames(counts))

  # A comma-delimited file must be readable via `sep`, which the original
  # implementation ignored.
  csv <- tempfile(fileext = ".csv")
  utils::write.table(
    cbind(Cell_id = rownames(counts), counts),
    csv, sep = ",", quote = FALSE, row.names = FALSE
  )
  backCsv <- readInputTable(csv, sep = ",")
  expect_equal(dim(backCsv), dim(counts))
})

test_that("readInputTable errors on a missing file", {
  expect_error(readInputTable(file.path(tempdir(), "nope.tsv")),
               "does not exist")
})

test_that("saveSummaryStats writes a table and returns the statistics", {
  refs <- make_test_references()
  session <- local_copyscat_session(refs)
  scCNVCaller$meanReadsPerCell <- 1234

  stats <- saveSummaryStats("_stats")
  expect_equal(stats$meanReadsPerCell, 1234)
  expect_true(file.exists(
    file.path(session$outputDir, "testsample_stats.tsv")
  ))
})
