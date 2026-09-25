# Shared fixtures for the CopyscAT test suite.

#' Build a small, self-consistent set of reference files in a temp directory.
#'
#' Returns the three paths plus the bin table they were generated from, so
#' tests can construct matrices whose row order matches the references.
make_test_references <- function(dir = tempfile("copyscat_refs"),
                                 chroms = c("chr1", "chr2", "chrX", "chrY"),
                                 binsPerChrom = 10L,
                                 binSize = 1e6) {
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)

  bins <- do.call(rbind, lapply(chroms, function(cc) {
    data.frame(
      chrom = cc,
      start = seq(0, by = binSize, length.out = binsPerChrom),
      end   = seq(binSize, by = binSize, length.out = binsPerChrom),
      stringsAsFactors = FALSE
    )
  }))

  # Arm labels: first half p, centromere in the middle, rest q.
  arm <- rep("p", nrow(bins))
  for (cc in chroms) {
    idx <- which(bins$chrom == cc)
    mid <- idx[ceiling(length(idx) / 2)]
    arm[idx[idx < mid]] <- "p"
    arm[mid] <- "cen"
    arm[idx[idx > mid]] <- "q"
  }
  bins$arm <- arm

  # Deterministic CpG densities, comfortably above the usual minCPG cutoffs.
  bins$cpg <- 400L + (seq_len(nrow(bins)) %% 7L) * 25L

  chromSizeFile <- file.path(dir, "test_chrom_sizes.tsv")
  cpgFile       <- file.path(dir, "test_cpg_densities.tsv")
  cytobandFile  <- file.path(dir, "test_cytoband_densities.tsv")

  utils::write.table(
    data.frame(chrom = chroms, length = binsPerChrom * binSize),
    chromSizeFile, sep = "\t", quote = FALSE,
    col.names = FALSE, row.names = FALSE
  )
  utils::write.table(
    bins[, c("chrom", "start", "end", "cpg")],
    cpgFile, sep = "\t", quote = FALSE,
    col.names = FALSE, row.names = FALSE
  )
  utils::write.table(
    bins[, c("chrom", "start", "end", "arm")],
    cytobandFile, sep = "\t", quote = FALSE,
    col.names = FALSE, row.names = FALSE
  )

  list(
    dir = dir, bins = bins, binSize = binSize,
    genomeFile = chromSizeFile, cpgFile = cpgFile,
    cytobandFile = cytobandFile
  )
}

#' Initialise a full CopyscAT session against the temp references.
local_copyscat_session <- function(refs = make_test_references(),
                                   outputDir = tempfile("copyscat_out"),
                                   cellSuffix = "-1",
                                   env = parent.frame()) {
  dir.create(outputDir, recursive = TRUE, showWarnings = FALSE)
  initialiseEnvironment(
    genomeFile   = refs$genomeFile,
    cytobandFile = refs$cytobandFile,
    cpgFile      = refs$cpgFile,
    binSize      = refs$binSize,
    minFrags     = 100,
    cellSuffix   = cellSuffix,
    lowerTrim    = 0.5,
    upperTrim    = 0.8
  )
  setOutputFile(outputDir, "testsample")
  withr::defer(reset_copyscat_session(), envir = env)
  c(refs, list(outputDir = outputDir))
}

#' Empty the session environment between tests.
reset_copyscat_session <- function() {
  rm(list = ls(envir = CopyscAT::scCNVCaller, all.names = TRUE),
     envir = CopyscAT::scCNVCaller)
}

#' Simulate a raw count matrix (cells x bins) with optional CNVs.
#'
#' @param refs Output of make_test_references().
#' @param nCells Number of cells.
#' @param baseCount Mean counts per bin in a diploid cell.
#' @param gains Named list mapping chromosome -> multiplier, applied to the
#'   first `propAltered` of cells.
#' @param propAltered Fraction of cells carrying `gains`.
simulate_counts <- function(refs, nCells = 60L, baseCount = 300,
                            gains = list(), propAltered = 0.5,
                            cellSuffix = "-1", seed = 42) {
  set.seed(seed)
  bins <- refs$bins
  binNames <- paste(bins$chrom, bins$start, sep = "_")
  cellNames <- paste0("CELL", sprintf("%03d", seq_len(nCells)), cellSuffix)

  mat <- matrix(
    stats::rpois(nCells * nrow(bins), lambda = baseCount),
    nrow = nCells, ncol = nrow(bins),
    dimnames = list(cellNames, binNames)
  )

  if (length(gains) > 0L) {
    altered <- seq_len(max(1L, floor(nCells * propAltered)))
    for (cc in names(gains)) {
      cols <- which(bins$chrom == cc)
      mat[altered, cols] <- round(mat[altered, cols] * gains[[cc]])
    }
  }
  as.data.frame(mat)
}
