#' Load genome references and analysis parameters for a session
#'
#' Reads the three reference tables produced by [generateReferences()] and
#' stores them, together with the binning parameters, in the [scCNVCaller]
#' session environment. Run this once per R session before any other step.
#'
#' @param genomeFile Path to a two-column chromosome-size table
#'   (chromosome, length), tab separated and without a header.
#' @param cytobandFile Path to the binned cytoband table produced by
#'   [generateReferences()].
#' @param cpgFile Path to the binned CpG density table produced by
#'   [generateReferences()]. Four columns: chromosome, start, end, density.
#' @param binSize Genomic bin size in base pairs. Must match the `tileWidth`
#'   used to build the reference files. Default `1e6`.
#' @param minFrags Minimum number of fragments for a cell to be considered.
#' @param cellSuffix Character vector of barcode suffixes present in the input
#'   matrix, for example `c("-1", "-2")`.
#' @param lowerTrim Lower quantile for the trimmed mean used by [cutAverage()].
#' @param upperTrim Upper quantile for the trimmed mean used by [cutAverage()].
#'
#' @return Invisibly `TRUE`. Called for the side effect of populating
#'   [scCNVCaller].
#' @family session setup
#' @export
#' @examples
#' \dontrun{
#' initialiseEnvironment(
#'   genomeFile   = "hg38_chrom_sizes.tsv",
#'   cytobandFile = "hg38_1e+06_cytoband_densities_granges.tsv",
#'   cpgFile      = "hg38_1e+06_cpg_densities.tsv",
#'   binSize      = 1e6,
#'   cellSuffix   = c("-1")
#' )
#' }
initialiseEnvironment <- function(genomeFile,
                                  cytobandFile,
                                  cpgFile,
                                  binSize = 1e6,
                                  minFrags = 1e4,
                                  cellSuffix = c("-1"),
                                  lowerTrim = 0.5,
                                  upperTrim = 0.8) {
  files <- c(genomeFile = genomeFile, cytobandFile = cytobandFile,
             cpgFile = cpgFile)
  for (nm in names(files)) {
    if (!is.character(files[[nm]]) || length(files[[nm]]) != 1L) {
      stop(sprintf("`%s` must be a single file path.", nm), call. = FALSE)
    }
    if (!file.exists(files[[nm]])) {
      stop(sprintf("`%s` does not exist: %s", nm, files[[nm]]), call. = FALSE)
    }
  }
  if (!is.numeric(binSize) || length(binSize) != 1L || binSize <= 0) {
    stop("`binSize` must be a single positive number.", call. = FALSE)
  }
  if (!is.character(cellSuffix) || length(cellSuffix) < 1L) {
    stop("`cellSuffix` must be a non-empty character vector.", call. = FALSE)
  }
  if (!is.numeric(lowerTrim) || !is.numeric(upperTrim) ||
      lowerTrim < 0 || upperTrim > 1 || lowerTrim >= upperTrim) {
    stop("`lowerTrim` and `upperTrim` must satisfy 0 <= lowerTrim < upperTrim <= 1.",
         call. = FALSE)
  }

  scCNVCaller$genomeFile <- genomeFile
  scCNVCaller$cytobandFile <- cytobandFile
  scCNVCaller$cpgFile <- cpgFile

  scCNVCaller$chrom_sizes <- read_chrom_sizes(genomeFile)
  scCNVCaller$cpg_data <- read_cpg_data(cpgFile)
  scCNVCaller$cytoband_data <- read_cytoband_data(cytobandFile)

  if (nrow(scCNVCaller$cpg_data) != nrow(scCNVCaller$cytoband_data)) {
    warning(
      sprintf(
        paste0("CpG file has %d bins but the cytoband file has %d. Both are ",
               "indexed positionally against the input matrix, so they must ",
               "describe the same bins at the same bin size."),
        nrow(scCNVCaller$cpg_data), nrow(scCNVCaller$cytoband_data)
      ),
      call. = FALSE
    )
  }

  scCNVCaller$binSize <- binSize
  scCNVCaller$minFrags <- minFrags
  scCNVCaller$cellSuffix <- cellSuffix
  scCNVCaller$lowerTrim <- lowerTrim
  scCNVCaller$upperTrim <- upperTrim

  scCNVCaller$meanReadsPerCell <- 0
  scCNVCaller$startingCellCount <- 0
  scCNVCaller$cellsPassingFilter <- 0
  scCNVCaller$blacklistCount <- 0
  scCNVCaller$finalFilterCells <- 0
  scCNVCaller$isMale <- FALSE

  invisible(TRUE)
}

#' Read a chromosome-size table
#' @param path File path.
#' @return A data frame with columns `chrom` and `length`.
#' @noRd
read_chrom_sizes <- function(path) {
  tbl <- utils::read.delim(path, stringsAsFactors = FALSE, header = FALSE)
  if (ncol(tbl) < 2L) {
    stop(sprintf("Chromosome size file needs 2 columns, found %d: %s",
                 ncol(tbl), path), call. = FALSE)
  }
  tbl <- tbl[, 1:2, drop = FALSE]
  colnames(tbl) <- c("chrom", "length")
  tbl
}

#' Read a binned CpG density table
#' @param path File path.
#' @return A data frame with columns `chrom`, `start`, `end`, `cpg_density`.
#' @noRd
read_cpg_data <- function(path) {
  tbl <- utils::read.table(path, stringsAsFactors = FALSE, header = FALSE)
  if (ncol(tbl) < 4L) {
    stop(sprintf("CpG density file needs 4 columns, found %d: %s",
                 ncol(tbl), path), call. = FALSE)
  }
  tbl <- tbl[, 1:4, drop = FALSE]
  colnames(tbl) <- c("chrom", "start", "end", "cpg_density")
  tbl
}

#' Read a binned cytoband table
#' @param path File path.
#' @return A data frame whose 4th column holds the chromosome arm label.
#' @noRd
read_cytoband_data <- function(path) {
  tbl <- utils::read.table(path, stringsAsFactors = FALSE, header = FALSE)
  if (ncol(tbl) < 4L) {
    stop(sprintf("Cytoband file needs at least 4 columns, found %d: %s",
                 ncol(tbl), path), call. = FALSE)
  }
  tbl
}

#' Set the output directory and file prefix
#'
#' All functions that write plots or tables prepend `outputDir` and
#' `outputFileSuffix` to their file names.
#'
#' @param outputDir Directory for output files. A trailing slash is added when
#'   missing. The directory must already exist.
#' @param outputFileSuffix Prefix used for every file written in this session,
#'   typically the sample name.
#'
#' @return Invisibly `TRUE`.
#' @family session setup
#' @export
#' @examples
#' setOutputFile(tempdir(), "demo_sample")
setOutputFile <- function(outputDir, outputFileSuffix) {
  if (!is.character(outputDir) || length(outputDir) != 1L) {
    stop("`outputDir` must be a single directory path.", call. = FALSE)
  }
  if (!is.character(outputFileSuffix) || length(outputFileSuffix) != 1L) {
    stop("`outputFileSuffix` must be a single string.", call. = FALSE)
  }
  if (!dir.exists(outputDir)) {
    stop(sprintf("`outputDir` does not exist: %s", outputDir), call. = FALSE)
  }

  prefix <- outputDir
  if (str_ends(prefix, "/", negate = TRUE)) {
    prefix <- str_c(prefix, "/", sep = "")
  }
  scCNVCaller$locPrefix <- prefix
  scCNVCaller$outPrefix <- outputFileSuffix
  invisible(TRUE)
}

#' Reload only the reference tables
#'
#' Replaces the genome references in [scCNVCaller] while leaving binning
#' parameters untouched. Use it to switch reference builds inside a session;
#' otherwise prefer [initialiseEnvironment()].
#'
#' @param chromSizeFile Path to a chromosome-size table.
#' @param cpgDataFile Path to a binned CpG density table.
#' @param cytobandFile Path to a binned cytoband table.
#'
#' @return Invisibly `TRUE`.
#' @family session setup
#' @export
initializeStandards <- function(chromSizeFile, cpgDataFile, cytobandFile) {
  for (path in c(chromSizeFile, cpgDataFile, cytobandFile)) {
    if (!file.exists(path)) {
      stop(sprintf("Reference file does not exist: %s", path), call. = FALSE)
    }
  }
  scCNVCaller$chrom_sizes <- read_chrom_sizes(chromSizeFile)
  scCNVCaller$cpg_data <- read_cpg_data(cpgDataFile)
  scCNVCaller$cytoband_data <- read_cytoband_data(cytobandFile)
  invisible(TRUE)
}

#' Read a preprocessed single-cell fragment matrix
#'
#' Reads the count matrix written by `process_fragment_file.py`: cells in rows,
#' genomic bins in columns, with the first column holding the barcode.
#'
#' @param inputFile Path to the count matrix.
#' @param sep Field separator. Default tab.
#'
#' @return A data frame of counts with barcodes as row names.
#' @family data import
#' @export
#' @examples
#' \dontrun{
#' scData <- readInputTable("my_sample_matrix.tsv")
#' }
readInputTable <- function(inputFile, sep = "\t") {
  if (!is.character(inputFile) || length(inputFile) != 1L) {
    stop("`inputFile` must be a single file path.", call. = FALSE)
  }
  if (!file.exists(inputFile)) {
    stop(sprintf("`inputFile` does not exist: %s", inputFile), call. = FALSE)
  }

  scData <- fread(inputFile, header = TRUE, stringsAsFactors = FALSE,
                  sep = sep, data.table = FALSE)
  if (ncol(scData) < 2L) {
    stop(
      sprintf(
        "Read only %d column from %s. Check that `sep` matches the file.",
        ncol(scData), inputFile
      ),
      call. = FALSE
    )
  }

  colnames(scData)[1] <- "Cell_id"
  if (anyDuplicated(scData$Cell_id) > 0L) {
    warning("Duplicate barcodes in the first column; row names were made unique.",
            call. = FALSE)
    rownames(scData) <- make.unique(as.character(scData$Cell_id))
  } else {
    rownames(scData) <- scData$Cell_id
  }
  scData$Cell_id <- NULL
  scData
}

#' Write per-session summary statistics
#'
#' @param outputSuffix Suffix appended to the output file name.
#'
#' @return A named list of the statistics, invisibly written to
#'   `<locPrefix><outPrefix><outputSuffix>.tsv`.
#' @family reporting
#' @export
saveSummaryStats <- function(outputSuffix = "_stats") {
  assert_initialised(c("locPrefix", "outPrefix"), "saveSummaryStats")
  summaryStats <- list(
    meanReadsPerCell   = scCNVCaller$meanReadsPerCell,
    startingCellCount  = scCNVCaller$startingCellCount,
    cellsPassingFilter = scCNVCaller$cellsPassingFilter,
    blacklistCount     = scCNVCaller$blacklistCount,
    finalFilterCells   = scCNVCaller$finalFilterCells
  )
  write.table(
    summaryStats,
    file = str_c(scCNVCaller$locPrefix, scCNVCaller$outPrefix, outputSuffix,
                 ".tsv"),
    quote = FALSE
  )
  summaryStats
}
