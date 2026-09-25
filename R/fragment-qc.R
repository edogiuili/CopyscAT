#' Per-cell signal from one or more CopyscAT matrices
#'
#' Summarises how much signal each barcode carries, so that sequencing depth
#' can be compared across samples before running the CNV workflow.
#'
#' The values are **summed fragment lengths in base pairs**, not fragment
#' counts. `process_fragment_file.py` accumulates `end - start - 1` for every
#' fragment, so the row sum of a CopyscAT matrix is the total genomic span
#' covered by a cell. Use [fragmentsPerCell()] when you need actual fragment
#' counts.
#'
#' @param matrices Either a named list or character vector of matrix file
#'   paths, or a named list of already-loaded data frames. Names are used as
#'   sample labels.
#' @param sep Field separator, when reading from file. Default tab.
#'
#' @return A data frame with one row per cell and columns `sample`, `barcode`,
#'   `signal_bp` (total base pairs covered) and `nonzero_bins` (bins with any
#'   signal).
#' @family quality control
#' @seealso [fragmentsPerCell()] for true fragment counts,
#'   [plotCellDistribution()] and [plotCellKnee()] to visualise the result.
#' @export
#' @examples
#' data(scDataSamp)
#' signal <- signalPerCell(list(demo = scDataSamp))
#' summary(signal$signal_bp)
signalPerCell <- function(matrices, sep = "\t") {
  if (length(matrices) == 0L) {
    stop("`matrices` is empty.", call. = FALSE)
  }
  if (is.null(names(matrices)) || any(!nzchar(names(matrices)))) {
    stop("`matrices` must be named; the names are used as sample labels.",
         call. = FALSE)
  }

  perSample <- lapply(names(matrices), function(sampleName) {
    entry <- matrices[[sampleName]]

    if (is.character(entry)) {
      if (length(entry) != 1L) {
        stop(sprintf("`matrices[['%s']]` must be a single file path.",
                     sampleName), call. = FALSE)
      }
      if (!file.exists(entry)) {
        stop(sprintf("Matrix file for sample '%s' does not exist: %s",
                     sampleName, entry), call. = FALSE)
      }
      counts <- readInputTable(entry, sep = sep)
    } else if (is.data.frame(entry) || is.matrix(entry)) {
      counts <- as.data.frame(entry)
    } else {
      stop(sprintf(
        "`matrices[['%s']]` must be a file path, data frame or matrix, not %s.",
        sampleName, class(entry)[1]
      ), call. = FALSE)
    }

    numericCols <- vapply(counts, is.numeric, logical(1))
    if (!any(numericCols)) {
      stop(sprintf("Sample '%s' has no numeric columns to summarise.",
                   sampleName), call. = FALSE)
    }
    counts <- counts[, numericCols, drop = FALSE]

    barcodes <- rownames(counts)
    if (is.null(barcodes)) barcodes <- as.character(seq_len(nrow(counts)))

    data.frame(
      sample       = sampleName,
      barcode      = barcodes,
      signal_bp    = rowSums(counts, na.rm = TRUE),
      nonzero_bins = rowSums(counts > 0, na.rm = TRUE),
      stringsAsFactors = FALSE,
      row.names    = NULL
    )
  })

  do.call(rbind, perSample)
}

#' Fragment counts per cell from fragment files
#'
#' Counts the fragments recorded for each barcode. Unlike [signalPerCell()],
#' which reports base pairs covered, this is the number of fragments and so is
#' directly comparable with the `minFrags` argument of
#' [initialiseEnvironment()] and the `-f` option of
#' `process_fragment_file.py`.
#'
#' Reads either a 10x `fragments.tsv.gz` (barcode in the fourth column,
#' `#` comment lines skipped) or a CellRanger `singlecell.csv`. For a
#' `singlecell.csv` the `passed_filters` column is used, which counts only
#' fragments passing CellRanger's own filters and so will not exactly match a
#' direct count of the fragment file.
#'
#' @param files Named character vector of paths. Names are used as sample
#'   labels. A path ending in `.csv` is read as a CellRanger
#'   `singlecell.csv`; anything else is treated as a fragment file.
#' @param minFragments Drop barcodes with fewer fragments than this. Default
#'   `0`, which keeps everything.
#'
#' @return A data frame with columns `sample`, `barcode` and `n_fragments`.
#' @family quality control
#' @seealso [signalPerCell()] for the base-pair signal already in a CopyscAT
#'   matrix.
#' @export
#' @examples
#' \dontrun{
#' frags <- fragmentsPerCell(
#'   c(tumour1 = "tumour1/fragments.tsv.gz"),
#'   minFragments = 1000
#' )
#' plotCellDistribution(frags, "n_fragments", cutoff = 1e4)
#' }
fragmentsPerCell <- function(files, minFragments = 0) {
  if (length(files) == 0L) {
    stop("`files` is empty.", call. = FALSE)
  }
  if (is.null(names(files)) || any(!nzchar(names(files)))) {
    stop("`files` must be named; the names are used as sample labels.",
         call. = FALSE)
  }
  absent <- files[!file.exists(files)]
  if (length(absent) > 0L) {
    stop(sprintf("These fragment files do not exist: %s",
                 paste(absent, collapse = ", ")), call. = FALSE)
  }

  perSample <- lapply(names(files), function(sampleName) {
    path <- files[[sampleName]]

    if (grepl("\\.csv$", path, ignore.case = TRUE)) {
      cellranger <- fread(path, data.table = FALSE)
      needed <- c("barcode", "passed_filters")
      if (!all(needed %in% colnames(cellranger))) {
        stop(sprintf(
          paste0("'%s' looks like a CellRanger singlecell.csv but lacks the ",
                 "columns %s."),
          path, paste(setdiff(needed, colnames(cellranger)), collapse = " and ")
        ), call. = FALSE)
      }
      counted <- data.frame(
        sample = sampleName,
        barcode = as.character(cellranger$barcode),
        n_fragments = as.numeric(cellranger$passed_filters),
        stringsAsFactors = FALSE
      )
    } else {
      fragments <- fread(
        file = path, header = FALSE, select = 4L,
        col.names = "barcode", showProgress = FALSE,
        data.table = FALSE
      )
      # cellranger-atac 2.x writes '#' comment lines above the data.
      fragments <- fragments[!startsWith(fragments$barcode, "#"), ,
                             drop = FALSE]
      tallied <- table(fragments$barcode)
      counted <- data.frame(
        sample = sampleName,
        barcode = names(tallied),
        n_fragments = as.numeric(tallied),
        stringsAsFactors = FALSE
      )
    }

    counted <- counted[counted$n_fragments >= minFragments, , drop = FALSE]
    rownames(counted) <- NULL
    counted
  })

  do.call(rbind, perSample)
}

#' Summarise per-cell depth by sample
#'
#' @param cellData Output of [signalPerCell()] or [fragmentsPerCell()].
#' @param valueColumn Column to summarise. Defaults to whichever of
#'   `signal_bp` or `n_fragments` is present.
#'
#' @return A data frame with one row per sample: cell count, minimum, quartiles,
#'   median and maximum.
#' @family quality control
#' @export
#' @examples
#' data(scDataSamp)
#' summariseCellDepth(signalPerCell(list(demo = scDataSamp)))
summariseCellDepth <- function(cellData, valueColumn = NULL) {
  valueColumn <- resolve_value_column(cellData, valueColumn)

  parts <- split(cellData[[valueColumn]], cellData$sample)
  summarised <- lapply(names(parts), function(sampleName) {
    values <- parts[[sampleName]]
    data.frame(
      sample = sampleName,
      cells  = length(values),
      min    = min(values, na.rm = TRUE),
      q25    = unname(quantile(values, 0.25, na.rm = TRUE)),
      median = median(values, na.rm = TRUE),
      q75    = unname(quantile(values, 0.75, na.rm = TRUE)),
      max    = max(values, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, summarised)
}

#' Plot the per-cell depth distribution
#'
#' Draws one histogram per sample on a shared log scale, which makes the split
#' between real cells and ambient barcodes easy to see. A bimodal shape is
#' normal: the left mode is empty or ambient barcodes, the right mode cells.
#'
#' @param cellData Output of [signalPerCell()] or [fragmentsPerCell()].
#' @param valueColumn Column to plot. Defaults to whichever of `signal_bp` or
#'   `n_fragments` is present.
#' @param cutoff Optional threshold to mark with a dashed line, for example the
#'   `minFrags` you intend to use.
#' @param bins Number of histogram bins.
#' @param logScale Use a log10 x axis. Default `TRUE`.
#'
#' @return A [ggplot2::ggplot()] object.
#' @family quality control
#' @export
#' @examples
#' data(scDataSamp)
#' plotCellDistribution(signalPerCell(list(demo = scDataSamp)))
plotCellDistribution <- function(cellData, valueColumn = NULL, cutoff = NULL,
                                 bins = 60, logScale = TRUE) {
  valueColumn <- resolve_value_column(cellData, valueColumn)

  plotted <- cellData
  if (logScale) {
    nonPositive <- sum(plotted[[valueColumn]] <= 0, na.rm = TRUE)
    if (nonPositive > 0L) {
      warning(sprintf(
        "%d cells with a value of zero or less were dropped for the log scale.",
        nonPositive
      ), call. = FALSE)
      plotted <- plotted[plotted[[valueColumn]] > 0, , drop = FALSE]
    }
  }

  p <- ggplot(plotted, aes(x = .data[[valueColumn]])) +
    ggplot2::geom_histogram(bins = bins, fill = "#2c5f7c", colour = NA,
                            alpha = 0.85) +
    facet_wrap(~ sample, scales = "free_y") +
    ggplot2::labs(
      x = if (logScale) paste0(valueColumn, " (log10)") else valueColumn,
      y = "Cells"
    ) +
    ggplot2::theme_bw() +
    theme(strip.background = element_rect(fill = "grey92"))

  if (logScale) p <- p + ggplot2::scale_x_log10()
  if (!is.null(cutoff)) {
    p <- p + ggplot2::geom_vline(xintercept = cutoff, linetype = 2,
                                 colour = "#c0392b")
  }
  p
}

#' Plot ranked per-cell depth as a knee plot
#'
#' Barcodes are ranked by depth within each sample and drawn on log-log axes.
#' The elbow marks the transition from cells to background, and samples whose
#' curve sits below the others are shallower and may need different
#' `minFrags` or `blacklistCutoff` values.
#'
#' @param cellData Output of [signalPerCell()] or [fragmentsPerCell()].
#' @param valueColumn Column to plot. Defaults to whichever of `signal_bp` or
#'   `n_fragments` is present.
#' @param cutoff Optional threshold to mark with a dashed horizontal line.
#'
#' @return A [ggplot2::ggplot()] object.
#' @family quality control
#' @export
#' @examples
#' data(scDataSamp)
#' plotCellKnee(signalPerCell(list(demo = scDataSamp)))
plotCellKnee <- function(cellData, valueColumn = NULL, cutoff = NULL) {
  valueColumn <- resolve_value_column(cellData, valueColumn)

  ranked <- do.call(rbind, lapply(split(cellData, cellData$sample), function(x) {
    x <- x[order(-x[[valueColumn]]), , drop = FALSE]
    x$rank <- seq_len(nrow(x))
    x
  }))
  ranked <- ranked[ranked[[valueColumn]] > 0, , drop = FALSE]

  p <- ggplot(ranked, aes(x = .data$rank, y = .data[[valueColumn]],
                          colour = .data$sample)) +
    ggplot2::geom_line(linewidth = 0.7) +
    ggplot2::scale_x_log10() +
    ggplot2::scale_y_log10() +
    ggplot2::labs(x = "Barcode rank", y = valueColumn, colour = "Sample") +
    ggplot2::theme_bw()

  if (!is.null(cutoff)) {
    p <- p + ggplot2::geom_hline(yintercept = cutoff, linetype = 2,
                                 colour = "#c0392b")
  }
  p
}

#' Pick the depth column to use
#'
#' @param cellData A data frame from [signalPerCell()] or [fragmentsPerCell()].
#' @param valueColumn Explicit column name, or `NULL` to detect one.
#' @return The column name, as a string.
#' @noRd
resolve_value_column <- function(cellData, valueColumn = NULL) {
  if (!is.data.frame(cellData) || nrow(cellData) == 0L) {
    stop("`cellData` must be a non-empty data frame from signalPerCell() or fragmentsPerCell().",
         call. = FALSE)
  }
  if (!"sample" %in% colnames(cellData)) {
    stop("`cellData` needs a `sample` column.", call. = FALSE)
  }

  if (is.null(valueColumn)) {
    candidates <- intersect(c("signal_bp", "n_fragments"), colnames(cellData))
    if (length(candidates) == 0L) {
      stop(
        paste0("Could not find a depth column. Expected `signal_bp` or ",
               "`n_fragments`; pass `valueColumn` explicitly."),
        call. = FALSE
      )
    }
    return(candidates[1])
  }

  if (!valueColumn %in% colnames(cellData)) {
    stop(sprintf("`%s` is not a column of `cellData`. Available: %s",
                 valueColumn, paste(colnames(cellData), collapse = ", ")),
         call. = FALSE)
  }
  valueColumn
}
