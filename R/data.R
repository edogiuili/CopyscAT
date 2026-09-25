#' Example processed fragment matrix
#'
#' Binned ATAC fragment counts for 1371 adult glioblastoma cells, in the format
#' produced by `process_fragment_file.py`. Used by the examples and tests.
#'
#' @format A data frame with 1371 rows (cells) and 3103 columns (1 Mb genomic
#'   bins). Row names are cell barcodes; column names are `chrom_position`.
#' @source Adult glioblastoma single-cell ATAC sequencing.
#' @family example data
"scDataSamp"

#' Example cluster assignments
#'
#' Cluster labels for the cells in [scDataSamp], as consumed by
#' [smoothClusters()].
#'
#' @format A data frame with one row per barcode: the barcode in the first
#'   column and its cluster label in the second.
#' @source Adult glioblastoma single-cell ATAC sequencing.
#' @family example data
"scDataSampClusters"
