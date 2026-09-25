#' @keywords internal
"_PACKAGE"

## usethis namespace: start
#' @importFrom changepoint cpt.meanvar cpts param.est
#' @importFrom data.table data.table fread setkey setnames transpose :=
#' @importFrom dplyr arrange bind_rows desc filter group_by if_else inner_join
#'   left_join mutate mutate_all mutate_at mutate_if n rename row_number select
#'   select_if slice_head summarise summarise_at summarise_if summarize_at
#'   summarize_if between lag
#' @importFrom dplyr %>%
#' @importFrom edgeR cpm
#' @importFrom ggplot2 aes element_line element_rect element_text facet_wrap
#'   geom_histogram geom_hline geom_line geom_smooth geom_violin geom_vline
#'   ggplot ggtitle labs scale_x_log10 scale_y_log10 theme theme_bw xlab
#' @importFrom gplots heatmap.2
#' @importFrom grDevices colorRampPalette dev.off pdf
#' @importFrom graphics abline hist legend plot segments title
#' @importFrom mclust Mclust mclustBIC priorControl
#' @importFrom rlang .data sym
#' @importFrom stats IQR as.dendrogram dist cutree median pnorm qnorm quantile
#'   rnorm sd var
#'   predict
#' @importFrom stringr fixed str_c str_detect str_ends str_remove str_split
#' @importFrom tibble as_tibble column_to_rownames rownames_to_column
#' @importFrom tidyr gather separate spread
#' @importFrom utils head setTxtProgressBar tail txtProgressBar write.table
#'   read.table
#' @importFrom viridis viridis
## usethis namespace: end
NULL

# Quiet R CMD check notes about columns referenced via non-standard evaluation
# inside dplyr/data.table pipelines.
utils::globalVariables(c(
  ".", "..chromXName", "..chromYName", "..zero_list", ".SD", "1", "2",
  "Alteration", "alteration", "arm", "Average", "Barcode", "barcode",
  "blacklist", "Cell", "Cell_id", "cellName", "Cells", "Chrom", "chrom",
  "chromArm", "chromMatch", "clust", "Clust", "clust_diff", "clust_text",
  "cluster", "count", "counts", "cpg", "cpgNum", "Density", "Diff", "dm",
  "end", "endCoord", "genes", "hgnc_symbol", "init", "init2", "interval",
  "keeper", "Loc", "loc1", "maxClust", "Mean", "mean", "medianNorm", "min",
  "name", "neoplastic", "neoplasticFALSE", "neoplasticTRUE", "newIndex",
  "normCluster", "num", "Point", "Pos", "pos", "pos_b", "ratio",
  "rank", "raw_medians", "rowname", "segment", "Start", "start", "startChrom",
  "startCoord", "tiles.end", "tiles.start", "touching", "transcript_length",
  "V1", "V2", "V3", "V4", "V5", "value", "Value", "var", "zoffset", "zscore"
))
