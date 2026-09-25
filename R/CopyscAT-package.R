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
#'   geom_smooth geom_violin ggplot ggtitle theme xlab
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
  ".", ".SD", "Average", "Barcode", "Cell", "Cell_id", "Cells", "Chrom",
  "Density", "Diff", "Loc", "Mean", "Point", "Pos", "V1", "V2", "V4",
  "alteration", "arm", "barcode", "blacklist", "chrom", "chromArm",
  "chromMatch", "clust", "clust_diff", "clust_text", "cluster", "cpg",
  "cpgNum", "count",
  "counts", "end", "genes", "hgnc_symbol", "init", "init2", "interval",
  "keeper", "loc1", "maxClust", "mean", "min", "name", "neoplastic",
  "neoplasticFALSE", "neoplasticTRUE", "normCluster", "num", "pos", "pos_b",
  "ratio", "raw_medians", "rowname", "segment", "start", "tiles.end",
  "tiles.start", "transcript_length", "value", "var", "zoffset", "zscore"
,
  "..chromXName", "..chromYName", "..zero_list", "Alteration", "Clust", "Start", "Value", "cellName", "endCoord", "medianNorm", "startChrom", "startCoord", "touching", "V3", "V5", "dm", "newIndex"
))
