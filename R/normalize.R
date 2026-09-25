#' Normalise a raw single-cell fragment matrix
#'
#' First analysis step. Drops outlier cells above `upperFilterQuantile`,
#' flags uninformative ("blacklist") bins that are near-empty across most
#' cells, removes cells with too many empty bins, then scales to counts per
#' million with [edgeR::cpm()].
#'
#' @param inputMatrix Raw count matrix from [readInputTable()], cells in rows
#'   and bins in columns.
#' @param logNorm Log-transform during CPM normalisation. Default `FALSE`.
#' @param maxZero Maximum number of empty bins tolerated per cell, after
#'   uninformative bins are excluded.
#' @param imputeZeros Replace zero counts with the per-bin median. Off by
#'   default and not recommended.
#' @param blacklistProp Proportion of cells that must fall below
#'   `blacklistCutoff` for a bin to be treated as uninformative. Lower is
#'   stricter.
#' @param priorCount Prior count added before log normalisation. Only used
#'   when `logNorm = TRUE`.
#' @param blacklistCutoff Signal below which a bin counts as empty in a cell.
#' @param dividingFactor Deprecated. Scaling divisor applied after CPM; set to
#'   `1` to disable.
#' @param upperFilterQuantile Cells with total signal above this quantile are
#'   dropped as likely doublets.
#'
#' @return A data frame with `chrom` and `pos` columns, a logical `blacklist`
#'   column, a `raw_medians` column, and one column per surviving cell.
#' @family normalisation
#' @export
#' @examples
#' \dontrun{
#' scData_k_norm <- normalizeMatrixN(
#'   scDataSamp,
#'   blacklistProp = 0.8, blacklistCutoff = 125, dividingFactor = 1
#' )
#' }
normalizeMatrixN <- function(inputMatrix,logNorm=FALSE,maxZero=2000,imputeZeros=FALSE,blacklistProp=0.8, priorCount=1,blacklistCutoff=100,dividingFactor=1e6,upperFilterQuantile=0.95)
{
  assert_initialised(c("locPrefix", "outPrefix", "cellSuffix"),
                     "normalizeMatrixN")
  if (!is.data.frame(inputMatrix) && !is.matrix(inputMatrix)) {
    stop("`inputMatrix` must be a data frame or matrix.", call. = FALSE)
  }
  if (nrow(inputMatrix) == 0L || ncol(inputMatrix) == 0L) {
    stop("`inputMatrix` is empty.", call. = FALSE)
  }

  sc_t<-data.table(t(inputMatrix))
  #sc_t
  cellReads<-data.table::transpose(sc_t[,lapply(.SD,sum)],keep.names="Cell")
  pdf(str_c(scCNVCaller$locPrefix,scCNVCaller$outPrefix,"_signal_distribution.pdf"),width=6,height=4)
  hist(cellReads$V1, breaks=50,main=scCNVCaller$outPrefix,xlab = "Signal")
  abline(v=quantile(cellReads$V1,upperFilterQuantile),col=c("red"),lty=2)
  dev.off()
  
  readsCells=cellReads[,mean(V1)]
  #apply quantile filter
  sc_t[,cellReads[cellReads[,V1>quantile(cellReads$V1,upperFilterQuantile)]]$Cell:=NULL]
  nCells<-nrow(cellReads)
  scCNVCaller$startingCellCount<-nCells
  message("Total number of starting cells: ", nCells,
          " Average reads per cell: ", mean(readsCells))
  scCNVCaller$meanReadsPerCell<-readsCells
  blacklistPropCutoff=blacklistProp*nCells
  #good
  #find bad columns
  sc_lines<-data.table(inputMatrix)
  #blacklistCutoff = 500
  sc_lines<-sc_lines[,lapply(.SD,function(x) x<blacklistCutoff)][,lapply(.SD,sum)]
  sc_pos<-data.table::transpose(sc_lines,keep.names = "Pos")
  blacklistRegions<-sc_pos[which(sc_pos[,V1>=blacklistPropCutoff]),]$Pos
  message(sprintf("Blacklisted regions: %d", length(blacklistRegions)))
  message(sprintf("Total bins: %d", nrow(sc_pos)))
  scCNVCaller$blacklistCount <- length(blacklistRegions)
  if (length(blacklistRegions) == nrow(sc_pos))
  {
    stop(
      sprintf(
        paste0("Every one of the %d bins was blacklisted, so no signal is ",
               "left to normalise. Lower `blacklistCutoff` (currently %s) or ",
               "raise `blacklistProp` (currently %s)."),
        nrow(sc_pos), format(blacklistCutoff), format(blacklistProp)
      ),
      call. = FALSE
    )
  }
  sc_t2<-sc_t[which(sc_pos[,V1<blacklistPropCutoff]),lapply(.SD,function(x) as.numeric(x<blacklistCutoff))][,lapply(.SD,sum)]
  
  sc_zeros<-data.table::transpose(sc_t2,keep.names="Cells")
  #zero_cutoff=zero_cutoff
  zero_list<-sc_zeros[which(sc_zeros[,V1<maxZero])]$Cells
  scCNVCaller$cellsPassingFilter<-length(zero_list)
  if (length(zero_list) == 0L) {
    stop(
      sprintf(
        paste0("No cell passed the empty-bin filter (`maxZero` = %s). Raise ",
               "`maxZero` or lower `blacklistCutoff`."),
        format(maxZero)
      ),
      call. = FALSE
    )
  }
  message(scCNVCaller$cellsPassingFilter, " cells passing filter")
  #get low confidence regions
  #message(blacklistRegions$Chrom)
  #hist(zeros$Value,breaks=100)
  tmp1 <- as.data.frame(sc_t,stringsAsFactors=FALSE) %>% dplyr::select(dplyr::all_of(zero_list))
  raw_medians<-t(data.table::transpose(sc_t[,..zero_list])[,lapply(.SD,median)])
  #tmp1
  tmp1<-cbind.data.frame(tmp1,raw_medians,stringsAsFactors=FALSE)
  #impute zeros in cells passing filter
  if (imputeZeros==TRUE)
  {
    tmp3 <- tmp1 %>% mutate_at(dplyr::vars(dplyr::ends_with(scCNVCaller$cellSuffix)), ~ if_else(. == 0, raw_medians, .)) %>% dplyr::select(-raw_medians)
  }
  else
  {
    tmp3 <- tmp1  # %>% dplyr::select(-raw_medians)
  }
  # cpm() derives library sizes from column sums and rejects any column that
  # sums to zero, reporting only "library sizes should be greater than zero"
  # with no indication of which cell or which parameter is responsible. A
  # barcode with no signal in the retained bins reaches cpm() whenever
  # `maxZero` exceeds the bin count, which disables the empty-bin filter.
  librarySizes <- colSums(tmp3, na.rm = TRUE)
  emptyCells <- names(librarySizes)[librarySizes <= 0]
  if (length(emptyCells) > 0L) {
    shown <- utils::head(emptyCells, 5)
    filterDisabled <- maxZero >= nrow(tmp3)
    stop(
      sprintf(
        paste0(
          "%d of %d barcodes have no signal left after filtering, and ",
          "edgeR::cpm() cannot normalise a cell whose counts sum to zero.\n",
          "  Affected barcodes: %s%s\n",
          "  Retained bins: %d   blacklisted: %d of %d\n",
          "  maxZero = %s%s\n",
          "To fix, do one of:\n",
          "  * lower `maxZero` below the bin count (%d) so empty barcodes ",
          "are filtered out;\n",
          "  * lower `blacklistCutoff` (currently %s) if too many bins were ",
          "blacklisted;\n",
          "  * raise `minFrags` upstream, or drop empty barcodes with ",
          "inputMatrix[rowSums(inputMatrix) > 0, ]."
        ),
        length(emptyCells), length(librarySizes),
        paste(shown, collapse = ", "),
        if (length(emptyCells) > length(shown)) {
          sprintf(" (and %d more)", length(emptyCells) - length(shown))
        } else "",
        nrow(tmp3), length(blacklistRegions), nrow(sc_pos),
        format(maxZero),
        if (filterDisabled) {
          sprintf(" -- at or above the bin count (%d), so the empty-bin filter is disabled",
                  nrow(tmp3))
        } else "",
        nrow(tmp3), format(blacklistCutoff)
      ),
      call. = FALSE
    )
  }

  #now normalize quantiles to account for differences in coverage
  scData_n<-cpm(tmp3,log = logNorm,prior.count = priorCount)
  
  rownames(scData_n)<-colnames(inputMatrix)
  colnames(scData_n)<-colnames(tmp3)
  scData_nc_split <- rownames_to_column(as.data.frame(scData_n),var = "Loc") %>% mutate(blacklist=(Loc %in% blacklistRegions)) %>% separate(col=Loc,into=c("chrom","pos"),sep="_",extra="merge",fill="right")
  scData_k<- scData_nc_split %>% mutate_at(dplyr::vars(dplyr::ends_with(scCNVCaller$cellSuffix)), list(~ (. / dividingFactor)))
  return(scData_k)
}
