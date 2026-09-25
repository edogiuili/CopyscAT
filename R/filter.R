#' Filter low-quality cells
#'
#' Drops cells covering too few chromosome arms, then drops cells whose total
#' signal lies more than `signalSDcut` standard deviations from the mean, which
#' removes likely doublets and debris.
#'
#' @param inputMatrix Collapsed matrix from [collapseChrom3N()].
#' @param minimumSegments Minimum number of arms with signal above
#'   `minDensity` for a cell to be kept.
#' @param minDensity Signal above which an arm counts as covered.
#' @param signalSDcut Standard deviations from mean total signal beyond which
#'   a cell is treated as an outlier.
#'
#' @return The input matrix with failing cell columns removed.
#' @family quality control
#' @export
#' @examples
#' \dontrun{
#' scData_collapse <- filterCells(
#'   scData_collapse, minimumSegments = 40, minDensity = 0.1
#' )
#' }
filterCells <- function(inputMatrix,minimumSegments=40,minDensity=0,signalSDcut=2)
{
  assert_initialised("cellSuffix", "filterCells")
  tmp<-inputMatrix
  initialCells<-ncol(tmp)-1
  cellQuality<-tmp %>% 
    gather(Cell,Density,dplyr::ends_with(scCNVCaller$cellSuffix)) %>% dplyr::select(chrom,Cell,Density) %>%
    spread(chrom,Density) 
  cellQuality <- cellQuality %>% mutate(count=rowSums(.[2:ncol(cellQuality)]>minDensity))
  barcodesPassing <- (cellQuality %>% dplyr::filter(count>=minimumSegments))$Cell
  message(sprintf("%d of %d barcodes passed the segment-coverage filter",
                  length(barcodesPassing), initialCells))
  if (length(barcodesPassing) == 0L) {
    stop(
      sprintf(
        paste0("No cell covered at least %d segments above density %s. Lower ",
               "`minimumSegments` or `minDensity`."),
        minimumSegments, format(minDensity)
      ),
      call. = FALSE
    )
  }

  #filter scData_prechrom for quality checks
  tmp <- tmp %>% dplyr::select(c("chrom", dplyr::all_of(barcodesPassing)))
  
  #average total signal per cell
  #remove outlier cells
  aveSignal <- tmp %>% summarize_at(dplyr::vars(dplyr::ends_with(scCNVCaller$cellSuffix)),list(sum)) %>% gather(Cell,Average)
  sd_aveSignal <- sd(aveSignal$Average)
  mean_aveSignal <- mean(aveSignal$Average)
  #sd_aveSignal
  #mean_aveSignal - signalSDcut*sd_aveSignal
  #select cells
  cellsPassingB <- aveSignal %>% dplyr::filter(between(Average,mean_aveSignal - signalSDcut*sd_aveSignal,mean_aveSignal + signalSDcut*sd_aveSignal)) %>% dplyr::select(Cell)
  message(sprintf("%d of %d barcodes passed the signal-outlier filter",
                  length(cellsPassingB$Cell), length(barcodesPassing)))
  if (length(cellsPassingB$Cell) == 0L) {
    stop(
      sprintf(
        "No cell survived the signal-outlier filter; try raising `signalSDcut` (currently %s).",
        format(signalSDcut)
      ),
      call. = FALSE
    )
  }
  #filter scData_prechrom for quality checks
  scData_chrom_filtered <- tmp %>% dplyr::select(c("chrom", dplyr::all_of(cellsPassingB$Cell)))
  #average total signal per cell
  scCNVCaller$finalFilterCells<-length(cellsPassingB$Cell)
  return(scData_chrom_filtered)
}

#recompute range values after normalization and filtering
#' Compute per-arm centre and spread
#'
#' Returns the summary value and interquartile range of each chromosome arm
#' across cells. Feed the result to [scaleMatrix()] and
#' [identifyCNVClusters()].
#'
#' @param inputMatrix Collapsed, filtered matrix.
#' @param summaryFunction Function used for the centre. Default [cutAverage()].
#'
#' @return A list of two data frames, each with columns `chrom` and `Value`:
#'   the centres and the IQRs.
#' @family normalisation
#' @export
#' @examples
#' \dontrun{
#' median_iqr <- computeCenters(scData_collapse, summaryFunction = cutAverage)
#' }
computeCenters <- function(inputMatrix,summaryFunction=cutAverage)
{
  DT1<-data.table(Chrom=inputMatrix$chrom,Cells=inputMatrix %>% dplyr::select(dplyr::ends_with(scCNVCaller$cellSuffix)))
  DT1b<-data.table::transpose(DT1,make.names = "Chrom")
  DT1b2<-DT1b[1:(nrow(DT1b)),lapply(.SD,summaryFunction)]
  DT1b_IQR<-DT1b[1:(nrow(DT1b)),lapply(.SD,IQR)]

  return(list(
    rownames_to_column(data.frame(Value=t(DT1b2)),"chrom"),
    rownames_to_column(data.frame(Value=t(DT1b_IQR)),"chrom"))
    )
}
#' Scale arm signal to Z-scores
#'
#' Centres and scales each arm using the centres and IQRs from
#' [computeCenters()], and drops `chrYp`.
#'
#' @param inputMatrix Collapsed, filtered matrix.
#' @param median_iqr_list List of centres and IQRs from [computeCenters()].
#' @param spread Return long format (`Cell`, `Density`) instead of wide.
#'
#' @return The scaled matrix, wide or long depending on `spread`.
#' @family normalisation
#' @export
scaleMatrix <- function(inputMatrix,median_iqr_list,spread=FALSE)
{
  #FINAL FILTER: remove Y chromosome
  tmp<-inputMatrix
  median_list<-median_iqr_list[[1]]
  iqr_list<-median_iqr_list[[2]]
  tmp[2:ncol(tmp)]=t(scale(t(tmp[2:ncol(tmp)]),center=median_list$Value,scale=iqr_list$Value))
  tmp <- tmp %>% dplyr::filter(chrom!="chrYp") 
  if (spread==TRUE)
  {
    scData_chrom_spread <- tmp %>% gather(Cell,Density,2:ncol(tmp))
    return(scData_chrom_spread)
  }
  else
  {
    return(tmp)
  }
}
