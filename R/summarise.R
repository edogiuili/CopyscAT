#' Trimmed mean between the configured quantiles
#'
#' Default summary function throughout CopyscAT. Averages the values lying
#' between the `lowerTrim` and `upperTrim` quantiles set by
#' [initialiseEnvironment()], which resists the heavy tails typical of
#' single-cell ATAC signal.
#'
#' @param inputVector Numeric vector.
#'
#' @return A single numeric value, or `NA_real_` when the quantiles are
#'   undefined.
#' @family summarisation
#' @export
#' @examples
#' \dontrun{
#' cutAverage(rnorm(100))
#' }
cutAverage <- function(inputVector)
{
  #remove top and bottom 5%
  #other option is to use denstiy
  #return(tmp_val_dens$x[which.max(tmp_val_dens$y)])
  
  #OPTION 5, use average of quantiles; 20 and 80
  assert_initialised(c("lowerTrim", "upperTrim"), "cutAverage")
  quants_tmp <- quantile(inputVector, c(scCNVCaller$lowerTrim, scCNVCaller$upperTrim),
                         na.rm = TRUE, type = 8)
  if (is.na(quants_tmp[1])) {
    return(NA_real_)
  }
  # The quantiles are computed with na.rm = TRUE, so drop NAs here too;
  # otherwise `inputVector[boolVal]` keeps an NA and mean() returns NA.
  boolVal <- !is.na(inputVector) &
    inputVector >= quants_tmp[1] & inputVector <= quants_tmp[2]
  if (!any(boolVal)) {
    return(NA_real_)
  }
  mean(inputVector[boolVal])
}
#' Merge adjacent genomic bins
#'
#' @param inputMatrix Matrix with `chrom` and a position column.
#' @param binSizeRatio Factor by which to widen bins. `2` doubles bin width.
#' @param inColumn Name of the position column to rebin.
#'
#' @return The matrix regrouped on the widened bins, with numeric columns
#'   summed.
#' @family summarisation
#' @export
scaleMatrixBins <- function(inputMatrix, binSizeRatio,inColumn)
{
  newMatrix<-inputMatrix %>% mutate(pos_b=floor(as.numeric(!!rlang::sym(inColumn))/(scCNVCaller$binSize*binSizeRatio))*(scCNVCaller$binSize*binSizeRatio))
  newMatrix<-newMatrix %>% dplyr::select(-inColumn) %>% group_by(pos_b,chrom) %>% summarise_if(is.numeric,list(sum)) %>% arrange(chrom,pos_b)
  
  return(newMatrix)
}
