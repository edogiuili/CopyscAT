#' Collapse binned signal to chromosome arms
#'
#' Removes centromeric and blacklisted bins, divides signal by a power
#' transform of CpG density to correct the accessibility bias, then summarises
#' each chromosome arm per cell. Also infers sample sex from the X/Y signal
#' difference and records it in `scCNVCaller$isMale`.
#'
#' @param inputMatrix Normalised matrix from [normalizeMatrixN()].
#' @param minimumSegments Unused; kept for backward compatibility.
#' @param summaryFunction Function used to summarise each arm. Default
#'   [cutAverage()].
#' @param logTrans Apply a log transform to the CpG scaling factor.
#' @param binExpand Merge this many adjacent bins before collapsing. `1`
#'   leaves bins untouched.
#' @param minimumChromValue Drop arms whose median signal falls below this.
#' @param tssEnrich Weight applied to the CpG scaling factor.
#' @param logBase Log base used when `logTrans = TRUE`.
#' @param minCPG Drop bins with CpG density below this.
#' @param powVal Exponent of the CpG correction. `0.70`-`0.75` suits hg19 and
#'   hg38.
#'
#' @return A data frame with a `chrom` column of arm labels, one column per
#'   cell, and a `medianNorm` column.
#' @family normalisation
#' @export
#' @examples
#' \dontrun{
#' scData_collapse <- collapseChrom3N(
#'   scData_k_norm, summaryFunction = cutAverage,
#'   minimumChromValue = 100, minCPG = 300, powVal = 0.73
#' )
#' }
collapseChrom3N<-function(inputMatrix,minimumSegments=40,summaryFunction=cutAverage,logTrans=FALSE,binExpand=1,minimumChromValue=2,tssEnrich=5,logBase=2,minCPG=300,powVal=0.73)
{
  assert_initialised(c("cytoband_data", "cpg_data"), "collapseChrom3N")
  if (nrow(inputMatrix) != nrow(scCNVCaller$cytoband_data) ||
      nrow(inputMatrix) != nrow(scCNVCaller$cpg_data)) {
    stop(
      sprintf(
        paste0("`inputMatrix` has %d bins but the references have %d ",
               "(cytoband) and %d (CpG). All three are matched by position, ",
               "so they must come from the same bin size."),
        nrow(inputMatrix), nrow(scCNVCaller$cytoband_data),
        nrow(scCNVCaller$cpg_data)
      ),
      call. = FALSE
    )
  }

  #remove centromeric bins
  #these are ordered in the same fashion so should be accurate
  #divide signal by cpg density + 1 (To avoid dividing by zero)
  #remove blacklist regions
  scData_k <- inputMatrix %>% mutate(chromArm=scCNVCaller$cytoband_data$V4,cpg=scCNVCaller$cpg_data$cpg_density) %>% mutate(chrom=str_c(chrom,chromArm)) %>% dplyr::filter(cpg>minCPG)
  if (binExpand>1)
  {
    #collapse bins
    scData_k <- scaleMatrixBins(scData_k,binSizeRatio = binExpand,"pos")
    #rename column
    scData_k <- scData_k %>% rename(pos = pos_b)
  }
  #gender determination
  sckn<-data.table(scData_k)
  total_cpg<-sckn[,summaryFunction(cpg),by="chrom"]
  sckn[,c("blacklist","raw_medians","chromArm","cpg","pos"):=NULL]
  #tail(colnames(sckn))
  setkey(sckn,chrom)
  chromXName="chrXq"
  chromYName="chrYq"
  x_tmp<-(tail(scCNVCaller$cytoband_data[which(scCNVCaller$cytoband_data$V1=="chrX"),c(1,4)],n=1))
  y_tmp<-(tail(scCNVCaller$cytoband_data[which(scCNVCaller$cytoband_data$V1=="chrY"),c(1,4)],n=1))
  # Sex inference needs both sex chromosomes in the reference and in the data.
  # Genomes lacking them (or subsets built for testing) simply skip the step.
  canInferSex <- nrow(x_tmp) > 0L && nrow(y_tmp) > 0L
  if (canInferSex) {
    chromXName<-paste(x_tmp$V1,x_tmp$V4,sep="")
    chromYName<-paste(y_tmp$V1,y_tmp$V4,sep="")
    canInferSex <- all(c(chromXName, chromYName) %in% sckn$chrom)
  }
  if (!canInferSex) {
    message("chrX/chrY not both present; skipping sex inference ",
            "(scCNVCaller$isMale left as ", scCNVCaller$isMale, ").")
  }
  if (canInferSex) {
  sckn<-sckn[c(chromXName,chromYName),lapply(.SD,quantile,probs=0.8),by="chrom"]
  #nrow(sckn)
  #read cpg values
  sckn[,cpg:=total_cpg[chrom %in% c(chromXName,chromYName)]$V1]
  xy_signal<-data.table::transpose(sckn[,lapply(.SD,"/",1+cpg),by="chrom"][,cpg:=NULL],make.names = "chrom")[,lapply(.SD,quantile,0.8)]
  
  #xy_signal
  # Sex inference is always run on the untransformed signal: the 0.25 cutoff
  # below was calibrated on that scale. Kept deliberately independent of the
  # `logTrans` argument, which controls the CpG correction further down.
  sexCutoff=5e-7
  sexUsesLogTrans=FALSE
  if (sexUsesLogTrans)
  {
    xy_signal<-data.table::transpose(sckn[,lapply(.SD,"/",1+log(cpg,base=2)),by="chrom"][,cpg:=NULL],make.names = "chrom")[,lapply(.SD,quantile,0.8)]
    
  } else
  {
    xy_signal<-data.table::transpose(sckn[,lapply(.SD,"/",1+cpg),by="chrom"][,cpg:=NULL],make.names = "chrom")[,lapply(.SD,quantile,0.8)]
    
    sexCutoff=0.25
  }
  diffXY<-xy_signal[,..chromXName]-xy_signal[,..chromYName]
  
  if (length(diffXY) > 0L && !is.na(diffXY) && abs(diffXY)<sexCutoff)
  {
    scCNVCaller$isMale=TRUE
  }
  }
  #top part works awesome
  
  #now to convert rest
  sckn<-data.table(scData_k)
  sckn <- sckn[chromArm!="cen"][blacklist==FALSE]
  #total_cpg
  #setnames(total_cpg,"V1","totalcpg")
  #remove raw_medians,blacklist,chromArm,cpg
  sckn<-sckn[,c("blacklist","raw_medians","chromArm","pos"):=NULL]
  #mm I see the issue - the cutAverage upper quantile busts / overcorrects the cpg - use a median instead
  cpg_density<-sckn[,cpg]
  median_chromosome_density<-sckn[1:(nrow(sckn)),lapply(.SD,summaryFunction),by="chrom"]
  med_temp<-data.table::transpose(median_chromosome_density[,c("cpg"):=NULL],make.names="chrom")
  median_chrom_signal<-data.table::transpose(med_temp[,lapply(.SD,summaryFunction)],keep.names="chrom")
  
  
  #now append 
  #summaryFunction
  
  
  
  #TODO: cpg of 200 and medianNorm cutoffs only work for 1e6 bins; need to change for others
  if (logTrans==TRUE)
  {
    #5.1e-6 old 
    #5.1e-06, totalcpg 130 for 200k
    #THIS ONE
    sckn<-sckn[,lapply(.SD,"/",1+log(tssEnrich*cpg^powVal,base=logBase)),by="chrom"][,cpg:=NULL][,lapply(.SD,summaryFunction),by="chrom"]
    #+log(totalcpg,base=2)
    #was 200
  }
  else
  {
    #THIS ONE
    sckn<-sckn[,lapply(.SD,"/",1+tssEnrich*cpg^powVal),by="chrom"][,cpg:=NULL][,lapply(.SD,summaryFunction),by="chrom"]
  }
  
  sckn<-sckn[,medianNorm:=median_chrom_signal$V1]
  sckn<-sckn[which(sckn[,medianNorm>minimumChromValue])]
  #return the appropriate matrices
  return(as.data.frame(sckn))
}
