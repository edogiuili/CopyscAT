#' Smooth copy number calls within cell clusters
#'
#' Replaces per-cell copy number estimates with a per-cluster consensus, which
#' suppresses the noise inherent in single-cell ATAC coverage. A cluster is
#' called altered on an arm when more than `percentPositive` of its cells are.
#'
#' @param inputClusters Two-column data frame of barcode and cluster label.
#' @param inputCNVList Per-cell copy number estimates, normally
#'   `final_cnv_list[[3]]` from [annotateCNV4()].
#' @param inputCNVClusterFile Optional path to a CSV of copy number calls, used
#'   instead of `inputCNVList`.
#' @param percentPositive Fraction of a cluster that must carry an alteration
#'   for the cluster to be called altered.
#' @param removeEmpty Drop arms where every cluster has the same call.
#'
#' @return A data frame of barcodes with their cluster-level copy number calls.
#' @family CNV annotation
#' @export
#' @examples
#' \dontrun{
#' smoothedCNVList <- smoothClusters(
#'   scDataSampClusters, inputCNVList = final_cnv_list[[3]],
#'   percentPositive = 0.4
#' )
#' }
smoothClusters <- function(inputClusters,inputCNVList,inputCNVClusterFile="",percentPositive=0.5,removeEmpty=TRUE) 
{
  if (inputCNVClusterFile != "" && !file.exists(inputCNVClusterFile)) {
    stop(sprintf("`inputCNVClusterFile` does not exist: %s", inputCNVClusterFile),
         call. = FALSE)
  }
  if (ncol(inputClusters) < 2L) {
    stop("`inputClusters` needs at least 2 columns: barcode and cluster.",
         call. = FALSE)
  }
  #inputClusters<-read.table(inputClusterFile,stringsAsFactors=FALSE,header=TRUE,sep=",")
  colnames(inputClusters)[2]<-"clust"
  colnames(inputClusters)[1]<-"Barcode"
  inputCNV<-""
  if (inputCNVClusterFile!="")
  {
    inputCNV<-read.table(inputCNVClusterFile,stringsAsFactors=FALSE,header=TRUE,sep=",")
    inputCNV <- inputCNV %>% mutate_at(dplyr::vars(dplyr::starts_with("chr")),list(as.double))
    
  }
  else
  {
    #input from CNV calls
    inputCNV<-inputCNVList
    #colnames(inputCNV)[1]<-"Barcode"
  }
  colnames(inputCNV)[1]<-"Barcode"
  
  #colnames(inputClusters)[1]<-"Barcode"
  b1<-left_join(inputClusters,inputCNV,by="Barcode")
  #SMOOTH
  b1c<-b1 %>% mutate(clust=inputClusters[,2])
  specRound <- function(x,boost=0.1)
  {
    return (2+sign(x-2)*round(abs(x-2)+boost,digits=0))
  }
  b1c_round<-b1c %>% mutate_at(dplyr::vars(dplyr::starts_with("chr")),list(~ if_else(is.na(.), 2, .))) %>% 
    group_by(clust) %>% summarise_at(dplyr::vars(dplyr::starts_with("chr")),list(mean))  %>% 
    mutate_at(dplyr::vars(dplyr::starts_with("chr")),list(~ specRound(., 0.5 - percentPositive)))
  if (removeEmpty)
  {
  isMultiple<-function(x){
    x2<-as.numeric(x)
    if (max(x2)==min(x2)){
      return(TRUE)
    }
    else
    {
      return(FALSE)
    }}
  bb<-data.table(b1c_round)
  bb[,.SD,.SDcols=colnames(bb)[which(bb[,lapply(.SD,isMultiple),.SDcols=-c("clust")]==FALSE)]]
  b1c_round<-bb
  }
  cluster_clean<-left_join(inputClusters,b1c_round,by="clust")
  return(cluster_clean)
}
