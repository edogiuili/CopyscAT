#' Annotate copy number calls as loss, neutral or gain
#'
#' Labels each arm relative to the inferred normal cluster and writes the
#' per-cell calls to disk. Superseded by [annotateCNV4()], which reports
#' estimated absolute copy number instead of three-level calls.
#'
#' @param cnvResults Output of [clusterCNV()].
#' @param saveOutput Write results to the session output path.
#' @param maxClust2 Maximum clusters used upstream.
#' @param outputSuffix Suffix for the output files.
#'
#' @return A list of three elements: ternary calls per cell, the candidate CNV
#'   table, and binarised calls.
#' @family CNV annotation
#' @seealso [annotateCNV4()] for the current workflow.
#' @export
annotateCNV3 <- function(cnvResults,saveOutput=TRUE,maxClust2=4,outputSuffix="_1")
{
  if (saveOutput) {
    assert_initialised(c("locPrefix", "outPrefix"), "annotateCNV3")
  }
  cell_assignments<-cnvResults[[1]]
  cell_assignments<-column_to_rownames(cell_assignments,var = "Cells")
  chrom_clusters_final<-cnvResults[[2]]
  colnames(chrom_clusters_final)=c("Chrom",seq(1:(maxClust2-2)),"0")
  
  #identify possible CNVs - can use bins to test
  possCNV <- cell_assignments %>% summarise_if(is.numeric,list(max)) %>% gather(Chrom,maxClust) %>% dplyr::filter(maxClust>1)
  possCNV <- possCNV %>% mutate(Type="Unknown")
  #output CNV list by cluster
  #WITHOUT BLANKS
  
  #WITH BLANKS
  consensus_CNV_clusters<-rownames_to_column(cell_assignments) %>% dplyr::filter(str_detect(rowname,"X",negate=TRUE)) %>% gather(Chrom,clust,2:(ncol(cell_assignments)+1)) %>% dplyr::filter(Chrom %in% possCNV$Chrom) %>% spread(Chrom,clust) %>% arrange(rowname)
  
  if(saveOutput==TRUE)
  {
    write.table(x=consensus_CNV_clusters,file=str_c(scCNVCaller$locPrefix,scCNVCaller$outPrefix,outputSuffix,"_cnv.csv"),quote=FALSE,row.names = FALSE,sep=",")
  }
  
  colnames(chrom_clusters_final)<-c("Chrom","1","2","0","T1","T2","T3")
  possCNV <- possCNV %>% mutate(normCluster=0)
  
  for (suspCNV in possCNV$Chrom)
  {
    clusterList<-unique(chrom_clusters_final %>% dplyr::filter(Chrom==suspCNV) %>% dplyr::select(-Chrom))
    if (possCNV$maxClust[possCNV$Chrom==suspCNV]==2)
    {
      abn<-min(abs(clusterList[1:possCNV$maxClust[possCNV$Chrom==suspCNV]]))
      #first item is normal one
      abn_index<-which(abs(clusterList[1:possCNV$maxClust[possCNV$Chrom==suspCNV]])==abn)
      #compare to median autosomal value
      #if value is less than autosomal, suspect a loss, otherwise suspect a gain to be abnormal
     # possCNV
      possCNV$normCluster[possCNV$Chrom==suspCNV] =abn_index
    }
    else
    {
      #assume middle cluster is normal
      possCNV$normCluster[possCNV$Chrom==suspCNV] =2
    }
  }
  #merge cell assignments and annotate
  v1<-rownames_to_column(cell_assignments) %>% dplyr::filter(str_detect(rowname,"X",negate=TRUE)) %>% gather(Chrom,clust,2:(ncol(cell_assignments)+1)) %>% dplyr::filter(Chrom %in% possCNV$Chrom)
  v2<-left_join(v1,possCNV,by="Chrom") %>% mutate(clust_text=if_else(clust>normCluster,3,if_else(clust<normCluster,1,2))) %>% mutate(clust_text=if_else(clust==0,2,clust_text))
  v2 <- v2 %>% dplyr::select(rowname,Chrom,clust_text) %>% spread(Chrom,clust_text) %>% arrange(rowname)
  if(saveOutput==TRUE)
  {
    write.table(x=v2,file=str_c(scCNVCaller$locPrefix,scCNVCaller$outPrefix,outputSuffix,"_cnv_anno_binary.csv"),quote=FALSE,row.names = FALSE,sep=",")
  }
  colnames(possCNV)[1]<-"chrom"
  
  final_cnv_binarized<-inner_join(consensus_CNV_clusters %>% gather(chrom,clust,2:ncol(consensus_CNV_clusters)),possCNV,by="chrom") %>% mutate(clust=if_else(clust==normCluster,0,1)) %>% dplyr::select(rowname,chrom,clust) %>% spread(chrom,clust)
  return(list(v2,possCNV,final_cnv_binarized))
}
#' Estimate absolute copy number per cell
#'
#' Converts per-arm cluster Z-scores into estimated absolute copy number,
#' centred so that the modal state is 2. Arms whose copy number range is
#' narrower than `filterRange`, or whose smallest altered group has fewer than
#' `minAlteredCells` cells, are dropped.
#'
#' @param cnvResults Output of [clusterCNV()].
#' @param saveOutput Write results to the session output path.
#' @param maxClust2 Maximum clusters used upstream.
#' @param outputSuffix Suffix for the output files.
#' @param sdCNV Standard deviation of the copy number prior.
#' @param filterResults Drop arms that look unaltered.
#' @param filterRange Minimum copy number range for an arm to be retained.
#' @param minAlteredCells Minimum cells in the smallest group of an arm.
#'
#' @return A list of three elements: per-arm cluster copy numbers, the
#'   candidate CNV table, and per-cell copy number estimates.
#' @family CNV annotation
#' @export
#' @examples
#' \dontrun{
#' final_cnv_list <- annotateCNV4(
#'   candidate_cnvs_clean, saveOutput = TRUE,
#'   outputSuffix = "clean_cnv", sdCNV = 0.5, filterRange = 0.8
#' )
#' }
annotateCNV4 <- function(cnvResults,saveOutput=TRUE,maxClust2=4,outputSuffix="_1",sdCNV=0.6,filterResults=TRUE,filterRange=0.8,minAlteredCells=40)
{
  if (saveOutput) {
    assert_initialised(c("locPrefix", "outPrefix"), "annotateCNV4")
  }
  cell_assignments<-cnvResults[[1]]
  cell_assignments<-column_to_rownames(cell_assignments,var = "Cells")
  chrom_clusters_final<-cnvResults[[2]]
  shift_val=0
  if (length(which(chrom_clusters_final$V2==0))>1)
  {
    shift_val = mean(chrom_clusters_final$V1[which(chrom_clusters_final$V2==0)])
  }
  #identify possible CNVs - can use bins to test
  possCNV <- cell_assignments %>% summarise_if(is.numeric,list(max)) %>% gather(Chrom,maxClust) %>% dplyr::filter(maxClust>1)
  possCNV <- possCNV %>% mutate(Type="Unknown")
  #label by Z scores
  chrom_clusters_final<-chrom_clusters_final %>% gather("cluster","zscore",dplyr::starts_with("V")) %>% mutate(cluster=str_remove(cluster,"V")) 
  chrom_clusters_final$zscore<-sapply(sapply(chrom_clusters_final$zscore-shift_val,pnorm,log.p=TRUE),qnorm,2,sdCNV,log.p=TRUE)
  #replace each column
  trimmedCNV<-vector()
  thresholdVal=filterRange
  for (chrom in possCNV$Chrom)
  {
    zscores<-chrom_clusters_final %>% dplyr::filter(Chrom==chrom) %>% dplyr::select(zscore)
    cell_assignments[,chrom]<-as.numeric(as.character(factor(cell_assignments[,chrom],levels=seq(from=0,to=6),labels = c("2",zscores$zscore))))
    if (diff(range(zscores))>=thresholdVal)
    {
      trimmedCNV<-append(trimmedCNV,chrom)
    }
  }
  #output CNV list by cluster
  #WITHOUT BLANKS
  
  #WITH BLANKS
  if (!filterResults)
  {
    consensus_CNV_clusters<-rownames_to_column(cell_assignments) %>% dplyr::filter(str_detect(rowname,"X",negate=TRUE)) %>% gather(Chrom,clust,2:(ncol(cell_assignments)+1)) %>% dplyr::filter(Chrom %in% possCNV$Chrom) %>% spread(Chrom,clust) %>% arrange(rowname)
  }
  else
  {
    #identify chromosomes not used
    unused = cnvResults[[1]] %>% dplyr::filter(str_detect(Cells,"X",negate=TRUE)) %>% gather(chrom, cluster, dplyr::starts_with("chr")) %>% group_by(chrom,cluster) %>% dplyr::filter(cluster!=0) %>% summarise(num=n()) %>% group_by(chrom) %>% summarise(min=min(num)) %>% dplyr::filter(min<minAlteredCells)
    trimmedCNV<-trimmedCNV[!(trimmedCNV %in% unused$chrom)]
    consensus_CNV_clusters<-rownames_to_column(cell_assignments) %>% dplyr::filter(str_detect(rowname,"X",negate=TRUE)) %>% gather(Chrom,clust,2:(ncol(cell_assignments)+1)) %>% dplyr::filter(Chrom %in% trimmedCNV) %>% spread(Chrom,clust) %>% arrange(rowname)
  }
  if (saveOutput==TRUE)
  {
    write.table(x=consensus_CNV_clusters,file=str_c(scCNVCaller$locPrefix,scCNVCaller$outPrefix,outputSuffix,"_cnv_scores.csv"),quote=FALSE,row.names = FALSE,sep=",")
  }
  return(list(chrom_clusters_final,possCNV,consensus_CNV_clusters))
}
#' Estimate absolute copy number against known normal cells
#'
#' Like [annotateCNV4()], but anchors copy number 2 to a supplied set of
#' non-neoplastic barcodes rather than to the modal state. More accurate when
#' the normal population is reliable, and more prone to false positives when
#' the data are noisy.
#'
#' @param cnvResults Output of [clusterCNV()].
#' @param expectedNormals Barcodes of known or inferred normal cells, for
#'   example `nmf_results$normalBarcodes`.
#' @param saveOutput Write results to the session output path.
#' @param maxClust2 Maximum clusters used upstream.
#' @param outputSuffix Suffix for the output files.
#' @param sdCNV Standard deviation of the copy number prior.
#' @param filterResults Drop arms that look unaltered.
#' @param filterRange Minimum copy number range for an arm to be retained.
#' @param minAlteredCellProp Proportion of the normal population expected to
#'   be unaltered on a retained arm.
#'
#' @return A list of three elements, as for [annotateCNV4()].
#' @family CNV annotation
#' @export
annotateCNV4B <- function(cnvResults,expectedNormals,saveOutput=TRUE,maxClust2=4,outputSuffix="_1",sdCNV=0.6,filterResults=TRUE,filterRange=0.8,minAlteredCellProp=0.75)
{
  if (saveOutput) {
    assert_initialised(c("locPrefix", "outPrefix"), "annotateCNV4B")
  }
  if (length(expectedNormals) == 0L) {
    stop("`expectedNormals` is empty; annotateCNV4B needs normal barcodes to anchor copy number 2.",
         call. = FALSE)
  }
  knownNormals <- intersect(expectedNormals, cnvResults[[1]]$Cells)
  if (length(knownNormals) == 0L) {
    stop("None of `expectedNormals` appear in the cluster assignments.",
         call. = FALSE)
  }
  if (length(knownNormals) < length(expectedNormals)) {
    warning(
      sprintf("%d of %d `expectedNormals` are absent from the cluster assignments and were ignored.",
              length(expectedNormals) - length(knownNormals), length(expectedNormals)),
      call. = FALSE
    )
  }
  cell_assignments<-cnvResults[[1]]
  normal_clusters<-cnvResults[[1]] %>% dplyr::filter(Cells %in% expectedNormals) %>% summarise_if(is.numeric,mean) %>% mutate_all(round)
  cell_assignments<-column_to_rownames(cell_assignments,var = "Cells")
  chrom_clusters_final<-cnvResults[[2]]
  chrom_clusters_final<-chrom_clusters_final %>% mutate(norm=t(normal_clusters)) %>% mutate(zoffset=if_else(norm=="1",V1,V2))
  #no offset for X or Y
  chrom_clusters_final$zoffset[chrom_clusters_final$Chrom %in% c("chrXp","chrXq","chrYp","chrYq")]<-0
  shift_val=0
  if (length(which(chrom_clusters_final$V2==0))>1)
  {
    shift_val = mean(chrom_clusters_final$V1[which(chrom_clusters_final$V2==0)])
  }
  #identify possible CNVs - can use bins to test
  possCNV <- cell_assignments %>% summarise_if(is.numeric,list(max)) %>% gather(Chrom,maxClust) %>% dplyr::filter(maxClust>1)
  possCNV <- possCNV %>% mutate(Type="Unknown")
  #label by Z scores
  chrom_clusters_final<-chrom_clusters_final %>% gather("cluster","zscore",dplyr::starts_with("V")) %>% mutate(cluster=str_remove(cluster,"V")) 
  chrom_clusters_final$zscore<-sapply(sapply(chrom_clusters_final$zscore-chrom_clusters_final$zoffset,pnorm,log.p=TRUE),qnorm,2,sdCNV,log.p=TRUE)
  #replace each column
  trimmedCNV<-vector()
  thresholdVal=filterRange
  for (chrom in possCNV$Chrom)
  {
    zscores<-chrom_clusters_final %>% dplyr::filter(Chrom==chrom) %>% dplyr::select(zscore)
    cell_assignments[,chrom]<-as.numeric(as.character(factor(cell_assignments[,chrom],levels=seq(from=0,to=6),labels = c("2",zscores$zscore))))
    if (diff(range(zscores))>=thresholdVal)
    {
      trimmedCNV<-append(trimmedCNV,chrom)
    }
  }
  #output CNV list by cluster
  #WITHOUT BLANKS
  
  #WITH BLANKS
  if (!filterResults)
  {
    consensus_CNV_clusters<-rownames_to_column(cell_assignments) %>% dplyr::filter(str_detect(rowname,"X",negate=TRUE)) %>% gather(Chrom,clust,2:(ncol(cell_assignments)+1)) %>% dplyr::filter(Chrom %in% possCNV$Chrom) %>% spread(Chrom,clust) %>% arrange(rowname)
  }
  else
  {
    #identify chromosomes not used
    #filter those that have fewer altered than expected # of normal cells
    unused = cnvResults[[1]] %>% dplyr::filter(str_detect(Cells,"X",negate=TRUE)) %>% gather(chrom, cluster, dplyr::starts_with("chr")) %>% group_by(chrom,cluster) %>% dplyr::filter(cluster!=0) %>% summarise(num=n()) %>% group_by(chrom) %>% summarise(min=min(num)) %>% dplyr::filter(min<(minAlteredCellProp*length(expectedNormals)))
    trimmedCNV<-trimmedCNV[!(trimmedCNV %in% unused$chrom)]
    consensus_CNV_clusters<-rownames_to_column(cell_assignments) %>% dplyr::filter(str_detect(rowname,"X",negate=TRUE)) %>% gather(Chrom,clust,2:(ncol(cell_assignments)+1)) %>% dplyr::filter(Chrom %in% trimmedCNV) %>% spread(Chrom,clust) %>% arrange(rowname)
  }
  if (saveOutput==TRUE)
  {
    write.table(x=consensus_CNV_clusters,file=str_c(scCNVCaller$locPrefix,scCNVCaller$outPrefix,outputSuffix,"_cnv_scores.csv"),quote=FALSE,row.names = FALSE,sep=",")
  }
  return(list(chrom_clusters_final,possCNV,consensus_CNV_clusters))
}
