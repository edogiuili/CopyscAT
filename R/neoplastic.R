#' Separate neoplastic from non-neoplastic cells
#'
#' Factorises the collapsed arm matrix with \pkg{NMF}, clusters cells on the
#' factor loadings, and treats the cluster with the lowest variance across arms
#' as the non-neoplastic population. Use the returned barcodes as the baseline
#' for [identifyCNVClusters()] or [annotateCNV4B()].
#'
#' Inspect the heatmap it writes and adjust `cutHeight` and `nmfComponents`
#' when the clusters do not separate cleanly.
#'
#' @param inputMatrix Collapsed, filtered matrix from [filterCells()].
#' @param estimatedCellularity Expected tumour cellularity. Works best below
#'   about 0.9.
#' @param nmfComponents Number of NMF components.
#' @param outputHeatmap Write a heatmap and violin plot to the output path.
#' @param cutHeight Dendrogram cut height, as a fraction of maximum height.
#' @param methodHclust Linkage method passed to [fastcluster::hclust()].
#'
#' @return A list with `cellAssigns` (cluster per barcode), `normalBarcodes`
#'   and `clusterNormal` (the index of the non-neoplastic cluster).
#' @family neoplastic classification
#' @export
#' @examples
#' \dontrun{
#' nmf_results <- identifyNonNeoplastic(
#'   scData_collapse, methodHclust = "ward.D", cutHeight = 0.4
#' )
#' }
identifyNonNeoplastic <- function(inputMatrix,estimatedCellularity=0.8,nmfComponents=5,outputHeatmap=TRUE,cutHeight=0.6,methodHclust="ward.D")
{
  #uses NMF and fast hclust packages
  message("Running NMF")
  res <- NMF::nmf(column_to_rownames(inputMatrix,var="chrom"), c(nmfComponents), "brunet", seed="nndsvd", .stop=NMF::nmf.stop.threshold(0.1), maxIter=2500)
  message("Computing clusters")
  #dist<-dist(t(coef(res)),method="euclidean")
  tscale<-scale(x=t(NMF::coef(res)),center=TRUE)
  dist=dist(tscale,method="euclidean")
  
  clusts<-fastcluster::hclust(dist,method = methodHclust)
  #average or ward
  #clusts<-agnes(t(coef(res)),diss = FALSE,metric="euclidean",method="ward")
  #YEAH THIS WORKS AWESOME
  #some tricksy samples may need 4
  
  #add PDF
  cell_assigns<-cutree(clusts,h=max(clusts$height)*cutHeight)
  
  if (outputHeatmap==TRUE)
  {
    pdf(file=str_c(scCNVCaller$locPrefix,scCNVCaller$outPrefix,"_nmf_heatmap.pdf"),width=6,height=6)
    heatmap.2(NMF::coef(res),Rowv=FALSE,Colv=as.dendrogram(clusts),dendrogram="column",density.info="none",trace="none",scale="none",labCol=FALSE,col=colorRampPalette(viridis(5)),symkey=FALSE,useRaster=TRUE,ColSideColors=viridis(n=length(unique(as.character(cell_assigns))))[cell_assigns])
    legend("topright",fill=viridis(n=3),x.intersp = 0.8,y.intersp=0.8,legend=unique(as.character(cell_assigns)),horiz = FALSE,cex = 0.9,border=TRUE,bty="n")
    dev.off()
  }
  by_cluster<-left_join(rownames_to_column(data.frame(t(column_to_rownames(inputMatrix,"chrom"))),var="barcode"),
                        rownames_to_column(data.frame(cluster=cell_assigns),var="barcode"),by="barcode")
  #now estimate the residuals etc for each cluster
  if (outputHeatmap==TRUE)
  {
  pdf(file=str_c(scCNVCaller$locPrefix,scCNVCaller$outPrefix,"_nmf_violinplot.pdf"),width=9,height=6)
  print(ggplot(by_cluster %>% gather(dplyr::starts_with("chr"),key="chrom",value="value"),aes(chrom,value)) + geom_violin() + theme(axis.text.x=element_text(angle=-90,vjust = 0.5,hjust = 0, color = "#000000",size = 6)) + facet_wrap(~cluster))
  dev.off()
  }
  target_cluster_1=data.frame(t(column_to_rownames(by_cluster %>% group_by(cluster) %>% summarise_if(is.numeric,list(mean)),var="cluster"))) %>% summarise_if(is.numeric,list(var)) %>% gather(key="cluster",value="var") %>% mutate(cluster=as.numeric(str_remove(cluster,"X"))) %>% arrange(var) %>% slice_head(n=1) %>% dplyr::select(cluster)
  
  tumor_cell_ids=names(which(cell_assigns!=target_cluster_1))
  normalCluster=unlist(target_cluster_1)
  normal_cell_ids=names(which(cell_assigns==unlist(normalCluster)))
  return(list(cellAssigns=cell_assigns,normalBarcodes=normal_cell_ids,clusterNormal=normalCluster))
}
#' Estimate the fraction of cycling cells
#'
#' Fits a Gaussian mixture to the total signal on one chromosome to separate
#' 2N from 4N cells. Chromosome X is the default because it is often unaltered;
#' pick another with `chromToUse` when X is affected in your samples.
#'
#' @param inputMatrix Raw matrix, before normalisation.
#' @param sampName Sample name used in plot titles.
#' @param cutoff Minimum signal for a cell to be considered.
#' @param maxG Maximum mixture components.
#' @param logTrans Log-transform the signal first.
#' @param modelName \pkg{mclust} model name. Default `"E"` (equal variance).
#' @param chromToUse Chromosome used for the estimate.
#'
#' @return A named integer vector of mixture assignments per barcode.
#' @family neoplastic classification
#' @export
estimateCellCycleFraction<-function(inputMatrix,sampName,cutoff=10000,maxG=4,logTrans=FALSE,modelName="E",chromToUse="chrX")
{
  dt<-data.table(t(inputMatrix))
  chrXRows<-which(str_detect(colnames(inputMatrix),chromToUse))
  chrXsignal<-t(dt[chrXRows,lapply(.SD,sum)])
  #print histogram of signal
  print(hist(chrXsignal[chrXsignal>cutoff],breaks=500,main = sprintf("%s_X_signal",sampName)))
  #cutoff to remove near-zero/high-noise cells
  chrXsignal<-chrXsignal[chrXsignal>cutoff,]
  #option to log transform - may work better for some samples
  if (logTrans)
  {
    chrXsignal<-log(chrXsignal+1,base=2)
  }
  #enable gaussian noise - here I use an equal variance model, because I find it works better
  fit<-Mclust(chrXsignal,G=1:maxG,modelNames = c(modelName),initialization = list(noise=TRUE))
  #print results of clustering
  print(plot(fit,what="classification"))
  print(summary(fit))
  print(fit$parameters)
  return(fit$classification)
}
