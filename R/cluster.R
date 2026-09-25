#' Simulate pseudodiploid control cells
#'
#' Draws normally distributed values around `x` to act as a synthetic diploid
#' baseline when a sample contains too few non-neoplastic cells.
#'
#' @param x Mean value to simulate around.
#' @param num Number of cells to simulate.
#' @param sd_fact Standard deviation as a fraction of `x`.
#'
#' @return A numeric vector of length `num`.
#' @family CNV calling
#' @export
makeFakeCells<-function(x,num=300,sd_fact=0.1){
  a<-rnorm(n=num,mean=x,sd=sd_fact*x)
  return(a)
}
#' Call copy number states by Gaussian mixture decomposition
#'
#' Fits a [mclust::Mclust()] mixture per chromosome arm and assigns each cell
#' to a component. Components that are too close, too small or too uncertain
#' are merged or discarded. Optionally adds simulated pseudodiploid cells
#' ([makeFakeCells()]) so that a baseline exists in highly aneuploid samples.
#'
#' Writes a multi-page PDF of per-arm classification plots to the session
#' output path.
#'
#' @param inputMatrix Collapsed, filtered matrix.
#' @param median_iqr Centres and IQRs from [computeCenters()].
#' @param useDummyCells Add simulated pseudodiploid cells.
#' @param propDummy Simulated cells as a proportion of real cells.
#' @param deltaMean Minimum separation between component means.
#' @param subsetSize Cells sampled for mixture initialisation. Must not exceed
#'   the number of available cells.
#' @param fakeCellSD Standard deviation of simulated cells, as a fraction of
#'   the signal range. `0.10`-`0.20` works well.
#' @param summaryFunction Summary function for arms. Default [cutAverage()].
#' @param minDiff Minimum Z-score gap between retained clusters.
#' @param deltaBIC2 Minimum BIC gain to accept more than two components.
#' @param bicMinimum Minimum BIC gain to accept two components over one.
#' @param minMix Minimum mixing proportion for the smaller component.
#' @param uncertaintyCutoff Cells above this classification uncertainty are
#'   left unassigned (cluster `0`).
#' @param mergeCutoff Mean difference below which wide components are split
#'   rather than merged.
#' @param summarySuffix Suffix for the diagnostic PDF.
#' @param shrinky Shrinkage passed to [mclust::priorControl()].
#' @param maxClust Maximum number of components per arm.
#' @param IQRCutoff Quantile of per-arm IQRs used as the common scale.
#' @param medianQuantileCutoff Quantile of per-arm centres used as the common
#'   centre. Set to `-1` to centre on `normalCells` instead.
#' @param normalCells Barcodes of known normal cells. Required when
#'   `medianQuantileCutoff = -1`.
#' @param verbose Print the per-arm mixture summary to the console. The
#'   classification plots are written to the PDF either way.
#'
#' @return A list of three elements: per-cell cluster assignments, per-arm
#'   cluster means, and the centres/IQRs actually used.
#' @family CNV calling
#' @export
#' @examples
#' \dontrun{
#' candidate_cnvs <- identifyCNVClusters(
#'   scData_collapse, median_iqr,
#'   useDummyCells = TRUE, subsetSize = 600, maxClust = 4
#' )
#' }
identifyCNVClusters <- function(inputMatrix, median_iqr, useDummyCells = TRUE,propDummy=0.25, deltaMean=0.03, subsetSize=500, fakeCellSD=0.1, summaryFunction=cutAverage, minDiff=0.25,deltaBIC2=0.25, bicMinimum=0.1,minMix=0.3,uncertaintyCutoff=0.55, mergeCutoff=3,summarySuffix="",shrinky=0,maxClust=4,IQRCutoff=0.25,medianQuantileCutoff=0.3,normalCells=NULL,verbose=FALSE)
{
  assert_initialised(c("locPrefix", "outPrefix", "cellSuffix"),
                     "identifyCNVClusters")
  if (medianQuantileCutoff == -1 && is.null(normalCells)) {
    stop(
      paste0("`medianQuantileCutoff = -1` centres the data on known normal ",
             "cells, so `normalCells` must be supplied."),
      call. = FALSE
    )
  }
  if (!is.null(normalCells)) {
    absent <- setdiff(normalCells, colnames(inputMatrix))
    if (length(absent) > 0L) {
      stop(
        sprintf("%d of %d `normalCells` are not columns of `inputMatrix`, e.g. %s",
                length(absent), length(normalCells),
                paste(utils::head(absent, 3), collapse = ", ")),
        call. = FALSE
      )
    }
  }

  pdf(file=str_c(scCNVCaller$locPrefix,scCNVCaller$outPrefix,"_cnv_summary",summarySuffix,".pdf"),width=8,height=6)
  on.exit(grDevices::dev.off(), add = TRUE)
  
  scData_chrom<-inputMatrix
  median_chrom_signal<-median_iqr[[1]]
  median_chrom_signal_filter <- median_chrom_signal %>% dplyr::filter(chrom %in% scData_chrom$chrom) %>% dplyr::filter(chrom!="chrX",chrom!="chrY")
  IQR_chrom_signal<-median_iqr[[2]]
  IQR_chrom_signal_filter <- IQR_chrom_signal %>% dplyr::filter(chrom %in% median_chrom_signal_filter$chrom)
  if (medianQuantileCutoff!=-1)
  {
     median_chrom_signal_filter$Value<-median(median_chrom_signal_filter$Value)
  }
  #IQR_chrom_signal_filter$Value<-quantile(IQR_chrom_signal_filter$Value,0.25)
  #IQR_chrom_signal_filter$Value<-1/binSize
  #generate simulated cells
  fakeCellNum=floor(ncol(scData_chrom)*propDummy)
  fakeCells<-sapply(median_chrom_signal_filter$Value,makeFakeCells,fakeCellNum,fakeCellSD)
  fakeCellsSorted<-cbind.data.frame(chrom=median_chrom_signal_filter$chrom,t(fakeCells),stringsAsFactors=FALSE)
  colnames(fakeCellsSorted)[2:ncol(fakeCellsSorted)]=str_c("X",colnames(fakeCellsSorted)[2:ncol(fakeCellsSorted)])
  
  #generate scaffold with simulated cells
  if (useDummyCells)
  {
    scData_chrom2 <- left_join(inputMatrix,fakeCellsSorted,by="chrom")
    scData_chrom2 <- scData_chrom2 %>% dplyr::filter(chrom!="chrYp")
  }
  else
  {
    scData_chrom2 <- scData_chrom
  }
  #summaryFunction
  
  med_iqr2<-list()
  if (medianQuantileCutoff!=-1)
  {
     med_iqr2<-computeCenters(inputMatrix = scData_chrom2,summaryFunction = summaryFunction)
  }
  else
  {
    med_iqr2<-computeCenters(inputMatrix = scData_chrom2 %>% dplyr::select(chrom, dplyr::all_of(normalCells)),summaryFunction = summaryFunction) 
  }
  median_chrom_signal_filter <- med_iqr2[[1]] %>% dplyr::filter(Value>0)
  IQR_chrom_signal_filter <- med_iqr2[[2]] %>% dplyr::filter(chrom %in% median_chrom_signal_filter$chrom)
  #clean up again
  #changed from median
  if (medianQuantileCutoff!=-1)
  {
    median_chrom_signal_filter$Value<-quantile(median_chrom_signal_filter$Value,medianQuantileCutoff)
    IQR_chrom_signal_filter$Value<-quantile(IQR_chrom_signal_filter$Value,IQRCutoff)
  }
  #add simulated cells if needed
  
  # scData_chrom2 already includes the simulated cells when useDummyCells is
  # TRUE, so scaling is identical either way.
  scData_chrom_spread <- scaleMatrix(scData_chrom2,median_iqr_list = list(median_chrom_signal_filter,IQR_chrom_signal_filter),spread = TRUE)
  #NOW DO GAUSSIAN CLUSTERING
  print(ggplot(scData_chrom_spread,aes(chrom,Density))+geom_violin(scale="width") + theme(axis.text.x=element_text(angle=-90,hjust = 0,vjust=0.5),
                                                                                          
                                                                                          axis.text.y=element_text(color="#000000"),
                                                                                          axis.line.x.bottom = element_line(colour="#000000"),
                                                                                          axis.line.y.left =  element_line(colour="#000000"),
                                                                                          panel.background = element_rect(fill="#ffffff"),
                                                                                          plot.background = element_rect(color="#ffffff",fill="#ffffff")) +
          xlab("Segment")) 
  
  
  
  #initialize variables for clustering
  cell_ids=unique(scData_chrom_spread$Cell)
  chrom_clusters<-data.frame(Chrom=median_chrom_signal_filter$chrom,V1=0,V2=0,V3=0,V4=0,V5=0,V6=0,stringsAsFactors=FALSE)
  cell_assignments<-data.frame(Cells=cell_ids,stringsAsFactors = FALSE)
  
  #initialization
  set.seed(0)
  nAvailableCells <- length(unique(scData_chrom_spread$Cell))
  if (subsetSize > nAvailableCells) {
    warning(
      sprintf(
        "`subsetSize` (%d) exceeds the %d available cells; using all of them.",
        subsetSize, nAvailableCells
      ),
      call. = FALSE
    )
    subsetSize <- nAvailableCells
  }
  defSubset=sample(seq_len(nAvailableCells),subsetSize)
  print(ggplot(scData_chrom_spread %>% dplyr::filter(Cell %in% unique(scData_chrom_spread$Cell)[defSubset]),aes(chrom,Density))+geom_violin(scale="width") + theme(axis.text.x=element_text(angle=-90,hjust = 0,vjust=0.5),
                                                                                          
                                                                                          axis.text.y=element_text(color="#000000"),
                                                                                          axis.line.x.bottom = element_line(colour="#000000"),
                                                                                          axis.line.y.left =  element_line(colour="#000000"),
                                                                                          panel.background = element_rect(fill="#ffffff"),
                                                                                          plot.background = element_rect(color="#ffffff",fill="#ffffff")) +
          xlab("Segment")) 
  
  #can cut down if the mixing proportion minimum is less than X% (e.g. 20)
  for (m5 in seq_along(median_chrom_signal_filter$chrom))
  {
    set.seed(0)
    
    #hc doesn't wo
    dens1<-scData_chrom_spread %>% dplyr::filter(chrom==median_chrom_signal_filter$chrom[m5]) %>% dplyr::select(Density)
    initParams=list(noise=TRUE,subset=defSubset) #hcPairs=randomPairs(dens1$Density,modelName = "V")) #,noise=sample(1:length(dens1$Density),10)) 
    if (verbose) message(median_chrom_signal_filter$chrom[m5])
    fit = Mclust(dens1$Density,G=1:(maxClust-1),modelNames=c("V"), prior = priorControl(shrinkage=shrinky), initialization = initParams) 
    #retype if greater than two components
    if (fit$G>2)
    {
      if ((fit$BIC[fit$G]-fit$BIC[fit$G-1])<deltaBIC2)
      {
        fit = Mclust(dens1$Density,G=1:2,prior = priorControl(),modelNames="V",initialization = initParams)
      }
      else
      {
        range_lists<-list()
        #now do the matching for 3
        overlapping_clusters<-list()
        clust_list1<-c()
        for (c1 in seq_len(fit$G))
        {
          if (length(which(fit$classification==c1))>0)
          {
            range_lists[[c1]]<-range(fit$data[fit$classification==c1])
          #ignore empty clusters
            clust_list1<-append(clust_list1,c1)
          }
        }
        for (c1 in clust_list1)
        {
          for (c2 in clust_list1)
          {
            if (c1!=c2)
            {
              isOverlap = (range_lists[[c1]][1]<range_lists[[c2]][1]) & (range_lists[[c1]][2]>range_lists[[c2]][2])
              if (isOverlap==TRUE)
              {
                #with c1 being the mama cluster
                overlapping_clusters<-c(overlapping_clusters,list(c(c1,c2)))
              }
            }
          }
        }
        #process overlaps
        if (length(overlapping_clusters)>0)
        {
        for (k0 in seq_along(overlapping_clusters))
        {
          overlap_lists<-overlapping_clusters[[k0]]
          if (overlap_lists[1]!=overlap_lists[2])
          {
          #test difference in means
          diff_mean<-fit$parameters$mean[overlap_lists[2]]-fit$parameters$mean[overlap_lists[1]]
          if (abs(diff_mean) > deltaMean)
          {
            #test quantiles of larger list
            #if meets criteria then divide as previously and recompute means
            if (diff_mean<0)
            {
              #transfer small 2 to cluster 1
              fit$classification[fit$classification==overlap_lists[1] & fit$data<=fit$parameters$mean[overlap_lists[2]]]<-overlap_lists[2]
              
            }
            else
            {
              #2 is shifted right of one so right fill the cluster
              fit$classification[fit$classification==overlap_lists[1] & fit$data>=fit$parameters$mean[overlap_lists[2]]]<-overlap_lists[2]
              #relabel others
            }
            #recompute means
            fit$parameters$mean[overlap_lists[1]]<-mean(fit$data[fit$classification==overlap_lists[1]])
            fit$parameters$mean[overlap_lists[2]]<-mean(fit$data[fit$classification==overlap_lists[2]])
          }
          else
          {
            #merge small into big
            if (verbose) message(sprintf("transfer %d into %d",overlap_lists[2],overlap_lists[1]))
            fit$classification[fit$classification==overlap_lists[2]]<-overlap_lists[1]
            for (k1 in seq_along(overlapping_clusters))
                  {
                    if (verbose) message(sprintf("overlap %d; list %d",overlapping_clusters[[k1]][1],overlap_lists[1]))
                    if (overlapping_clusters[[k1]][1]==overlap_lists[2])
                    {
                      overlapping_clusters[[k1]][1]=overlap_lists[1]
                    }
                  }
          }
          }
          
        }
        }
        #check if two 
        if (length(which(fit$classification==2))==0)
        {
          #push pop 3 into 2
          fit$classification[fit$classification==3]<-2
          fit$parameters$mean[2]<-fit$parameters$mean[3]
        }
        #range1<-range(fit$data[fit$classification==1])
        #range2<-range(fit$data[fit$classification==2])
        #range3<-range(fit$data[fit$classification==3])
      }
    }
    if (fit$G==2)
    {
      #if delta BIC less than X do it
      deltaBIC = fit$BIC[2]-fit$BIC[1]
      #10-15 BIC minimum
      #or smaller pop less than predefined percentage
      arrangedMix <- fit$parameters$pro[order(fit$parameters$pro)]
      propCutoff=minMix
      minProp = 0
      if (0 %in% fit$classification)
      {
        minProp = arrangedMix[2]
      }
      else
      {
        minProp = arrangedMix[3]
      }
      if ((deltaBIC < bicMinimum) | (minProp < propCutoff))
      {
        
        #redo as single gaussian
        fit = Mclust(dens1$Density,G=1,prior = priorControl(),modelNames="V",initialization = initParams)
      }
      else
      {
        #retype if difference between means is less than cutoff
        diffMean <- (fit$parameters$mean[2]-fit$parameters$mean[1])
        
        varChange <- 0
        if (fit$modelName=="V")
        {
          arrangeVar <- fit$parameters$variance$scale[order(fit$parameters$variance$scale)]
          varChange <- (arrangeVar[2]/arrangeVar[1])
          
        }
        #if the larger slops greater than smaller then we need to fix it
       # if ((diffMean < deltaMean)) # | ((varChange > 4) & (diffMean < deltaMean * 4)))
       # {
          
      #  }
       # else
       # {
          if ((varChange > 1.8) & (diffMean<mergeCutoff))
          {
            #cut the distribution
            largerSpread <- order(fit$parameters$variance$scale)[2]
            if (largerSpread==2)
            {
              #transfer small 2 to cluster 1
              fit$classification[fit$classification==2 & fit$data<fit$parameters$mean[1]]<-1
            }
            else
            {
              #transfer larger 1 to cluster 2
              fit$classification[fit$classification==1 & fit$data>fit$parameters$mean[2]]<-2
            }
            #recompute means
            fit$parameters$mean[1]<-mean(fit$data[fit$classification==1])
            fit$parameters$mean[2]<-mean(fit$data[fit$classification==2])
          }
     #   }
      }
    }
    #test uncertainty
    fit$classification[fit$uncertainty>uncertaintyCutoff]<-0
    #print plot
    print(plot(fit,what="classification")+title(median_chrom_signal_filter$chrom[m5]))
    if (verbose) print(summary(fit))
    #can do adaptations here to check for overfitting
    chrom_clusters[chrom_clusters$Chrom==median_chrom_signal_filter$chrom[m5],2:(2+length(fit$parameters$mean)-1)]<-fit$parameters$mean
    cell_assignments <- cbind.data.frame(cell_assignments,fit$classification,stringsAsFactors=FALSE)
    
    names(cell_assignments)[m5+1]=median_chrom_signal_filter$chrom[m5]
    #TODO: compute LRT statistic
  }
  names(cell_assignments)[2:(length(median_chrom_signal_filter$chrom)+1)]<-median_chrom_signal_filter$chrom
  #try to assign groups based on imputed clusters (device closed via on.exit)
  return(list(cell_assignments,chrom_clusters,med_iqr2))
}

#OPTION: make fake cells to use as a normalizing factor (if low non-tumour content)
# #PARAMETERS for clustering
# useDummyCells=FALSE
# #maximum G value
# maxClust=4
# #minimum difference between cluster means in fit
# deltaMean=0.10
# #as above but for merging clusters
# minDiff=0.25
# #minimums for BIC of G>2 and G=2 clusters
# deltaBIC2 = 50
# bicMinimum = 5

#cluster the CNVs
#' Merge and renumber adjacent copy number clusters
#'
#' Post-processes [identifyCNVClusters()] output: orders clusters by mean,
#' merges pairs separated by less than `minDiff`, and renumbers so that
#' cluster indices are contiguous.
#'
#' @param initialResultList Output of [identifyCNVClusters()].
#' @param medianIQR Centres and IQRs, normally `initialResultList[[3]]`.
#' @param maxClust Maximum clusters used in [identifyCNVClusters()].
#' @param minDiff Minimum Z-score separation between retained clusters.
#'
#' @return A list of two elements: cleaned per-cell assignments and cleaned
#'   per-arm cluster means.
#' @family CNV calling
#' @export
#' @examples
#' \dontrun{
#' candidate_cnvs_clean <- clusterCNV(
#'   candidate_cnvs, medianIQR = candidate_cnvs[[3]], minDiff = 1.5
#' )
#' }
clusterCNV<-function(initialResultList,medianIQR,maxClust=4,minDiff=0.25){
  #return a list of CNV
  #MERGE AND COLLAPSE CLUSTERS
  chrom_clusters<-initialResultList[[2]]
  cell_assignments<-initialResultList[[1]]
  distinctClusts=1
  chrom_clusters_final<-chrom_clusters
  
  for (m5 in seq_along(medianIQR[[1]]$chrom)){
    clust_list=chrom_clusters %>% dplyr::filter(Chrom==medianIQR[[1]]$chrom[m5])
    
    zeroClusters<-length(which(clust_list[2:maxClust]==0))
    numClusts<-maxClust-zeroClusters-1
    #now order vector
    activeClusters<-clust_list[2:(2+numClusts-1)]
    hasZero<-length(which(cell_assignments[,m5+1]==0))>0
    if (numClusts>1)
    {
      #reorder the keys and values
      #reorder keys
      #copy to final list
      #only reorder if needed to save resources
      if (!is_sorted_ascending(as.numeric(activeClusters[1,])))
      {
        #reorder and reassign
        chrom_clusters_final[m5,2:(2+numClusts-1)]<-as.matrix(activeClusters)[base::order(as.numeric(activeClusters[1,]))]
        #now reorder the values by the vector using the factor levels trick
        cell_assignments[,m5+1]<-factor(cell_assignments[,m5+1])
        newOrder<-base::order(as.numeric(activeClusters[1,]))
        if (hasZero) {
          newOrder<-c(0,newOrder)
        }
        #newOrder
        #TODO: what about zeros
        levels(cell_assignments[,m5+1])<-newOrder
        #tmp_fact<-factor(dm_list_test2$newIndex)
        cell_assignments[,m5+1]<-as.numeric(levels(cell_assignments[,m5+1])[cell_assignments[,m5+1]])
      }
    }
    minDiffA=minDiff-0.00001
    while (minDiffA>0 & minDiffA<minDiff)
    {
      tmpa<-data.frame(init=seq(2,1+numClusts))
      tmpa<-tmpa %>% mutate(init2=lag(init,default=2)) %>% mutate(clust_diff=(as.double(chrom_clusters_final[m5,init])-as.double(chrom_clusters_final[m5,init2]))) %>% dplyr::filter(clust_diff!=0) %>% arrange(clust_diff)
      if (nrow(tmpa)==0)
      {
        minDiffA=0
      }
      else
      {
        minDiffA=min(tmpa$clust_diff)
        #skip looping if clusters far apart
        if (minDiffA<minDiff)
        {
          #now iterate
          for (j1 in seq_len(nrow(tmpa)))
          {
            if (tmpa$clust_diff[j1]<minDiff)
            {
              #merge these indices
              #weighted average
              number_a<-length(which(cell_assignments[,m5+1]==(tmpa$init[j1]-1)))
              number_b<-length(which(cell_assignments[,m5+1]==(tmpa$init2[j1]-1)))
              if ((number_a + number_b)!=0)
              {
                
                new_mean<-(chrom_clusters_final[m5,tmpa$init[j1]]*number_a+chrom_clusters_final[m5,tmpa$init2[j1]]*number_b)/(number_a+number_b)
                chrom_clusters_final[m5,tmpa$init[j1]]<-new_mean
                chrom_clusters_final[m5,tmpa$init2[j1]]<-new_mean
                cell_assignments[cell_assignments[,m5+1]==(tmpa$init2[j1]-1),m5+1]=tmpa$init[j1]-1
                break()
              }
              else
              {
                #cluster is empty thus not assignable
                minDiffA = 0
              }
            }
          }
        }
      }
    }
    #refactor the list - or could use dense_rank - todo later
    minCellAssign=1
    if (hasZero){
      minCellAssign=0
    }
    cell_assignments[,m5+1]<-factor(cell_assignments[,m5+1])
    uniq_clusts<-unique(as.numeric(chrom_clusters_final[m5,2:ncol(chrom_clusters_final)]))
    levels(cell_assignments[,m5+1])<-seq(minCellAssign,length(uniq_clusts))
    #tmp_fact<-factor(dm_list_test2$newIndex)
    
    
    chrom_clusters_final[m5,2:ncol(chrom_clusters_final)]<-0
    chrom_clusters_final[m5,2:(1+length(uniq_clusts))]<-uniq_clusts
    cell_assignments[,m5+1]<-as.numeric(levels(cell_assignments[,m5+1])[cell_assignments[,m5+1]])
    #TODO: merge cluster values and medians into a list which we can merge (e.g. chr0_0,0.5))
    
  }
  # chrom_clusters_final
  return(list(cell_assignments,chrom_clusters_final))
}
