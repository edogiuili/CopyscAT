#' Call regions of loss of heterozygosity
#'
#' Runs changepoint analysis per chromosome on the normalised matrix to find
#' contiguous stretches of reduced signal, then keeps those recurring in at
#' least `lossCutoffCells` cells.
#'
#' @param inputMatrixIn Normalised matrix from [normalizeMatrixN()].
#' @param lossCutoff Maximum Z-score for a bin to count as lost.
#' @param uncertaintyCutLoss Maximum mixture uncertainty for a cell call.
#' @param diffThreshold Difference threshold between adjacent bins.
#' @param minLength Minimum region length in base pairs.
#' @param minSeg Minimum number of bins per changepoint segment.
#' @param lossCutoffCells Minimum cells carrying a region to keep it.
#' @param targetFun Spread function used to score regions. Default [stats::IQR()].
#' @param quantileLimit Quantile of the per-cell signal used to call losses.
#' @param cpgCutoff Drop bins with CpG density at or below this.
#' @param meanThreshold Minimum mean signal for a bin to be considered.
#' @param dummyQuantile,dummyPercentile,dummySd Parameters of the simulated
#'   diploid baseline used to score candidate regions.
#'
#' @return A list whose second element holds the called LOH intervals, in the
#'   form consumed by [annotateLosses()].
#' @family LOH
#' @export
getLOHRegions <- function(inputMatrixIn,lossCutoff=(-0.25), uncertaintyCutLoss=0.5, diffThreshold=0.9, minLength=3e6, minSeg=3, lossCutoffCells=100,targetFun=IQR,quantileLimit=0.3,cpgCutoff=0,meanThreshold=4,dummyQuantile=0.5,dummyPercentile=0.2,dummySd=0.2)
{
  #c("E","V")
  assert_initialised(c("locPrefix", "outPrefix", "cpg_data", "cytoband_data"),
                     "getLOHRegions")
  pdf(str_c(scCNVCaller$locPrefix,scCNVCaller$outPrefix,"_LOH.pdf"),width=6,height=4)
  on.exit(grDevices::dev.off(), add = TRUE)
  inputMatrix<-data.table(inputMatrixIn)
  inputMatrix<-inputMatrix[,cpg:=scCNVCaller$cpg_data$cpg_density]
  inputMatrix<-inputMatrix[,arm:=scCNVCaller$cytoband_data$V4]
  inputMatrix<-inputMatrix[which(inputMatrix[,arm!="cen" & blacklist==0 & cpg>cpgCutoff & chrom!="chrY"])]
  inputMatrix<-inputMatrix[,c("raw_medians","blacklist","cpg","arm"):=NULL]
  #inputMatrix<-input
  #check if raw_medians
  #mutate
  chromList<-unique(inputMatrix$chrom)
  # chromList<-c("chr9","chr10")
  dm_per_cell_vals<-data.frame(cellName=colnames(inputMatrix  %>% dplyr::select(dplyr::ends_with(scCNVCaller$cellSuffix))),stringsAsFactors=FALSE)
  alteration_list=c()
  last_coords=c(0,0)
  alteration_delta=c()
  sliceList<-function(fit)
  {
    #reassign items
    #assume there are only two
    #recompute means
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
     mean1a<-mean(fit$data[fit$classification==1])
     mean2b<-mean(fit$data[fit$classification==2])
    
    
    return(fit$classification)
  }
  #may not need IQR, may work with mean/median instead
  for (targetChrom in chromList)
  {
    inputMatrixK<-inputMatrix[which(inputMatrix[,chrom==targetChrom])]
    tmp_coords<-inputMatrixK[,c("chrom","pos")] 
    IQRv<-data.table::transpose(inputMatrixK[,chrom:=NULL],make.names = "pos")[,lapply(.SD,quantile,quantileLimit,na.rm=TRUE)]
    IQRs<-scale(t(IQRv),center=TRUE,scale=TRUE)
    if (is.nan(as.vector(IQRs)[1])) {
      next
    }
    expectedSignalN<-data.table::transpose(inputMatrixK[,lapply(.SD,max),by="pos"],make.names="pos")[,lapply(.SD,mean)]
    cm<-cpt.meanvar(data=as.vector(IQRs),test.stat="Normal", penalty="AIC",method = "PELT",minseglen = minSeg)
    cptlist<-t(rbind(cm@param.est$mean,cm@cpts))
    colnames(cptlist)<-c("Mean","Point")
    cptlist<-as_tibble(cptlist) %>% mutate(Diff = Point - lag(Point))
    plot(cm@data.set,xlab="Chromosome bin",ylim=c(-5,5),ylab="Z-score")
    p1<-as_tibble(cptlist) %>% mutate(Start=lag(Point))
    p1$Start[1]<-0
    p1<-p1 %>% dplyr::select(Start,Mean,Point,Mean)
    for (a in seq_len(nrow(p1)))
    {
      segments(p1$Start[a],p1$Mean[a],p1$Point[a],p1$Mean[a],col="red")
    }
    
    if (is.na(cptlist$Diff[1]))
    {
      cptlist$Diff[1]=cptlist$Point[1]-1
    }
    #set up list of double minutes; cutoff 0.002
    #1.05e-4
    coord_list<-vector(mode = "list")
    d_loss<-cptlist %>% dplyr::filter(Mean<(lossCutoff))
    d_lossb<-d_loss
    if (nrow(d_loss)>0)
    {
      d_lossb<-d_loss %>% mutate(startCoord=as.numeric(tmp_coords$pos[Point-Diff]),endCoord=as.numeric(tmp_coords$pos[Point]),startChrom=tmp_coords$chrom[Point])
      tmpb_merged<-d_lossb %>% mutate(touching=(startCoord==lag(endCoord)))
      tmpb_merged$touching[1]<-FALSE #otherwise NA
      tmpb_merged <- tmpb_merged %>% mutate(startCoord=if_else(touching==TRUE,lag(startCoord),startCoord),Diff=if_else(touching==TRUE,Diff+lag(Diff),Diff))
      if (length(which(tmpb_merged$touching==TRUE)) > 0)
      {
        tmp_merged_final<-tmpb_merged %>% dplyr::filter(!(row_number() %in% (which(tmpb_merged$touching==TRUE)-1))) %>% dplyr::select(-touching)
        d_lossb<-tmp_merged_final
      }
      d_lossb <- d_lossb %>% arrange(startChrom,as.numeric(startCoord))
    }
    #coord_list<-append(coord_list,list(c(last_coords[1],last_coords[2])))
    #change this to loop
    if (nrow(d_lossb)>=1)
    {
      for (i in seq_len(nrow(d_lossb)))
      {
        set.seed(12354)
        last_coords[1]<-d_lossb$startCoord[i]
        last_coords[2]<-d_lossb$endCoord[i]
        alterationName<-str_c(targetChrom,sprintf(fmt="%d",last_coords[1]),sprintf(fmt="%d",last_coords[2]),sep="_")
        
        posList=seq(from=last_coords[1],to=last_coords[2],by = 1e6)
        posList<-sprintf(fmt="%d",posList)
        #posList
        #or max
        t2d<-inputMatrixK[which(inputMatrixK[,pos %in% posList])][,lapply(.SD,max),.SDcols=-c("pos")] # %>% select(-cpg,-pos,-chrom) %>% summarise_if(is.numeric,max) # %>% select(-cpg)
        #generate random vals for zeros
        deadVal=unlist(quantile(expectedSignalN,dummyQuantile))
        if (length(t2d)==length(which(t2d==0)))
        {
          next
        }
        t2d<-t(t2d)
        #ignore the zeros, just make fake cells here, since we generate this already
        
        
        nCells=nrow(t2d)
        valsToFill<-data.frame(vals=rnorm(n=floor(dummyPercentile * nCells),mean=deadVal,sd = dummySd*deadVal))
        rownames(valsToFill)<-paste("X",seq(from=1,to=floor(dummyPercentile * nCells)),"-1",sep="")
        #t2d[which(t2d==0)]<-valsToFill
        t2d<-rbind(t2d,as.matrix(valsToFill))
        t2d_trans<-as.vector(t(sqrt(t2d+deadVal)))
        
        hist(t2d_trans,main=alterationName,breaks = 50)
        fit1<-Mclust(t2d_trans,G=2:3,modelNames=c("V"),initialization = list(noise=TRUE))
        if (length(fit1)==0)
        {
          fit1<-Mclust(t2d_trans,G=1:2,modelNames=c("V"),initialization = list(noise=TRUE))
          if (length(fit1)==0)
          {
            next()
          }
        }
        #plot(fit1)
        #if (!is.na(fit1))
        #{
        if (fit1$G>2)
        {
          #collapse clusters
          #which diff is bigger
          plot(fit1,what="classification")
          diff_2=abs(fit1$parameters$mean[3]-fit1$parameters$mean[2])
          diff_1=abs(fit1$parameters$mean[2]-fit1$parameters$mean[1])
          clust1_mean<-mean(fit1$parameters$mean[1:2])
          clust2_mean<-fit1$parameters$mean[3]
          if (diff_2>2.5*diff_1)
          {
            fit1$classification[fit1$classification==2]<-1
            fit1$classification[fit1$classification==3]<-2
          }
          else
          {
            #fit1$classification[fit1$classification==2]<-1
            fit1$classification[fit1$classification==3]<-2
            clust1_mean<-fit1$parameters$mean[1]
            clust2_mean<-mean(fit1$parameters$mean[2:3])
          }
          dVal2<-unlist(quantile(expectedSignalN,0.6))
          plot(fit1,what="classification")
          #recompute means
          fit1$parameters$mean[1]<-mean(fit1$data[fit1$classification==1])
          fit1$parameters$mean[2]<-mean(fit1$data[fit1$classification==2])
          fit1$classification<-sliceList(fit1)
          plot(fit1,what="classification")
          #signalDiff=2^(-1*clust2_mean)-dVal2
          signalDiff=dVal2-clust1_mean^2
          
          delta_mean = clust2_mean-clust1_mean
          #      delta_mean
          if (abs(delta_mean)>diffThreshold)
          {
            alteration_list<-cbind(alteration_list,alterationName)
            alteration_delta<-cbind(alteration_delta,signalDiff)
            fit1$classification[fit1$uncertainty>uncertaintyCutLoss]<-0
            dm_per_cell_vals<-cbind.data.frame(dm_per_cell_vals,fit1$classification[1:nCells])
          }
        }
        if (fit1$G==2)
        {
          plot(fit1,what="classification")
          
          #check for and collapse cluster assignments
          #cells in 2 with val < 1 and vice versa
          #fit1$classification[fit1$data<fit1$parameters$mean[1] & fit1$classification==2]<-1
          #fit1$classification[fit1$data>fit1$parameters$mean[2] & fit1$classification==1]<-2
          fit1$classification<-sliceList(fit1)
          plot(fit1,what="classification")
          delta_mean = fit1$parameters$mean[2]-fit1$parameters$mean[1]
          dVal2<-unlist(quantile(expectedSignalN,0.6))
          
          signalDiff=dVal2-fit1$parameters$mean[1]^2
          if (abs(delta_mean)>diffThreshold)
          {
            alteration_list<-cbind(alteration_list,alterationName)
            alteration_delta<-cbind(alteration_delta,signalDiff)
            #add delta mean
            #uncertainty less than 0.2
            fit1$classification[fit1$uncertainty>uncertaintyCutLoss]<-0
            dm_per_cell_vals<-cbind.data.frame(dm_per_cell_vals,fit1$classification[1:nCells])
          }
          # }
          #  }
        }
      }
    }
    #cleanup
    
  }
  #final cleanup; the PDF device is closed via on.exit()
  if (length(alteration_list)>0)
  {
    colnames(dm_per_cell_vals)[2:ncol(dm_per_cell_vals)]<-alteration_list
    
    #cutoff for minimum # of cells
    #count # of cells in each cluster
    #rowSums(.[2:ncol(cellQuality)]>0)
    loss_cluster_counts<-dm_per_cell_vals %>% gather(Alteration,Clust,2:ncol(dm_per_cell_vals)) %>% group_by(Alteration) %>% count(Alteration,Clust)
    min_loss<-loss_cluster_counts %>% spread(Clust,n) %>% dplyr::select(-'0') %>% mutate(min=min(`1`,`2`))
    min_loss[is.na(min_loss)]<-0
    #cut appropriately
    dm_per_cell_vals <- dm_per_cell_vals %>% dplyr::select(cellName, dplyr::all_of(min_loss$Alteration[min_loss$min>lossCutoffCells])) %>% mutate_at(dplyr::vars(dplyr::starts_with('chr')), list(~ if_else(. == 1, -1, .))) %>% mutate_at(dplyr::vars(dplyr::starts_with('chr')), list(~ if_else(. == 2, 1, .)))
    return(list(dm_per_cell_vals,alteration_list,alteration_delta))
  }
  return(list())
}
#' Annotate LOH regions with overlapping genes
#'
#' Queries Ensembl through \pkg{biomaRt} for genes overlapping each called
#' region and writes the annotated table to the session output path. Requires
#' network access.
#'
#' @param lossFile Output of [getLOHRegions()].
#'
#' @return A tibble of regions with a comma-separated `genes` column.
#' @family LOH
#' @export
annotateLosses<-function(lossFile)
{
  if (!requireNamespace("biomaRt", quietly = TRUE)) {
    stop(
      paste0("annotateLosses() needs the 'biomaRt' package. Install it with ",
             'BiocManager::install("biomaRt").'),
      call. = FALSE
    )
  }
  assert_initialised(c("locPrefix", "outPrefix"), "annotateLosses")
  coords<-t(as.data.frame(sapply(lossFile[[2]],str_split,"[._]"),stringsAsFactors = FALSE))
  colnames(coords)<-c("chrom","start","end")
  coords_final<-as_tibble(coords) %>% dplyr::select(chrom,start,end) %>% mutate(genes="")
  #annotate
  ensembl <- biomaRt::useMart("ensembl")
  ensembl = biomaRt::useDataset("hsapiens_gene_ensembl",mart=ensembl)
  filts=c('chromosome_name','start','end')
  for (i4 in seq_len(nrow(coords_final)))
  {
    chrname=str_remove(coords_final$chrom[i4],"chr")
    vals=list(chromosome_name=chrname,start=coords_final$start[i4],end=coords_final$end[i4])
    genes_found<-biomaRt::getBM(attributes=c('hgnc_symbol','transcript_length'), 
                       filters = filts, 
                       values = vals, 
                       mart = ensembl)
    coords_final$genes[i4]<-str_c(unique((genes_found %>% dplyr::filter(transcript_length>400) %>% arrange(desc(transcript_length)))$hgnc_symbol),collapse=",")
  }
  write.table(x=coords_final,file=str_c(scCNVCaller$locPrefix,scCNVCaller$outPrefix,"_loss_annotations.csv"),quote=FALSE,row.names = FALSE,sep=",")
  return(coords_final)
  
}
