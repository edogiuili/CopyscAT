#' Call sub-chromosomal gains and losses against a normal population
#'
#' Fits a negative binomial model per bin comparing neoplastic with
#' non-neoplastic cells, segments the resulting fold-change track by
#' changepoint analysis, and assigns per-cell copy number to each segment by
#' linear discriminant analysis. Requires the normal population identified by
#' [identifyNonNeoplastic()].
#'
#' @param inputMatrix Normalised matrix from [normalizeMatrixN()].
#' @param clusterResults Output of [identifyNonNeoplastic()].
#' @param minSeglen Minimum bins per changepoint segment.
#' @param lossThreshold Fold change below which a segment is a loss.
#' @param gainThreshold Fold change above which a segment is a gain.
#' @param minCutoff Minimum cells in the smallest group of a retained segment.
#' @param rangeThreshold Minimum fold-change range for a retained segment.
#'
#' @return A data frame of barcodes with estimated copy number per called
#'   segment, named `chrom:start-end`.
#' @family CNV calling
#' @export
#' @examples
#' \dontrun{
#' altered_segments <- getAlteredSegments(
#'   scData_k_norm, nmf_results, lossThreshold = 0.6
#' )
#' }
getAlteredSegments <- function(inputMatrix,clusterResults,minSeglen=3,lossThreshold=0.5,gainThreshold=1.6,minCutoff=50, rangeThreshold=0.4)
{
  assert_initialised(c("locPrefix", "outPrefix", "cellSuffix"),
                     "getAlteredSegments")
  if (!all(c("cellAssigns", "clusterNormal") %in% names(clusterResults))) {
    stop(
      "`clusterResults` must be the list returned by identifyNonNeoplastic().",
      call. = FALSE
    )
  }
  coefs_list=data.frame()
  
  curChrom = ""
  
  #iterate overall chromosomes
  chromList<-unique(inputMatrix$chrom)
  for (curChrom in chromList)
  {
    
    test1<-inputMatrix %>% dplyr::filter(chrom==curChrom)
    
    for (i in seq_along(unique(test1$pos)))
    {
      test_counts<-t(test1[i,] %>% dplyr::select(dplyr::ends_with(scCNVCaller$cellSuffix)))
      merged_by_cluster<-inner_join(rownames_to_column(data.frame(counts=test_counts[,1]),var="barcode"),
                                    rownames_to_column(data.frame(cluster=clusterResults$cellAssigns),var="barcode"),
                                    by="barcode")
      merged_by_cluster<-merged_by_cluster %>% mutate(neoplastic=(cluster!=clusterResults$clusterNormal)) %>%  mutate(cluster=factor(cluster),neoplastic=factor(neoplastic))
      merged_by_cluster %>% group_by(neoplastic) %>% summarise(mean=median(counts),counts=n())
      
      if (mean(merged_by_cluster$counts)>5)
      {
        fit=MASS::glm.nb(counts~neoplastic+0,data=merged_by_cluster,link=identity,init.theta=0.5,start=c(100,100))
        
        if (nrow(coefs_list)>0)
        {
          coefs_list=rbind(coefs_list,data.frame(t(fit$coefficients),segment=test1[i,"pos"],chrom=curChrom))
        }
        else
        {
          coefs_list=data.frame(t(fit$coefficients),segment=test1[i,"pos"],chrom=curChrom)
        }
      }
    }
  }
  #iterate through chromosomes
  allCpts=list()
  #print outputs of each chromosome
  pdf(str_c(scCNVCaller$locPrefix,scCNVCaller$outPrefix,"_segment_scatterplots.pdf"),width=5,height=3)
  on.exit(grDevices::dev.off(), add = TRUE)
  for (curChrom in chromList)
  {
    if (curChrom!="chrY" & curChrom  != "chrM")
    {
      #compute log ratio
      withRatio<-coefs_list %>% dplyr::filter(chrom==curChrom) %>% mutate(ratio=log(neoplasticTRUE+1,base=2) - log(neoplasticFALSE+1,base=2)) %>% dplyr::select(-chrom)
      withRatio$ratio = 2^withRatio$ratio
      
      print(ggplot(withRatio %>% mutate(segment=as.numeric(segment)),aes(segment,ratio)) + geom_smooth(level=0.95,span=0.3) + ggtitle(curChrom))
      cpts=changepoint::cpt.meanvar(as.numeric(unlist(withRatio %>% dplyr::select(ratio))),penalty = "AIC",method = "BinSeg",minseglen = minSeglen,Q=20) # can change to 3  for fewer bits
      cpt_list<-data.frame(end=cpts@cpts, start=lag(cpts@cpts,default=0)+1,mean=cpts@param.est$mean) %>% dplyr::filter(mean>gainThreshold | mean<lossThreshold) %>% mutate(chrom=curChrom)
      #append start and end based on loc
      cpt_list$end=withRatio$segment[cpt_list$end]
      cpt_list$start=withRatio$segment[cpt_list$start]
      allCpts=append(allCpts,list(cpt_list))
    }
  }
  all_roi<-bind_rows(allCpts)
  #iterate through regions
  allVals=list()
  for (j in seq_len(nrow(all_roi)))
  {
    curStart=all_roi$start[j]
    curEnd=all_roi$end[j]
    cellSubset=inputMatrix %>% dplyr::filter(chrom==all_roi$chrom[j])
    locStart=which(cellSubset$pos==curStart)
    locEnd=which(cellSubset$pos==curEnd)
    cellSubset=cellSubset[locStart:locEnd,]
    merged_by_cluster<-inner_join(rownames_to_column(data.frame(mean=colMeans(cellSubset %>% dplyr::select(dplyr::ends_with(scCNVCaller$cellSuffix)))),var="barcode"),
                                  rownames_to_column(data.frame(cluster=clusterResults$cellAssigns),var="barcode"),
                                    by="barcode")
    merged_by_cluster<-merged_by_cluster %>% mutate(neoplastic=(cluster!=clusterResults$clusterNormal)) %>%  mutate(cluster=factor(cluster),neoplastic=factor(neoplastic))
    merged_by_cluster %>% group_by(neoplastic) %>% summarise(mean=median(mean),counts=n())
    #can  we feed residuals into  the decomposition?
    #estimated ratio
    testdata=(merged_by_cluster %>% dplyr::select(-barcode,-neoplastic))
    discrims<-MASS::lda(cluster~.,data=testdata)
    preds = discrims %>% predict(testdata)
    preds$class
    estimated_losses=dplyr::bind_cols(left_join(data.frame(cluster=preds$class),
                                         rownames_to_column(data.frame(discrims$means),var="cluster")) %>% mutate(mean=mean/discrims$means[clusterResults$clusterNormal]),barcode=merged_by_cluster$barcode) %>% mutate(alteration=j)
    allVals=append(allVals,list(estimated_losses))
  }
  #if range less than 0.4 remove column
  valsTable=bind_rows(allVals)
  alterationsTable=valsTable %>% dplyr::select(mean,barcode,alteration) # %>% spread(value=mean,key=barcode)
  #name alterations
  all_roi<-all_roi %>% mutate(name=sprintf("%s:%s-%s",chrom,start,end),alteration=row_number())
  final_alterations<-left_join(alterationsTable,all_roi %>% dplyr::select(alteration,name),by="alteration")
  final_alterations<-final_alterations %>% dplyr::select(barcode, mean,name) %>% spread(key=name,value=mean)
  all_ranges<-apply(final_alterations %>% dplyr::select(-barcode),2,max)-apply(final_alterations %>% dplyr::select(-barcode),2,min)
  which(all_ranges>rangeThreshold)
  final_alterations<-final_alterations[c("barcode",names(which(all_ranges>rangeThreshold)))]
  #then remove ones with max of one group < X (e.g.  50)
  clean_list=unlist(final_alterations %>% gather("alteration","value",dplyr::starts_with("chr")) %>% group_by(alteration,value) %>% summarise(min=n()) %>% arrange(alteration,min) %>% slice_head(n=1) %>% dplyr::filter(min>minCutoff) %>% dplyr::select(alteration))
  clean_list<-as.vector(clean_list)
  power2<-function(x)
  {
    return(2^x)
  }
  #report as powers of 2 -  estimated CN
  final_alterations<-final_alterations %>% dplyr::select(c("barcode",clean_list)) %>% mutate_if(is.numeric,list(power2))
  
  return(final_alterations)
}
