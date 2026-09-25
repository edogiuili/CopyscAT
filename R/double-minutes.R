#' Detect focal amplifications in a single cell
#'
#' Scales one cell's signal on one chromosome to a robust Z-score, then runs
#' changepoint analysis ([changepoint::cpt.meanvar()]) to find segments whose
#' mean exceeds `peakCutoff` (amplifications) or falls below `lossCutoff`
#' (losses). Used per cell by [identifyDoubleMinutes()].
#'
#' @param inputMatrix Single-chromosome matrix with a `loc1` column of bin
#'   labels and one column per cell.
#' @param targetCell Column index of the cell to analyse.
#' @param doPlot Draw the changepoint fit. Default `FALSE`.
#' @param penalty_type Penalty passed to [changepoint::cpt.meanvar()].
#' @param doLosses Return losses below `lossCutoff` instead of gains.
#' @param peakCutoff Minimum segment Z-score to call an amplification.
#' @param lossCutoff Maximum segment Z-score to call a loss.
#' @param minThreshold Skip the cell unless its maximum Z-score exceeds this.
#'   Raising it speeds up large datasets at the cost of sensitivity.
#'
#' @return A character vector of interval labels (`chrom_start.chrom_end`), or
#'   `NULL` when the cell has no qualifying segment.
#' @family focal amplification
#' @export
getDoubleMinutes = function(inputMatrix,targetCell,doPlot=FALSE,penalty_type="SIC",doLosses=FALSE,peakCutoff=5,lossCutoff=-1,minThreshold=4)
{
  medianDensity<-median(inputMatrix[,targetCell])
  if (medianDensity>0){
    tmp<-inputMatrix[,targetCell]
    #tmp[tmp>4e-6]<-medianDensity
    #consider using IQR Here
    tmp<-as.numeric(scale(tmp,center=median(tmp),scale=IQR(tmp)),stringsAsFactors=FALSE)
    #adjusted minseglength 
    if (max(tmp)>minThreshold)
    {
    cm<-cpt.meanvar(data=tmp,test.stat="Normal", penalty=penalty_type,method = "PELT")
    #,minseglen = 3)  
    #TODO: need to document double minutes in some way (e.g. presence of each per cell) - this works
    #output cptlist
    if (doPlot==TRUE)
    {
      print(plot(cm,xlab="Chromosome bin",ylim=c(0,50),ylab="Chromosome Z-score"))
    }
    cptlist<-t(rbind(cm@param.est$mean,cm@cpts))
    colnames(cptlist)<-c("Mean","Point")
    cptlist<-as_tibble(cptlist) %>% mutate(Diff = Point - lag(Point))
    #set up list of double minutes; cutoff 0.002
    d_minutes<-cptlist %>% dplyr::filter(Mean>peakCutoff)
    d_losses<-cptlist %>% dplyr::filter(Mean<lossCutoff)
    if (!doLosses)
    {
      if (nrow(d_minutes)>0)
      {
        #for each row - to extend for more than one spike
        d_min_coords<-vector(mode="character",length=nrow(d_minutes))
        
        for (k in seq_len(nrow(d_minutes)))
        {
          if (is.na(d_minutes$Diff[k]))
          {
            #assign if first value to a start of 1
            d_minutes$Diff[k]=d_minutes$Point[k]-1
          }
          #just do column file here
          coords_temp<-inputMatrix[(d_minutes$Point[k] - d_minutes$Diff[k]):d_minutes$Point[k],]$loc1
          tmp_coords<-sprintf(fmt="%s.%s",coords_temp[1],coords_temp[length(coords_temp)-1])
          d_min_coords[k]=tmp_coords
        }
        return(d_min_coords)
      }
    }
    else {
      #look for losses
      if (nrow(d_losses)>0)
      {
        #for each row - to extend for more than one spike
        d_loss_coords<-vector(mode="character",length=nrow(d_losses))
        
        for (k in seq_len(nrow(d_losses)))
        {
          if (is.na(d_losses$Diff[k]))
          {
            #assign if first value to a start of 1
            d_losses$Diff[k]=d_losses$Point[k]-1
          }
          #just do column file here
          coords_temp<-inputMatrix[(d_losses$Point[k] - d_losses$Diff[k]):d_losses$Point[k],]$loc1
          tmp_coords<-sprintf(fmt="%s.%s",coords_temp[1],coords_temp[length(coords_temp)-1])
          d_loss_coords[k]=tmp_coords
        }
        #return list of double minutes for target cell
        #d_min_coords<-inputMatrix[(d_minutes$Point - d_minutes$Diff):d_minutes$Point,]$loc1 
        #d_min_coords
        
        return(d_loss_coords)
      } 
    }
    }
    else  {
      return(NULL)
    }
  }
  else {
    return(NULL)
  }
}
#' Test whether two bin intervals overlap
#'
#' Intervals are encoded as `chrom_start.chrom_end`, the format produced by
#' [getDoubleMinutes()].
#'
#' @param interval.1,interval.2 Interval labels.
#'
#' @return A character vector of length 3: merged left end, merged right end,
#'   and `"TRUE"`/`"FALSE"` for overlap. Character, not numeric, because the
#'   ends are bin labels.
#' @family focal amplification
#' @export
testOverlap <- function(interval.1,interval.2){
  isOverlap=FALSE
  leftEnd="0"
  rightEnd="0"
  interval.1.sep<-sapply(str_split(interval.1,fixed(".")),str_split,fixed("_"))
  interval.2.sep<-sapply(str_split(interval.2,fixed(".")),str_split,fixed("_"))
  
  if (as.numeric(interval.1.sep[[1]][2])<=as.numeric(interval.2.sep[[1]][2]))
  {
    leftEnd=interval.1.sep[[1]][2]
    #next check other end; if bigger then overlap is true
    if (as.numeric(interval.1.sep[[2]][2])>=as.numeric(interval.2.sep[[1]][2]))
    {
      isOverlap=TRUE
      rightEnd=interval.1.sep[[2]][2]
      if (as.numeric(interval.2.sep[[2]][2])>=as.numeric(interval.1.sep[[2]][2]))
      {
        rightEnd=interval.2.sep[[2]][2]
      }
    }
  } else {
    #right interval starts first  
    leftEnd=interval.2.sep[[1]][2]
    #next check other end; if bigger then overlap is true
    if (as.numeric(interval.2.sep[[2]][2])>=as.numeric(interval.1.sep[[1]][2]))
    {
      isOverlap=TRUE
      rightEnd=interval.2.sep[[2]][2]
      if (as.numeric(interval.2.sep[[2]][2])<=as.numeric(interval.1.sep[[2]][2]))
      {
        rightEnd=interval.1.sep[[2]][2]
      }
      
    }
  }
  return(c(leftEnd,rightEnd,isOverlap))
}

#PURPOSE: identify double minutes on all chromosomes using changepoint analysis
#' Detect recurrent focal amplifications across all cells
#'
#' Runs [getDoubleMinutes()] on every cell and chromosome, merges overlapping
#' intervals, and keeps only events recurring in enough cells. Slow on large
#' datasets; raise `minThreshold` to speed it up.
#'
#' @param inputMatrix Normalised matrix from [normalizeMatrixN()].
#' @param minCells Minimum cells carrying an event to keep it (first pass).
#' @param qualityCutoff2 Minimum cells carrying an event to keep it (final
#'   pass).
#' @param peakCutoff Minimum segment Z-score to call an amplification.
#' @param lossCutoff Maximum segment Z-score to call a loss.
#' @param doPlots Write a QC plot every `imageNumber` cells.
#' @param imageNumber Plot interval, in cells.
#' @param logTrans Treat the input as log-transformed when scaling by CpG.
#' @param cpgTransform Divide signal by CpG density before calling.
#' @param doLosses Call losses instead of gains.
#' @param minThreshold Z-score below which a cell is skipped.
#'
#' @return A data frame with one row per cell and one logical column per
#'   retained event, or `NULL` when nothing passes the filters.
#' @family focal amplification
#' @export
#' @examples
#' \dontrun{
#' dm_candidates <- identifyDoubleMinutes(
#'   scData_k_norm, minCells = 100, qualityCutoff2 = 100, minThreshold = 4
#' )
#' }
identifyDoubleMinutes<-function(inputMatrix,minCells=100,qualityCutoff2=100,peakCutoff=5,lossCutoff=-1,doPlots=FALSE,imageNumber=1000,logTrans=FALSE,cpgTransform=FALSE,doLosses=FALSE,minThreshold=4)
{
  assert_initialised(c("locPrefix", "outPrefix", "cellSuffix", "cpg_data",
                       "cytoband_data", "chrom_sizes"),
                     "identifyDoubleMinutes")
  if (nrow(inputMatrix) != nrow(scCNVCaller$cpg_data)) {
    stop(
      sprintf(
        paste0("`inputMatrix` has %d bins but the loaded CpG reference has ",
               "%d. They are matched by position, so the reference bin size ",
               "must match the matrix."),
        nrow(inputMatrix), nrow(scCNVCaller$cpg_data)
      ),
      call. = FALSE
    )
  }
  dm_per_cell<-data.frame(cellName=colnames(inputMatrix %>% dplyr::select(dplyr::ends_with(scCNVCaller$cellSuffix))),stringsAsFactors=FALSE)
  dm_per_cell[,seq(from=2,to=5)]<-FALSE
  
  #initialize repeat index
  firstRepeatIndex=2
  #MINIMUM # of cells with alteration - 50 is defualt
  qualityCutoff=minCells
  #test - can remove centromeres
  #cytoband_data$V4
  scData_k_cpg<-inputMatrix %>% mutate(cpg=scCNVCaller$cpg_data$cpg_density) %>% mutate(arm=scCNVCaller$cytoband_data$V4) %>% dplyr::filter(arm!="cen", blacklist==0) %>% dplyr::select(-arm,-blacklist) #%>%  #cpg+
  if (cpgTransform==TRUE)
  {
    if (logTrans==FALSE)
    {
      scData_k_cpg <- scData_k_cpg %>% mutate_at(dplyr::vars(dplyr::ends_with(scCNVCaller$cellSuffix)), list(~ . / (cpg + 1)))
    }
    else
    {
      scData_k_cpg <- scData_k_cpg %>% mutate_at(dplyr::vars(dplyr::ends_with(scCNVCaller$cellSuffix)), list(~ . / (log(cpg, base = 2) + 1)))
    }
  }
  for (m in 1:(nrow(scCNVCaller$chrom_sizes)-3))
  {
  
    #loop through chromosomes
    currentChrom = scCNVCaller$chrom_sizes$chrom[m]
    message(str_c("processing ",currentChrom))
    #NEW:
    scData_1h <- data.frame(scData_k_cpg %>% dplyr::filter(chrom==currentChrom ) %>% mutate(loc1=paste(chrom,pos,sep="_")) %>% dplyr::select(-chrom,-pos,-cpg))
    #seleect peak cutoff
    
    #scData_1h[,2:ncol(scData_1h)]<-scale(scData_1h[,2:ncol(scData_1h)],center=TRUE,scale=TRUE)
    dm_list=data.frame(dm="null",stringsAsFactors=FALSE)
    if (m>1)
    {
      #rm(dm_per_cell_clean)
      rm(recurrentEvents)
    }
    #initialize temporary matrix
    dm_per_cell_temp<-as.data.frame(matrix(nrow = nrow(dm_per_cell),ncol=5))
    colnames(dm_per_cell_temp)<-c("cellName","V2","V3","V4","V5")
    dm_per_cell_temp[,1]<-dm_per_cell$cellName
    #blank out temporary file
    #need to be smarter with this bit here but I can't be bothered RN
    dm_per_cell_temp[,seq(from=2,to=5)]<-FALSE
    
    #START TESTING HERE
    #temp_dm<-getDoubleMinutes(scData_1h,17,doPlot=TRUE)
    #dm_index<-which(dm_list$dm==temp_dm_boundaries)
    #dm_index
    plottableCell=FALSE
    plotEachTime=imageNumber
    pb<-txtProgressBar(min=0,max=nrow(dm_per_cell)-1,style=3)
    for (cellNum in 1:(nrow(dm_per_cell)-1))
    {
      setTxtProgressBar(pb,cellNum)
      # medianDensity<-median(scData_1h[,5])
      # if (medianDensity>0){
      #   #replace telomere values
      #   tmp[tmp>4e-6]<-medianDensity
      #   plot(tmp)
      #   #consider using IQR Here
      if ((cellNum %% plotEachTime)==0)
      {
        plottableCell=TRUE
      }
      if (doPlots & plottableCell)
      {
        pdf(str_c(scCNVCaller$locPrefix,scCNVCaller$outPrefix,"_cell",cellNum,"_",currentChrom,".pdf"),width=6,height=4)
        temp_dm<-getDoubleMinutes(scData_1h,cellNum,doLosses=doLosses,peakCutoff = peakCutoff,lossCutoff = lossCutoff,doPlot=TRUE,minThreshold=minThreshold)
        dev.off()
        plottableCell=FALSE
      }
      else
      {
        temp_dm<-getDoubleMinutes(scData_1h,cellNum,doLosses=doLosses,peakCutoff = peakCutoff,lossCutoff = lossCutoff,doPlot=FALSE,minThreshold=minThreshold)
      }
      
      if (!is.null(temp_dm))
      {
        #iterate through each DM identified
        for (j in seq_along(temp_dm))
        {
          #LENGTH CHECK on temp_dm
          #get boundaries - this will be our identifier for this location
          #OBSOLETE: temp_dm_boundaries<-sprintf(fmt="%s.%s",temp_dm[j,1],temp_dm[j,ncol(temp_dm)])
          #this works to ID if this vector already exists
          #,arr.ind=TRUE
          dm_index<-which(dm_list==temp_dm[j])
          if (length(dm_index)==0)
          {
            dm_list<-rbind.data.frame(dm_list,data.frame(dm=temp_dm[j],stringsAsFactors = FALSE),stringsAsFactors=FALSE)
            dm_per_cell_temp[cellNum,nrow(dm_list)]=TRUE
            #TODO: get way to name this after the interval (name column name)
          } else {
            
            # dm_index2<-which(dm_list==temp_dm_boundaries,arr.ind=TRUE)
            #because cellName is a column, and first value of dm is always blank, so don't add one
            dm_per_cell_temp[cellNum,dm_index]=TRUE
          }
        }
      }
    }
    close(pb)
    names(dm_per_cell_temp)[2:nrow(dm_list)]<-dm_list$dm[2:nrow(dm_list)]
    #dm_per_cell
    #test for overlaps - start index 2
    #start one smaller
    
    dm_list_test2=dm_list %>% mutate(newIndex=1)
    
    chrom_name=currentChrom
    #now loop through all intervals - if detected greater than one
    if (nrow(dm_list_test2)>1)
    {
      for (j1 in 2:(nrow(dm_list_test2)))
      {
        if(dm_list_test2$newIndex[j1]==1)
        {
          dm_list_test2$newIndex[j1]=j1
        }
        #(j1+1)
        for (j2 in 2:(nrow(dm_list_test2)))
        {
          if (dm_list_test2$newIndex[j2]!=j1)
          {
            interval.1<-dm_list_test2$dm[j1]
            interval.2<-dm_list_test2$dm[j2]
            overlapResult<-testOverlap(interval.1,interval.2)
            if (overlapResult[3]=="TRUE")
            {
              #update table
              dm_list_test2$dm[j1]=str_c(chrom_name,"_",overlapResult[1],".",chrom_name,"_",overlapResult[2])
              dm_list_test2$dm[j2]=str_c(chrom_name,"_",overlapResult[1],".",chrom_name,"_",overlapResult[2])
              dm_list_test2$newIndex[j2]=dm_list_test2$newIndex[j1]
            }
          }
        }
      }
      #NOTE: still some minor kinks in the overlap function; it works but the intervals aren't coming up perfect
      #TODO: fudge indices to empty blanks
      tmp_fact<-factor(dm_list_test2$newIndex)
      levels(tmp_fact)<-seq(1,length(levels(tmp_fact)))
      dm_list_test2$newIndex<-as.numeric(levels(tmp_fact)[tmp_fact])
      #option: can clean up dataframe here and replace NA with FALSE, and then proceed instead of making a giant one up front
      dm_per_cell_temp[is.na(dm_per_cell_temp)]<-FALSE
      #OK that worked, now to amalgamate cell-wide data
      for (j3 in 2:nrow(dm_list_test2))
      {
        #combine and rename columns appropriately
        dm_per_cell_temp[,dm_list_test2$newIndex[j3]]<-(dm_per_cell_temp[,j3] | dm_per_cell_temp[,dm_list_test2$newIndex[j3]])
      }
      #now rename colnames
      for (j3 in 2:nrow(dm_list_test2))
      {
        colnames(dm_per_cell_temp)[dm_list_test2$newIndex[j3]]<-dm_list_test2$dm[j3]
      }
      dm_per_cell_temp<-dm_per_cell_temp[,1:(max(dm_list_test2$newIndex))]
      
      
      recurrentEvents <- dm_per_cell_temp  %>% summarize_if(is.logical,sum) %>% select_if(function(x) any(x>qualityCutoff))
      #select only these columns that meet cutoff
      #check that recurrent events have been identified
      if (dim(recurrentEvents)[2]!=0)
      {
        message("copying recurrent events")
        dm_per_cell_clean <- dm_per_cell_temp %>% dplyr::select(cellName,colnames(recurrentEvents))
        #copy over columns and column names and increment column counter - ncol-2
        dm_per_cell[,firstRepeatIndex:(firstRepeatIndex+ncol(dm_per_cell_clean)-2)]=dm_per_cell_clean[,2:(ncol(dm_per_cell_clean))]
        #check indexes in following line
        colnames(dm_per_cell)[firstRepeatIndex:(firstRepeatIndex+ncol(dm_per_cell_clean)-2)]<-colnames(dm_per_cell_clean[2:(ncol(dm_per_cell_clean))])
        firstRepeatIndex=firstRepeatIndex + ncol(dm_per_cell_clean)-1
      }
    }
    #beautiful
  }
  #now clean this
  if (firstRepeatIndex>2)
  {
  dm_per_cell_clean<-dm_per_cell[,1:(firstRepeatIndex-1)]
  
  #can cutoff here too (def 50-75)
  secondThreshold=qualityCutoff2
  keepAlterations <- dm_per_cell_clean %>% summarize_if(is.logical,sum) %>% gather(Chrom,Value) %>% dplyr::filter(Value>secondThreshold)
  dm_per_cell_clean <- dm_per_cell_clean %>% dplyr::select(cellName, keepAlterations$Chrom)
  return(dm_per_cell_clean)
  }
  message("No recurrent focal amplifications passed the filters.")
  NULL
}
