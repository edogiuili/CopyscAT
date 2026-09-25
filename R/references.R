#' Build CopyscAT reference files for a genome
#'
#' Tiles a \pkg{BSgenome} object into fixed-width bins and writes the three
#' reference tables that [initialiseEnvironment()] expects: CpG density,
#' cytoband arms and chromosome sizes. Cytoband and CpG tracks are fetched from
#' UCSC, so this needs network access. Run once per genome and bin size.
#'
#' @param genomeObject A \pkg{BSgenome} object, for example
#'   `BSgenome.Hsapiens.UCSC.hg38`.
#' @param genomeText UCSC genome name, for example `"hg38"`.
#' @param tileWidth Bin width in base pairs. Must match the `binSize` later
#'   passed to [initialiseEnvironment()].
#' @param outputDir Directory for the three output files.
#'
#' @return Invisibly a character vector of the files written.
#' @family session setup
#' @export
#' @examples
#' \dontrun{
#' library(BSgenome.Hsapiens.UCSC.hg38)
#' generateReferences(
#'   BSgenome.Hsapiens.UCSC.hg38, genomeText = "hg38",
#'   tileWidth = 1e6, outputDir = "."
#' )
#' }
generateReferences <- function(genomeObject,genomeText="hg38",tileWidth=1e6,outputDir="~")
{
  needed <- c("GenomicRanges", "GenomeInfoDb", "rtracklayer", "jsonlite",
              "S4Vectors", "IRanges")
  absent <- needed[!vapply(needed, requireNamespace, logical(1), quietly = TRUE)]
  if (length(absent) > 0L) {
    stop(
      sprintf(
        "generateReferences() needs these packages: %s. Install them with BiocManager::install().",
        paste(absent, collapse = ", ")
      ),
      call. = FALSE
    )
  }
  if (!dir.exists(outputDir)) {
    stop(sprintf("`outputDir` does not exist: %s", outputDir), call. = FALSE)
  }
  fileSuffixes=str_c(outputDir,"/",genomeText,"_",tileWidth,sep="")
  message("Output to ", fileSuffixes)
  #chrom sizes - 
  chrom_sizes<-rownames_to_column(data.frame(length=GenomeInfoDb::seqlengths(genomeObject)),var="chrom") %>% mutate(keeper=str_detect(chrom,pattern="alt|random|chrUn|fix|_|MT",negate=TRUE)) %>% dplyr::filter(keeper==TRUE) %>% dplyr::select(-keeper)
  chroms<-GenomicRanges::GRanges(seqnames=genomeObject@seqinfo)
  
  #https://bioconductor.org/packages/devel/bioc/vignettes/rtracklayer/inst/doc/rtracklayer.pdf
  #a<-tile(chroms,width=1e6)
  #tile the genome
  tiles<-GenomicRanges::tileGenome(GenomeInfoDb::seqlengths(chroms),tilewidth=tileWidth,cut.last.tile.in.chrom=TRUE)
  #remove excess chroms
  tiles<-GenomeInfoDb::keepStandardChromosomes(tiles,pruning.mode = "tidy")
  tbl.cytobands<-GenomicRanges::GRanges()
  tbl.cpgIslandExt<-GenomicRanges::GRanges()
  seqInfo<-rtracklayer::SeqinfoForUCSCGenome(genomeText)
  if (is.null(seqInfo))
  {
    print("fallback to UCSC REST API (manual table retrieval)")
    cpgResults = jsonlite::fromJSON(sprintf("https://api.genome.ucsc.edu/getData/track?genome=%s;track=cpgIslandExtUnmasked",genomeText),simplifyMatrix=TRUE)
    cytoBandResults = jsonlite::fromJSON(sprintf("https://api.genome.ucsc.edu/getData/track?genome=%s;track=cytoBand",genomeText),simplifyMatrix=TRUE)
    #select specifically
    cytoBandRanges<-GenomicRanges::GRanges(cytoBandResults$cytoBand[[1]])
    cpgResultRanges=GenomicRanges::GRanges()
    for (i in names(cpgResults$cpgIslandExtUnmasked))
    {
      if (length(cpgResults$cpgIslandExtUnmasked[[i]])>0)
      {
        cpgResultRanges=c(cpgResultRanges,GenomicRanges::GRanges(cpgResults$cpgIslandExtUnmasked[[i]]))
      }
    }
    tbl.cytobands<-cytoBandRanges
    tbl.cpgIslandExt<-cpgResultRanges
  }
  else
  {
    mySession = rtracklayer::browserSession("UCSC")
    rtracklayer::genome(mySession) <- genomeText
    tbl.cytobands <- GenomeInfoDb::keepStandardChromosomes(rtracklayer::track(
      rtracklayer::ucscTableQuery(mySession, table="cytoBand")),pruning.mode = "tidy")
    tbl.cpgIslandExt <- GenomeInfoDb::keepStandardChromosomes(rtracklayer::track(
      rtracklayer::ucscTableQuery(mySession, table="cpgIslandExtUnmasked")),pruning.mode = "tidy")
  }
  #sort seqlevels on cytobands
  GenomeInfoDb::seqlevels(tbl.cytobands) <- GenomeInfoDb::sortSeqlevels(GenomeInfoDb::seqlevels(tbl.cytobands))
  tbl.cytobands <- GenomicRanges::sort(tbl.cytobands)
  tbl.cpgIslandExt <- GenomicRanges::sort(tbl.cpgIslandExt)
  #sort cytoband list
  
  #now intersect these to generate the reference files
  #http://web.mit.edu/~r/current/arch/i386_linux26/lib/R/library/GenomicRanges/html/findOverlaps-methods.html
  #MATCHES FOR Cytoband
  centromeres<-which(S4Vectors::mcols(tbl.cytobands)$gieStain=="acen")
  S4Vectors::mcols(tbl.cytobands)$name[centromeres]<-"cen"
  S4Vectors::mcols(tbl.cytobands)$gieStain<-NULL
  tbl.cytobands$name<-as.character(tbl.cytobands$name)
  #pad names
  tbl.cytobands$name[which(tbl.cytobands$name=="")]="-"
  tbl.cytobands$name<-str_remove(tbl.cytobands$name,"[0-9.]+")
  matches<-GenomicRanges::findOverlaps(tbl.cytobands,tiles,select="all",ignore.strand=TRUE)
  #matches
  a3<-tiles[S4Vectors::subjectHits(matches)]
  # S4Vectors::subjectHits(matches)
  #sort matches
  matches2 = matches[order(S4Vectors::queryHits(matches))]
  #tbl.cytobands
  S4Vectors::mcols(a3) <- cbind.data.frame(
    S4Vectors::mcols(a3),
    S4Vectors::mcols(tbl.cytobands[S4Vectors::queryHits(matches)]))
  #label Cytobands
  #https://web.mit.edu/~r/current/arch/i386_linux26/lib/R/library/GenomicRanges/html/GRanges-class.html
  #ALL CLEAN :)
  #S4Vectors::mcols(empties)<-cbind.data.frame(S4Vectors::mcols(empties),name="")
  #empties
  empties<-tiles[-unique(S4Vectors::queryHits(matches))]
  S4Vectors::mcols(empties)<-cbind.data.frame(S4Vectors::mcols(empties),name="p")
  a3<-GenomicRanges::sort(append(GenomicRanges::GRanges(a3),empties))
  a3<-unique(a3)
  #add empties to zero (cpgNum = 0)
  #a3<-GenomicRanges::sort(append(GenomicRanges::GRanges(a3),empties))
  
  #overlap stuff again to remove doubles
  
  #label CPGs
  tbl.cpgIslandExt<-GenomicRanges::sort(tbl.cpgIslandExt)
  matches<-GenomicRanges::findOverlaps(tiles,tbl.cpgIslandExt,select="all",ignore.strand=TRUE)
  #sum up the stuff
  #https://www.rdocumentation.org/packages/S4Vectors/versions/0.4.0/topics/Hits-class
  b1<-tbl.cpgIslandExt[S4Vectors::subjectHits(matches)]
  
  b1_df <- cbind.data.frame(
    data.frame(S4Vectors::mcols(b1)),
    data.frame(tiles=IRanges::ranges(tiles[S4Vectors::queryHits(matches)])),
    chromMatch=GenomeInfoDb::seqnames(tiles[S4Vectors::queryHits(matches)]))
  #NEED TO FIND TILES THAT AREN'T MATCHES
  length(unique(S4Vectors::queryHits(matches)))
  b1t<-as_tibble(b1_df)
  #subjectHits vs queryHits
  b1t<-b1t %>% mutate(interval=str_c(chromMatch,tiles.start,tiles.end,sep="-"))
  b1t$interval
  b1t2<-b1t %>% group_by(interval) %>% summarise_at(dplyr::vars(cpgNum),sum) %>% separate(interval,into=c("chrom","start","end"),sep="-")
  #missing a few here
  #GET EMPTIES
  empties<-tiles[-unique(S4Vectors::queryHits(matches))]
  S4Vectors::mcols(empties)<-cbind.data.frame(S4Vectors::mcols(empties),cpgNum=0)
  #add empties to zero (cpgNum = 0)
  
  cpg_densities<-GenomicRanges::sort(append(GenomicRanges::GRanges(b1t2),empties))
  #TODO: fill with missing zeros
  write.table(as_tibble(GenomicRanges::sort(cpg_densities))[,c(1:3,6)],str_c(fileSuffixes,"_cpg_densities.tsv",sep=""),quote=FALSE,sep="\t",col.names = FALSE,row.names=FALSE)
  write.table(as_tibble(a3)[,c(1:3,6)],str_c(fileSuffixes,"_cytoband_densities_granges.tsv",sep=""),quote=FALSE,sep="\t",col.names = FALSE,row.names=FALSE)
  write.table(as_tibble(chrom_sizes),str_c(outputDir,"/",genomeText,"_chrom_sizes.tsv",sep=""),quote=FALSE,sep="\t",col.names = FALSE,row.names=FALSE)

  invisible(c(
    str_c(fileSuffixes, "_cpg_densities.tsv", sep = ""),
    str_c(fileSuffixes, "_cytoband_densities_granges.tsv", sep = ""),
    str_c(outputDir, "/", genomeText, "_chrom_sizes.tsv", sep = "")
  ))
}
