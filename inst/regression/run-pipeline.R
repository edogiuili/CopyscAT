# Runs the same pipeline against whichever CopyscAT is on .libPaths(),
# and saves every intermediate to an RDS for comparison.
args <- commandArgs(trailingOnly=TRUE)
lib <- args[1]; tag <- args[2]; outrds <- args[3]
.libPaths(c(lib, .libPaths()))
suppressPackageStartupMessages(library(CopyscAT))
# the original needs these attached; the new one doesn't, but attaching
# them keeps the two runs on identical footing.
suppressPackageStartupMessages({library(dplyr); library(tidyr); library(tibble)
  library(data.table); library(stringr); library(mclust); library(edgeR)})

cat(tag, "version:", as.character(packageVersion("CopyscAT")), "\n")

SEED <- 7
set.seed(SEED)

# ---- build references (same code for both runs) ----
dir <- file.path(tempdir(), paste0("refs_", tag)); dir.create(dir, showWarnings=FALSE)
chroms <- paste0("chr", c(1:6, "X", "Y")); binsPerChrom <- 25L; binSize <- 1e6
bins <- do.call(rbind, lapply(chroms, function(cc) data.frame(
  chrom=cc, start=seq(0,by=binSize,length.out=binsPerChrom),
  end=seq(binSize,by=binSize,length.out=binsPerChrom), stringsAsFactors=FALSE)))
arm <- rep("p", nrow(bins))
for (cc in chroms) { idx <- which(bins$chrom==cc); mid <- idx[ceiling(length(idx)/2)]
  arm[idx[idx<mid]] <- "p"; arm[mid] <- "cen"; arm[idx[idx>mid]] <- "q" }
bins$arm <- arm
bins$cpg <- 400L + (seq_len(nrow(bins)) %% 7L)*25L
gf <- file.path(dir,"cs.tsv"); cf <- file.path(dir,"cpg.tsv"); yf <- file.path(dir,"cyto.tsv")
write.table(data.frame(chrom=chroms,length=binsPerChrom*binSize), gf, sep="\t",
            quote=FALSE, col.names=FALSE, row.names=FALSE)
write.table(bins[,c("chrom","start","end","cpg")], cf, sep="\t", quote=FALSE,
            col.names=FALSE, row.names=FALSE)
write.table(bins[,c("chrom","start","end","arm")], yf, sep="\t", quote=FALSE,
            col.names=FALSE, row.names=FALSE)

out <- file.path(tempdir(), paste0("out_", tag)); dir.create(out, showWarnings=FALSE)
initialiseEnvironment(genomeFile=gf, cytobandFile=yf, cpgFile=cf, binSize=binSize,
                      minFrags=100, cellSuffix="-1", lowerTrim=0.5, upperTrim=0.8)
setOutputFile(out, "cmp")

# ---- identical simulated counts ----
set.seed(SEED)
nCells <- 150L
binNames <- paste(bins$chrom, bins$start, sep="_")
cellNames <- paste0("CELL", sprintf("%03d", seq_len(nCells)), "-1")
mat <- matrix(rpois(nCells*nrow(bins), lambda=500), nrow=nCells,
              dimnames=list(cellNames, binNames))
altered <- seq_len(nCells/2)
mat[altered, which(bins$chrom=="chr4")] <- round(mat[altered, which(bins$chrom=="chr4")]*2)
counts <- as.data.frame(mat)

res <- list()
res$norm <- normalizeMatrixN(counts, logNorm=FALSE, maxZero=2000, imputeZeros=FALSE,
                blacklistProp=0.8, blacklistCutoff=10, dividingFactor=1,
                upperFilterQuantile=0.99)
res$collapse <- collapseChrom3N(res$norm, summaryFunction=cutAverage, binExpand=1,
                minimumChromValue=0, logTrans=FALSE, tssEnrich=1, logBase=2,
                minCPG=100, powVal=0.73)
res$isMale <- scCNVCaller$isMale
res$filter <- filterCells(res$collapse, minimumSegments=1, minDensity=0, signalSDcut=3)
res$centers <- computeCenters(res$filter, summaryFunction=cutAverage)
res$scaled <- scaleMatrix(res$filter, res$centers, spread=FALSE)
set.seed(SEED)
res$cand <- identifyCNVClusters(res$filter, res$centers, useDummyCells=TRUE,
                propDummy=0.25, minMix=0.01, deltaMean=0.03, deltaBIC2=0.25,
                bicMinimum=0.1, subsetSize=80, fakeCellSD=0.08,
                uncertaintyCutoff=0.55, summaryFunction=cutAverage, maxClust=4,
                mergeCutoff=3, IQRCutoff=0.2, medianQuantileCutoff=0.4)
res$clean <- clusterCNV(initialResultList=res$cand, medianIQR=res$cand[[3]], minDiff=1.0)
res$final <- annotateCNV4(res$clean, saveOutput=FALSE, outputSuffix="_c",
                sdCNV=0.5, filterResults=TRUE, filterRange=0.8)
res$summary <- list(mean=scCNVCaller$meanReadsPerCell,
                    start=scCNVCaller$startingCellCount,
                    passing=scCNVCaller$cellsPassingFilter,
                    blacklist=scCNVCaller$blacklistCount,
                    finalFilter=scCNVCaller$finalFilterCells)
saveRDS(res, outrds)
cat(tag, "DONE\n")
