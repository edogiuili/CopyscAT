# Calling copy number variants from single-cell ATAC data

CopyscAT infers copy number alterations from single-cell ATAC sequencing
without a matched normal control. It bins fragment counts across the
genome, corrects the accessibility bias using CpG density, and
decomposes the per-chromosome-arm signal into copy number states with
Gaussian mixture models.

The chunks below are not evaluated when the vignette is built, because a
full run needs reference files and a fragment matrix that are too large
to ship.

## Installation

``` r

# install.packages("remotes")
remotes::install_github("edogiuili/CopyscAT")
```

Some dependencies come from Bioconductor:

``` r

# install.packages("BiocManager")
BiocManager::install(c("edgeR", "biomaRt", "rtracklayer", "GenomicRanges"))
```

## Preparing input

Fragment matrices are produced from a `fragments.tsv.gz` file with the
`process_fragment_file.py` script shipped with the package:

``` sh
python3 process_fragment_file.py -i fragments.tsv.gz -o sample_matrix.tsv \
    -b 1000000 -f 1000 -g hg38_chrom_sizes.tsv
```

## One-time reference generation

Reference files depend on the genome build and the bin size. Generate
them once and reuse them:

``` r

library(CopyscAT)
library(BSgenome.Hsapiens.UCSC.hg38)

generateReferences(
  BSgenome.Hsapiens.UCSC.hg38,
  genomeText = "hg38",
  tileWidth  = 1e6,
  outputDir  = "."
)
```

This writes a chromosome-size table, a binned cytoband table and a
binned CpG density table. The package also ships prebuilt references for
hg19 and hg38 at 1 Mb resolution in `hg19_references/` and
`hg38_references/`.

## Initialising a session

CopyscAT keeps references and parameters in a session environment, so
this must be run once per R session. `binSize` must match the
`tileWidth` used above.

``` r

initialiseEnvironment(
  genomeFile   = "hg38_chrom_sizes.tsv",
  cytobandFile = "hg38_1e+06_cytoband_densities_granges.tsv",
  cpgFile      = "hg38_1e+06_cpg_densities.tsv",
  binSize      = 1e6,
  minFrags     = 1e4,
  cellSuffix   = c("-1"),
  lowerTrim    = 0.5,
  upperTrim    = 0.8
)

setOutputFile(".", "samp_dataset")
```

`cellSuffix` must match the barcode suffix in your matrix, otherwise the
column selection throughout the pipeline silently matches nothing.

## Checking depth before you start

Sequencing depth drives most of the parameter choices below, so it is
worth looking at it before committing to thresholds — especially when
several samples will be analysed together.

``` r

depth <- signalPerCell(c(
  tumour1 = "tumour1_matrix.tsv",
  tumour2 = "tumour2_matrix.tsv"
))

summariseCellDepth(depth)
plotCellDistribution(depth, cutoff = 1e6)
plotCellKnee(depth)
```

[`plotCellDistribution()`](https://edogiuili.github.io/CopyscAT/reference/plotCellDistribution.md)
draws one histogram per sample on a log scale. A bimodal shape is
expected: the left mode is ambient or empty barcodes, the right mode
real cells.
[`plotCellKnee()`](https://edogiuili.github.io/CopyscAT/reference/plotCellKnee.md)
overlays the samples as ranked curves, whose elbow marks the transition
between the two. A sample whose curve sits well below the others is
shallower and may need its own `minFrags` and `blacklistCutoff`.

Note the units. `process_fragment_file.py` accumulates fragment
*lengths*, so the row sum of a CopyscAT matrix is the number of base
pairs a cell covers, not the number of fragments — which is why the
column is called `signal_bp`. For fragment counts, read the source files
instead:

``` r

frags <- fragmentsPerCell(c(tumour1 = "tumour1/fragments.tsv.gz"))

# Or, much faster, from CellRanger output:
frags <- fragmentsPerCell(c(tumour1 = "tumour1/singlecell.csv"))

plotCellDistribution(frags, cutoff = 1e4)
```

These counts are directly comparable with the `minFrags` argument of
[`initialiseEnvironment()`](https://edogiuili.github.io/CopyscAT/reference/initialiseEnvironment.md)
and the `-f` option of `process_fragment_file.py`.

## Normalisation and quality control

``` r

scData <- readInputTable("sample_matrix.tsv")

scData_k_norm <- normalizeMatrixN(
  scData,
  logNorm             = FALSE,
  maxZero             = 2000,
  imputeZeros         = FALSE,
  blacklistProp       = 0.8,
  blacklistCutoff     = 125,
  dividingFactor      = 1,
  upperFilterQuantile = 0.95
)

scData_collapse <- collapseChrom3N(
  scData_k_norm,
  summaryFunction   = cutAverage,
  binExpand         = 1,
  minimumChromValue = 100,
  logTrans          = FALSE,
  tssEnrich         = 1,
  logBase           = 2,
  minCPG            = 300,
  powVal            = 0.73
)

scData_collapse <- filterCells(
  scData_collapse, minimumSegments = 40, minDensity = 0.1
)

graphCNVDistribution(scData_collapse, outputSuffix = "violins")
median_iqr <- computeCenters(scData_collapse, summaryFunction = cutAverage)
```

## Calling chromosome-arm CNVs

### Using all cells as the baseline

Suitable when a substantial diploid population is present.

``` r

candidate_cnvs <- identifyCNVClusters(
  scData_collapse, median_iqr,
  useDummyCells        = TRUE,
  propDummy            = 0.25,
  minMix               = 0.01,
  deltaMean            = 0.03,
  deltaBIC2            = 0.25,
  bicMinimum           = 0.1,
  subsetSize           = 600,
  fakeCellSD           = 0.08,
  uncertaintyCutoff    = 0.55,
  summaryFunction      = cutAverage,
  maxClust             = 4,
  mergeCutoff          = 3,
  IQRCutoff            = 0.2,
  medianQuantileCutoff = 0.4
)

candidate_cnvs_clean <- clusterCNV(
  initialResultList = candidate_cnvs,
  medianIQR         = candidate_cnvs[[3]],
  minDiff           = 1.5
)

final_cnv_list <- annotateCNV4(
  candidate_cnvs_clean,
  saveOutput    = TRUE,
  outputSuffix  = "clean_cnv",
  sdCNV         = 0.5,
  filterResults = TRUE,
  filterRange   = 0.8
)
```

`subsetSize` must not exceed the number of cells that survive filtering;
it is clamped with a warning if it does.

### Using inferred non-neoplastic cells as the baseline

Use this when the calls from the first approach look mis-centred. It
works best below about 90% tumour cellularity.

``` r

nmf_results <- identifyNonNeoplastic(
  scData_collapse, methodHclust = "ward.D", cutHeight = 0.4
)

print(paste("Normal cluster is:", nmf_results$clusterNormal))

write.table(
  rownames_to_column(data.frame(nmf_results$cellAssigns), var = "Barcode"),
  file = str_c(scCNVCaller$locPrefix, scCNVCaller$outPrefix,
               "_nmf_clusters.csv"),
  quote = FALSE, row.names = FALSE, sep = ","
)

median_iqr <- computeCenters(
  scData_collapse %>% dplyr::select(chrom, nmf_results$normalBarcodes),
  summaryFunction = cutAverage
)

candidate_cnvs <- identifyCNVClusters(
  scData_collapse, median_iqr,
  useDummyCells        = TRUE,
  medianQuantileCutoff = -1,
  normalCells          = nmf_results$normalBarcodes,
  subsetSize           = 800,
  fakeCellSD           = 0.09,
  uncertaintyCutoff    = 0.65,
  maxClust             = 4,
  IQRCutoff            = 0.25
)
```

Setting `medianQuantileCutoff = -1` requires `normalCells`; the call
errors if it is missing.

To anchor copy number 2 to the normal population explicitly, use
[`annotateCNV4B()`](https://edogiuili.github.io/CopyscAT/reference/annotateCNV4B.md).
It is more accurate with a reliable normal set and more prone to false
positives on noisy data.

``` r

final_cnv_list <- annotateCNV4B(
  candidate_cnvs_clean, nmf_results$normalBarcodes,
  saveOutput = TRUE, outputSuffix = "clean_cnv_b2",
  sdCNV = 0.6, filterResults = TRUE, filterRange = 0.4,
  minAlteredCellProp = 0.5
)
```

## Smoothing calls across clusters

``` r

smoothedCNVList <- smoothClusters(
  scDataSampClusters,
  inputCNVList    = final_cnv_list[[3]],
  percentPositive = 0.4,
  removeEmpty     = FALSE
)
```

## Focal amplifications

This step is slow. Raising `minThreshold` skips quiet cells and speeds
it up at the cost of sensitivity.

``` r

dm_candidates <- identifyDoubleMinutes(
  scData_k_norm, minCells = 100, qualityCutoff2 = 100, minThreshold = 4
)

write.table(
  dm_candidates,
  file = str_c(scCNVCaller$locPrefix, scCNVCaller$outPrefix, "_dm.csv"),
  quote = FALSE, row.names = FALSE, sep = ","
)
```

It returns `NULL` and reports it when nothing passes the filters.

Byte-compiling the function can shave time off a large run:

``` r

dmRead <- compiler::cmpfun(identifyDoubleMinutes)
dm_candidates <- dmRead(scData_k_norm, minCells = 100, minThreshold = 4)
```

## Segmental gains and losses

Requires the normal population from
[`identifyNonNeoplastic()`](https://edogiuili.github.io/CopyscAT/reference/identifyNonNeoplastic.md).

``` r

altered_segments <- getAlteredSegments(
  scData_k_norm, nmf_results, lossThreshold = 0.6
)

write.table(
  altered_segments,
  file = str_c(scCNVCaller$locPrefix, scCNVCaller$outPrefix, "_segments.csv"),
  quote = FALSE, row.names = FALSE, sep = ","
)
```

## Cycling cells

``` r

barcodeCycling <- estimateCellCycleFraction(
  scData, sampName = "sample", cutoff = 1000
)

# The largest component is taken as the 2N population.
write.table(
  barcodeCycling[order(names(barcodeCycling))] == max(barcodeCycling),
  file = str_c(scCNVCaller$locPrefix, scCNVCaller$outPrefix,
               "_cycling_cells.tsv"),
  sep = "\t", quote = FALSE, row.names = TRUE, col.names = FALSE
)
```

Chromosome X is the default because it is usually unaltered. Pick
another with `chromToUse` if X is affected in your samples.

## Session summary

``` r

saveSummaryStats()
```
