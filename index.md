# CopyscAT

Copy number variant inference from single-cell ATAC sequencing, without
a matched normal control.

CopyscAT bins fragment counts across the genome, corrects the chromatin
accessibility bias using CpG density, and decomposes the
per-chromosome-arm signal with Gaussian mixture models to assign copy
number states to individual cells. It also detects double minutes and
other focal amplifications by changepoint analysis, calls regions of
loss of heterozygosity, and separates neoplastic from non-neoplastic
cells by non-negative matrix factorisation.

This is a fork of [spcdot/CopyscAT](https://github.com/spcdot/CopyscAT).

## Installation

``` r

# install.packages("remotes")
remotes::install_github("edogiuili/CopyscAT")
```

Several dependencies come from Bioconductor and are best installed
first:

``` r

# install.packages("BiocManager")
BiocManager::install(c("edgeR", "biomaRt", "rtracklayer", "GenomicRanges"))
```

[`generateReferences()`](https://edogiuili.github.io/CopyscAT/reference/generateReferences.md)
additionally needs a `BSgenome` package for your genome, for example
`BSgenome.Hsapiens.UCSC.hg38`.

## Quick start

``` r

library(CopyscAT)

# Once per genome build and bin size.
library(BSgenome.Hsapiens.UCSC.hg38)
generateReferences(BSgenome.Hsapiens.UCSC.hg38, genomeText = "hg38",
                   tileWidth = 1e6, outputDir = ".")

# Once per R session.
initialiseEnvironment(
  genomeFile   = "hg38_chrom_sizes.tsv",
  cytobandFile = "hg38_1e+06_cytoband_densities_granges.tsv",
  cpgFile      = "hg38_1e+06_cpg_densities.tsv",
  binSize      = 1e6,
  minFrags     = 1e4,
  cellSuffix   = c("-1")
)
setOutputFile(".", "samp_dataset")

# Normalise, collapse to chromosome arms, filter.
scData          <- readInputTable("sample_matrix.tsv")
scData_k_norm   <- normalizeMatrixN(scData, blacklistCutoff = 125,
                                    dividingFactor = 1)
scData_collapse <- collapseChrom3N(scData_k_norm, minimumChromValue = 100,
                                   tssEnrich = 1, minCPG = 300)
scData_collapse <- filterCells(scData_collapse, minimumSegments = 40,
                               minDensity = 0.1)

# Call CNVs.
median_iqr           <- computeCenters(scData_collapse)
candidate_cnvs       <- identifyCNVClusters(scData_collapse, median_iqr)
candidate_cnvs_clean <- clusterCNV(candidate_cnvs, candidate_cnvs[[3]],
                                   minDiff = 1.5)
final_cnv_list       <- annotateCNV4(candidate_cnvs_clean, saveOutput = TRUE,
                                     outputSuffix = "clean_cnv")
```

See
[`vignette("copyscat")`](https://edogiuili.github.io/CopyscAT/articles/copyscat.md)
for the full workflow, including focal amplifications, LOH, and using
inferred normal cells as the baseline.

## Preparing input matrices

Fragment matrices are generated from a `fragments.tsv.gz` file:

``` sh
python3 process_fragment_file.py -i fragments.tsv.gz -o sample_matrix.tsv \
    -b 1000000 -f 1000 -g hg38_chrom_sizes.tsv
```

Run `python3 process_fragment_file.py --help` for the full parameter
list.

Prebuilt 1 Mb reference files for hg19 and hg38 are included in
`hg19_references/` and `hg38_references/`.

## Checking sequencing depth

Depth drives most of the parameter choices, so inspect it before setting
thresholds — particularly when comparing several samples:

``` r

depth <- signalPerCell(c(tumour1 = "t1_matrix.tsv", tumour2 = "t2_matrix.tsv"))

summariseCellDepth(depth)             # per-sample quartiles
plotCellDistribution(depth, cutoff = 1e6)   # histograms, one panel per sample
plotCellKnee(depth)                   # ranked barcode curves
```

[`signalPerCell()`](https://edogiuili.github.io/CopyscAT/reference/signalPerCell.md)
reports **base pairs**, not fragments: `process_fragment_file.py` sums
fragment lengths, so a matrix row sum is the genomic span a cell covers.
For fragment counts, read the source files:

``` r

frags <- fragmentsPerCell(c(tumour1 = "tumour1/fragments.tsv.gz"))
# or, much faster, from CellRanger:
frags <- fragmentsPerCell(c(tumour1 = "tumour1/singlecell.csv"))
```

These are directly comparable with `minFrags` and the `-f` option of
`process_fragment_file.py`.

## Notes on parameters

- `cellSuffix` must match the barcode suffix in your matrix (`-1`, `-2`,
  …). If it does not, column selection silently matches nothing.
- `binSize` must match the `tileWidth` used to build the reference
  files.
- `subsetSize` in
  [`identifyCNVClusters()`](https://edogiuili.github.io/CopyscAT/reference/identifyCNVClusters.md)
  is clamped to the available cell count, with a warning.
- [`identifyNonNeoplastic()`](https://edogiuili.github.io/CopyscAT/reference/identifyNonNeoplastic.md)
  works best below about 90% tumour cellularity.
- `maxZero` caps how many low-signal bins a cell may have, so a value
  above your bin count disables the filter and lets empty barcodes
  through to [`edgeR::cpm()`](https://rdrr.io/pkg/edgeR/man/cpm.html).

## Development

``` r

devtools::load_all(".")   # load without installing
devtools::test()          # run the test suite
devtools::document()      # regenerate NAMESPACE and man/
devtools::check()         # full R CMD check
```

## Citation

Nikolic A, et al. *Copy-scAT: Deconvoluting single-cell chromatin
accessibility of genetic subclones in cancer.* Science Advances (2021).
<https://doi.org/10.1126/sciadv.abg6045>

## License

GPL-3. Copyright (C) 2020 University of Calgary.
