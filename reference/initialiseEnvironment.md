# Load genome references and analysis parameters for a session

Reads the three reference tables produced by
[`generateReferences()`](https://edogiuili.github.io/CopyscAT/reference/generateReferences.md)
and stores them, together with the binning parameters, in the
[scCNVCaller](https://edogiuili.github.io/CopyscAT/reference/scCNVCaller.md)
session environment. Run this once per R session before any other step.

## Usage

``` r
initialiseEnvironment(
  genomeFile,
  cytobandFile,
  cpgFile,
  binSize = 1e+06,
  minFrags = 10000,
  cellSuffix = c("-1"),
  lowerTrim = 0.5,
  upperTrim = 0.8
)
```

## Arguments

- genomeFile:

  Path to a two-column chromosome-size table (chromosome, length), tab
  separated and without a header.

- cytobandFile:

  Path to the binned cytoband table produced by
  [`generateReferences()`](https://edogiuili.github.io/CopyscAT/reference/generateReferences.md).

- cpgFile:

  Path to the binned CpG density table produced by
  [`generateReferences()`](https://edogiuili.github.io/CopyscAT/reference/generateReferences.md).
  Four columns: chromosome, start, end, density.

- binSize:

  Genomic bin size in base pairs. Must match the `tileWidth` used to
  build the reference files. Default `1e6`.

- minFrags:

  Minimum number of fragments for a cell to be considered.

- cellSuffix:

  Character vector of barcode suffixes present in the input matrix, for
  example `c("-1", "-2")`.

- lowerTrim:

  Lower quantile for the trimmed mean used by
  [`cutAverage()`](https://edogiuili.github.io/CopyscAT/reference/cutAverage.md).

- upperTrim:

  Upper quantile for the trimmed mean used by
  [`cutAverage()`](https://edogiuili.github.io/CopyscAT/reference/cutAverage.md).

## Value

Invisibly `TRUE`. Called for the side effect of populating
[scCNVCaller](https://edogiuili.github.io/CopyscAT/reference/scCNVCaller.md).

## See also

Other session setup:
[`generateReferences()`](https://edogiuili.github.io/CopyscAT/reference/generateReferences.md),
[`initializeStandards()`](https://edogiuili.github.io/CopyscAT/reference/initializeStandards.md),
[`setOutputFile()`](https://edogiuili.github.io/CopyscAT/reference/setOutputFile.md)

## Examples

``` r
if (FALSE) { # \dontrun{
initialiseEnvironment(
  genomeFile   = "hg38_chrom_sizes.tsv",
  cytobandFile = "hg38_1e+06_cytoband_densities_granges.tsv",
  cpgFile      = "hg38_1e+06_cpg_densities.tsv",
  binSize      = 1e6,
  cellSuffix   = c("-1")
)
} # }
```
