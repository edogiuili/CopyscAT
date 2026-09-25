# Build CopyscAT reference files for a genome

Tiles a BSgenome object into fixed-width bins and writes the three
reference tables that
[`initialiseEnvironment()`](https://edogiuili.github.io/CopyscAT/reference/initialiseEnvironment.md)
expects: CpG density, cytoband arms and chromosome sizes. Cytoband and
CpG tracks are fetched from UCSC, so this needs network access. Run once
per genome and bin size.

## Usage

``` r
generateReferences(
  genomeObject,
  genomeText = "hg38",
  tileWidth = 1e+06,
  outputDir = "~"
)
```

## Arguments

- genomeObject:

  A BSgenome object, for example `BSgenome.Hsapiens.UCSC.hg38`.

- genomeText:

  UCSC genome name, for example `"hg38"`.

- tileWidth:

  Bin width in base pairs. Must match the `binSize` later passed to
  [`initialiseEnvironment()`](https://edogiuili.github.io/CopyscAT/reference/initialiseEnvironment.md).

- outputDir:

  Directory for the three output files.

## Value

Invisibly a character vector of the files written.

## See also

Other session setup:
[`initialiseEnvironment()`](https://edogiuili.github.io/CopyscAT/reference/initialiseEnvironment.md),
[`initializeStandards()`](https://edogiuili.github.io/CopyscAT/reference/initializeStandards.md),
[`setOutputFile()`](https://edogiuili.github.io/CopyscAT/reference/setOutputFile.md)

## Examples

``` r
if (FALSE) { # \dontrun{
library(BSgenome.Hsapiens.UCSC.hg38)
generateReferences(
  BSgenome.Hsapiens.UCSC.hg38, genomeText = "hg38",
  tileWidth = 1e6, outputDir = "."
)
} # }
```
