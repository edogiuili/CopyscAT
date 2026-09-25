# Plot the distribution of collapsed signal

Writes a violin plot of per-arm signal across cells to the session
output path. Useful as a QC check before and after scaling.

## Usage

``` r
graphCNVDistribution(inputMatrix, outputSuffix = "_all")
```

## Arguments

- inputMatrix:

  Collapsed matrix from
  [`collapseChrom3N()`](https://edogiuili.github.io/CopyscAT/reference/collapseChrom3N.md).

- outputSuffix:

  Suffix for the PDF file name.

## Value

Invisibly `NULL`. Called for the file it writes.

## See also

Other reporting:
[`saveSummaryStats()`](https://edogiuili.github.io/CopyscAT/reference/saveSummaryStats.md)

## Examples

``` r
if (FALSE) { # \dontrun{
graphCNVDistribution(scData_collapse, outputSuffix = "violins")
} # }
```
