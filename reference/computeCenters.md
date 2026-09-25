# Compute per-arm centre and spread

Returns the summary value and interquartile range of each chromosome arm
across cells. Feed the result to
[`scaleMatrix()`](https://edogiuili.github.io/CopyscAT/reference/scaleMatrix.md)
and
[`identifyCNVClusters()`](https://edogiuili.github.io/CopyscAT/reference/identifyCNVClusters.md).

## Usage

``` r
computeCenters(inputMatrix, summaryFunction = cutAverage)
```

## Arguments

- inputMatrix:

  Collapsed, filtered matrix.

- summaryFunction:

  Function used for the centre. Default
  [`cutAverage()`](https://edogiuili.github.io/CopyscAT/reference/cutAverage.md).

## Value

A list of two data frames, each with columns `chrom` and `Value`: the
centres and the IQRs.

## See also

Other normalisation:
[`collapseChrom3N()`](https://edogiuili.github.io/CopyscAT/reference/collapseChrom3N.md),
[`normalizeMatrixN()`](https://edogiuili.github.io/CopyscAT/reference/normalizeMatrixN.md),
[`scaleMatrix()`](https://edogiuili.github.io/CopyscAT/reference/scaleMatrix.md)

## Examples

``` r
if (FALSE) { # \dontrun{
median_iqr <- computeCenters(scData_collapse, summaryFunction = cutAverage)
} # }
```
