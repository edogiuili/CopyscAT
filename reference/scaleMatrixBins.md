# Merge adjacent genomic bins

Merge adjacent genomic bins

## Usage

``` r
scaleMatrixBins(inputMatrix, binSizeRatio, inColumn)
```

## Arguments

- inputMatrix:

  Matrix with `chrom` and a position column.

- binSizeRatio:

  Factor by which to widen bins. `2` doubles bin width.

- inColumn:

  Name of the position column to rebin.

## Value

The matrix regrouped on the widened bins, with numeric columns summed.

## See also

Other summarisation:
[`cutAverage()`](https://edogiuili.github.io/CopyscAT/reference/cutAverage.md)
