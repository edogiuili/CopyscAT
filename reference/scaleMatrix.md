# Scale arm signal to Z-scores

Centres and scales each arm using the centres and IQRs from
[`computeCenters()`](https://edogiuili.github.io/CopyscAT/reference/computeCenters.md),
and drops `chrYp`.

## Usage

``` r
scaleMatrix(inputMatrix, median_iqr_list, spread = FALSE)
```

## Arguments

- inputMatrix:

  Collapsed, filtered matrix.

- median_iqr_list:

  List of centres and IQRs from
  [`computeCenters()`](https://edogiuili.github.io/CopyscAT/reference/computeCenters.md).

- spread:

  Return long format (`Cell`, `Density`) instead of wide.

## Value

The scaled matrix, wide or long depending on `spread`.

## See also

Other normalisation:
[`collapseChrom3N()`](https://edogiuili.github.io/CopyscAT/reference/collapseChrom3N.md),
[`computeCenters()`](https://edogiuili.github.io/CopyscAT/reference/computeCenters.md),
[`normalizeMatrixN()`](https://edogiuili.github.io/CopyscAT/reference/normalizeMatrixN.md)
