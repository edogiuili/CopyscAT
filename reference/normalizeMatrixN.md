# Normalise a raw single-cell fragment matrix

First analysis step. Drops outlier cells above `upperFilterQuantile`,
flags uninformative ("blacklist") bins that are near-empty across most
cells, removes cells with too many empty bins, then scales to counts per
million with [`edgeR::cpm()`](https://rdrr.io/pkg/edgeR/man/cpm.html).

## Usage

``` r
normalizeMatrixN(
  inputMatrix,
  logNorm = FALSE,
  maxZero = 2000,
  imputeZeros = FALSE,
  blacklistProp = 0.8,
  priorCount = 1,
  blacklistCutoff = 100,
  dividingFactor = 1e+06,
  upperFilterQuantile = 0.95
)
```

## Arguments

- inputMatrix:

  Raw count matrix from
  [`readInputTable()`](https://edogiuili.github.io/CopyscAT/reference/readInputTable.md),
  cells in rows and bins in columns.

- logNorm:

  Log-transform during CPM normalisation. Default `FALSE`.

- maxZero:

  Maximum number of empty bins tolerated per cell, after uninformative
  bins are excluded.

- imputeZeros:

  Replace zero counts with the per-bin median. Off by default and not
  recommended.

- blacklistProp:

  Proportion of cells that must fall below `blacklistCutoff` for a bin
  to be treated as uninformative. Lower is stricter.

- priorCount:

  Prior count added before log normalisation. Only used when
  `logNorm = TRUE`.

- blacklistCutoff:

  Signal below which a bin counts as empty in a cell.

- dividingFactor:

  Deprecated. Scaling divisor applied after CPM; set to `1` to disable.

- upperFilterQuantile:

  Cells with total signal above this quantile are dropped as likely
  doublets.

## Value

A data frame with `chrom` and `pos` columns, a logical `blacklist`
column, a `raw_medians` column, and one column per surviving cell.

## See also

Other normalisation:
[`collapseChrom3N()`](https://edogiuili.github.io/CopyscAT/reference/collapseChrom3N.md),
[`computeCenters()`](https://edogiuili.github.io/CopyscAT/reference/computeCenters.md),
[`scaleMatrix()`](https://edogiuili.github.io/CopyscAT/reference/scaleMatrix.md)

## Examples

``` r
if (FALSE) { # \dontrun{
scData_k_norm <- normalizeMatrixN(
  scDataSamp,
  blacklistProp = 0.8, blacklistCutoff = 125, dividingFactor = 1
)
} # }
```
