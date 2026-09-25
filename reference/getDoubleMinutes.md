# Detect focal amplifications in a single cell

Scales one cell's signal on one chromosome to a robust Z-score, then
runs changepoint analysis
([`changepoint::cpt.meanvar()`](https://rdrr.io/pkg/changepoint/man/cpt.meanvar.html))
to find segments whose mean exceeds `peakCutoff` (amplifications) or
falls below `lossCutoff` (losses). Used per cell by
[`identifyDoubleMinutes()`](https://edogiuili.github.io/CopyscAT/reference/identifyDoubleMinutes.md).

## Usage

``` r
getDoubleMinutes(
  inputMatrix,
  targetCell,
  doPlot = FALSE,
  penalty_type = "SIC",
  doLosses = FALSE,
  peakCutoff = 5,
  lossCutoff = -1,
  minThreshold = 4
)
```

## Arguments

- inputMatrix:

  Single-chromosome matrix with a `loc1` column of bin labels and one
  column per cell.

- targetCell:

  Column index of the cell to analyse.

- doPlot:

  Draw the changepoint fit. Default `FALSE`.

- penalty_type:

  Penalty passed to
  [`changepoint::cpt.meanvar()`](https://rdrr.io/pkg/changepoint/man/cpt.meanvar.html).

- doLosses:

  Return losses below `lossCutoff` instead of gains.

- peakCutoff:

  Minimum segment Z-score to call an amplification.

- lossCutoff:

  Maximum segment Z-score to call a loss.

- minThreshold:

  Skip the cell unless its maximum Z-score exceeds this. Raising it
  speeds up large datasets at the cost of sensitivity.

## Value

A character vector of interval labels (`chrom_start.chrom_end`), or
`NULL` when the cell has no qualifying segment.

## See also

Other focal amplification:
[`identifyDoubleMinutes()`](https://edogiuili.github.io/CopyscAT/reference/identifyDoubleMinutes.md),
[`testOverlap()`](https://edogiuili.github.io/CopyscAT/reference/testOverlap.md)
