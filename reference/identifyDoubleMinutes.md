# Detect recurrent focal amplifications across all cells

Runs
[`getDoubleMinutes()`](https://edogiuili.github.io/CopyscAT/reference/getDoubleMinutes.md)
on every cell and chromosome, merges overlapping intervals, and keeps
only events recurring in enough cells. Slow on large datasets; raise
`minThreshold` to speed it up.

## Usage

``` r
identifyDoubleMinutes(
  inputMatrix,
  minCells = 100,
  qualityCutoff2 = 100,
  peakCutoff = 5,
  lossCutoff = -1,
  doPlots = FALSE,
  imageNumber = 1000,
  logTrans = FALSE,
  cpgTransform = FALSE,
  doLosses = FALSE,
  minThreshold = 4
)
```

## Arguments

- inputMatrix:

  Normalised matrix from
  [`normalizeMatrixN()`](https://edogiuili.github.io/CopyscAT/reference/normalizeMatrixN.md).

- minCells:

  Minimum cells carrying an event to keep it (first pass).

- qualityCutoff2:

  Minimum cells carrying an event to keep it (final pass).

- peakCutoff:

  Minimum segment Z-score to call an amplification.

- lossCutoff:

  Maximum segment Z-score to call a loss.

- doPlots:

  Write a QC plot every `imageNumber` cells.

- imageNumber:

  Plot interval, in cells.

- logTrans:

  Treat the input as log-transformed when scaling by CpG.

- cpgTransform:

  Divide signal by CpG density before calling.

- doLosses:

  Call losses instead of gains.

- minThreshold:

  Z-score below which a cell is skipped.

## Value

A data frame with one row per cell and one logical column per retained
event, or `NULL` when nothing passes the filters.

## See also

Other focal amplification:
[`getDoubleMinutes()`](https://edogiuili.github.io/CopyscAT/reference/getDoubleMinutes.md),
[`testOverlap()`](https://edogiuili.github.io/CopyscAT/reference/testOverlap.md)

## Examples

``` r
if (FALSE) { # \dontrun{
dm_candidates <- identifyDoubleMinutes(
  scData_k_norm, minCells = 100, qualityCutoff2 = 100, minThreshold = 4
)
} # }
```
