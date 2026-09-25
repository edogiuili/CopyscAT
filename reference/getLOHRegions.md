# Call regions of loss of heterozygosity

Runs changepoint analysis per chromosome on the normalised matrix to
find contiguous stretches of reduced signal, then keeps those recurring
in at least `lossCutoffCells` cells.

## Usage

``` r
getLOHRegions(
  inputMatrixIn,
  lossCutoff = (-0.25),
  uncertaintyCutLoss = 0.5,
  diffThreshold = 0.9,
  minLength = 3e+06,
  minSeg = 3,
  lossCutoffCells = 100,
  targetFun = IQR,
  quantileLimit = 0.3,
  cpgCutoff = 0,
  meanThreshold = 4,
  dummyQuantile = 0.5,
  dummyPercentile = 0.2,
  dummySd = 0.2
)
```

## Arguments

- inputMatrixIn:

  Normalised matrix from
  [`normalizeMatrixN()`](https://edogiuili.github.io/CopyscAT/reference/normalizeMatrixN.md).

- lossCutoff:

  Maximum Z-score for a bin to count as lost.

- uncertaintyCutLoss:

  Maximum mixture uncertainty for a cell call.

- diffThreshold:

  Difference threshold between adjacent bins.

- minLength:

  Minimum region length in base pairs.

- minSeg:

  Minimum number of bins per changepoint segment.

- lossCutoffCells:

  Minimum cells carrying a region to keep it.

- targetFun:

  Spread function used to score regions. Default
  [`stats::IQR()`](https://rdrr.io/r/stats/IQR.html).

- quantileLimit:

  Quantile of the per-cell signal used to call losses.

- cpgCutoff:

  Drop bins with CpG density at or below this.

- meanThreshold:

  Minimum mean signal for a bin to be considered.

- dummyQuantile, dummyPercentile, dummySd:

  Parameters of the simulated diploid baseline used to score candidate
  regions.

## Value

A list whose second element holds the called LOH intervals, in the form
consumed by
[`annotateLosses()`](https://edogiuili.github.io/CopyscAT/reference/annotateLosses.md).

## See also

Other LOH:
[`annotateLosses()`](https://edogiuili.github.io/CopyscAT/reference/annotateLosses.md)
