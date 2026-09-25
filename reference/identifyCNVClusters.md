# Call copy number states by Gaussian mixture decomposition

Fits a
[`mclust::Mclust()`](https://mclust-org.github.io/mclust/reference/Mclust.html)
mixture per chromosome arm and assigns each cell to a component.
Components that are too close, too small or too uncertain are merged or
discarded. Optionally adds simulated pseudodiploid cells
([`makeFakeCells()`](https://edogiuili.github.io/CopyscAT/reference/makeFakeCells.md))
so that a baseline exists in highly aneuploid samples.

## Usage

``` r
identifyCNVClusters(
  inputMatrix,
  median_iqr,
  useDummyCells = TRUE,
  propDummy = 0.25,
  deltaMean = 0.03,
  subsetSize = 500,
  fakeCellSD = 0.1,
  summaryFunction = cutAverage,
  minDiff = 0.25,
  deltaBIC2 = 0.25,
  bicMinimum = 0.1,
  minMix = 0.3,
  uncertaintyCutoff = 0.55,
  mergeCutoff = 3,
  summarySuffix = "",
  shrinky = 0,
  maxClust = 4,
  IQRCutoff = 0.25,
  medianQuantileCutoff = 0.3,
  normalCells = NULL,
  verbose = FALSE
)
```

## Arguments

- inputMatrix:

  Collapsed, filtered matrix.

- median_iqr:

  Centres and IQRs from
  [`computeCenters()`](https://edogiuili.github.io/CopyscAT/reference/computeCenters.md).

- useDummyCells:

  Add simulated pseudodiploid cells.

- propDummy:

  Simulated cells as a proportion of real cells.

- deltaMean:

  Minimum separation between component means.

- subsetSize:

  Cells sampled for mixture initialisation. Must not exceed the number
  of available cells.

- fakeCellSD:

  Standard deviation of simulated cells, as a fraction of the signal
  range. `0.10`-`0.20` works well.

- summaryFunction:

  Summary function for arms. Default
  [`cutAverage()`](https://edogiuili.github.io/CopyscAT/reference/cutAverage.md).

- minDiff:

  Minimum Z-score gap between retained clusters.

- deltaBIC2:

  Minimum BIC gain to accept more than two components.

- bicMinimum:

  Minimum BIC gain to accept two components over one.

- minMix:

  Minimum mixing proportion for the smaller component.

- uncertaintyCutoff:

  Cells above this classification uncertainty are left unassigned
  (cluster `0`).

- mergeCutoff:

  Mean difference below which wide components are split rather than
  merged.

- summarySuffix:

  Suffix for the diagnostic PDF.

- shrinky:

  Shrinkage passed to
  [`mclust::priorControl()`](https://mclust-org.github.io/mclust/reference/priorControl.html).

- maxClust:

  Maximum number of components per arm.

- IQRCutoff:

  Quantile of per-arm IQRs used as the common scale.

- medianQuantileCutoff:

  Quantile of per-arm centres used as the common centre. Set to `-1` to
  centre on `normalCells` instead.

- normalCells:

  Barcodes of known normal cells. Required when
  `medianQuantileCutoff = -1`.

- verbose:

  Print the per-arm mixture summary to the console. The classification
  plots are written to the PDF either way.

## Value

A list of three elements: per-cell cluster assignments, per-arm cluster
means, and the centres/IQRs actually used.

## Details

Writes a multi-page PDF of per-arm classification plots to the session
output path.

## See also

Other CNV calling:
[`clusterCNV()`](https://edogiuili.github.io/CopyscAT/reference/clusterCNV.md),
[`getAlteredSegments()`](https://edogiuili.github.io/CopyscAT/reference/getAlteredSegments.md),
[`makeFakeCells()`](https://edogiuili.github.io/CopyscAT/reference/makeFakeCells.md)

## Examples

``` r
if (FALSE) { # \dontrun{
candidate_cnvs <- identifyCNVClusters(
  scData_collapse, median_iqr,
  useDummyCells = TRUE, subsetSize = 600, maxClust = 4
)
} # }
```
