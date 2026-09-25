# Plot the per-cell depth distribution

Draws one histogram per sample on a shared log scale, which makes the
split between real cells and ambient barcodes easy to see. A bimodal
shape is normal: the left mode is empty or ambient barcodes, the right
mode cells.

## Usage

``` r
plotCellDistribution(
  cellData,
  valueColumn = NULL,
  cutoff = NULL,
  bins = 60,
  logScale = TRUE
)
```

## Arguments

- cellData:

  Output of
  [`signalPerCell()`](https://edogiuili.github.io/CopyscAT/reference/signalPerCell.md)
  or
  [`fragmentsPerCell()`](https://edogiuili.github.io/CopyscAT/reference/fragmentsPerCell.md).

- valueColumn:

  Column to plot. Defaults to whichever of `signal_bp` or `n_fragments`
  is present.

- cutoff:

  Optional threshold to mark with a dashed line, for example the
  `minFrags` you intend to use.

- bins:

  Number of histogram bins.

- logScale:

  Use a log10 x axis. Default `TRUE`.

## Value

A
[`ggplot2::ggplot()`](https://ggplot2.tidyverse.org/reference/ggplot.html)
object.

## See also

Other quality control:
[`filterCells()`](https://edogiuili.github.io/CopyscAT/reference/filterCells.md),
[`fragmentsPerCell()`](https://edogiuili.github.io/CopyscAT/reference/fragmentsPerCell.md),
[`plotCellKnee()`](https://edogiuili.github.io/CopyscAT/reference/plotCellKnee.md),
[`signalPerCell()`](https://edogiuili.github.io/CopyscAT/reference/signalPerCell.md),
[`summariseCellDepth()`](https://edogiuili.github.io/CopyscAT/reference/summariseCellDepth.md)

## Examples

``` r
data(scDataSamp)
plotCellDistribution(signalPerCell(list(demo = scDataSamp)))
```
