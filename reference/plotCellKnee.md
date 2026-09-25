# Plot ranked per-cell depth as a knee plot

Barcodes are ranked by depth within each sample and drawn on log-log
axes. The elbow marks the transition from cells to background, and
samples whose curve sits below the others are shallower and may need
different `minFrags` or `blacklistCutoff` values.

## Usage

``` r
plotCellKnee(cellData, valueColumn = NULL, cutoff = NULL)
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

  Optional threshold to mark with a dashed horizontal line.

## Value

A
[`ggplot2::ggplot()`](https://ggplot2.tidyverse.org/reference/ggplot.html)
object.

## See also

Other quality control:
[`filterCells()`](https://edogiuili.github.io/CopyscAT/reference/filterCells.md),
[`fragmentsPerCell()`](https://edogiuili.github.io/CopyscAT/reference/fragmentsPerCell.md),
[`plotCellDistribution()`](https://edogiuili.github.io/CopyscAT/reference/plotCellDistribution.md),
[`signalPerCell()`](https://edogiuili.github.io/CopyscAT/reference/signalPerCell.md),
[`summariseCellDepth()`](https://edogiuili.github.io/CopyscAT/reference/summariseCellDepth.md)

## Examples

``` r
data(scDataSamp)
plotCellKnee(signalPerCell(list(demo = scDataSamp)))
```
