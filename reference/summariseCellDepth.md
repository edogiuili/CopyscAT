# Summarise per-cell depth by sample

Summarise per-cell depth by sample

## Usage

``` r
summariseCellDepth(cellData, valueColumn = NULL)
```

## Arguments

- cellData:

  Output of
  [`signalPerCell()`](https://edogiuili.github.io/CopyscAT/reference/signalPerCell.md)
  or
  [`fragmentsPerCell()`](https://edogiuili.github.io/CopyscAT/reference/fragmentsPerCell.md).

- valueColumn:

  Column to summarise. Defaults to whichever of `signal_bp` or
  `n_fragments` is present.

## Value

A data frame with one row per sample: cell count, minimum, quartiles,
median and maximum.

## See also

Other quality control:
[`filterCells()`](https://edogiuili.github.io/CopyscAT/reference/filterCells.md),
[`fragmentsPerCell()`](https://edogiuili.github.io/CopyscAT/reference/fragmentsPerCell.md),
[`plotCellDistribution()`](https://edogiuili.github.io/CopyscAT/reference/plotCellDistribution.md),
[`plotCellKnee()`](https://edogiuili.github.io/CopyscAT/reference/plotCellKnee.md),
[`signalPerCell()`](https://edogiuili.github.io/CopyscAT/reference/signalPerCell.md)

## Examples

``` r
data(scDataSamp)
summariseCellDepth(signalPerCell(list(demo = scDataSamp)))
#>   sample cells    min     q25  median      q75      max
#> 1   demo  1371 559574 2952958 7806915 13142928 98118506
```
