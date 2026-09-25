# Filter low-quality cells

Drops cells covering too few chromosome arms, then drops cells whose
total signal lies more than `signalSDcut` standard deviations from the
mean, which removes likely doublets and debris.

## Usage

``` r
filterCells(inputMatrix, minimumSegments = 40, minDensity = 0, signalSDcut = 2)
```

## Arguments

- inputMatrix:

  Collapsed matrix from
  [`collapseChrom3N()`](https://edogiuili.github.io/CopyscAT/reference/collapseChrom3N.md).

- minimumSegments:

  Minimum number of arms with signal above `minDensity` for a cell to be
  kept.

- minDensity:

  Signal above which an arm counts as covered.

- signalSDcut:

  Standard deviations from mean total signal beyond which a cell is
  treated as an outlier.

## Value

The input matrix with failing cell columns removed.

## Examples

``` r
if (FALSE) { # \dontrun{
scData_collapse <- filterCells(
  scData_collapse, minimumSegments = 40, minDensity = 0.1
)
} # }
```
