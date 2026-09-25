# Annotate LOH regions with overlapping genes

Queries Ensembl through biomaRt for genes overlapping each called region
and writes the annotated table to the session output path. Requires
network access.

## Usage

``` r
annotateLosses(lossFile)
```

## Arguments

- lossFile:

  Output of
  [`getLOHRegions()`](https://edogiuili.github.io/CopyscAT/reference/getLOHRegions.md).

## Value

A tibble of regions with a comma-separated `genes` column.

## See also

Other LOH:
[`getLOHRegions()`](https://edogiuili.github.io/CopyscAT/reference/getLOHRegions.md)
