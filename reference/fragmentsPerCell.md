# Fragment counts per cell from fragment files

Counts the fragments recorded for each barcode. Unlike
[`signalPerCell()`](https://edogiuili.github.io/CopyscAT/reference/signalPerCell.md),
which reports base pairs covered, this is the number of fragments and so
is directly comparable with the `minFrags` argument of
[`initialiseEnvironment()`](https://edogiuili.github.io/CopyscAT/reference/initialiseEnvironment.md)
and the `-f` option of `process_fragment_file.py`.

## Usage

``` r
fragmentsPerCell(files, minFragments = 0)
```

## Arguments

- files:

  Named character vector of paths. Names are used as sample labels. A
  path ending in `.csv` is read as a CellRanger `singlecell.csv`;
  anything else is treated as a fragment file.

- minFragments:

  Drop barcodes with fewer fragments than this. Default `0`, which keeps
  everything.

## Value

A data frame with columns `sample`, `barcode` and `n_fragments`.

## Details

Reads either a 10x `fragments.tsv.gz` (barcode in the fourth column, `#`
comment lines skipped) or a CellRanger `singlecell.csv`. For a
`singlecell.csv` the `passed_filters` column is used, which counts only
fragments passing CellRanger's own filters and so will not exactly match
a direct count of the fragment file.

## See also

[`signalPerCell()`](https://edogiuili.github.io/CopyscAT/reference/signalPerCell.md)
for the base-pair signal already in a CopyscAT matrix.

Other quality control:
[`filterCells()`](https://edogiuili.github.io/CopyscAT/reference/filterCells.md),
[`plotCellDistribution()`](https://edogiuili.github.io/CopyscAT/reference/plotCellDistribution.md),
[`plotCellKnee()`](https://edogiuili.github.io/CopyscAT/reference/plotCellKnee.md),
[`signalPerCell()`](https://edogiuili.github.io/CopyscAT/reference/signalPerCell.md),
[`summariseCellDepth()`](https://edogiuili.github.io/CopyscAT/reference/summariseCellDepth.md)

## Examples

``` r
if (FALSE) { # \dontrun{
frags <- fragmentsPerCell(
  c(tumour1 = "tumour1/fragments.tsv.gz"),
  minFragments = 1000
)
plotCellDistribution(frags, "n_fragments", cutoff = 1e4)
} # }
```
