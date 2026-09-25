# Per-cell signal from one or more CopyscAT matrices

Summarises how much signal each barcode carries, so that sequencing
depth can be compared across samples before running the CNV workflow.

## Usage

``` r
signalPerCell(matrices, sep = "\t")
```

## Arguments

- matrices:

  Either a named list or character vector of matrix file paths, or a
  named list of already-loaded data frames. Names are used as sample
  labels.

- sep:

  Field separator, when reading from file. Default tab.

## Value

A data frame with one row per cell and columns `sample`, `barcode`,
`signal_bp` (total base pairs covered) and `nonzero_bins` (bins with any
signal).

## Details

The values are **summed fragment lengths in base pairs**, not fragment
counts. `process_fragment_file.py` accumulates `end - start - 1` for
every fragment, so the row sum of a CopyscAT matrix is the total genomic
span covered by a cell. Use
[`fragmentsPerCell()`](https://edogiuili.github.io/CopyscAT/reference/fragmentsPerCell.md)
when you need actual fragment counts.

## See also

[`fragmentsPerCell()`](https://edogiuili.github.io/CopyscAT/reference/fragmentsPerCell.md)
for true fragment counts,
[`plotCellDistribution()`](https://edogiuili.github.io/CopyscAT/reference/plotCellDistribution.md)
and
[`plotCellKnee()`](https://edogiuili.github.io/CopyscAT/reference/plotCellKnee.md)
to visualise the result.

Other quality control:
[`filterCells()`](https://edogiuili.github.io/CopyscAT/reference/filterCells.md),
[`fragmentsPerCell()`](https://edogiuili.github.io/CopyscAT/reference/fragmentsPerCell.md),
[`plotCellDistribution()`](https://edogiuili.github.io/CopyscAT/reference/plotCellDistribution.md),
[`plotCellKnee()`](https://edogiuili.github.io/CopyscAT/reference/plotCellKnee.md),
[`summariseCellDepth()`](https://edogiuili.github.io/CopyscAT/reference/summariseCellDepth.md)

## Examples

``` r
data(scDataSamp)
signal <- signalPerCell(list(demo = scDataSamp))
summary(signal$signal_bp)
#>     Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
#>   559574  2952958  7806915  9518312 13142928 98118506 
```
