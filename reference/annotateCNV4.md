# Estimate absolute copy number per cell

Converts per-arm cluster Z-scores into estimated absolute copy number,
centred so that the modal state is 2. Arms whose copy number range is
narrower than `filterRange`, or whose smallest altered group has fewer
than `minAlteredCells` cells, are dropped.

## Usage

``` r
annotateCNV4(
  cnvResults,
  saveOutput = TRUE,
  maxClust2 = 4,
  outputSuffix = "_1",
  sdCNV = 0.6,
  filterResults = TRUE,
  filterRange = 0.8,
  minAlteredCells = 40
)
```

## Arguments

- cnvResults:

  Output of
  [`clusterCNV()`](https://edogiuili.github.io/CopyscAT/reference/clusterCNV.md).

- saveOutput:

  Write results to the session output path.

- maxClust2:

  Maximum clusters used upstream.

- outputSuffix:

  Suffix for the output files.

- sdCNV:

  Standard deviation of the copy number prior.

- filterResults:

  Drop arms that look unaltered.

- filterRange:

  Minimum copy number range for an arm to be retained.

- minAlteredCells:

  Minimum cells in the smallest group of an arm.

## Value

A list of three elements: per-arm cluster copy numbers, the candidate
CNV table, and per-cell copy number estimates.

## See also

Other CNV annotation:
[`annotateCNV3()`](https://edogiuili.github.io/CopyscAT/reference/annotateCNV3.md),
[`annotateCNV4B()`](https://edogiuili.github.io/CopyscAT/reference/annotateCNV4B.md),
[`smoothClusters()`](https://edogiuili.github.io/CopyscAT/reference/smoothClusters.md)

## Examples

``` r
if (FALSE) { # \dontrun{
final_cnv_list <- annotateCNV4(
  candidate_cnvs_clean, saveOutput = TRUE,
  outputSuffix = "clean_cnv", sdCNV = 0.5, filterRange = 0.8
)
} # }
```
