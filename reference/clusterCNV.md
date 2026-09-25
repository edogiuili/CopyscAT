# Merge and renumber adjacent copy number clusters

Post-processes
[`identifyCNVClusters()`](https://edogiuili.github.io/CopyscAT/reference/identifyCNVClusters.md)
output: orders clusters by mean, merges pairs separated by less than
`minDiff`, and renumbers so that cluster indices are contiguous.

## Usage

``` r
clusterCNV(initialResultList, medianIQR, maxClust = 4, minDiff = 0.25)
```

## Arguments

- initialResultList:

  Output of
  [`identifyCNVClusters()`](https://edogiuili.github.io/CopyscAT/reference/identifyCNVClusters.md).

- medianIQR:

  Centres and IQRs, normally `initialResultList[[3]]`.

- maxClust:

  Maximum clusters used in
  [`identifyCNVClusters()`](https://edogiuili.github.io/CopyscAT/reference/identifyCNVClusters.md).

- minDiff:

  Minimum Z-score separation between retained clusters.

## Value

A list of two elements: cleaned per-cell assignments and cleaned per-arm
cluster means.

## See also

Other CNV calling:
[`getAlteredSegments()`](https://edogiuili.github.io/CopyscAT/reference/getAlteredSegments.md),
[`identifyCNVClusters()`](https://edogiuili.github.io/CopyscAT/reference/identifyCNVClusters.md),
[`makeFakeCells()`](https://edogiuili.github.io/CopyscAT/reference/makeFakeCells.md)

## Examples

``` r
if (FALSE) { # \dontrun{
candidate_cnvs_clean <- clusterCNV(
  candidate_cnvs, medianIQR = candidate_cnvs[[3]], minDiff = 1.5
)
} # }
```
