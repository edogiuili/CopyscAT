# Separate neoplastic from non-neoplastic cells

Factorises the collapsed arm matrix with NMF, clusters cells on the
factor loadings, and treats the cluster with the lowest variance across
arms as the non-neoplastic population. Use the returned barcodes as the
baseline for
[`identifyCNVClusters()`](https://edogiuili.github.io/CopyscAT/reference/identifyCNVClusters.md)
or
[`annotateCNV4B()`](https://edogiuili.github.io/CopyscAT/reference/annotateCNV4B.md).

## Usage

``` r
identifyNonNeoplastic(
  inputMatrix,
  estimatedCellularity = 0.8,
  nmfComponents = 5,
  outputHeatmap = TRUE,
  cutHeight = 0.6,
  methodHclust = "ward.D"
)
```

## Arguments

- inputMatrix:

  Collapsed, filtered matrix from
  [`filterCells()`](https://edogiuili.github.io/CopyscAT/reference/filterCells.md).

- estimatedCellularity:

  Expected tumour cellularity. Works best below about 0.9.

- nmfComponents:

  Number of NMF components.

- outputHeatmap:

  Write a heatmap and violin plot to the output path.

- cutHeight:

  Dendrogram cut height, as a fraction of maximum height.

- methodHclust:

  Linkage method passed to
  [`fastcluster::hclust()`](https://rdrr.io/pkg/fastcluster/man/hclust.html).

## Value

A list with `cellAssigns` (cluster per barcode), `normalBarcodes` and
`clusterNormal` (the index of the non-neoplastic cluster).

## Details

Inspect the heatmap it writes and adjust `cutHeight` and `nmfComponents`
when the clusters do not separate cleanly.

## See also

Other neoplastic classification:
[`estimateCellCycleFraction()`](https://edogiuili.github.io/CopyscAT/reference/estimateCellCycleFraction.md)

## Examples

``` r
if (FALSE) { # \dontrun{
nmf_results <- identifyNonNeoplastic(
  scData_collapse, methodHclust = "ward.D", cutHeight = 0.4
)
} # }
```
