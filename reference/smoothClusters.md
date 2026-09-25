# Smooth copy number calls within cell clusters

Replaces per-cell copy number estimates with a per-cluster consensus,
which suppresses the noise inherent in single-cell ATAC coverage. A
cluster is called altered on an arm when more than `percentPositive` of
its cells are.

## Usage

``` r
smoothClusters(
  inputClusters,
  inputCNVList,
  inputCNVClusterFile = "",
  percentPositive = 0.5,
  removeEmpty = TRUE
)
```

## Arguments

- inputClusters:

  Two-column data frame of barcode and cluster label.

- inputCNVList:

  Per-cell copy number estimates, normally `final_cnv_list[[3]]` from
  [`annotateCNV4()`](https://edogiuili.github.io/CopyscAT/reference/annotateCNV4.md).

- inputCNVClusterFile:

  Optional path to a CSV of copy number calls, used instead of
  `inputCNVList`.

- percentPositive:

  Fraction of a cluster that must carry an alteration for the cluster to
  be called altered.

- removeEmpty:

  Drop arms where every cluster has the same call.

## Value

A data frame of barcodes with their cluster-level copy number calls.

## See also

Other CNV annotation:
[`annotateCNV3()`](https://edogiuili.github.io/CopyscAT/reference/annotateCNV3.md),
[`annotateCNV4()`](https://edogiuili.github.io/CopyscAT/reference/annotateCNV4.md),
[`annotateCNV4B()`](https://edogiuili.github.io/CopyscAT/reference/annotateCNV4B.md)

## Examples

``` r
if (FALSE) { # \dontrun{
smoothedCNVList <- smoothClusters(
  scDataSampClusters, inputCNVList = final_cnv_list[[3]],
  percentPositive = 0.4
)
} # }
```
