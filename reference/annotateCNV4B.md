# Estimate absolute copy number against known normal cells

Like
[`annotateCNV4()`](https://edogiuili.github.io/CopyscAT/reference/annotateCNV4.md),
but anchors copy number 2 to a supplied set of non-neoplastic barcodes
rather than to the modal state. More accurate when the normal population
is reliable, and more prone to false positives when the data are noisy.

## Usage

``` r
annotateCNV4B(
  cnvResults,
  expectedNormals,
  saveOutput = TRUE,
  maxClust2 = 4,
  outputSuffix = "_1",
  sdCNV = 0.6,
  filterResults = TRUE,
  filterRange = 0.8,
  minAlteredCellProp = 0.75
)
```

## Arguments

- cnvResults:

  Output of
  [`clusterCNV()`](https://edogiuili.github.io/CopyscAT/reference/clusterCNV.md).

- expectedNormals:

  Barcodes of known or inferred normal cells, for example
  `nmf_results$normalBarcodes`.

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

- minAlteredCellProp:

  Proportion of the normal population expected to be unaltered on a
  retained arm.

## Value

A list of three elements, as for
[`annotateCNV4()`](https://edogiuili.github.io/CopyscAT/reference/annotateCNV4.md).

## See also

Other CNV annotation:
[`annotateCNV3()`](https://edogiuili.github.io/CopyscAT/reference/annotateCNV3.md),
[`annotateCNV4()`](https://edogiuili.github.io/CopyscAT/reference/annotateCNV4.md),
[`smoothClusters()`](https://edogiuili.github.io/CopyscAT/reference/smoothClusters.md)
