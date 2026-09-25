# Annotate copy number calls as loss, neutral or gain

Labels each arm relative to the inferred normal cluster and writes the
per-cell calls to disk. Superseded by
[`annotateCNV4()`](https://edogiuili.github.io/CopyscAT/reference/annotateCNV4.md),
which reports estimated absolute copy number instead of three-level
calls.

## Usage

``` r
annotateCNV3(cnvResults, saveOutput = TRUE, maxClust2 = 4, outputSuffix = "_1")
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

## Value

A list of three elements: ternary calls per cell, the candidate CNV
table, and binarised calls.

## See also

[`annotateCNV4()`](https://edogiuili.github.io/CopyscAT/reference/annotateCNV4.md)
for the current workflow.

Other CNV annotation:
[`annotateCNV4()`](https://edogiuili.github.io/CopyscAT/reference/annotateCNV4.md),
[`annotateCNV4B()`](https://edogiuili.github.io/CopyscAT/reference/annotateCNV4B.md),
[`smoothClusters()`](https://edogiuili.github.io/CopyscAT/reference/smoothClusters.md)
