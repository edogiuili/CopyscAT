# Simulate pseudodiploid control cells

Draws normally distributed values around `x` to act as a synthetic
diploid baseline when a sample contains too few non-neoplastic cells.

## Usage

``` r
makeFakeCells(x, num = 300, sd_fact = 0.1)
```

## Arguments

- x:

  Mean value to simulate around.

- num:

  Number of cells to simulate.

- sd_fact:

  Standard deviation as a fraction of `x`.

## Value

A numeric vector of length `num`.

## See also

Other CNV calling:
[`clusterCNV()`](https://edogiuili.github.io/CopyscAT/reference/clusterCNV.md),
[`getAlteredSegments()`](https://edogiuili.github.io/CopyscAT/reference/getAlteredSegments.md),
[`identifyCNVClusters()`](https://edogiuili.github.io/CopyscAT/reference/identifyCNVClusters.md)
