# Call sub-chromosomal gains and losses against a normal population

Fits a negative binomial model per bin comparing neoplastic with
non-neoplastic cells, segments the resulting fold-change track by
changepoint analysis, and assigns per-cell copy number to each segment
by linear discriminant analysis. Requires the normal population
identified by
[`identifyNonNeoplastic()`](https://edogiuili.github.io/CopyscAT/reference/identifyNonNeoplastic.md).

## Usage

``` r
getAlteredSegments(
  inputMatrix,
  clusterResults,
  minSeglen = 3,
  lossThreshold = 0.5,
  gainThreshold = 1.6,
  minCutoff = 50,
  rangeThreshold = 0.4
)
```

## Arguments

- inputMatrix:

  Normalised matrix from
  [`normalizeMatrixN()`](https://edogiuili.github.io/CopyscAT/reference/normalizeMatrixN.md).

- clusterResults:

  Output of
  [`identifyNonNeoplastic()`](https://edogiuili.github.io/CopyscAT/reference/identifyNonNeoplastic.md).

- minSeglen:

  Minimum bins per changepoint segment.

- lossThreshold:

  Fold change below which a segment is a loss.

- gainThreshold:

  Fold change above which a segment is a gain.

- minCutoff:

  Minimum cells in the smallest group of a retained segment.

- rangeThreshold:

  Minimum fold-change range for a retained segment.

## Value

A data frame of barcodes with estimated copy number per called segment,
named `chrom:start-end`.

## See also

Other CNV calling:
[`clusterCNV()`](https://edogiuili.github.io/CopyscAT/reference/clusterCNV.md),
[`identifyCNVClusters()`](https://edogiuili.github.io/CopyscAT/reference/identifyCNVClusters.md),
[`makeFakeCells()`](https://edogiuili.github.io/CopyscAT/reference/makeFakeCells.md)

## Examples

``` r
if (FALSE) { # \dontrun{
altered_segments <- getAlteredSegments(
  scData_k_norm, nmf_results, lossThreshold = 0.6
)
} # }
```
