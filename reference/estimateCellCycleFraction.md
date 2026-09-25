# Estimate the fraction of cycling cells

Fits a Gaussian mixture to the total signal on one chromosome to
separate 2N from 4N cells. Chromosome X is the default because it is
often unaltered; pick another with `chromToUse` when X is affected in
your samples.

## Usage

``` r
estimateCellCycleFraction(
  inputMatrix,
  sampName,
  cutoff = 10000,
  maxG = 4,
  logTrans = FALSE,
  modelName = "E",
  chromToUse = "chrX"
)
```

## Arguments

- inputMatrix:

  Raw matrix, before normalisation.

- sampName:

  Sample name used in plot titles.

- cutoff:

  Minimum signal for a cell to be considered.

- maxG:

  Maximum mixture components.

- logTrans:

  Log-transform the signal first.

- modelName:

  mclust model name. Default `"E"` (equal variance).

- chromToUse:

  Chromosome used for the estimate.

## Value

A named integer vector of mixture assignments per barcode.

## See also

Other neoplastic classification:
[`identifyNonNeoplastic()`](https://edogiuili.github.io/CopyscAT/reference/identifyNonNeoplastic.md)
