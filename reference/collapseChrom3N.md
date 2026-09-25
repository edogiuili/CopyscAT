# Collapse binned signal to chromosome arms

Removes centromeric and blacklisted bins, divides signal by a power
transform of CpG density to correct the accessibility bias, then
summarises each chromosome arm per cell. Also infers sample sex from the
X/Y signal difference and records it in `scCNVCaller$isMale`.

## Usage

``` r
collapseChrom3N(
  inputMatrix,
  minimumSegments = 40,
  summaryFunction = cutAverage,
  logTrans = FALSE,
  binExpand = 1,
  minimumChromValue = 2,
  tssEnrich = 5,
  logBase = 2,
  minCPG = 300,
  powVal = 0.73
)
```

## Arguments

- inputMatrix:

  Normalised matrix from
  [`normalizeMatrixN()`](https://edogiuili.github.io/CopyscAT/reference/normalizeMatrixN.md).

- minimumSegments:

  Unused; kept for backward compatibility.

- summaryFunction:

  Function used to summarise each arm. Default
  [`cutAverage()`](https://edogiuili.github.io/CopyscAT/reference/cutAverage.md).

- logTrans:

  Apply a log transform to the CpG scaling factor.

- binExpand:

  Merge this many adjacent bins before collapsing. `1` leaves bins
  untouched.

- minimumChromValue:

  Drop arms whose median signal falls below this.

- tssEnrich:

  Weight applied to the CpG scaling factor.

- logBase:

  Log base used when `logTrans = TRUE`.

- minCPG:

  Drop bins with CpG density below this.

- powVal:

  Exponent of the CpG correction. `0.70`-`0.75` suits hg19 and hg38.

## Value

A data frame with a `chrom` column of arm labels, one column per cell,
and a `medianNorm` column.

## See also

Other normalisation:
[`computeCenters()`](https://edogiuili.github.io/CopyscAT/reference/computeCenters.md),
[`normalizeMatrixN()`](https://edogiuili.github.io/CopyscAT/reference/normalizeMatrixN.md),
[`scaleMatrix()`](https://edogiuili.github.io/CopyscAT/reference/scaleMatrix.md)

## Examples

``` r
if (FALSE) { # \dontrun{
scData_collapse <- collapseChrom3N(
  scData_k_norm, summaryFunction = cutAverage,
  minimumChromValue = 100, minCPG = 300, powVal = 0.73
)
} # }
```
