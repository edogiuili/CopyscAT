# Trimmed mean between the configured quantiles

Default summary function throughout CopyscAT. Averages the values lying
between the `lowerTrim` and `upperTrim` quantiles set by
[`initialiseEnvironment()`](https://edogiuili.github.io/CopyscAT/reference/initialiseEnvironment.md),
which resists the heavy tails typical of single-cell ATAC signal.

## Usage

``` r
cutAverage(inputVector)
```

## Arguments

- inputVector:

  Numeric vector.

## Value

A single numeric value, or `NA_real_` when the quantiles are undefined.

## See also

Other summarisation:
[`scaleMatrixBins()`](https://edogiuili.github.io/CopyscAT/reference/scaleMatrixBins.md)

## Examples

``` r
if (FALSE) { # \dontrun{
cutAverage(rnorm(100))
} # }
```
