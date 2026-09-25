# Write per-session summary statistics

Write per-session summary statistics

## Usage

``` r
saveSummaryStats(outputSuffix = "_stats")
```

## Arguments

- outputSuffix:

  Suffix appended to the output file name.

## Value

A named list of the statistics, invisibly written to
`<locPrefix><outPrefix><outputSuffix>.tsv`.

## See also

Other reporting:
[`graphCNVDistribution()`](https://edogiuili.github.io/CopyscAT/reference/graphCNVDistribution.md)
