# Test whether two bin intervals overlap

Intervals are encoded as `chrom_start.chrom_end`, the format produced by
[`getDoubleMinutes()`](https://edogiuili.github.io/CopyscAT/reference/getDoubleMinutes.md).

## Usage

``` r
testOverlap(interval.1, interval.2)
```

## Arguments

- interval.1, interval.2:

  Interval labels.

## Value

A character vector of length 3: merged left end, merged right end, and
`"TRUE"`/`"FALSE"` for overlap. Character, not numeric, because the ends
are bin labels.

## See also

Other focal amplification:
[`getDoubleMinutes()`](https://edogiuili.github.io/CopyscAT/reference/getDoubleMinutes.md),
[`identifyDoubleMinutes()`](https://edogiuili.github.io/CopyscAT/reference/identifyDoubleMinutes.md)
