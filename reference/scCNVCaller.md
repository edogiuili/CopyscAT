# Package-internal state for a CopyscAT session

CopyscAT keeps genome references, binning parameters and output paths in
a dedicated environment rather than passing them through every call.
Populate it with
[`initialiseEnvironment()`](https://edogiuili.github.io/CopyscAT/reference/initialiseEnvironment.md)
and
[`setOutputFile()`](https://edogiuili.github.io/CopyscAT/reference/setOutputFile.md)
before running any analysis step.

## Usage

``` r
scCNVCaller
```

## Format

An [environment](https://rdrr.io/r/base/environment.html).

## Details

The environment is exported because the original tutorial workflow reads
`scCNVCaller$locPrefix` and `scCNVCaller$outPrefix` directly when
writing user-side output files.
