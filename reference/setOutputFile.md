# Set the output directory and file prefix

All functions that write plots or tables prepend `outputDir` and
`outputFileSuffix` to their file names.

## Usage

``` r
setOutputFile(outputDir, outputFileSuffix)
```

## Arguments

- outputDir:

  Directory for output files. A trailing slash is added when missing.
  The directory must already exist.

- outputFileSuffix:

  Prefix used for every file written in this session, typically the
  sample name.

## Value

Invisibly `TRUE`.

## See also

Other session setup:
[`generateReferences()`](https://edogiuili.github.io/CopyscAT/reference/generateReferences.md),
[`initialiseEnvironment()`](https://edogiuili.github.io/CopyscAT/reference/initialiseEnvironment.md),
[`initializeStandards()`](https://edogiuili.github.io/CopyscAT/reference/initializeStandards.md)

## Examples

``` r
setOutputFile(tempdir(), "demo_sample")
```
