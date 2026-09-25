# Reload only the reference tables

Replaces the genome references in
[scCNVCaller](https://edogiuili.github.io/CopyscAT/reference/scCNVCaller.md)
while leaving binning parameters untouched. Use it to switch reference
builds inside a session; otherwise prefer
[`initialiseEnvironment()`](https://edogiuili.github.io/CopyscAT/reference/initialiseEnvironment.md).

## Usage

``` r
initializeStandards(chromSizeFile, cpgDataFile, cytobandFile)
```

## Arguments

- chromSizeFile:

  Path to a chromosome-size table.

- cpgDataFile:

  Path to a binned CpG density table.

- cytobandFile:

  Path to a binned cytoband table.

## Value

Invisibly `TRUE`.

## See also

Other session setup:
[`generateReferences()`](https://edogiuili.github.io/CopyscAT/reference/generateReferences.md),
[`initialiseEnvironment()`](https://edogiuili.github.io/CopyscAT/reference/initialiseEnvironment.md),
[`setOutputFile()`](https://edogiuili.github.io/CopyscAT/reference/setOutputFile.md)
