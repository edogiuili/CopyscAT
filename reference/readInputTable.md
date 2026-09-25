# Read a preprocessed single-cell fragment matrix

Reads the count matrix written by `process_fragment_file.py`: cells in
rows, genomic bins in columns, with the first column holding the
barcode.

## Usage

``` r
readInputTable(inputFile, sep = "\t")
```

## Arguments

- inputFile:

  Path to the count matrix.

- sep:

  Field separator. Default tab.

## Value

A data frame of counts with barcodes as row names.

## Examples

``` r
if (FALSE) { # \dontrun{
scData <- readInputTable("my_sample_matrix.tsv")
} # }
```
