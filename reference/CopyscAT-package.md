# CopyscAT: Copy Number Variant Inference from Single-Cell ATAC Sequencing

Identifies large-scale and focal copy number alterations from
single-cell ATAC sequencing data without requiring a matched normal
control. Signal is binned across the genome, normalised for CpG density
and sequencing depth, then decomposed with Gaussian mixture models to
assign per-cell copy number states. Also detects double minutes and
other focal amplifications via changepoint analysis, estimates regions
of loss of heterozygosity, and separates neoplastic from non-neoplastic
cells by non-negative matrix factorisation.

## See also

Useful links:

- <https://github.com/edogiuili/CopyscAT>

- Report bugs at <https://github.com/edogiuili/CopyscAT/issues>

## Author

**Maintainer**: Ana Nikolic <nikolica@ucalgary.ca>

Other contributors:

- Edoardo Giuili <edoardo.giuili@ugent.be> \[contributor\]
