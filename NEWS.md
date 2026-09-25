# CopyscAT 1.0.0

Restructuring release. Existing scripts should continue to work.

CNV calls are unchanged. This was verified by installing v0.40 and v1.0.0 side
by side and running both over the same simulated dataset with a fixed seed:
the outputs of `collapseChrom3N()`, `filterCells()`, `computeCenters()`,
`scaleMatrix()`, `identifyCNVClusters()`, `clusterCNV()` and `annotateCNV4()`
are numerically identical.

Two outputs do differ, both because v0.40 was wrong:

* The `pos` column of `normalizeMatrixN()` now holds full coordinates. v0.40
  split bin names on any non-alphanumeric character, so a position written in
  scientific notation was truncated at the exponent sign — `1e+06` became
  `"1e"`, and every bin on a chromosome collapsed to the same handful of
  labels. Downstream stages key on `chrom`, which is why the CNV calls were
  unaffected. Counts were never altered.
* `saveSummaryStats()` reports accurate cell counts. On the verification run
  v0.40 reported `finalFilterCells = 148` and `startingCellCount = 0` where
  the true values were 147 and 150. Both versions produced the same 147-cell
  matrix; only the reported numbers were wrong.

## New features

* Added five functions for inspecting per-cell sequencing depth before running
  the CNV workflow:
  * `signalPerCell()` summarises the signal each barcode carries, from matrix
    files or in-memory data frames, across any number of samples.
  * `fragmentsPerCell()` counts fragments per barcode from a 10x
    `fragments.tsv.gz` or a CellRanger `singlecell.csv`.
  * `summariseCellDepth()` returns per-sample quartiles.
  * `plotCellDistribution()` draws faceted log-scale histograms, with an
    optional line at a proposed cutoff.
  * `plotCellKnee()` draws ranked barcode curves, whose elbow marks the
    transition from cells to background.

  Note that `signalPerCell()` reports base pairs, not fragments:
  `process_fragment_file.py` accumulates fragment *lengths*, so the row sum of
  a CopyscAT matrix is the genomic span a cell covers. The two functions
  therefore return differently named columns, `signal_bp` and `n_fragments`,
  and the plotting functions use whichever is present.

## Package structure

* Split the single 2504-line `cnv_caller_functions_packaged.R` into topical
  files (`setup.R`, `normalize.R`, `collapse.R`, `cluster.R`, `annotate.R`,
  `loh.R`, and others).
* Moved 23 packages from `Depends` to `Imports`, so loading CopyscAT no longer
  attaches them to the user's search path. `NAMESPACE` now declares explicit
  imports instead of relying on attachment.
* Dropped `Rtsne`, `FNN`, `igraph`, `zoo`, `limma` and `sp`, none of which
  were used.
* Moved `biomaRt`, `rtracklayer`, `GenomicRanges` and `jsonlite` to `Suggests`.
  `annotateLosses()` and `generateReferences()` now report clearly when a
  suggested package is missing.
* Added `.gitignore` and `.Rbuildignore`; untracked `.DS_Store`, `.Rhistory`
  and a stale duplicate of the source file.

## Bug fixes

* `readInputTable()` ignored its `sep` argument and always read tab-separated
  files.
* `cutAverage()` returned `NA` whenever the input contained any `NA`, despite
  computing its quantiles with `na.rm = TRUE`.
* `filterCells()` recorded the pre-outlier-filter cell count in
  `finalFilterCells`, so `saveSummaryStats()` overstated surviving cells. Its
  two progress messages also reported the same denominator.
* `collapseChrom3N()` aborted with `argument is of length zero` on genomes
  without both chrX and chrY; sex inference is now skipped with a message.
* `collapseChrom3N()` reassigned `logTrans <- FALSE` mid-function, making the
  `logTrans = TRUE` branch unreachable. Sex inference keeps its original
  untransformed behaviour via a separate local variable.
* `clusterCNV()` called `S4Vectors::isSorted()`, available only because a
  Bioconductor package happened to be attached. Replaced with a local helper.
* `identifyCNVClusters()` failed with an opaque `sample()` error when
  `subsetSize` exceeded the cell count; it is now clamped with a warning.
* `mclust::mclustBIC` was not imported, so `identifyCNVClusters()` failed with
  `could not find function "mclustBIC"` unless mclust was attached.
* `getAlteredSegments()` hardcoded the barcode suffix `-1` instead of using
  `scCNVCaller$cellSuffix`.
* `graphics::segments()` and `stats::as.dendrogram()` were used without being
  imported.
* `normalizeMatrixN()` split bin names with an unanchored `separate()`, which
  truncated coordinates containing extra separators.
* PDF devices are now closed with `on.exit()`, so an error mid-analysis no
  longer leaves a device open and the file truncated.
* `annotateCNV4B()` failed with "`condition` must be a logical vector, not a
  logical matrix" on current dplyr. `t()` on the one-row normal-cluster frame
  returns an n-by-1 matrix; older dplyr collapsed it inside `mutate()`, newer
  versions do not. The value is now converted with `as.vector()`, and a guard
  reports a mismatch between the number of cluster indices and chromosome arms.
* `normalizeMatrixN()` let a barcode with no signal reach `edgeR::cpm()`, which
  rejected it with "library sizes should be greater than zero" — naming neither
  the barcode nor the parameter at fault. This happens whenever `maxZero` is at
  or above the bin count, which disables the empty-bin filter. The new check
  lists the affected barcodes, the retained and blacklisted bin counts, whether
  the filter is disabled, and the three ways to resolve it.

## Error handling

* `initialiseEnvironment()` and `setOutputFile()` validate paths, bin sizes,
  trim quantiles and barcode suffixes, and warn when the CpG and cytoband
  references describe different numbers of bins.
* Functions that need an initialised session now say so instead of failing on
  a `NULL` lookup further downstream.
* Conditions that previously printed a message and returned `NULL` — every bin
  blacklisted, no cell passing a filter — now raise errors naming the
  parameter to adjust.
* `annotateCNV4B()` validates `expectedNormals` against the cluster
  assignments.

## Documentation and infrastructure

* Rewrote all roxygen documentation: every parameter documented, `@return`
  sections added, functions grouped into families for the reference index.
* Added a `copyscat` vignette covering the full workflow.
* Added a pkgdown site configuration and a rewritten README.
* Added a testthat suite of 122 tests, including end-to-end integration tests
  that verify a simulated copy number gain is recovered by the pipeline.
* Added GitHub Actions for R CMD check (Linux/macOS/Windows), test coverage,
  linting and pkgdown deployment.
* `identifyCNVClusters()` gained a `verbose` argument; per-arm mixture
  summaries are no longer printed unconditionally.
