# Contributing to CopyscAT

## Development setup

``` r

install.packages(c("devtools", "roxygen2", "testthat", "pkgdown", "lintr"))
BiocManager::install(c("edgeR", "biomaRt", "rtracklayer", "GenomicRanges"))
```

``` r

devtools::load_all(".")   # load without installing
devtools::test()          # run the test suite
devtools::document()      # regenerate NAMESPACE and man/ after editing roxygen
devtools::check()         # full R CMD check
lintr::lint_package()     # style and correctness lints
```

## Layout

Source files are grouped by analysis stage, and each holds the functions
for one step of the workflow:

| File | Contents |
|----|----|
| `R/zzz.R` | Session environment and the initialisation guard |
| `R/setup.R` | [`initialiseEnvironment()`](https://edogiuili.github.io/CopyscAT/reference/initialiseEnvironment.md), [`setOutputFile()`](https://edogiuili.github.io/CopyscAT/reference/setOutputFile.md), input reading |
| `R/normalize.R` | [`normalizeMatrixN()`](https://edogiuili.github.io/CopyscAT/reference/normalizeMatrixN.md) |
| `R/collapse.R` | [`collapseChrom3N()`](https://edogiuili.github.io/CopyscAT/reference/collapseChrom3N.md) |
| `R/filter.R` | [`filterCells()`](https://edogiuili.github.io/CopyscAT/reference/filterCells.md), [`computeCenters()`](https://edogiuili.github.io/CopyscAT/reference/computeCenters.md), [`scaleMatrix()`](https://edogiuili.github.io/CopyscAT/reference/scaleMatrix.md) |
| `R/summarise.R` | [`cutAverage()`](https://edogiuili.github.io/CopyscAT/reference/cutAverage.md), [`scaleMatrixBins()`](https://edogiuili.github.io/CopyscAT/reference/scaleMatrixBins.md) |
| `R/cluster.R` | [`identifyCNVClusters()`](https://edogiuili.github.io/CopyscAT/reference/identifyCNVClusters.md), [`clusterCNV()`](https://edogiuili.github.io/CopyscAT/reference/clusterCNV.md) |
| `R/annotate.R` | [`annotateCNV3()`](https://edogiuili.github.io/CopyscAT/reference/annotateCNV3.md), [`annotateCNV4()`](https://edogiuili.github.io/CopyscAT/reference/annotateCNV4.md), [`annotateCNV4B()`](https://edogiuili.github.io/CopyscAT/reference/annotateCNV4B.md) |
| `R/double-minutes.R` | Focal amplification detection |
| `R/loh.R` | LOH calling and gene annotation |
| `R/neoplastic.R` | NMF-based normal/tumour separation, cell cycle |
| `R/segments.R` | Sub-chromosomal gains and losses |
| `R/references.R` | [`generateReferences()`](https://edogiuili.github.io/CopyscAT/reference/generateReferences.md) |
| `R/plots.R`, `R/smooth.R` | Plotting and cluster smoothing |

Put a new function in the file matching its stage rather than creating a
new one, unless it starts a genuinely new stage.

## Conventions

- Declare every external function in `R/CopyscAT-package.R` with
  `@importFrom`, or qualify it as `pkg::fun()`. Do not rely on a package
  being attached.
- Put a package in `Suggests` when only one entry point needs it, and
  guard the call with
  [`requireNamespace()`](https://rdrr.io/r/base/ns-load.html).
- Validate arguments at the top of exported functions, and call
  `assert_initialised()` when the function reads the session
  environment.
- Report failures with [`stop()`](https://rdrr.io/r/base/stop.html) and
  recoverable problems with
  [`warning()`](https://rdrr.io/r/base/warning.html) or
  [`message()`](https://rdrr.io/r/base/message.html). Do not print a
  message and return `NULL`; that surfaces as a confusing error further
  down the pipeline.
- Name the parameter to adjust in error messages, for example “Lower
  `blacklistCutoff` (currently 125)”.
- Wrap [`pdf()`](https://rdrr.io/r/grDevices/pdf.html) with
  `on.exit(dev.off(), add = TRUE)` so errors do not leave a truncated
  file.
- Use [`seq_len()`](https://rdrr.io/r/base/seq.html) and
  [`seq_along()`](https://rdrr.io/r/base/seq.html) rather than `1:n`,
  which iterates once when `n` is zero.
- Add columns referenced by non-standard evaluation to the
  [`globalVariables()`](https://rdrr.io/r/utils/globalVariables.html)
  call in `R/CopyscAT-package.R`.

## Tests

Tests live in `tests/testthat/`. `tests/testthat/helper-copyscat.R`
provides `make_test_references()`, `local_copyscat_session()` and
`simulate_counts()`, which build a small synthetic genome and a count
matrix with known copy number alterations. Use them rather than the
shipped example data, which is too large and too slow for a test.

- Unit tests: `test-summarise.R`, `test-overlap.R`, `test-setup.R`
- Error paths: `test-error-handling.R`
- Integration: `test-integration-pipeline.R`,
  `test-integration-cnv-calling.R`

A bug fix should come with a test that fails before it.

## Style

The upstream code predates tidyverse style, and `.lintr` disables the
purely cosmetic linters so that real findings stay visible. Match the
surrounding style when editing existing code; write new code in the
house style above.
