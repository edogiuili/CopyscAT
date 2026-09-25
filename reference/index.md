# Package index

## Session setup

Build reference files and load them into the session. Run these before
any analysis step.

- [`generateReferences()`](https://edogiuili.github.io/CopyscAT/reference/generateReferences.md)
  : Build CopyscAT reference files for a genome
- [`initialiseEnvironment()`](https://edogiuili.github.io/CopyscAT/reference/initialiseEnvironment.md)
  : Load genome references and analysis parameters for a session
- [`initializeStandards()`](https://edogiuili.github.io/CopyscAT/reference/initializeStandards.md)
  : Reload only the reference tables
- [`setOutputFile()`](https://edogiuili.github.io/CopyscAT/reference/setOutputFile.md)
  : Set the output directory and file prefix

## Data import

Read a preprocessed fragment matrix.

- [`readInputTable()`](https://edogiuili.github.io/CopyscAT/reference/readInputTable.md)
  : Read a preprocessed single-cell fragment matrix

## Normalisation

Normalise raw counts, correct for CpG density and collapse binned signal
to chromosome arms.

- [`collapseChrom3N()`](https://edogiuili.github.io/CopyscAT/reference/collapseChrom3N.md)
  : Collapse binned signal to chromosome arms
- [`computeCenters()`](https://edogiuili.github.io/CopyscAT/reference/computeCenters.md)
  : Compute per-arm centre and spread
- [`normalizeMatrixN()`](https://edogiuili.github.io/CopyscAT/reference/normalizeMatrixN.md)
  : Normalise a raw single-cell fragment matrix
- [`scaleMatrix()`](https://edogiuili.github.io/CopyscAT/reference/scaleMatrix.md)
  : Scale arm signal to Z-scores

## Quality control

Inspect per-cell sequencing depth across samples, and remove
low-coverage cells and likely doublets.

- [`filterCells()`](https://edogiuili.github.io/CopyscAT/reference/filterCells.md)
  : Filter low-quality cells
- [`fragmentsPerCell()`](https://edogiuili.github.io/CopyscAT/reference/fragmentsPerCell.md)
  : Fragment counts per cell from fragment files
- [`plotCellDistribution()`](https://edogiuili.github.io/CopyscAT/reference/plotCellDistribution.md)
  : Plot the per-cell depth distribution
- [`plotCellKnee()`](https://edogiuili.github.io/CopyscAT/reference/plotCellKnee.md)
  : Plot ranked per-cell depth as a knee plot
- [`signalPerCell()`](https://edogiuili.github.io/CopyscAT/reference/signalPerCell.md)
  : Per-cell signal from one or more CopyscAT matrices
- [`summariseCellDepth()`](https://edogiuili.github.io/CopyscAT/reference/summariseCellDepth.md)
  : Summarise per-cell depth by sample

## Summarisation

Summary functions used to reduce per-bin signal to a single value, and
helpers for merging bins.

- [`cutAverage()`](https://edogiuili.github.io/CopyscAT/reference/cutAverage.md)
  : Trimmed mean between the configured quantiles
- [`scaleMatrixBins()`](https://edogiuili.github.io/CopyscAT/reference/scaleMatrixBins.md)
  : Merge adjacent genomic bins

## CNV calling

Decompose per-arm signal into copy number states and clean up the
resulting clusters.

- [`clusterCNV()`](https://edogiuili.github.io/CopyscAT/reference/clusterCNV.md)
  : Merge and renumber adjacent copy number clusters
- [`getAlteredSegments()`](https://edogiuili.github.io/CopyscAT/reference/getAlteredSegments.md)
  : Call sub-chromosomal gains and losses against a normal population
- [`identifyCNVClusters()`](https://edogiuili.github.io/CopyscAT/reference/identifyCNVClusters.md)
  : Call copy number states by Gaussian mixture decomposition
- [`makeFakeCells()`](https://edogiuili.github.io/CopyscAT/reference/makeFakeCells.md)
  : Simulate pseudodiploid control cells

## CNV annotation

Convert cluster assignments into copy number estimates.

- [`annotateCNV3()`](https://edogiuili.github.io/CopyscAT/reference/annotateCNV3.md)
  : Annotate copy number calls as loss, neutral or gain
- [`annotateCNV4()`](https://edogiuili.github.io/CopyscAT/reference/annotateCNV4.md)
  : Estimate absolute copy number per cell
- [`annotateCNV4B()`](https://edogiuili.github.io/CopyscAT/reference/annotateCNV4B.md)
  : Estimate absolute copy number against known normal cells
- [`smoothClusters()`](https://edogiuili.github.io/CopyscAT/reference/smoothClusters.md)
  : Smooth copy number calls within cell clusters

## Focal amplifications

Detect double minutes and other focal amplifications.

- [`getDoubleMinutes()`](https://edogiuili.github.io/CopyscAT/reference/getDoubleMinutes.md)
  : Detect focal amplifications in a single cell
- [`identifyDoubleMinutes()`](https://edogiuili.github.io/CopyscAT/reference/identifyDoubleMinutes.md)
  : Detect recurrent focal amplifications across all cells
- [`testOverlap()`](https://edogiuili.github.io/CopyscAT/reference/testOverlap.md)
  : Test whether two bin intervals overlap

## Loss of heterozygosity

Call LOH regions and annotate them with overlapping genes.

- [`annotateLosses()`](https://edogiuili.github.io/CopyscAT/reference/annotateLosses.md)
  : Annotate LOH regions with overlapping genes
- [`getLOHRegions()`](https://edogiuili.github.io/CopyscAT/reference/getLOHRegions.md)
  : Call regions of loss of heterozygosity

## Neoplastic classification

Separate tumour from normal cells and estimate cycling fraction.

- [`estimateCellCycleFraction()`](https://edogiuili.github.io/CopyscAT/reference/estimateCellCycleFraction.md)
  : Estimate the fraction of cycling cells
- [`identifyNonNeoplastic()`](https://edogiuili.github.io/CopyscAT/reference/identifyNonNeoplastic.md)
  : Separate neoplastic from non-neoplastic cells

## Reporting

Plots and summary statistics.

- [`graphCNVDistribution()`](https://edogiuili.github.io/CopyscAT/reference/graphCNVDistribution.md)
  : Plot the distribution of collapsed signal
- [`saveSummaryStats()`](https://edogiuili.github.io/CopyscAT/reference/saveSummaryStats.md)
  : Write per-session summary statistics

## Example data

- [`scDataSamp`](https://edogiuili.github.io/CopyscAT/reference/scDataSamp.md)
  : Example processed fragment matrix
- [`scDataSampClusters`](https://edogiuili.github.io/CopyscAT/reference/scDataSampClusters.md)
  : Example cluster assignments

## Session state

- [`scCNVCaller`](https://edogiuili.github.io/CopyscAT/reference/scCNVCaller.md)
  : Package-internal state for a CopyscAT session
- [`CopyscAT`](https://edogiuili.github.io/CopyscAT/reference/CopyscAT-package.md)
  [`CopyscAT-package`](https://edogiuili.github.io/CopyscAT/reference/CopyscAT-package.md)
  : CopyscAT: Copy Number Variant Inference from Single-Cell ATAC
  Sequencing
