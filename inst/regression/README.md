# Cross-version regression check

Verifies that a change to CopyscAT leaves CNV calls untouched, by running the
same pipeline under two installed versions and comparing every intermediate
numerically.

`run-pipeline.R` builds a small synthetic genome, simulates counts with a fixed
seed and a known gain on chr4, runs the full workflow, and saves every stage to
an RDS. `compare-runs.R` diffs two such files with `all.equal()`.

## Usage

Install the reference version into its own library:

```sh
mkdir -p /tmp/origlib
git worktree add /tmp/copyscat-ref <ref>       # tag, branch or commit
R CMD INSTALL -l /tmp/origlib --no-docs /tmp/copyscat-ref
```

Install the working version the same way, then run both and compare:

```sh
Rscript inst/regression/run-pipeline.R /tmp/origlib REF /tmp/res_REF.rds
Rscript inst/regression/run-pipeline.R /tmp/newlib  NEW /tmp/res_NEW.rds
Rscript inst/regression/compare-runs.R /tmp/res_REF.rds /tmp/res_NEW.rds
```

## Interpreting the result

Every stage should report `IDENTICAL`. A `DIFFERS` line is not automatically a
regression, but it must be explained before release — see the v1.0.0 entry in
`NEWS.md`, where two stages differ because the older version was wrong.

Both runs set the same seed before `identifyCNVClusters()`, which samples cells
for mixture initialisation. Comparing across machines or R versions can still
shift results, because the mixture fits are sensitive to BLAS differences; run
both halves on the same machine.
