av <- commandArgs(TRUE)
a <- readRDS(av[1]); b <- readRDS(av[2])
cmp <- function(name, x, y) {
  r <- all.equal(x, y, tolerance = 1e-9)
  if (isTRUE(r)) { cat(sprintf("  %-26s IDENTICAL\n", name)); TRUE }
  else { cat(sprintf("  %-26s DIFFERS\n", name))
         cat(paste0("       ", utils::head(r, 6), collapse="\n"), "\n"); FALSE }
}
cat("=== stage-by-stage (ORIG v0.40 vs NEW v1.0.0) ===\n")
ok <- c(
 cmp("normalizeMatrixN",      a$norm,       b$norm),
 cmp("collapseChrom3N",       a$collapse,   b$collapse),
 cmp("sex inference (isMale)",a$isMale,     b$isMale),
 cmp("filterCells",           a$filter,     b$filter),
 cmp("computeCenters",        a$centers,    b$centers),
 cmp("scaleMatrix",           a$scaled,     b$scaled),
 cmp("identifyCNV [assign]",  a$cand[[1]],  b$cand[[1]]),
 cmp("identifyCNV [means]",   a$cand[[2]],  b$cand[[2]]),
 cmp("clusterCNV [assign]",   a$clean[[1]], b$clean[[1]]),
 cmp("clusterCNV [means]",    a$clean[[2]], b$clean[[2]]),
 cmp("annotateCNV4 [scores]", a$final[[3]], b$final[[3]])
)
cat("\n=== summary statistics ===\n")
for (k in names(a$summary)) {
  same <- isTRUE(all.equal(a$summary[[k]], b$summary[[k]]))
  cat(sprintf("  %-12s orig=%-10s new=%-10s %s\n", k, format(a$summary[[k]]),
      format(b$summary[[k]]), if (same) "same" else "<-- DIFFERS"))
}
cat("\n=== verdict ===\n")
cat(sum(ok), "of", length(ok), "stages numerically identical\n")
