# ============================================================================
# 04_empirical_squamates.R
#
# Empirical application of AncCond: Does body size predict transitions
# between oviparity and viviparity in squamate reptiles?
#
# Data:
#    - Tonini et al. 2016 posterior trees (100 trees, pruned via VertLife)
#    - SquamBase (Meiri 2024) trait data: max SVL + reproductive mode
#
# Outputs:
#    - Console: AncCond results summary
#    - ../results/empirical_squamate_results.csv
#    - ../figures/figure6_empirical.pdf
#
# Dependencies: ape, phytools, parallel
# ============================================================================

library(ape)
library(phytools)
library(parallel)

# ── Detect cores ────────────────────────────────────────────────────────────
n.cores <- 13
cat("Using", n.cores, "cores\n")

# Initialize the cluster
cl <- makeCluster(n.cores)

# ── Paths (relative to manuscript/scripts/) ──
source("AncCond.R")
tree.file  <- "../data/output.nex"
trait.file <- "../data/squamate_traits.csv"
out.csv    <- "../results/empirical_squamate_results.csv"
out.fig    <- "../figures/figure6_empirical.pdf"

# ── 1. Read trees ──
cat("Reading posterior trees...\n")
trees <- read.nexus(tree.file)
n.trees <- length(trees)
cat("  ", n.trees, "posterior trees,", length(trees[[1]]$tip.label), "tips each\n")

# ── 2. Read trait data ──
dat <- read.csv(trait.file, stringsAsFactors = FALSE)
cat("Trait data (raw):", nrow(dat), "rows\n")

# Deduplicate: keep mean log_svl per species (some species have multiple entries)
dups <- sum(duplicated(dat$species))
if (dups > 0) cat("  Deduplicating", dups, "duplicate species (using mean SVL)...\n")
dat <- aggregate(log_svl_mm ~ species + repro_mode, data = dat, FUN = mean)
cat("Trait data (unique):", nrow(dat), "species\n")
cat("  Oviparous:", sum(dat$repro_mode == "Oviparous"), "\n")
cat("  Viviparous:", sum(dat$repro_mode == "Viviparous"), "\n")

# ── 3. Match to first tree (all trees have same tip set) ──
tree1 <- trees[[1]]
shared <- intersect(tree1$tip.label, dat$species)
cat("Species in tree AND traits:", length(shared), "\n")

rownames(dat) <- dat$species
dat.m <- dat[shared, ]

# Continuous trait: log(max SVL in mm)
x <- setNames(dat.m$log_svl_mm, dat.m$species)

# Discrete trait: 0 = oviparous (ancestral), 1 = viviparous (derived)
y <- setNames(as.integer(dat.m$repro_mode == "Viviparous"), dat.m$species)

cat("  State 0 (oviparous):", sum(y == 0), "\n")
cat("  State 1 (viviparous):", sum(y == 1), "\n")

# ── Quick sanity check: run tree 1 serially first ──
n.maps <- 50
n.sims <- 500

cat("\n── SANITY CHECK: Tree 1 (serial) ──\n")
t0 <- proc.time()
tr1 <- drop.tip(trees[[1]], setdiff(trees[[1]]$tip.label, shared))
cat("  Tree pruned:", length(tr1$tip.label), "tips\n")
cat("  Running AncCond (n.maps=", n.maps, ", n.sims=", n.sims, ")...\n")
res1 <- AncCond(tr1, x, y, n.maps = n.maps, n.sims = n.sims, hypothesis = "none")
t1 <- proc.time() - t0
cat("  Elapsed:", round(t1[3], 1), "sec\n")
cat("  Mk model:", res1$model, "\n")
cat("  q01 =", round(res1$Q[1,2], 4), "   q10 =", round(res1$Q[2,1], 4), "\n")
cat("  T.obs.01 =", round(res1$T.obs.01, 4),
    " (SVL =", round(exp(res1$T.obs.01), 1), "mm)\n")
cat("  T.obs.10 =", round(res1$T.obs.10, 4),
    " (SVL =", round(exp(res1$T.obs.10), 1), "mm)\n")
cat("  p.01 =", round(res1$p.01, 4), "\n")
cat("  p.10 =", round(res1$p.10, 4), "\n")
cat("  p.combined =", round(res1$p.combined, 4), "\n")
cat("  Estimated total time:", round(t1[3] * n.trees / n.cores / 60, 1), "min on",
    n.cores, "cores\n")

# Store tree 1 nulls for figure
null.01.tree1 <- res1$T.null.01
null.10.tree1 <- res1$T.null.10
obs.01.tree1  <- res1$T.obs.01
obs.10.tree1  <- res1$T.obs.10

cat("\n── Tree 1 looks OK? Launching parallel run on remaining", n.trees - 1, "trees ──\n\n")

# ── 4. Parallel Setup ──

# Worker function — Define BEFORE clusterExport
run.one.tree <- function(i, trees, shared, x, y, n.maps, n.sims) {
  library(ape)
  library(phytools)
  tr <- drop.tip(trees[[i]], setdiff(trees[[i]]$tip.label, shared))
  res <- AncCond(tr, x, y, n.maps = n.maps, n.sims = n.sims, hypothesis = "none")
  list(
    tree       = i,
    T.obs.01   = res$T.obs.01,
    T.obs.10   = res$T.obs.10,
    p.01       = res$p.01,
    p.10       = res$p.10,
    p.upper.01 = res$p.upper.01,
    p.lower.01 = res$p.lower.01,
    p.upper.10 = res$p.upper.10,
    p.lower.10 = res$p.lower.10,
    p.combined = res$p.combined,
    q01        = res$Q[1, 2],
    q10        = res$Q[2, 1],
    model      = res$model
  )
}

# Export environment and source the required file on workers
clusterEvalQ(cl, {
  library(ape)
  library(phytools)
  source("AncCond.R")
})

# Export data AND the function to the cluster
clusterExport(cl, c("trees", "shared", "x", "y", "n.maps", "n.sims", "run.one.tree"))

# ── 5. Parallel run across remaining posterior trees ──
t0.all <- proc.time()

# Replaced mclapply with parLapply
par.results <- parLapply(cl, 2:n.trees, function(i) {
  run.one.tree(i, trees, shared, x, y, n.maps, n.sims)
})

t1.all <- proc.time() - t0.all
cat("Parallel run complete:", round(t1.all[3] / 60, 1), "min\n\n")

# Close cluster
stopCluster(cl)

# ── 6. Combine results ──
results <- data.frame(
  tree        = integer(n.trees),
  T.obs.01    = numeric(n.trees),
  T.obs.10    = numeric(n.trees),
  p.01        = numeric(n.trees),
  p.10        = numeric(n.trees),
  p.upper.01  = numeric(n.trees),
  p.lower.01  = numeric(n.trees),
  p.upper.10  = numeric(n.trees),
  p.lower.10  = numeric(n.trees),
  p.combined  = numeric(n.trees),
  q01         = numeric(n.trees),
  q10         = numeric(n.trees),
  model       = character(n.trees),
  stringsAsFactors = FALSE
)

# Fill tree 1
results$tree[1]        <- 1
results$T.obs.01[1]    <- res1$T.obs.01
results$T.obs.10[1]    <- res1$T.obs.10
results$p.01[1]        <- res1$p.01
results$p.10[1]        <- res1$p.10
results$p.upper.01[1]  <- res1$p.upper.01
results$p.lower.01[1]  <- res1$p.lower.01
results$p.upper.10[1]  <- res1$p.upper.10
results$p.lower.10[1]  <- res1$p.lower.10
results$p.combined[1]  <- res1$p.combined
results$q01[1]         <- res1$Q[1, 2]
results$q10[1]         <- res1$Q[2, 1]
results$model[1]       <- res1$model

# Check for errors in parallel results
n.errors <- 0
for (j in seq_along(par.results)) {
  r <- par.results[[j]]
  i <- j + 1  # tree index
  if (inherits(r, "try-error") || is.null(r)) {
    cat("  WARNING: Tree", i, "failed.\n")
    n.errors <- n.errors + 1
    next
  }
  results$tree[i]        <- r$tree
  results$T.obs.01[i]    <- r$T.obs.01
  results$T.obs.10[i]    <- r$T.obs.10
  results$p.01[i]        <- r$p.01
  results$p.10[i]        <- r$p.10
  results$p.upper.01[i]  <- r$p.upper.01
  results$p.lower.01[i]  <- r$p.lower.01
  results$p.upper.10[i]  <- r$p.upper.10
  results$p.lower.10[i]  <- r$p.lower.10
  results$p.combined[i]  <- r$p.combined
  results$q01[i]         <- r$q01
  results$q10[i]         <- r$q10
  results$model[i]       <- r$model
}

if (n.errors > 0) cat("\n", n.errors, "trees failed!\n")

# ── 7. Summary and Save ──
cat("\n══════════════════════════════════════════════\n")
cat("SUMMARY ACROSS", n.trees, "POSTERIOR TREES\n")
cat("══════════════════════════════════════════════\n\n")

cat("Trees with p.combined < 0.05:", sum(results$p.combined < 0.05, na.rm=TRUE), "/", n.trees, "\n\n")

write.csv(results, out.csv, row.names = FALSE)
cat("Results saved to:", out.csv, "\n")

# ── 8. Figure 6: null distribution + observed value ──
pdf(out.fig, width = 7, height = 4)
par(mfrow = c(1, 2), mar = c(4.5, 4.5, 2.5, 1), cex.lab = 1.2, cex.axis = 1)

# Panel A: 0→1 transitions
valid.01 <- null.01.tree1[is.finite(null.01.tree1)]
if (length(valid.01) > 10) {
  hist(valid.01, breaks = 30, col = "grey80", border = "grey50",
       main = expression("Gains of viviparity (0" %->% "1)"),
       xlab = "Ancestral log(SVL mm)", ylab = "Frequency",
       xlim = range(c(valid.01, obs.01.tree1), na.rm = TRUE))
  abline(v = obs.01.tree1, col = "firebrick", lwd = 2.5, lty = 1)
  mtext("A", side = 3, adj = -0.1, line = 0.8, font = 2, cex = 1.3)
}

# Panel B: 1→0 transitions
valid.10 <- null.10.tree1[is.finite(null.10.tree1)]
if (length(valid.10) > 10) {
  hist(valid.10, breaks = 30, col = "grey80", border = "grey50",
       main = expression("Losses of viviparity (1" %->% "0)"),
       xlab = "Ancestral log(SVL mm)", ylab = "Frequency",
       xlim = range(c(valid.10, obs.10.tree1), na.rm = TRUE))
  abline(v = obs.10.tree1, col = "steelblue", lwd = 2.5, lty = 1)
  mtext("B", side = 3, adj = -0.1, line = 0.8, font = 2, cex = 1.3)
}

dev.off()
cat("Figure saved to:", out.fig, "\n")
cat("\nDone.\n")