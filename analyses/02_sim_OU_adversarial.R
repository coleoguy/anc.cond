# ============================================================================
# R3-SIM-1 + R1-SIM-1: OU Adversarial Scenario (Causation Reversal Test)
#
# Tests whether AncCond produces false positives when the DISCRETE trait
# drives the CONTINUOUS trait (opposite of what AncCond assumes).
#
# Setup:
#    1. Discrete trait evolves under Mk (standard forward/reverse rates).
#    2. Continuous trait evolves under two-optima OU where the optimum
#       switches depending on the current discrete state along the branch.
#    3. AncCond is run on the resulting tip data.
#
# Parallelised with parallel socket clusters (parLapply) for cross-platform support.
#
# Dependencies: ape, phytools, parallel, R.utils
# ============================================================================

library(ape)
library(phytools)
library(parallel)
library(R.utils)        # for withTimeout
source("AncCond.R")

# ── Detect cores (leave 2 free for the OS) ────────────────────────────────
n.cores <- max(1L, detectCores() - 2L)
cat("Using", n.cores, "cores\n")

# ── SETUP SOCKET CLUSTER (Cross-Platform) ──────────────────────────────────
cl <- makeCluster(n.cores)

# Initialize workers: load libraries and source AncCond.R on every node
clusterEvalQ(cl, {
  library(ape)
  library(phytools)
  library(parallel)
  library(R.utils)
  source("AncCond.R")
})

# ── Simulate continuous trait under state-dependent OU ────────────────────
.sim.OU.mapped <- function(tree, maps, theta, alpha, sigma, x0 = NULL) {
  n.tips  <- length(tree$tip.label)
  n.total <- n.tips + tree$Nnode
  root    <- n.tips + 1L
  edges   <- tree$edge
  
  vals <- numeric(n.total)
  
  # root state: stationary distribution around the starting optimum
  if (is.null(x0)) {
    root.edges <- which(edges[, 1] == root)
    root.state <- names(maps[[root.edges[1]]])[1]
    stat.var   <- sigma^2 / (2 * alpha)
    x0 <- rnorm(1, mean = theta[root.state], sd = sqrt(stat.var))
  }
  vals[root] <- x0
  
  for (k in seq_len(nrow(edges))) {
    pa <- edges[k, 1]
    ch <- edges[k, 2]
    x  <- vals[pa]
    
    seg    <- maps[[k]]
    states <- names(seg)
    for (j in seq_along(seg)) {
      dt <- seg[j]
      st <- states[j]
      th <- theta[st]
      e.at <- exp(-alpha * dt)
      mu   <- th + (x - th) * e.at
      v    <- (sigma^2 / (2 * alpha)) * (1 - exp(-2 * alpha * dt))
      x    <- rnorm(1, mean = mu, sd = sqrt(max(v, 0)))
    }
    vals[ch] <- x
  }
  
  out <- vals[1:n.tips]
  names(out) <- tree$tip.label
  out
}


# ── Simulation parameters ─────────────────────────────────────────────────
n.tips    <- 200
n.reps    <- 1000
n.maps.ac <- 100
n.sims.ac <- 500

# Mk rates for the discrete trait (symmetric, moderate)
q.rate <- 0.5
Q.mk   <- matrix(c(-q.rate, q.rate, q.rate, -q.rate), 2, 2,
                 dimnames = list(c("0","1"), c("0","1")))

# Parameter sweep: vary alpha (pull-back) and delta-theta
param.grid <- expand.grid(
  alpha    = c(0.5, 2, 8),
  delta.th = c(2, 4, 8),
  stringsAsFactors = FALSE
)
param.grid$sigma <- 1.0


# ── Worker function: one replicate ────────────────────────────────────────
.run.one.rep <- function(rep.id, a, dt, s, th) {
  set.seed(2026L * 1000L + rep.id)
  
  tree.i  <- pbtree(n = n.tips, scale = 1)
  h.i     <- sim.history(tree.i, Q.mk, nsim = 1)
  y.i     <- as.integer(h.i$states == "1")
  names(y.i) <- names(h.i$states)
  
  tree.cw <- reorder(tree.i, "cladewise")
  x.i     <- .sim.OU.mapped(tree.cw, h.i$maps, th, a, s)
  
  # Standard AncCond
  res.std <- tryCatch(
    withTimeout(
      AncCond(tree.i, x.i, y.i, n.maps = n.maps.ac, n.sims = n.sims.ac),
      timeout = 300
    ),
    TimeoutException = function(ex) {
      cat("    [TIMEOUT] AncCond std rep", rep.id, "\n")
      list(p.01 = NA_real_, p.combined = NA_real_)
    },
    error = function(e) list(p.01 = NA_real_, p.combined = NA_real_)
  )
  
  # Pruned AncCond (if feasible)
  n.zero <- sum(y.i == 0L)
  if (n.zero >= 10) {
    res.pr <- tryCatch(
      withTimeout(
        AncCond(tree.i, x.i, y.i, n.maps = n.maps.ac, n.sims = n.sims.ac,
                prune = TRUE),
        timeout = 300
      ),
      TimeoutException = function(ex) {
        cat("    [TIMEOUT] AncCond pruned rep", rep.id, "\n")
        list(p.01 = NA_real_, p.combined = NA_real_)
      },
      error = function(e) list(p.01 = NA_real_, p.combined = NA_real_)
    )
    data.frame(p.std       = res.std$p.01,
               p.pr        = res.pr$p.01,
               p.comb.std  = res.std$p.combined,
               p.comb.pr   = res.pr$p.combined,
               prune.ok    = TRUE)
  } else {
    data.frame(p.std       = res.std$p.01,
               p.pr        = NA_real_,
               p.comb.std  = res.std$p.combined,
               p.comb.pr   = NA_real_,
               prune.ok    = FALSE)
  }
}

# ── EXPORT TO CLUSTER ──────────────────────────────────────────────────────
# Send functions and parameters to worker nodes
clusterExport(cl, c(".sim.OU.mapped", ".run.one.rep", "n.tips", 
                    "n.maps.ac", "n.sims.ac", "Q.mk"))

# ── Chunked parallel runner ───────────────────────────────────────────────
chunk.size <- 25

.run.chunked <- function(rep.ids, ...) {
  # parLapply replaces mclapply for Windows/Mac/Linux support
  rep.list <- parLapply(cl, rep.ids, .run.one.rep, ...)
  
  bad <- vapply(rep.list, inherits, logical(1), "try-error")
  if (any(bad)) warning(sum(bad), " worker(s) failed in this chunk")
  do.call(rbind, rep.list[!bad])
}


# ── Main loop over parameter grid ────────────────────────────────────────
cat("========================================================\n")
cat("OU Adversarial Scenario: Discrete -> Continuous causation\n")
cat("========================================================\n")
cat("n.tips =", n.tips, "| n.reps =", n.reps, "\n")
cat("n.maps =", n.maps.ac, "| n.sims =", n.sims.ac, "\n")
cat("Cores:", n.cores, "| Chunk size:", chunk.size, "\n\n")

results <- vector("list", nrow(param.grid))

for (g in seq_len(nrow(param.grid))) {
  a  <- param.grid$alpha[g]
  dt <- param.grid$delta.th[g]
  s  <- param.grid$sigma[g]
  th <- c("0" = 0, "1" = dt)
  
  cat(sprintf("--- Grid %d/%d: alpha=%.1f, delta.theta=%.0f, sigma=%.1f ---\n",
              g, nrow(param.grid), a, dt, s))
  
  t.start <- proc.time()
  chunks  <- split(seq_len(n.reps), ceiling(seq_len(n.reps) / chunk.size))
  reps    <- NULL
  reps.raw <- NULL   # accumulate raw per-rep data for saving
  
  for (ci in seq_along(chunks)) {
    chunk.res <- .run.chunked(chunks[[ci]], a = a, dt = dt, s = s, th = th)
    reps <- rbind(reps, chunk.res)
    done <- nrow(reps)
    elapsed.so.far <- (proc.time() - t.start)["elapsed"]
    eta <- elapsed.so.far / done * (n.reps - done)
    cat(sprintf("    [%d/%d reps done | %.0fs elapsed | ~%.0fs remaining]\n",
                done, n.reps, elapsed.so.far, eta))
  }
  
  elapsed <- (proc.time() - t.start)["elapsed"]
  
  ok <- reps$prune.ok
  
  results[[g]] <- data.frame(
    alpha       = a,
    delta.theta = dt,
    sigma       = s,
    FP.standard = mean(reps$p.std < 0.05, na.rm = TRUE),
    FP.pruned   = mean(reps$p.pr[ok] < 0.05, na.rm = TRUE),
    n.prunable  = sum(ok),
    n.NA.std    = sum(is.na(reps$p.std)),
    mean.p.std  = mean(reps$p.std, na.rm = TRUE),
    mean.p.pr   = mean(reps$p.pr[ok], na.rm = TRUE)
  )
  
  # save raw per-rep data with condition labels
  reps$alpha    <- a
  reps$delta.th <- dt
  reps$sigma    <- s
  reps.raw <- rbind(reps.raw, reps)
  
  cat(sprintf("  FP (standard): %.3f | FP (pruned): %.3f  [%d/%d prunable, %d NA]\n",
              results[[g]]$FP.standard, results[[g]]$FP.pruned,
              sum(ok), n.reps, results[[g]]$n.NA.std))
  cat(sprintf("  Total: %.0f sec (%.1f sec/rep)\n\n", elapsed, elapsed / n.reps))
}

# Stop the cluster at the end of the script
stopCluster(cl)

res.df <- do.call(rbind, results)

cat("\n\n========================================\n")
cat("SUMMARY: False-positive rates (alpha=0.05, two-tailed)\n")
cat("========================================\n")
print(res.df, digits = 3, row.names = FALSE)

# ── Save results ──────────────────────────────────────────────────────────
saveRDS(res.df, "02_sim_OU_adversarial_results.rds")
write.csv(res.df, "02_sim_OU_adversarial_results.csv", row.names = FALSE)
write.csv(reps.raw, "02_sim_OU_adversarial_raw.csv", row.names = FALSE)
cat("\nResults saved to 02_sim_OU_adversarial_results.{rds,csv}\n")
cat("Raw per-rep data saved to 02_sim_OU_adversarial_raw.csv\n")