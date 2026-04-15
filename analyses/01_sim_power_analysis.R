# ============================================================================
# Main Power & False-Positive Simulation for AncCond
#
# Factorial design crossing:
#   - Tree size:        n.tips  = 50, 100, 200, 500
#   - Discrete rate:    q.rate  = 0.1, 0.5, 1.0, 2.0  (symmetric Mk)
#   - Signal strength:  s       = 1 (null), 2, 3, 5, 10
#
# Signal model (ancestral condition):
#   1. Simulate continuous trait x by BM on the original tree.
#   2. Reconstruct ancestral x at internal nodes; compute the average
#      x value along each branch = (parent_x + child_x) / 2.
#   3. Signal direction is randomly flipped per replicate (+1 or -1).
#      If +1: upper 25% branches get length * s, lower 25% get / s.
#      If -1: lower 25% branches get length * s, upper 25% get / s.
#      Middle 50% always unchanged.
#   4. Simulate discrete trait y under Mk on this stretched tree.
#   This creates signal: transitions concentrate where x is extreme
#   in the chosen direction. Random flipping ensures a fair comparison
#   to phyloglm/PGLS (inherently two-sided).
#   s = 1 gives the null (uniform scaling = no distortion).
#
# Rejection sampling: only datasets with at least min.state tips in each
# discrete state are accepted.  Discarded datasets are counted per
# condition and saved to 01_sim_discard_log.csv.
#
# Each condition: 200 accepted replicates.
# Reports:  power (prop. significant at alpha = 0.05) for s > 1
#           false-positive rate for s = 1
#
# Compares three methods on identical simulated datasets:
#   1. AncCond  — ancestral condition test (this package)
#   2. phyloglm — phylogenetic logistic regression (y ~ x)
#   3. PGLS     — phylogenetic GLS with BM correlation (x ~ y)
#
# Parallelised with parallel::mclapply (fork-safe on macOS/Linux).
# Output:  01_sim_power_raw.csv      (one row per replicate)
#          01_sim_power_results.csv   (one row per condition × method)
#          01_sim_discard_log.csv     (one row per condition: discard counts)
#
# Expected runtime: ~3-8 hours on 14 cores depending on n.tips = 500 runs.
#
# Dependencies: ape, phytools, phylolm, nlme, parallel
# ============================================================================

library(ape)
library(phytools)
library(phylolm)
library(nlme)
library(parallel)
source("AncCond.R")

# ── Detect cores ────────────────────────────────────────────────────────────
n.cores <- max(1L, detectCores() - 2L)
cat("Using", n.cores, "cores\n")

# ── Simulation parameters ───────────────────────────────────────────────────
n.tips.vec  <- c(50, 100, 200, 500)
q.rate.vec  <- c(0.1, 0.5, 1.0, 2.0)
s.vec       <- c(1, 2, 3, 5, 10)       # s = 1 is null
n.reps      <- 200
n.maps      <- 100
n.sims      <- 500
alpha       <- 0.05
min.state   <- 3L    # minimum tips in the minor state to accept a dataset
max.reject  <- 500L  # safety valve: max rejections per replicate before giving up

# Build parameter grid
grid <- expand.grid(n.tips = n.tips.vec,
                    q.rate = q.rate.vec,
                    s      = s.vec,
                    stringsAsFactors = FALSE)
cat(sprintf("Parameter grid: %d conditions x %d reps = %d total runs\n",
            nrow(grid), n.reps, nrow(grid) * n.reps))
cat(sprintf("Rejection threshold: min(n0, n1) >= %d\n", min.state))


# ── Helper: count true transitions from sim.history object ────────────────
# h$maps is a list (one per edge) of named numeric vectors.
# Names are states; a state change = a transition.
# Number of transitions on edge k = length(h$maps[[k]]) - 1.
.count.trans <- function(h) {
  sum(vapply(h$maps, function(m) length(m) - 1L, integer(1)))
}


# ── Helper: simulate one dataset (tree + x + y) with rejection ───────────
# Resamples until min(n0, n1) >= min.state.
# Returns list with tree, x, y, n0, n1, n.true.trans, n.discard, direction.
# Returns NULL if max.reject attempts exhausted (should not happen in practice).
.sim.dataset <- function(n.tips, q.rate, s, min.state, max.reject = 500L) {
  Q <- matrix(c(-q.rate, q.rate, q.rate, -q.rate), 2, 2)
  n.discard <- 0L

  repeat {
    if (n.discard >= max.reject) return(NULL)
    # 1. Pure-birth tree
    tree <- pbtree(n = n.tips, scale = 1)

    # 2. Continuous trait by BM on the original tree
    x <- fastBM(tree)

    # 3. Reconstruct ancestral x values at internal nodes
    anc.x <- fastAnc(tree, x)
    n.edges <- nrow(tree$edge)
    all.vals <- c(x[tree$tip.label], anc.x)
    names(all.vals)[seq_along(tree$tip.label)] <- seq_along(tree$tip.label)
    branch.avg <- numeric(n.edges)
    for (k in seq_len(n.edges)) {
      pa <- tree$edge[k, 1]
      ch <- tree$edge[k, 2]
      branch.avg[k] <- (all.vals[as.character(pa)] + all.vals[as.character(ch)]) / 2
    }

    # 4. Scale tree: randomly flip signal direction per replicate.
    tree.scaled <- tree
    direction <- sample(c(1L, -1L), 1)
    if (s != 1) {
      q25 <- quantile(branch.avg, 0.25)
      q75 <- quantile(branch.avg, 0.75)
      upper <- branch.avg >= q75
      lower <- branch.avg <= q25
      if (direction == 1L) {
        tree.scaled$edge.length[upper] <- tree$edge.length[upper] * s
        tree.scaled$edge.length[lower] <- tree$edge.length[lower] / s
      } else {
        tree.scaled$edge.length[lower] <- tree$edge.length[lower] * s
        tree.scaled$edge.length[upper] <- tree$edge.length[upper] / s
      }
    }

    # 5. Simulate discrete trait y under Mk on the stretched tree
    h <- sim.history(tree.scaled, Q, nsim = 1)
    y <- as.integer(h$states == "1")
    names(y) <- names(h$states)

    n0 <- sum(y == 0)
    n1 <- sum(y == 1)

    # Accept or reject
    if (min(n0, n1) >= min.state) {
      return(list(tree      = tree,
                  x         = x,
                  y         = y,
                  n0        = n0,
                  n1        = n1,
                  n.trans   = .count.trans(h),
                  n.discard = n.discard,
                  direction = direction,
                  Q         = Q))
    }
    n.discard <- n.discard + 1L
  }
}


# ── Helper: run one replicate ───────────────────────────────────────────────
# Returns a one-row data.frame with results.
run.one <- function(n.tips, q.rate, s, rep.id) {

  # Simulate dataset with rejection sampling
  dat <- .sim.dataset(n.tips, q.rate, s, min.state, max.reject)
  if (is.null(dat)) {
    # Safety valve: could not produce a valid dataset
    return(data.frame(
      n.tips = n.tips, q.rate = q.rate, s = s, direction = NA,
      rep = rep.id, n0 = NA, n1 = NA, n.true.trans = NA,
      n.discard = max.reject,
      T.obs.01 = NA, T.obs.10 = NA, p.anccond.01 = NA, p.anccond.10 = NA,
      p.upper.01 = NA, p.lower.01 = NA, p.upper.10 = NA, p.lower.10 = NA,
      p.phyloglm = NA, p.pgls = NA, stringsAsFactors = FALSE))
  }
  tree <- dat$tree;  x <- dat$x;  y <- dat$y;  Q <- dat$Q

  # Run AncCond on the ORIGINAL tree with x and y
  res.ac <- tryCatch(
    suppressWarnings(
      AncCond(tree, x, y,
              n.maps = n.maps,
              n.sims = n.sims,
              Q      = Q)
    ),
    error = function(e) list(T.obs.01 = NA, p.01 = NA, p.upper.01 = NA,
                             p.lower.01 = NA, T.obs.10 = NA, p.10 = NA,
                             p.upper.10 = NA, p.lower.10 = NA)
  )

  # Phylogenetic logistic regression: y ~ x
  p.phyloglm <- NA
  tryCatch({
    df.plg <- data.frame(y = y[tree$tip.label], x = x[tree$tip.label])
    fit.plg <- suppressWarnings(
      phyloglm(y ~ x, data = df.plg, phy = tree, method = "logistic_MPLE",
               boot = 0)
    )
    p.phyloglm <- coef(summary(fit.plg))["x", "p.value"]
  }, error = function(e) {})

  # PGLS: x ~ y  (continuous ~ discrete, BM correlation structure)
  p.pgls <- NA
  tryCatch({
    df.gls <- data.frame(x = x[tree$tip.label],
                         y = as.factor(y[tree$tip.label]),
                         row.names = tree$tip.label)
    fit.gls <- gls(x ~ y, data = df.gls,
                   correlation = corBrownian(phy = tree))
    p.pgls <- summary(fit.gls)$tTable["y1", "p-value"]
  }, error = function(e) {})

  data.frame(
    n.tips             = n.tips,
    q.rate             = q.rate,
    s                  = s,
    direction          = dat$direction,
    rep                = rep.id,
    n0                 = dat$n0,
    n1                 = dat$n1,
    n.true.trans       = dat$n.trans,
    n.discard          = dat$n.discard,
    T.obs.01           = res.ac$T.obs.01,
    T.obs.10           = res.ac$T.obs.10,
    p.anccond.01       = res.ac$p.01,
    p.anccond.10       = res.ac$p.10,
    p.upper.01         = res.ac$p.upper.01,
    p.lower.01         = res.ac$p.lower.01,
    p.upper.10         = res.ac$p.upper.10,
    p.lower.10         = res.ac$p.lower.10,
    p.phyloglm         = p.phyloglm,
    p.pgls             = p.pgls,
    stringsAsFactors   = FALSE
  )
}


# ── Main loop over conditions ───────────────────────────────────────────────
# Runs reps in batches of n.cores so progress prints between batches.
cat("\n=== Starting simulations ===\n")
all.raw     <- vector("list", nrow(grid))
discard.log <- vector("list", nrow(grid))

batch.size <- n.cores   # one batch = one mclapply call

t.total <- system.time({
  for (g in seq_len(nrow(grid))) {
    nt <- grid$n.tips[g]
    qr <- grid$q.rate[g]
    ss <- grid$s[g]

    cat(sprintf("\n[%3d/%d] n.tips=%3d  q=%.1f  s=%2d\n",
                g, nrow(grid), nt, qr, ss))

    cond.results <- vector("list", n.reps)
    done <- 0L
    cum.discard <- 0L

    t.cond <- system.time({
      while (done < n.reps) {
        this.batch <- min(batch.size, n.reps - done)
        ids <- (done + 1L):(done + this.batch)

        batch.out <- mclapply(ids, function(i) {
          run.one(nt, qr, ss, i)
        }, mc.cores = n.cores)

        for (j in seq_along(ids)) cond.results[[ids[j]]] <- batch.out[[j]]

        batch.df    <- do.call(rbind, batch.out)
        cum.discard <- cum.discard + sum(batch.df$n.discard, na.rm = TRUE)
        done        <- done + this.batch

        cat(sprintf("  %3d/%d done  |  rejected so far: %d\n",
                    done, n.reps, cum.discard))
      }
    })

    raw.df <- do.call(rbind, cond.results)
    all.raw[[g]] <- raw.df

    # discard summary for this condition
    total.discard <- sum(raw.df$n.discard, na.rm = TRUE)
    mean.trans    <- mean(raw.df$n.true.trans, na.rm = TRUE)
    discard.log[[g]] <- data.frame(
      n.tips         = nt,
      q.rate         = qr,
      s              = ss,
      n.reps         = n.reps,
      total.discard  = total.discard,
      mean.discard   = total.discard / n.reps,
      mean.true.trans = mean.trans,
      stringsAsFactors = FALSE
    )

    # quick summary (AncCond 0->1)
    p.01    <- raw.df$p.anccond.01
    valid   <- !is.na(p.01)
    n.valid <- sum(valid)
    n.sig   <- sum(p.01[valid] < alpha)

    cat(sprintf("  DONE  discard=%d  trans=%.1f  sig=%d/%d (%.1f%%)  [%.0fs]\n",
                total.discard, mean.trans, n.sig, n.valid,
                if (n.valid > 0) 100 * n.sig / n.valid else NA,
                t.cond["elapsed"]))
  }
})

cat(sprintf("\nTotal time: %.1f minutes\n", t.total["elapsed"] / 60))


# ── Combine raw results ────────────────────────────────────────────────────
raw.all     <- do.call(rbind, all.raw)
discard.all <- do.call(rbind, discard.log)


# ── Summarise by condition × method ──────────────────────────────────────────
methods <- c("anccond.01", "anccond.10", "phyloglm", "pgls")
p.cols  <- c("p.anccond.01", "p.anccond.10", "p.phyloglm", "p.pgls")

summary.list <- list()
for (g in seq_len(nrow(grid))) {
  sub <- raw.all[raw.all$n.tips == grid$n.tips[g] &
                  raw.all$q.rate == grid$q.rate[g] &
                  raw.all$s      == grid$s[g], ]
  for (m in seq_along(methods)) {
    pv      <- sub[[p.cols[m]]]
    valid   <- !is.na(pv)
    n.valid <- sum(valid)
    n.sig   <- sum(pv[valid] < alpha)

    summary.list[[length(summary.list) + 1]] <- data.frame(
      n.tips    = grid$n.tips[g],
      q.rate    = grid$q.rate[g],
      s         = grid$s[g],
      method    = methods[m],
      n.reps    = n.reps,
      n.valid   = n.valid,
      n.NA      = sum(!valid),
      n.sig     = n.sig,
      rate      = if (n.valid > 0) n.sig / n.valid else NA,
      stringsAsFactors = FALSE
    )
  }
}
summary.df <- do.call(rbind, summary.list)
summary.df$metric <- ifelse(summary.df$s == 1, "FPR", "power")


# ── Save outputs ────────────────────────────────────────────────────────────
write.csv(raw.all,     "01_sim_power_raw.csv",      row.names = FALSE)
write.csv(summary.df,  "01_sim_power_results.csv",  row.names = FALSE)
write.csv(discard.all, "01_sim_discard_log.csv",     row.names = FALSE)

cat("\n=== Results saved ===\n")
cat("  01_sim_power_raw.csv      — one row per replicate\n")
cat("  01_sim_power_results.csv  — one row per condition × method\n")
cat("  01_sim_discard_log.csv    — rejection counts per condition\n")


# ── Print summary table ────────────────────────────────────────────────────
cat("\n=== Summary (all methods) ===\n")
cat(sprintf("%-6s %-6s %-4s  %-10s %-6s %-4s %-6s  %s\n",
            "n.tips", "q.rate", "s", "method", "valid", "NA", "rate", "metric"))
cat(strrep("-", 65), "\n")
for (i in seq_len(nrow(summary.df))) {
  with(summary.df[i, ], {
    cat(sprintf("%-6d %-6.1f %-4d  %-10s %-6d %-4d %-6.3f  %s\n",
                n.tips, q.rate, s, method, n.valid, n.NA, rate, metric))
  })
}


# ── Discard summary ────────────────────────────────────────────────────────
cat("\n=== Rejection sampling: discarded datasets ===\n")
cat(sprintf("%-6s %-6s %-4s  %-8s %-10s %-10s\n",
            "n.tips", "q.rate", "s", "discard", "mean.disc", "mean.trans"))
cat(strrep("-", 55), "\n")
for (i in seq_len(nrow(discard.all))) {
  with(discard.all[i, ], {
    cat(sprintf("%-6d %-6.1f %-4d  %-8d %-10.1f %-10.1f\n",
                n.tips, q.rate, s, total.discard, mean.discard, mean.true.trans))
  })
}


# ── Quick sanity check: null FPR should be near alpha for all methods ───────
cat("\n=== Null (s=1) FPR by method ===\n")
for (m in methods) {
  null.m <- summary.df[summary.df$s == 1 & summary.df$method == m, ]
  cat(sprintf("  %-10s  mean=%.3f  range=%.3f-%.3f\n",
              m,
              mean(null.m$rate, na.rm = TRUE),
              min(null.m$rate, na.rm = TRUE),
              max(null.m$rate, na.rm = TRUE)))
}

cat("\n=== Done ===\n")
