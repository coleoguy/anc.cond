# ============================================================
# 3 - model_misspecification.R
#
# Robustness of AncCond to model misspecification.
#
# Two axes of misspecification:
#   A) Discrete trait model: generate SYM, analyse ARD (& vice versa)
#   B) Continuous trait model: generate multi-optima OU, analyse BM
#
# Factorial design (16 cells x 200 replicates = 3,200 AncCond runs):
#
#   cont_model     : BM | OU
#   signal         : none | strong
#     BM  -> none = sf 1 (no branch scaling), strong = sf 5
#     OU  -> none = equal theta, strong = theta = +/- 0.5
#   disc_gen       : SYM | ARD  (Q used to GENERATE discrete data)
#   disc_analysis  : SYM | ARD  (mat passed to AncCond for analysis)
#
# Expected behaviour under correct specification:
#   signal = none  -> FP rate ~5%
#   signal = strong -> power >> 5%
#
# Misspecification inflates FP or reduces power if models are
# sensitive to the mismatch.
# ============================================================

library(ape)
library(phytools)
library(geiger)
library(TreeSim)
library(foreach)
library(doParallel)

set.seed(42)

# ---- Source AncCond ----
source("../R/anc_cond.R")

# ---- Parameters ----
n_taxa <- 100
n_reps <- 200

# BM continuous
bm_sig2 <- 0.2

# OU continuous
ou_alpha <- 2.0
ou_sig2  <- 0.5
theta_none  <- c("1" = 0.0,  "2" = 0.0)   # same optimum (null: no OU signal)
theta_split <- c("1" = -0.5, "2" = 0.5)    # separated optima (OU signal)

# Discrete Q matrices for data generation
Q_SYM <- matrix(c(-0.6, 0.6,
                    0.6, -0.6), 2, 2, byrow = TRUE,
                dimnames = list(c("1","2"), c("1","2")))
Q_ARD <- matrix(c(-0.3, 0.3,
                    0.9, -0.9), 2, 2, byrow = TRUE,
                dimnames = list(c("1","2"), c("1","2")))

# AncCond mat vectors for analysis
mat_SYM <- c(0, 1, 1, 0)   # equal rates
mat_ARD <- c(0, 2, 1, 0)   # all rates different

# Branch scaling
sf_none   <- 1
sf_strong <- 5

# AncCond tuning (moderate settings for speed)
ac_nsim <- 50
ac_iter <- 100

# ============================================================
# Cell definitions
# ============================================================
cells <- expand.grid(
  cont_model    = c("BM", "OU"),
  signal        = c("none", "strong"),
  disc_gen      = c("SYM", "ARD"),
  disc_analysis = c("SYM", "ARD"),
  stringsAsFactors = FALSE
)
cat(nrow(cells), "cells x", n_reps, "reps =",
    nrow(cells) * n_reps, "total AncCond runs\n\n")

# ============================================================
# Simulate BD trees (shared across all cells)
# ============================================================
message("Simulating ", n_reps, " BD trees with ", n_taxa, " tips...")
bd_trees <- sim.bd.taxa(n = n_taxa, numbsim = n_reps,
                        lambda = 3, mu = 1, complete = FALSE)

# ============================================================
# Helper functions
# ============================================================

#' Simulate OU trait evolution on a simmap tree.
#'
#' Walks root-to-tip through the mapped edge segments.
#' Each segment's regime determines the OU optimum (theta).
#'
#' @param simmap_tree A simmap (phylo + $maps) from sim.history.
#' @param theta Named vector of optima per state, e.g. c("1"=-0.5, "2"=0.5).
#' @param alpha OU pull-back strength.
#' @param sigma2 BM-like variance rate.
#' @param x0 Root value.
#' @return Named numeric vector of tip values.
sim_ou_on_simmap <- function(simmap_tree, theta, alpha, sigma2, x0 = 0) {
  n_tips  <- length(simmap_tree$tip.label)
  n_total <- n_tips + simmap_tree$Nnode
  vals    <- rep(NA_real_, n_total)
  vals[n_tips + 1L] <- x0

  edge <- simmap_tree$edge

  # Build processing order: level-by-level from root
  # so that every parent is computed before its children
  order <- integer(nrow(edge))
  filled <- 0L
  ready  <- which(edge[, 1] == (n_tips + 1L))
  while (length(ready) > 0) {
    for (ei in ready) {
      filled <- filled + 1L
      order[filled] <- ei
    }
    child_nodes <- edge[ready, 2]
    ready <- which(edge[, 1] %in% child_nodes)
  }

  # Walk through edges in pre-order
  for (i in order) {
    pa <- edge[i, 1]
    ch <- edge[i, 2]
    v  <- vals[pa]
    map <- simmap_tree$maps[[i]]
    for (s in seq_along(map)) {
      st <- names(map)[s]
      t  <- map[s]
      th <- theta[st]
      mn <- th + (v - th) * exp(-alpha * t)
      vr <- (sigma2 / (2 * alpha)) * (1 - exp(-2 * alpha * t))
      v  <- rnorm(1, mean = mn, sd = sqrt(max(0, vr)))
    }
    vals[ch] <- v
  }

  tip_v <- vals[seq_len(n_tips)]
  names(tip_v) <- simmap_tree$tip.label
  tip_v
}

#' Extract tip discrete states from a simmap tree.
#'
#' Reads the last regime segment of each terminal edge.
get_tip_states <- function(simmap_tree) {
  n <- length(simmap_tree$tip.label)
  st <- character(n)
  for (i in seq_len(nrow(simmap_tree$edge))) {
    ch <- simmap_tree$edge[i, 2]
    if (ch <= n) {
      m <- simmap_tree$maps[[i]]
      st[ch] <- names(m)[length(m)]
    }
  }
  names(st) <- simmap_tree$tip.label
  st
}

#' Scale tree branches by continuous-trait ancestral estimates.
#'
#' Mirrors the procedure in 1 - data_simulation.R:
#'   branches with mean > 75th percentile are multiplied by sf,
#'   branches with mean < 25th percentile are divided by sf.
scale_tree <- function(tree, cont_trait, sf) {
  if (sf == 1) return(tree)
  est    <- anc.ML(tree, cont_trait, model = "BM")
  n_tips <- length(tree$tip.label)
  av     <- numeric(n_tips + tree$Nnode)
  av[seq_len(n_tips)] <- cont_trait[tree$tip.label]
  av[as.integer(names(est$ace))] <- est$ace
  bm    <- (av[tree$edge[, 1]] + av[tree$edge[, 2]]) / 2
  q_lo  <- quantile(bm, 0.25)
  q_hi  <- quantile(bm, 0.75)
  tr    <- tree
  el    <- tree$edge.length
  el[bm >= q_hi] <- el[bm >= q_hi] * sf
  el[bm <= q_lo] <- el[bm <= q_lo] / sf
  tr$edge.length <- el
  tr
}

#' Simulate non-degenerate discrete tip states.
#'
#' Retries until both states are present with at least 10 tips each.
sim_disc_tips <- function(tree, Q, max_tries = 500) {
  root <- sample(1:2, 1)
  for (a in seq_len(max_tries)) {
    out <- sim.char(tree, par = Q, model = "discrete",
                    root = root, nsim = 1)
    dv  <- as.integer(out[, , 1])
    names(dv) <- tree$tip.label
    tab <- table(dv)
    if (length(tab) == 2L && min(tab) >= 10) {
      return(as.character(dv))
    }
  }
  as.character(dv)
}

#' Simulate non-degenerate discrete trait history (simmap).
#'
#' Used for OU cells: we need the full regime painting.
sim_disc_history <- function(tree, Q, max_tries = 500) {
  for (a in seq_len(max_tries)) {
    sm   <- sim.history(tree, Q = Q, nsim = 1, message = FALSE)
    tips <- get_tip_states(sm)
    tab  <- table(tips)
    if (length(tab) == 2L && min(tab) >= 10) return(sm)
  }
  sm
}

# ============================================================
# Parallel setup
# ============================================================
n_cores <- max(1, detectCores() - 2)
cl <- makeCluster(n_cores)
registerDoParallel(cl)
message(sprintf("Parallel cluster: %d cores\n", n_cores))

# ============================================================
# Main loop: iterate over cells, parallelise reps within cell
# ============================================================
all_results <- data.frame()

for (ci in seq_len(nrow(cells))) {
  cell <- cells[ci, ]
  message(sprintf("[Cell %d/%d] cont=%s | signal=%s | gen=%s -> ana=%s",
                  ci, nrow(cells),
                  cell$cont_model, cell$signal,
                  cell$disc_gen, cell$disc_analysis))

  # -- Look up generation / analysis settings for this cell --
  Q_gen        <- if (cell$disc_gen == "SYM") Q_SYM else Q_ARD
  mat_analysis <- if (cell$disc_analysis == "SYM") mat_SYM else mat_ARD
  sf           <- if (cell$cont_model == "BM" && cell$signal == "strong") sf_strong else sf_none
  ou_theta     <- if (cell$signal == "strong") theta_split else theta_none

  # -- Generate all 200 datasets for this cell (serial; fast) --
  cell_tasks <- vector("list", n_reps)

  for (ri in seq_len(n_reps)) {
    tree <- bd_trees[[ri]]

    if (cell$cont_model == "BM") {
      # BM pathway: continuous first, scale, then discrete on scaled tree
      ct      <- fastBM(tree, sig2 = bm_sig2, a = 0)
      tr_s    <- scale_tree(tree, ct, sf)
      disc    <- sim_disc_tips(tr_s, Q_gen)
    } else {
      # OU pathway: discrete history first, then OU on simmap
      sm   <- sim_disc_history(tree, Q_gen)
      ct   <- sim_ou_on_simmap(sm, ou_theta, ou_alpha, ou_sig2, x0 = 0)
      disc <- get_tip_states(sm)
    }

    df <- data.frame(
      tip  = tree$tip.label,
      cont = ct[tree$tip.label],
      disc = disc[tree$tip.label],
      stringsAsFactors = FALSE
    )
    cell_tasks[[ri]] <- list(tree = tree, df = df, mat = mat_analysis)
  }

  # -- Run AncCond in parallel across 200 reps --
  cell_results <- foreach(
    t = cell_tasks,
    .packages = c("ape", "phytools", "diversitree"),
    .errorhandling = "pass"
  ) %dopar% {
    source("../R/anc_cond.R")
    res <- tryCatch(
      AncCond(tree  = t$tree,
              data  = t$df,
              mat   = t$mat,
              nsim  = 50,     # ac_nsim
              iter  = 100,    # ac_iter
              ncores = 1L),
      error = function(e) list(pvals = c(`12` = NA_real_, `21` = NA_real_))
    )
    c(p12 = unname(res$pvals["12"]),
      p21 = unname(res$pvals["21"]))
  }

  # -- Compile results for this cell --
  pvals <- do.call(rbind, cell_results)
  cell_df <- data.frame(
    cell_id       = ci,
    cont_model    = cell$cont_model,
    signal        = cell$signal,
    disc_gen      = cell$disc_gen,
    disc_analysis = cell$disc_analysis,
    rep           = seq_len(n_reps),
    p12           = pvals[, "p12"],
    p21           = pvals[, "p21"]
  )
  all_results <- rbind(all_results, cell_df)

  # -- Print cell summary --
  pct12 <- mean(cell_df$p12 < 0.05, na.rm = TRUE) * 100
  pct21 <- mean(cell_df$p21 < 0.05, na.rm = TRUE) * 100
  message(sprintf("  -> %%sig 1->2: %5.1f%%   2->1: %5.1f%%", pct12, pct21))

  # -- Checkpoint save after every cell --
  save(all_results, cells,
       file = "../results/misspec_checkpoint.RData")
}

stopCluster(cl)

# ============================================================
# Final summary
# ============================================================
cat("\n========== Model Misspecification Results ==========\n\n")
cat(sprintf("%-4s  %-5s  %-7s  %-5s  %-5s  %8s  %8s\n",
            "Cell", "Cont", "Signal", "Gen", "Ana", "%sig 12", "%sig 21"))
cat(strrep("-", 55), "\n")

for (ci in seq_len(nrow(cells))) {
  sub   <- all_results[all_results$cell_id == ci, ]
  pct12 <- mean(sub$p12 < 0.05, na.rm = TRUE) * 100
  pct21 <- mean(sub$p21 < 0.05, na.rm = TRUE) * 100
  cat(sprintf("%-4d  %-5s  %-7s  %-5s  %-5s  %7.1f%%  %7.1f%%\n",
              ci, cells$cont_model[ci], cells$signal[ci],
              cells$disc_gen[ci], cells$disc_analysis[ci],
              pct12, pct21))
}

# ---- Save final results ----
misspec_results <- all_results
save(misspec_results, cells,
     file = "../results/misspecification_results.RData")
message("\nFinal results saved to results/misspecification_results.RData")
