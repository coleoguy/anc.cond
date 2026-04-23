library(ape)
library(phytools)
library(parallel)
library(R.utils)
library(phylolm)

# ── 1. THE HYBRID SIMULATION ENGINE ───────────────────────────────────────
# State 0: Brownian Motion (Neutral Drift)
# State 1: Ornstein-Uhlenbeck (Selection toward theta_1)
.sim.hybrid.mapped <- function(tree, maps, theta_1, alpha, sigma, root_state) {
  n.tips  <- length(tree$tip.label)
  n.total <- n.tips + tree$Nnode
  root    <- n.tips + 1L
  edges   <- tree$edge
  vals    <- numeric(n.total)
  
  # Initialize root: Start at 0 (BM) or theta_1 (OU)
  vals[root] <- if(root_state == "1") theta_1 else 0
  
  for (k in seq_len(nrow(edges))) {
    pa <- edges[k, 1]; ch <- edges[k, 2]; x  <- vals[pa]
    
    seg <- maps[[k]]; states <- names(seg)
    
    for (j in seq_along(seg)) {
      dt <- seg[j]; st <- states[j]
      
      if (st == "1") {
        # OU transition: Pull toward theta_1
        e.at <- exp(-alpha * dt)
        mu   <- theta_1 + (x - theta_1) * e.at
        v    <- (sigma^2 / (2 * alpha)) * (1 - exp(-2 * alpha * dt))
        x    <- rnorm(1, mean = mu, sd = sqrt(max(v, 0)))
      } else {
        # BM transition: Random walk from current position
        x    <- rnorm(1, mean = x, sd = sqrt(sigma^2 * dt))
      }
    }
    vals[ch] <- x
  }
  
  out <- vals[1:n.tips]
  names(out) <- tree$tip.label
  return(out)
}

# ── 2. MASTER WORKER FUNCTION ─────────────────────────────────────────────
# This function executes a single replicate for a given parameter set
.run.simulation.cell <- function(idx, grid) {
  # Extract parameters
  p         <- grid[idx, ]
  n.tips    <- 200
  q.rate    <- 0.5
  n.maps.ac <- 100
  n.sims.ac <- 500
  sigma_val <- 1.0
  
  # Unique seed for this cell
  set.seed(2026L * 1000L + idx)
  
  # A. Generate Tree & History
  tree <- pbtree(n = n.tips, scale = 1)
  Q <- matrix(c(-q.rate, q.rate, q.rate, -q.rate), 2, 2, 
              dimnames = list(c("0","1"), c("0","1")))
  
  h <- sim.history(tree, Q, nsim = 1, message = FALSE, pi = c(0.5, 0.5))
  
  # B. Simulate Continuous Trait
  root_node <- n.tips + 1
  root_st <- names(h$maps[[which(tree$edge[,1] == root_node)[1]]])[1]
  x_raw <- .sim.hybrid.mapped(reorder(tree, "cladewise"), h$maps, 
                              theta_1 = p$delta.th, alpha = p$alpha, 
                              sigma = sigma_val, root_state = root_st)
  
  # C. DATA ALIGNMENT FIX
  # Explicitly index x and y by tree tip labels to prevent 'scrambled' data
  y_raw <- as.integer(h$states == "1")
  names(y_raw) <- names(h$states)
  
  y <- y_raw[tree$tip.label]
  x <- x_raw[tree$tip.label]
  
  # D. STATISTICAL TESTS
  
  # 1. Forward Test: phylolm (Trait ~ State)
  # Proves the discrete state impacts the continuous trait
  p.forward <- tryCatch({
    fit_f <- phylolm(x ~ as.factor(y), phy = tree, model = "BM")
    summary(fit_f)$coefficients["as.factor(y)1", "p.value"]
  }, error = function(e) NA_real_)
  
  # 2. Standard AncCond
  res.std <- tryCatch(
    withTimeout(AncCond(tree, x, y, n.maps = n.maps.ac, n.sims = n.sims.ac), timeout = 300),
    error = function(e) list(p.01 = NA_real_, p.combined = NA_real_)
  )
  
  # 3. Pruned AncCond (Condition: n.zero >= 10)
  n.zero <- sum(y == 0L)
  if (n.zero >= 10) {
    res.pr <- tryCatch(
      withTimeout(AncCond(tree, x, y, n.maps = n.maps.ac, n.sims = n.sims.ac, prune = TRUE), timeout = 300),
      error = function(e) list(p.01 = NA_real_, p.combined = NA_real_)
    )
    prune_ok <- TRUE
  } else {
    res.pr   <- list(p.01 = NA_real_, p.combined = NA_real_)
    prune_ok <- FALSE
  }
  
  # 4. phyloglm (State ~ Trait)
  p.plg <- tryCatch({
    fit_b <- phyloglm(y ~ x, data = data.frame(y=y, x=x), phy = tree, method = "logistic_MPLE")
    summary(fit_b)$coefficients["x", "p.value"]
  }, error = function(e) NA_real_)
  
  # E. Compile Results
  return(data.frame(
    alpha      = p$alpha, 
    delta.th   = p$delta.th, 
    rep        = p$rep,
    n.zero     = n.zero,
    p.forward  = p.forward,
    p.std      = res.std$p.01,
    p.pr       = res.pr$p.01,
    p.comb.std = res.std$p.combined,
    p.comb.pr  = res.pr$p.combined,
    p.phyloglm = p.plg,
    prune.ok   = prune_ok,
    stringsAsFactors = FALSE
  ))
}

# ── 3. EXECUTION SETUP ────────────────────────────────────────────────────
n.reps     <- 1000
param.grid <- expand.grid(
  alpha    = c(0.5, 2, 8),
  delta.th = c(2, 4, 8),
  rep      = 1:n.reps,
  stringsAsFactors = FALSE
)

# Detect Cores (Cross-platform)
n.cores <- max(1L, detectCores() - 2L)
cl <- makeCluster(n.cores)

# Export Functions and Global libraries to Workers
clusterExport(cl, c(".sim.hybrid.mapped", ".run.simulation.cell", "param.grid"))
clusterEvalQ(cl, {
  library(ape)
  library(phytools)
  library(R.utils)
  library(phylolm)
  source("AncCond.R") 
})

# ── 4. RUN FULL SIMULATION ────────────────────────────────────────────────
cat("Starting parallel grid search for", nrow(param.grid), "simulations...\n")
full_results_list <- parLapply(cl, 1:nrow(param.grid), .run.simulation.cell, grid = param.grid)
stopCluster(cl)

# ── 5. FINAL DATA PROCESSING ──────────────────────────────────────────────
final_results <- do.call(rbind, full_results_list)
save(final_results, file = "Hybrid_FullPipeline_Results_2026.RData")
cat("Simulation Complete. Results saved to Hybrid_FullPipeline_Results_2026.RData\n")

# Simple power summary for console verification
cat("\nSummary Power (Alpha = 0.05):\n")
final_results$sig_forward <- final_results$p.forward < 0.05
final_results$sig_std     <- final_results$p.std < 0.05
final_results$sig_plg     <- final_results$p.phyloglm < 0.05

aggregate(cbind(sig_forward, sig_std, sig_plg) ~ alpha + delta.th, 
          data = final_results, FUN = mean)