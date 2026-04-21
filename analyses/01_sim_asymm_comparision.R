# ============================================================================
# Asymmetric Power & False-Positive Simulation for AncCond
# Changes: min.state = 10% of taxa, Root state fixed to 0.
# ============================================================================

library(ape)
library(phytools)
library(phylolm)
library(nlme)
library(parallel)

# ── 1. Setup Parallel Environment ──────────────────────────────────────────
n.cores <- max(1L, detectCores() - 2L)
cat("Using", n.cores, "cores for Asymmetric Simulation\n")

cl <- makeCluster(n.cores)

clusterEvalQ(cl, {
  library(ape)
  library(phytools)
  library(phylolm)
  library(nlme)
  source("AncCond.R")
})

# ── 2. Parameters (Asymmetric Grid) ────────────────────────────────────────
n.tips.vec    <- c(50, 100, 200, 500)
q01.vec       <- c(0.1, 0.5, 1.0, 2.0) 
q10.vec       <- c(0.1) 
s.vec         <- c(1, 2, 3, 5, 10)     
n.reps        <- 200
n.maps        <- 100
n.sims        <- 500
alpha_level   <- 0.05
min.state.pct <- 0.10  # Changed to 10%
max.reject    <- 500L  

grid <- expand.grid(n.tips = n.tips.vec,
                    q01    = q01.vec,
                    q10    = q10.vec,
                    s      = s.vec,
                    stringsAsFactors = FALSE)

clusterExport(cl, c("n.maps", "n.sims", "min.state.pct", "max.reject"))
clusterSetRNGStream(cl, iseed = 2026)

# ── 3. Helper Functions ─────────────────────────────────────────────────────

.count.trans <- function(h) {
  sum(vapply(h$maps, function(m) length(m) - 1L, integer(1)))
}

.sim.dataset <- function(n.tips, q01, q10, s, min_pct, max.reject) {
  Q <- matrix(c(-q01, q01, q10, -q10), 2, 2, byrow = TRUE)
  rownames(Q) <- colnames(Q) <- c("0", "1")
  
  # Dynamic threshold based on percentage
  threshold <- ceiling(n.tips * min_pct)
  
  n.discard <- 0L
  repeat {
    if (n.discard >= max.reject) return(NULL)
    
    tree <- pbtree(n = n.tips, scale = 1)
    x <- fastBM(tree)
    anc.x <- fastAnc(tree, x)
    all.vals <- c(x[tree$tip.label], anc.x)
    names(all.vals)[seq_along(tree$tip.label)] <- seq_along(tree$tip.label)
    
    branch.avg <- numeric(nrow(tree$edge))
    for (k in seq_len(nrow(tree$edge))) {
      pa <- tree$edge[k, 1]; ch <- tree$edge[k, 2]
      branch.avg[k] <- (all.vals[as.character(pa)] + all.vals[as.character(ch)]) / 2
    }
    
    tree.scaled <- tree
    direction <- sample(c(1L, -1L), 1)
    if (s != 1) {
      q25 <- quantile(branch.avg, 0.25); q75 <- quantile(branch.avg, 0.75)
      if (direction == 1L) {
        tree.scaled$edge.length[branch.avg >= q75] <- tree.scaled$edge.length[branch.avg >= q75] * s
        tree.scaled$edge.length[branch.avg <= q25] <- tree.scaled$edge.length[branch.avg <= q25] / s
      } else {
        tree.scaled$edge.length[branch.avg <= q25] <- tree.scaled$edge.length[branch.avg <= q25] * s
        tree.scaled$edge.length[branch.avg >= q75] <- tree.scaled$edge.length[branch.avg >= q75] / s
      }
    }
    
    # ROOTED TO STATE ZERO: anc = "0"
    h <- sim.history(tree.scaled, Q, anc = "0", nsim = 1, message = FALSE)
    y <- as.integer(h$states == "1")
    names(y) <- names(h$states)
    
    if (min(sum(y == 0), sum(y == 1)) >= threshold) {
      return(list(tree=tree, x=x, y=y, n.trans=.count.trans(h), n.discard=n.discard, Q=Q))
    }
    n.discard <- n.discard + 1L
  }
}

run.one <- function(n.tips, q01, q10, s, rep.id, min_pct) {
  dat <- .sim.dataset(n.tips, q01, q10, s, min_pct, max.reject)
  if (is.null(dat)) return(NULL)
  
  res.ac <- tryCatch(
    suppressWarnings(AncCond(dat$tree, dat$x, dat$y, n.maps=n.maps, n.sims=n.sims, Q=dat$Q)),
    error = function(e) list(p.01=NA, p.10=NA)
  )
  
  p.plg <- tryCatch({
    fit <- phyloglm(y ~ x, data=data.frame(y=dat$y, x=dat$x), phy=dat$tree, method="logistic_MPLE")
    coef(summary(fit))["x", "p.value"]
  }, error = function(e) NA)
  
  p.pgls <- tryCatch({
    df <- data.frame(x=dat$x, y=as.factor(dat$y), row.names=dat$tree$tip.label)
    fit <- gls(x ~ y, data=df, correlation=corBrownian(phy=dat$tree))
    summary(fit)$tTable["y1", "p-value"]
  }, error = function(e) NA)
  
  data.frame(n.tips=n.tips, q01=q01, q10=q10, s=s, rep=rep.id,
             p.ac01=res.ac$p.01, p.ac10=res.ac$p.10, p.plg=p.plg, p.pgls=p.pgls,
             n.trans=dat$n.trans)
}

clusterExport(cl, c(".count.trans", ".sim.dataset", "run.one"))

# ── 4. Main Execution Loop ─────────────────────────────────────────────────
cat("\n=== Starting Asymmetric Run ===\n")
all.raw <- list()

for (g in seq_len(nrow(grid))) {
  nt <- grid$n.tips[g]; q01 <- grid$q01[g]; q10 <- grid$q10[g]; ss <- grid$s[g]
  cat(sprintf("[%d/%d] n=%d q01=%.1f q10=%.1f s=%d... ", g, nrow(grid), nt, q01, q10, ss))
  
  results_list <- parLapply(cl, 1:n.reps, function(i, nt_in, q01_in, q10_in, ss_in, mp_in) {
    run.one(nt_in, q01_in, q10_in, ss_in, i, mp_in)
  }, nt_in=nt, q01_in=q01, q10_in=q10, ss_in=ss, mp_in=min.state.pct)
  
  all.raw[[g]] <- do.call(rbind, results_list)
  cat("Done.\n")
}

stopCluster(cl)

# ── 5. Base R Summarization & File Output ──────────────────────────────────
raw_df <- do.call(rbind, all.raw)

# Calculate rejection rates grouped by conditions using aggregate()
summary_df <- aggregate(
  cbind(p.ac01, p.ac10, p.plg, p.pgls) ~ n.tips + q01 + q10 + s, 
  data = raw_df, 
  FUN = function(x) mean(x < alpha_level, na.rm = TRUE),
  na.action = na.pass
)

# Rename columns to match your preferred format
colnames(summary_df) <- c("n.tips", "q01", "q10", "s", "anccond.01", "anccond.10", "phyloglm", "pgls")

write.csv(raw_df, "01_sim_asymm_power_raw.csv", row.names = FALSE)
write.csv(summary_df, "01_sim_asymm_power_results.csv", row.names = FALSE)

cat("\nSimulation Complete. Files saved as 01_sim_asymm_power_*.csv\n")