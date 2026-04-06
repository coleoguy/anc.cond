# ============================================================================
# Main Power & False-Positive Simulation for AncCond
# ============================================================================

library(ape)
library(phytools)
library(phylolm)
library(nlme)
library(parallel)
library(foreach)
library(doParallel)
source("AncCondFast.R") # Sourcing the Fast version

# ── Detect environment and cores ────────────────────────────────────────────
n.cores <- max(1L, detectCores() - 2L)
is_windows <- .Platform$OS.type == "windows"
cat("Using", n.cores, "cores on", .Platform$OS.type, "\n")

# ── Simulation parameters ───────────────────────────────────────────────────
n.tips.vec  <- c(50, 100)
q.rate.vec  <- c(0.1, 0.5)
s.vec       <- c(1, 5)
n.reps      <- 50 
n.maps      <- 100 # Maps to n_maps
n.sims      <- 200 # Maps to n_iter (Bootstrap is slower, 200 is a good balance)
alpha       <- 0.05
min.state   <- 3L
max.reject  <- 500L

grid <- expand.grid(n.tips = n.tips.vec, q.rate = q.rate.vec, s = s.vec, stringsAsFactors = FALSE)

# ── Helper: Safe Extraction (Prevents crashes) ─────────────────────────────
safe_p <- function(res, name) {
  if (is.null(res) || !name %in% names(res)) return(NA)
  return(res[[name]])
}

.count.trans <- function(h) {
  sum(vapply(h$maps, function(m) length(m) - 1L, integer(1)))
}

.sim.dataset <- function(n.tips, q.rate, s, min.state, max.reject = 500L) {
  Q <- matrix(c(-q.rate, q.rate, q.rate, -q.rate), 2, 2)
  n.discard <- 0L
  repeat {
    if (n.discard >= max.reject) return(NULL)
    tree <- pbtree(n = n.tips, scale = 1)
    x <- fastBM(tree)
    anc.x <- fastAnc(tree, x)
    n.edges <- nrow(tree$edge)
    all.vals <- c(x[tree$tip.label], anc.x)
    names(all.vals)[seq_along(tree$tip.label)] <- seq_along(tree$tip.label)
    branch.avg <- numeric(n.edges)
    for (k in seq_len(n.edges)) {
      pa <- tree$edge[k, 1]; ch <- tree$edge[k, 2]
      branch.avg[k] <- (all.vals[as.character(pa)] + all.vals[as.character(ch)]) / 2
    }
    tree.scaled <- tree
    direction <- sample(c(1L, -1L), 1)
    if (s != 1) {
      q25 <- quantile(branch.avg, 0.25); q75 <- quantile(branch.avg, 0.75)
      upper <- branch.avg >= q75; lower <- branch.avg <= q25
      if (direction == 1L) {
        tree.scaled$edge.length[upper] <- tree.scaled$edge.length[upper] * s
        tree.scaled$edge.length[lower] <- tree.scaled$edge.length[lower] / s
      } else {
        tree.scaled$edge.length[lower] <- tree.scaled$edge.length[lower] * s
        tree.scaled$edge.length[upper] <- tree.scaled$edge.length[upper] / s
      }
    }
    h <- sim.history(tree.scaled, Q, nsim = 1)
    y <- as.integer(h$states == "1"); names(y) <- names(h$states)
    if (min(sum(y == 0), sum(y == 1)) >= min.state) {
      return(list(tree=tree, x=x, y=y, n0=sum(y==0), n1=sum(y==1), 
                  n.trans=.count.trans(h), n.discard=n.discard, direction=direction))
    }
    n.discard <- n.discard + 1L
  }
}

run.one <- function(n.tips, q.rate, s, rep.id) {
  dat <- .sim.dataset(n.tips, q.rate, s, min.state, max.reject)
  if (is.null(dat)) {
    return(data.frame(n.tips=n.tips, q.rate=q.rate, s=s, direction=NA, rep=rep.id, n0=NA, n1=NA, 
                      n.true.trans=NA, n.discard=max.reject, p.anccond.01=NA, p.anccond.10=NA, 
                      p.phyloglm=NA, p.pgls=NA, stringsAsFactors=FALSE))
  }
  
  # ── MINIMAL FIX: Format for AncCondFast (1/2 coding and 3-column DF) ──
  sim_data_01 <- data.frame(taxon = names(dat$x), cont = as.numeric(dat$x), disc = as.integer(dat$y + 1))
  
  # Call for 0->1
  res_01 <- tryCatch({
    AncCondFast(tree = dat$tree, data = sim_data_01, n_maps = n.maps, n_iter = n.sims, message = FALSE)
  }, error = function(e) NULL)
  
  # Call for 1->0 (Flip labels so 1 becomes state 1 and 0 becomes state 2)
  y_flipped <- ifelse(dat$y == 1, 1, 2)
  sim_data_10 <- data.frame(taxon = names(dat$x), cont = as.numeric(dat$x), disc = as.integer(y_flipped))
  
  res_10 <- tryCatch({
    AncCondFast(tree = dat$tree, data = sim_data_10, n_maps = n.maps, n_iter = n.sims, message = FALSE)
  }, error = function(e) NULL)
  
  # Standard methods
  p.phyloglm <- tryCatch({
    fit <- phyloglm(y ~ x, data = data.frame(y=dat$y, x=dat$x), phy = dat$tree, method = "logistic_MPLE")
    coef(summary(fit))["x", "p.value"]
  }, error = function(e) NA)
  
  p.pgls <- tryCatch({
    fit <- gls(x ~ as.factor(y), data = data.frame(x=dat$x, y=dat$y), correlation = corBrownian(phy = dat$tree))
    summary(fit)$tTable[2, 4]
  }, error = function(e) NA)
  
  data.frame(n.tips=n.tips, q.rate=q.rate, s=s, direction=dat$direction, rep=rep.id, 
             n0=dat$n0, n1=dat$n1, n.true.trans=dat$n.trans, n.discard=dat$n.discard,
             p.anccond.01=safe_p(res_01, "p_value"), p.anccond.10=safe_p(res_10, "p_value"),
             p.phyloglm=p.phyloglm, p.pgls=p.pgls, stringsAsFactors=FALSE)
}

# ── Parallel Setup ──────────────────────────────────────────────────────────
if (is_windows) {
  cl <- makeCluster(n.cores); registerDoParallel(cl)
  abs_path <- normalizePath("AncCondFast.R")
  clusterExport(cl, varlist = c("run.one", ".sim.dataset", ".count.trans", "safe_p", 
                                "min.state", "max.reject", "n.maps", "n.sims", "alpha", "grid", "abs_path"))
  clusterEvalQ(cl, { library(ape); library(phytools); library(phylolm); library(nlme); source(abs_path) })
}

# ── Main loop ───────────────────────────────────────────────────────────────
cat("\n=== Starting simulations ===\n")
all.raw <- list(); discard.log <- list()

for (g in seq_len(nrow(grid))) {
  nt <- grid$n.tips[g]; qr <- grid$q.rate[g]; ss <- grid$s[g]
  cat(sprintf("\n[%3d/%d] n.tips=%3d  q=%.1f  s=%2d\n", g, nrow(grid), nt, qr, ss))
  
  ids <- 1:n.reps
  batch.out <- if(is_windows) {
    foreach(i = ids) %dopar% { run.one(nt, qr, ss, i) }
  } else {
    mclapply(ids, function(i) run.one(nt, qr, ss, i), mc.cores = n.cores)
  }
  
  raw.df <- do.call(rbind, batch.out)
  all.raw[[g]] <- raw.df
  discard.log[[g]] <- data.frame(n.tips=nt, q.rate=qr, s=ss, total.discard=sum(raw.df$n.discard, na.rm=TRUE), 
                                 mean.true.trans=mean(raw.df$n.true.trans, na.rm=TRUE), stringsAsFactors=FALSE)
}

if (is_windows) stopCluster(cl)

# ── Final Summaries ─────────────────────────────────────────────────────────
raw.all <- do.call(rbind, all.raw)
methods <- c("anccond.01", "anccond.10", "phyloglm", "pgls")
p.cols  <- c("p.anccond.01", "p.anccond.10", "p.phyloglm", "p.pgls")

cat("\n=== Summary (all methods) ===\n")
for (g in seq_len(nrow(grid))) {
  sub <- raw.all[raw.all$n.tips == grid$n.tips[g] & raw.all$q.rate == grid$q.rate[g] & raw.all$s == grid$s[g], ]
  cat(sprintf("\nCond: n=%d q=%.1f s=%d\n", grid$n.tips[g], grid$q.rate[g], grid$s[g]))
  for (m in seq_along(methods)) {
    pv <- sub[[p.cols[m]]]; valid <- sum(!is.na(pv))
    rate <- if(valid > 0) mean(pv[!is.na(pv)] < alpha) else NA
    cat(sprintf("  %-10s: rate=%.3f (n=%d)\n", methods[m], rate, valid))
  }
}