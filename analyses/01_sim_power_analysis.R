# ============================================================================
# Main Power & False-Positive Simulation for AncCond (Fixed for Windows)
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
n.tips.vec    <- c(50, 100, 200, 500)
q.rate.vec    <- c(0.1, 0.5, 1.0, 2.0)
s.vec         <- c(1, 2, 3, 5, 10)       # s = 1 is null
n.reps        <- 200
n.maps        <- 100
n.sims        <- 500
alpha_level   <- 0.05
min.state.pct <- 0.10  
max.reject    <- 500L  

grid <- expand.grid(n.tips = n.tips.vec, q.rate = q.rate.vec, s = s.vec, stringsAsFactors = FALSE)

# ── Helper functions ──────────────────────────────────────────────────────
.count.trans <- function(h) {
  sum(vapply(h$maps, function(m) length(m) - 1L, integer(1)))
}

.sim.dataset <- function(n.tips, q.rate, s, min.state.pct, max.reject = 500L) {
  Q <- matrix(c(-q.rate, q.rate, q.rate, -q.rate), 2, 2)
  n.discard <- 0L
  threshold <- ceiling(n.tips * min.state.pct)
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
      if (direction == 1L) {
        tree.scaled$edge.length[branch.avg >= q75] <- tree.scaled$edge.length[branch.avg >= q75] * s
        tree.scaled$edge.length[branch.avg <= q25] <- tree.scaled$edge.length[branch.avg <= q25] / s
      } else {
        tree.scaled$edge.length[branch.avg <= q25] <- tree.scaled$edge.length[branch.avg <= q25] * s
        tree.scaled$edge.length[branch.avg >= q75] <- tree.scaled$edge.length[branch.avg >= q75] / s
      }
    }
    h <- sim.history(tree.scaled, Q, nsim = 1, message = FALSE)
    y <- as.integer(h$states == "1"); names(y) <- names(h$states)
    if (min(sum(y == 0), sum(y == 1)) >= threshold) {
      return(list(tree=tree, x=x, y=y, n.trans=.count.trans(h), n.discard=n.discard, Q=Q))
    }
    n.discard <- n.discard + 1L
  }
}

run.one <- function(n.tips, q.rate, s, rep.id, min.state.pct, n.maps, n.sims, max.reject) {
  dat <- .sim.dataset(n.tips, q.rate, s, min.state.pct, max.reject)
  if (is.null(dat)) return(data.frame(n.tips=n.tips, q.rate=q.rate, s=s, rep=rep.id, n.discard=max.reject, stringsAsFactors=F))
  res.ac <- tryCatch(suppressWarnings(AncCond(dat$tree, dat$x, dat$y, n.maps=n.maps, n.sims=n.sims, Q=dat$Q)),
                     error = function(e) list(p.01=NA, p.10=NA))
  p.plg <- tryCatch({ fit <- phyloglm(y ~ x, data=data.frame(y=dat$y, x=dat$x), phy=dat$tree, method="logistic_MPLE"); coef(summary(fit))["x", "p.value"] }, error = function(e) NA)
  p.pgls <- tryCatch({ df <- data.frame(x=dat$x, y=as.factor(dat$y), row.names=dat$tree$tip.label); fit <- gls(x ~ y, data=df, correlation=corBrownian(phy=dat$tree)); summary(fit)$tTable["y1", "p-value"] }, error = function(e) NA)
  data.frame(n.tips=n.tips, q.rate=q.rate, s=s, rep=rep.id, p.anccond.01=res.ac$p.01, p.anccond.10=res.ac$p.10, p.phyloglm=p.plg, p.pgls=p.pgls, n.true.trans=dat$n.trans, n.discard=dat$n.discard, stringsAsFactors=F)
}

# ── SETUP CLUSTER ──────────────────────────────────────────────────────────
cl <- makeCluster(n.cores)
clusterEvalQ(cl, { library(ape); library(phytools); library(phylolm); library(nlme); source("AncCond.R") })
clusterExport(cl, c(".sim.dataset", ".count.trans", "run.one", "min.state.pct", "n.maps", "n.sims", "max.reject"))

# ── Main loop ──────────────────────────────────────────────────────────────
cat("\n=== Starting simulations ===\n")
all.raw <- vector("list", nrow(grid))
discard.log <- vector("list", nrow(grid))
batch.size <- 25

t.total <- system.time({
  for (g in seq_len(nrow(grid))) {
    nt <- grid$n.tips[g]; qr <- grid$q.rate[g]; ss <- grid$s[g]
    cat(sprintf("\n[%3d/%d] n=%d q=%.1f s=%d... ", g, nrow(grid), nt, qr, ss))
    
    done <- 0L; cond.results <- list()
    while (done < n.reps) {
      this.batch <- min(batch.size, n.reps - done)
      ids <- (done + 1L):(done + this.batch)
      
      # Pass parameters explicitly to avoid function naming conflicts
      batch.out <- parLapply(cl, ids, function(i, nt_in, qr_in, ss_in) {
        run.one(nt_in, qr_in, ss_in, i, min.state.pct, n.maps, n.sims, max.reject)
      }, nt_in = nt, qr_in = qr, ss_in = ss)
      
      cond.results <- c(cond.results, batch.out)
      done <- done + this.batch
    }
    raw.cond <- do.call(rbind, cond.results)
    all.raw[[g]] <- raw.cond
    
    # Log discards for this condition
    discard.log[[g]] <- data.frame(n.tips=nt, q.rate=qr, s=ss, total.discard=sum(raw.cond$n.discard, na.rm=T))
    cat("Done.\n")
  }
})

stopCluster(cl)

library(dplyr)

# ── Combine results ───────────────────────────────────────────────────────
raw.all <- do.call(rbind, all.raw) 
discard.all <- do.call(rbind, discard.log)

summary.df <- aggregate(
  cbind(p.anccond.01, p.anccond.10, p.phyloglm, p.pgls) ~ n.tips + q.rate + s, 
  data = raw.all, 
  FUN = function(x) mean(x < alpha_level, na.rm = TRUE),
  na.action = na.pass
)

summary.df$metric <- ifelse(summary.df$s == 1, "FPR", "power")

write.csv(raw.all,     "01_sim_power_raw.csv",     row.names = FALSE)
write.csv(summary.df,  "01_sim_power_results.csv", row.names = FALSE)
write.csv(discard.all, "01_sim_discard_log.csv",   row.names = FALSE)

cat("\nAll files saved: 01_sim_power_raw.csv, 01_sim_power_results.csv, 01_sim_discard_log.csv\n")