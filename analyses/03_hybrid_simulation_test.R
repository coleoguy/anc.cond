library(ape)
library(phytools)

# ── 1. HYBRID ENGINE (SAME AS BEFORE) ────────────────────────────────────
.sim.hybrid.mapped <- function(tree, maps, theta_1, alpha, sigma, root_state) {
  n.tips <- length(tree$tip.label); n.total <- n.tips + tree$Nnode
  root <- n.tips + 1L; edges <- tree$edge; vals <- numeric(n.total)
  
  # Start at 0 for BM or theta_1 for OU
  vals[root] <- if(root_state == "1") theta_1 else 0
  
  for (k in seq_len(nrow(edges))) {
    pa <- edges[k, 1]; ch <- edges[k, 2]; x <- vals[pa]
    seg <- maps[[k]]; states <- names(seg)
    for (j in seq_along(seg)) {
      dt <- seg[j]; st <- states[j]
      if (st == "1") {
        e.at <- exp(-alpha * dt)
        mu <- theta_1 + (x - theta_1) * e.at
        v <- (sigma^2 / (2 * alpha)) * (1 - exp(-2 * alpha * dt))
        x <- rnorm(1, mean = mu, sd = sqrt(max(v, 0)))
      } else {
        x <- rnorm(1, mean = x, sd = sqrt(sigma^2 * dt))
      }
    }
    vals[ch] <- x
  }
  out <- vals[1:n.tips]; names(out) <- tree$tip.label
  return(out)
}

# ── 2. CALIBRATION RUNS ──────────────────────────────────────────────────
# Setup two test cases from your grid
test_cases <- list(
  weak   = list(alpha = 0.5, delta_theta = 2.0, title = "Weak: alpha=0.5, dTheta=2"),
  strong = list(alpha = 8.0, delta_theta = 8.0, title = "Strong: alpha=8, dTheta=8")
)

par(mfrow = c(2, 2)) # 2x2 grid for Phenograms and Boxplots
cols <- setNames(c("dodgerblue", "firebrick"), c("0", "1"))

for (case in test_cases) {
  # Generate Tree & History (using q.rate from block 2)
  tree <- pbtree(n = 200, scale = 1)
  Q <- matrix(c(-0.5, 0.5, 0.5, -0.5), 2, 2, dimnames = list(c("0","1"), c("0","1")))
  h <- sim.history(tree, Q, nsim = 1, message = FALSE, pi = c(0.5, 0.5))
  y <- as.integer(h$states == "1")
  
  # Simulate Trait
  root_st <- names(h$maps[[which(tree$edge[,1] == 201)[1]]])[1]
  x <- .sim.hybrid.mapped(reorder(tree, "cladewise"), h$maps, 
                          theta_1 = case$delta_theta, 
                          alpha = case$alpha, 
                          sigma = 1.0, 
                          root_state = root_st)
  
  # A. Plot Phenogram
  phenogram(h, x, colors = cols, main = case$title, fsize = 0.5)
  abline(h = case$delta_theta, col = "firebrick", lty = 2)
  
  # B. Plot Boxplot
  boxplot(x ~ y, col = cols, main = paste("Boxplot:", case$title), 
          names = c("BM", "OU"), ylab = "Value")
  abline(h = case$delta_theta, col = "firebrick", lty = 2)
}