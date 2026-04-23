# test_extreme_cases.R
# Smoke-test AncCond under two extremes:
#   1. Strong signal  -- discrete state perfectly tracks high/low cont. trait
#   2. Zero signal    -- discrete state is random w.r.t. the cont. trait
#
# Expected outcomes:
#   Strong signal  -> at least one p-value < 0.05
#   Zero signal    -> both p-values > 0.05 (most of the time)

library(ape)
library(phytools)

source("../R/anc_cond.R")

set.seed(123)

# ---- Shared tree --------------------------------------------------------
tree <- rtree(80)

# ---- Helper: build the 3-column data.frame AncCond expects -------------
make_data <- function(tree, cont, disc) {
  data.frame(
    species  = tree$tip.label,
    cont     = cont[tree$tip.label],
    disc     = disc[tree$tip.label],
    stringsAsFactors = FALSE
  )
}

# ========================================================================
# Case 1 — Strong signal
#
# Simulate a continuous trait (BM), then assign discrete state = 1 when
# the tip value is below the median and state = 2 when above.  This
# creates a strong deterministic link between trait value and state,
# so transitions from 1->2 should happen at high values and 2->1 at low
# values — AncCond should detect this.
# ========================================================================
cat("\n===== Case 1: Strong signal =====\n")

cont_strong <- fastBM(tree, sig2 = 1, a = 0)
med <- median(cont_strong)
disc_strong <- ifelse(cont_strong >= med, 2, 1)
names(disc_strong) <- names(cont_strong)

dat_strong <- make_data(tree, cont_strong, disc_strong)

res_strong <- AncCond(
  tree    = tree,
  data    = dat_strong,
  nsim    = 20,
  iter    = 100,
  message = FALSE
)
summary(res_strong)

strong_pass <- any(res_strong$pvals < 0.05, na.rm = TRUE)
cat("Strong-signal test:",
    if (strong_pass) "PASS (p < 0.05 detected)\n"
    else "WARN — no p < 0.05; may be a stochastic miss\n")

# ========================================================================
# Case 2 — Zero signal
#
# Same BM continuous trait, but the discrete state is drawn by simulating
# a Mk process on the *original* (unscaled) tree — completely independent
# of the continuous trait.  AncCond should find no significant departure
# from the null.
# ========================================================================
cat("\n===== Case 2: Zero signal =====\n")

cont_null <- fastBM(tree, sig2 = 1, a = 0)

# Simulate a discrete trait that is fully independent of cont_null.
Q_null <- matrix(c(-0.5, 0.5, 0.5, -0.5), 2, 2,
                 dimnames = list(c("1","2"), c("1","2")))
ok <- FALSE
attempts <- 0
while (!ok) {
  attempts <- attempts + 1
  sim_disc <- sim.history(tree, Q = Q_null, nsim = 1, message = FALSE)
  disc_null <- as.numeric(phytools::getStates(sim_disc, "tips"))
  names(disc_null) <- tree$tip.label
  tab <- table(disc_null)
  ok <- length(tab) == 2 && min(tab) >= 10
  if (attempts > 500) stop("Could not get a balanced discrete trait")
}

dat_null <- make_data(tree, cont_null, disc_null)

res_null <- AncCond(
  tree    = tree,
  data    = dat_null,
  nsim    = 20,
  iter    = 100,
  message = FALSE
)
summary(res_null)

null_pass <- all(is.na(res_null$pvals) | res_null$pvals > 0.05)
cat("Zero-signal test:",
    if (null_pass) "PASS (no p < 0.05)\n"
    else "WARN — false positive detected; can happen ~5% of the time\n")

# ---- Final verdict ------------------------------------------------------
cat("\n===== Summary =====\n")
cat("Strong signal p-values:", round(res_strong$pvals, 4), "\n")
cat("Zero signal   p-values:", round(res_null$pvals, 4), "\n")
if (strong_pass && null_pass) {
  cat("Both tests PASSED.\n")
} else {
  cat("One or both tests returned a warning — see details above.\n")
}
