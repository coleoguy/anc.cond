# ------------------------------------------------------------
# Birth–death (BD) tree simulation with FOR LOOPS (no helpers)
# Creates `bd_trees` in the Global Environment
# ------------------------------------------------------------

# install.packages(c("TreeSim","ape"))  # uncomment if needed
library(geiger)
library(TreeSim)
library(ape)

set.seed(123)

sizes  <- c(25, 50, 75, 100, 200)  # target tip counts
reps   <- 100                      # trees per size
lambda <- 1.0                      # birth/speciation rate
mu     <- 0.2                      # death/extinction rate

bd_trees <- setNames(vector("list", length(sizes)), as.character(sizes))

for (i in seq_along(sizes)) {
  n <- sizes[i]
  message(sprintf("Simulating %d trees with %d tips...", reps, n))
  
  trees_for_n <- vector("list", reps)  # preallocate
  k <- 0                               # how many accepted trees we have
  
  # Keep simulating until we have `reps` non-NULL trees
  while (k < reps) {
    batch <- sim.bd.taxa(n = n, numbsim = reps, lambda = lambda, mu = mu, complete = TRUE)
    
    # Loop through batch; accept non-NULL results
    for (j in seq_along(batch)) {
      tr <- batch[[j]]
      if (!is.null(tr)) {
        k <- k + 1
        trees_for_n[[k]] <- tr
        if (k == reps) break
      }
    }
  }
  
  bd_trees[[as.character(n)]] <- trees_for_n
}

# Sanity checks (optional):
#names(bd_trees)
#length(bd_trees[["25"]]); class(bd_trees[["25"]][[1]])

# ------------------------------------------------------------
# Simulate traits on BD trees using forward/reverse 2-state Mk
# - Discrete: Mk with forward_rate (1->2) and reverse_rate (2->1)
# - Continuous: Brownian motion with variance cont_sigma
# Produces: sim_results[["uni"]], sim_results[["bi"]]
#           sf_trees   [["uni"]], sf_trees   [["bi"]]
# ------------------------------------------------------------

library(ape)
library(phytools)

if (!exists("bd_trees")) stop("bd_trees not found in Global Environment.")
sizes <- as.integer(names(bd_trees))

# Scenarios
scenario_names <- c("uni", "bi")
forward_rates  <- c(0.2, 0.6)  # 1->2
reverse_rates  <- c(0.0, 0.6)  # 2->1

cont_sigma <- 0.2
verbose    <- TRUE
set.seed(42)

# Outputs
sim_results <- setNames(vector("list", length(scenario_names)), scenario_names)
sf_trees    <- setNames(vector("list", length(scenario_names)), scenario_names)

# ---- Loop over scenarios ----
for (s in seq_along(scenario_names)) {
  scen <- scenario_names[s]
  forward_rate <- forward_rates[s]
  reverse_rate <- reverse_rates[s]
  
  # Build matrices with different names per scenario
  if (scen == "uni") {
    q_matrix <- matrix(c(-forward_rate, forward_rate,
                         reverse_rate, -reverse_rate),
                       nrow = 2, byrow = TRUE)
    dimnames(q_matrix) <- list(c("1","2"), c("1","2"))
  } else {
    rate_matrix <- matrix(c(-forward_rate, forward_rate,
                            reverse_rate, -reverse_rate),
                          nrow = 2, byrow = TRUE)
    dimnames(rate_matrix) <- list(c("1","2"), c("1","2"))
  }
  
  # Stationary distribution (for BI root draws only)
  denom <- forward_rate + reverse_rate
  pi1 <- if (denom > 0) reverse_rate / denom else 1.0
  pi2 <- if (denom > 0) forward_rate / denom else 0.0
  
  # ---- Trait simulations ----
  per_scen_results <- setNames(vector("list", length(sizes)), as.character(sizes))
  
  for (i in seq_along(sizes)) {
    n <- sizes[i]
    trees_for_n <- bd_trees[[as.character(n)]]
    out_list    <- vector("list", length(trees_for_n))
    
    if (verbose) message(sprintf("[%s] Simulating traits for %d trees (n=%d)...", scen, length(trees_for_n), n))
    
    for (j in seq_along(trees_for_n)) {
      tr <- trees_for_n[[j]]
      
      # Continuous (BM)
      x_cont <- fastBM(tr, sig2 = cont_sigma, a = 0.0)
      
      # Root handling
      if (scen == "uni") {
        root_val <- "1"
      } else {
        root_val <- sample(c("1","2"), size = 1, prob = c(pi1, pi2))
      }
      
      # Discrete (Mk), using the scenario-specific matrix variable name
      if (scen == "uni") {
        x_disc <- rTraitDisc(phy = tr, model = q_matrix,   states = colnames(q_matrix),   root.value = root_val)
      } else {
        x_disc <- rTraitDisc(phy = tr, model = rate_matrix, states = colnames(rate_matrix), root.value = root_val)
      }
      x_disc <- setNames(as.character(x_disc), names(x_disc))
      
      # Store; keep matrix under distinct field names
      if (scen == "uni") {
        out_list[[j]] <- list(tree = tr, cont_trait = x_cont, disc_trait = x_disc,
                              root_state = root_val, Q = q_matrix, sigma2 = cont_sigma)
      } else {
        out_list[[j]] <- list(tree = tr, cont_trait = x_cont, disc_trait = x_disc,
                              root_state = root_val, rate_matrix = rate_matrix, sigma2 = cont_sigma)
      }
    }
    
    per_scen_results[[as.character(n)]] <- out_list
  }
  
  sim_results[[scen]] <- per_scen_results
  
  # ---- sf-scaled trees (sf = 1..10) using cont_trait we just simulated ----
  per_scen_sf <- setNames(vector("list", length(sizes)), as.character(sizes))
  
  for (si in seq_along(sizes)) {
    n <- sizes[si]
    tr_list <- bd_trees[[as.character(n)]]
    per_size_out <- vector("list", length(tr_list))
    
    if (verbose) message(sprintf("[%s] Scaling %d trees for n=%d ...", scen, length(tr_list), n))
    
    for (ti in seq_along(tr_list)) {
      tr <- tr_list[[ti]]
      
      x_cont <- sim_results[[scen]][[as.character(n)]][[ti]]$cont_trait
      if (is.null(x_cont)) stop(sprintf("[%s] Missing cont_trait for n=%d, tree %d.", scen, n, ti))
      
      anc_states <- fastAnc(tr, x_cont, vars = FALSE, CI = FALSE)
      
      Ntip  <- length(tr$tip.label)
      Nnode <- tr$Nnode
      total_nodes <- Ntip + Nnode
      all_states <- numeric(total_nodes)
      
      tip_vals <- x_cont[tr$tip.label]
      all_states[seq_len(Ntip)] <- as.numeric(tip_vals)
      
      int_ids <- as.integer(names(anc_states))
      all_states[int_ids] <- as.numeric(anc_states)
      
      E <- nrow(tr$edge)
      branch_means <- numeric(E)
      for (e in seq_len(E)) {
        parent <- tr$edge[e, 1]
        child  <- tr$edge[e, 2]
        branch_means[e] <- mean(c(all_states[parent], all_states[child]))
      }
      
      q_lo <- as.numeric(quantile(branch_means, probs = 0.25, type = 7, names = FALSE))
      q_hi <- as.numeric(quantile(branch_means, probs = 0.75, type = 7, names = FALSE))
      
      scaled_set <- setNames(vector("list", 10), paste0("sf", 1:10))
      for (sf in 1:10) {
        tr_scaled <- tr
        new_el <- tr$edge.length
        for (e in seq_len(E)) {
          bm <- branch_means[e]
          if (bm >= q_hi) {
            new_el[e] <- new_el[e] * sf
          } else if (bm <= q_lo) {
            new_el[e] <- new_el[e] / sf
          } else {
            new_el[e] <- new_el[e]
          }
        }
        tr_scaled$edge.length <- new_el
        scaled_set[[sf]] <- tr_scaled
      }
      
      per_size_out[[ti]] <- scaled_set
    }
    
    per_scen_sf[[as.character(n)]] <- per_size_out
  }
  
  sf_trees[[scen]] <- per_scen_sf
}

# Results:
# sim_results[["uni"]], sim_results[["bi"]]
# sf_trees   [["uni"]], sf_trees   [["bi"]]
# Example peeks:
# names(sim_results); names(sf_trees)
# table(sim_results[["bi"]][["25"]][[1]]$disc_trait)
# class(sf_trees[["bi"]][["100"]][[10]][["sf3"]])
