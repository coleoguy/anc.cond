# ------------------------------------------------------------
# sim_tree_data: BM tips -> scale (anc.ML) -> discrete per scaled tree
# - For loops only; minimal changes from your prior script
# - Stores per-sf discrete traits so analysis matches scaled trees
# ------------------------------------------------------------

library(ape)
library(phytools)
library(geiger)

set.seed(42)

# ============== 0) BD trees (reuse if already exist) ==============
if (!exists("bd_trees")) {
  sizes  <- c(25, 50, 75, 100, 200)
  reps   <- 100
  lambda <- 1.0
  mu     <- 0.2
  bd_trees <- setNames(vector("list", length(sizes)), as.character(sizes))
  
  for (i in seq_along(sizes)) {
    n <- sizes[i]
    message(sprintf("Simulating %d trees with %d tips...", reps, n))
    trees_for_n <- vector("list", reps)
    k <- 0
    while (k < reps) {
      batch <- TreeSim::sim.bd.taxa(n = n, numbsim = reps, lambda = lambda, mu = mu, complete = TRUE)
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
}
sizes <- as.integer(names(bd_trees))

# ============== 1) Scenarios & params ==============
scenario_names <- c("uni", "bi")
forward_rates  <- c(0.2, 0.6)   # 1->2
reverse_rates  <- c(0.0, 0.6)   # 2->1
cont_sigma     <- 0.2           # BM variance

# outputs
sim_results <- setNames(vector("list", length(scenario_names)), scenario_names)  # traits (now includes per-sf discrete)
sf_trees    <- setNames(vector("list", length(scenario_names)), scenario_names)  # scaled trees (sf1..sf10)

# ============== 2) Loop scenarios ==============
for (s in seq_along(scenario_names)) {
  scen <- scenario_names[s]
  forward_rate <- forward_rates[s]
  reverse_rate <- reverse_rates[s]
  
  # generator
  q_matrix <- matrix(c(-forward_rate, forward_rate,
                       reverse_rate, -reverse_rate),
                     nrow = 2, byrow = TRUE)
  dimnames(q_matrix) <- list(c("1","2"), c("1","2"))
  
  # stationary (for BI root draws)
  denom <- forward_rate + reverse_rate
  pi1 <- if (denom > 0) reverse_rate / denom else 1.0
  pi2 <- if (denom > 0) forward_rate / denom else 0.0
  
  per_scen_results <- setNames(vector("list", length(sizes)), as.character(sizes))
  per_scen_sf      <- setNames(vector("list", length(sizes)), as.character(sizes))
  
  # ---------- sizes ----------
  for (i in seq_along(sizes)) {
    n <- sizes[i]
    tr_list <- bd_trees[[as.character(n)]]
    
    trait_list <- vector("list", length(tr_list))  # per original tree
    sf_set_all <- vector("list", length(tr_list))  # per original tree: sf1..sf10 trees
    
    message(sprintf("[%s] Preparing %d trees (n=%d)...", scen, length(tr_list), n))
    
    # ---------- 100 trees ----------
    for (ti in seq_along(tr_list)) {
      tr <- tr_list[[ti]]
      
      # (A) Continuous tips (BM) on the ORIGINAL tree
      x_cont <- phytools::fastBM(tr, sig2 = cont_sigma, a = 0.0)  # named by tip
      
      # (B) Scale branches using anc.ML (BM ML reconstruction)
      fit_ml <- phytools::anc.ML(tr, x_cont, model = "BM")
      
      # assemble node values (tips + internals)
      Ntip  <- length(tr$tip.label)
      Nnode <- tr$Nnode
      total_nodes <- Ntip + Nnode
      all_states <- numeric(total_nodes)
      
      all_states[seq_len(Ntip)] <- as.numeric(x_cont[tr$tip.label])
      all_states[as.integer(names(fit_ml$ace))] <- as.numeric(fit_ml$ace)
      
      # branch means
      E <- nrow(tr$edge)
      branch_means <- numeric(E)
      for (e in seq_len(E)) {
        parent <- tr$edge[e, 1]
        child  <- tr$edge[e, 2]
        branch_means[e] <- (all_states[parent] + all_states[child]) / 2
      }
      q_lo <- as.numeric(quantile(branch_means, probs = 0.25, type = 7, names = FALSE))
      q_hi <- as.numeric(quantile(branch_means, probs = 0.75, type = 7, names = FALSE))
      
      # make 10 scaled versions: sf1..sf10
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
      sf_set_all[[ti]] <- scaled_set
      
      # (C) Discrete tips **on each SCALED tree** (keeps pipeline consistent)
      disc_by_sf <- setNames(vector("list", length(scaled_set)), names(scaled_set))
      root_by_sf <- setNames(vector("character", length(scaled_set)), names(scaled_set))
      
      for (sf_idx in seq_along(scaled_set)) {
        sf_name <- names(scaled_set)[sf_idx]
        tr_scaled <- scaled_set[[sf_idx]]
        
        # root per scenario
        if (scen == "uni") {
          root_state_num <- 1L
        } else {
          root_state_num <- sample(c(1L, 2L), size = 1, prob = c(pi1, pi2))
        }
        
        # ensure non-degenerate discrete tips per scaled tree
        max_attempts <- 1000L
        attempts <- 0L
        repeat {
          attempts <- attempts + 1L
          sim_out <- geiger::sim.char(tr_scaled,
                                      par   = q_matrix,
                                      model = "discrete",
                                      root  = root_state_num,
                                      nsim  = 1)
          disc_vec_num <- as.integer(sim_out[,,1])
          names(disc_vec_num) <- tr_scaled$tip.label
          
          # both states present & not near-fixed (5–95%)
          tab <- table(disc_vec_num)
          ok <- length(tab) == 2L && all(tab > 0L)
          if (ok) {
            p_min <- min(tab) / sum(tab)
            ok <- p_min > 0.05 && p_min < 0.95
          }
          if (ok || attempts >= max_attempts) break
        }
        
        disc_by_sf[[sf_idx]] <- setNames(as.character(disc_vec_num), tr_scaled$tip.label)
        root_by_sf[sf_idx]   <- as.character(root_state_num)
      }
      
      # (D) Store per-tree record
      if (scen == "uni") {
        # Keep backward-compat field `disc_trait` as sf1's vector
        trait_list[[ti]] <- list(
          tree         = tr,
          cont_trait   = x_cont,
          disc_trait   = disc_by_sf[["sf1"]],   # for false positive; analysis should use disc_by_sf
          disc_by_sf   = disc_by_sf,            # per-sf discrete vectors ("sf1".. "sf10")
          root_state   = "1",                   # uni root is always 1; per-sf kept in root_by_sf too
          root_by_sf   = root_by_sf,            # per-sf root state used for discrete sim
          Q            = q_matrix,
          sigma2       = cont_sigma
        )
      } else {
        trait_list[[ti]] <- list(
          tree         = tr,
          cont_trait   = x_cont,
          disc_trait   = disc_by_sf[["sf1"]],   # keep for false positive
          disc_by_sf   = disc_by_sf,            
          root_state   = NA_character_,         # varies per sf; see root_by_sf
          root_by_sf   = root_by_sf,            
          rate_matrix  = q_matrix,
          sigma2       = cont_sigma
        )
      }
    }
    
    per_scen_results[[as.character(n)]] <- trait_list
    per_scen_sf     [[as.character(n)]] <- sf_set_all
  }
  
  sim_results[[scen]] <- per_scen_results
  sf_trees   [[scen]] <- per_scen_sf
}

# ============== 3) Combine and save ==============
sim_tree_data <- list(
  bd_trees    = bd_trees,
  sim_results = sim_results,  # cont via fastBM; disc simulated per scaled tree (disc_by_sf)
  sf_trees    = sf_trees      # scaled trees (sf1..sf10)
)

save(sim_tree_data,
     file = "C:/Users/mcconnell.m.meghann/Documents/GitHub/anc.cond/results/sim_tree_data.RData")

message("Saved sim_tree_data to C:/Users/mcconnell.m.meghann/Documents/GitHub/anc.cond/results/sim_tree_data.RData")