# ------------------------------------------------------------
# sim_tree_data: BM tips -> scale (anc.ML) -> discrete per scaled tree
# - For loops only; minimal changes from your prior script
# - Stores per-sf discrete traits so analysis matches scaled trees
# ------------------------------------------------------------

library(ape)
library(phytools)
library(geiger)
library(TreeSim)

set.seed(42)

# ============== 0) BD trees (reuse if already exist) ==============
if (!exists("bd_trees")) {
  sizes  <- c(25, 50, 75, 100, 200)
  reps   <- 100
  lambda <- 3
  mu     <- 1
  bd_trees <- setNames(vector("list", 
                              length(sizes)), 
                       as.character(sizes))
  for (i in seq_along(sizes)) {
    n <- sizes[i]
    message(sprintf("Simulating %d trees with %d tips...", reps, n))
    bd_trees[[i]] <- TreeSim::sim.bd.taxa(n = n, numbsim = reps, 
                                          lambda = lambda, mu = mu, 
                                          complete = FALSE)
  }
}
rm(list=ls()[-1])

# ============== 1) Scenarios & params ==============
scenario_names <- c("uni", "bi")
# 1->2 first value is for uni second for bi
forward_rates  <- c(0.2, 0.6)
# 2->1 first value is for uni second for bi
reverse_rates  <- c(0.0, 0.6)   
cont_sigma     <- 0.2           # BM variance

# ======= 2) Continuous trait simulation ========
cont.traits <- setNames(vector("list", length(bd_trees)), 
                        names(bd_trees))  
for(i in seq_along(bd_trees)){
  for(j in seq_along(bd_trees[[i]]))
    cont.traits[[i]][[j]] <- phytools::fastBM(bd_trees[[i]][[j]], 
                                              sig2 = cont_sigma, a = 0.0)  
}
# ============== 3) Scale trees based on Cont. Trait ==============
# this will hold our scaled trees.
bd_trees_scaled <- vector("list", length(bd_trees))  # per original tree: sf1..sf10 trees
# here we loop through sizes first with i and then reps with j
for(i in seq_along(bd_trees)){
  cat(paste("Working on trees of size", names(bd_trees)[i]),"\n")
  for(j in seq_along(bd_trees[[i]])){
    tr <- bd_trees[[i]][[j]]
    ct <- cont.traits[[i]][[j]]
    ct_est <- phytools::anc.ML(tr, ct, model = "BM")
    Ntip  <- length(tr$tip.label)
    Nnode <- tr$Nnode
    total_nodes <- Ntip + Nnode
    all_states <- numeric(total_nodes)
    all_states[seq_len(Ntip)] <- as.numeric(ct[tr$tip.label])
    all_states[as.integer(names(ct_est$ace))] <- as.numeric(ct_est$ace)
    # branch means
    branch_means <- rep(0, nrow(tr$edge))
    for(j in seq_along(branch_means)) {
      parent <- tr$edge[j, 1]
      child  <- tr$edge[j, 2]
      branch_means[j] <- (all_states[parent] + all_states[child]) / 2
    }
    q_lo <- as.numeric(quantile(branch_means, probs = 0.25, 
                                type = 7, names = FALSE))
    q_hi <- as.numeric(quantile(branch_means, probs = 0.75, 
                                type = 7, names = FALSE))
    # make 10 scaled versions: sf1..sf10
    scaled_set <- setNames(vector("list", 10), paste0("sf", 1:10))
    for (sf in 1:10) {
      tr_scaled <- tr
      new_el <- tr$edge.length
      new_el[branch_means >= q_hi] <- new_el[branch_means >= q_hi] * sf
      new_el[branch_means <= q_lo] <- new_el[branch_means <= q_lo] / sf
      tr_scaled$edge.length <- new_el
      scaled_set[[sf]] <- tr_scaled
    }
    # TODO this line of code is not storing trees investigate
    # pick up here to get this loop working we want to leave
    # this loop with the below variable having 10 scaled trees for each
    # original tree.
    bd_trees_scaled[i][j] <- scaled_set
  }
}
# BM ML estimation





# outputs
sim_results <- setNames(vector("list", length(scenario_names)), 
                        scenario_names)  # traits (now includes per-sf discrete)
# sf_trees will hold trees that have been scaled
sf_trees    <- setNames(vector("list", length(scenario_names)), 
                        scenario_names)  # scaled trees (sf1..sf10)











# ============== 2) Loop scenarios ==============
for (s in seq_along(scenario_names)) {
  scen <- scenario_names[s]
  
  
  
  
  
  
  
  
  
  
  # build q matrix
  q_matrix <- matrix(c(-forward_rates[s], forward_rates[s],
                       reverse_rates[s], -reverse_rates[s]),
                     nrow = 2, byrow = TRUE)
  dimnames(q_matrix) <- list(c("1","2"), c("1","2"))
  # stationary (for BI root draws)
  ## dep denom <- forward_rates[s] + reverse_rates[s]
  ## dep pi1 <- if (denom > 0) reverse_rate / denom else 1.0
  ## dep pi2 <- if (denom > 0) forward_rate / denom else 0.0
  
  sizes <- as.integer(names(bd_trees))
  per_scen_results <- setNames(vector("list", length(sizes)), 
                               as.character(sizes))
  per_scen_sf      <- setNames(vector("list", length(sizes)), 
                               as.character(sizes))
  
  # ---------- sizes ----------
  for (i in seq_along(sizes)) {
    n <- sizes[i]
    tr_list <- bd_trees[[as.character(n)]]
    trait_list <- vector("list", length(tr_list))  # per original tree
    
    message(sprintf("[%s] Preparing %d trees (n=%d)...", scen, length(tr_list), n))
    
    # ---------- 100 trees ----------
    for (ti in seq_along(tr_list)) {
      tr <- tr_list[[ti]]
      
      
 
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
          # TODO fix to reflect manuscript
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