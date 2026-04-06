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

########## 0) BD trees (reuse if already exist) ###########
if (!exists("bd_trees")) {
  sizes  <- c(25, 50, 75, 100, 200)
  reps   <- 1000
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
# leave here with bd_trees (orig tree no scaling)
########## 0) BD trees (reuse if already exist) ###########

########## 1) Scenarios & params ###########
scenario_names <- c("uni", "bi")
# 1->2 first value is for uni second for bi
forward_rates  <- c(0.2, 0.6)
# 2->1 first value is for uni second for bi
reverse_rates  <- c(0.0, 0.6)   
cont_sigma     <- 0.2           # BM variance
########## 1) Scenarios & params ###########

######## 2) Continuous trait simulation ##########
cont.traits <- setNames(vector("list", length(bd_trees)), 
                        names(bd_trees))  
for(i in seq_along(bd_trees)){
  for(j in seq_along(bd_trees[[i]]))
    cont.traits[[i]][[j]] <- phytools::fastBM(bd_trees[[i]][[j]], 
                                              sig2 = cont_sigma, a = 0.0)  
}
rm(list=ls()[c(-1:-4,-7,-8)])
# leave with cont.traits
######## 2) End Continuous trait simulation ##########


########### 3) Scale trees based on Cont. Trait ########
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
    all_states <- Ntip + Nnode
    all_states[seq_len(Ntip)] <- as.numeric(ct[tr$tip.label])
    all_states[as.integer(names(ct_est$ace))] <- as.numeric(ct_est$ace)
    # branch means
    branch_means <- rep(0, nrow(tr$edge))
    for(k in seq_along(branch_means)) {
      parent <- tr$edge[k, 1]
      child  <- tr$edge[k, 2]
      branch_means[k] <- (all_states[parent] + all_states[child]) / 2
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
    bd_trees_scaled[[i]][[j]] <- scaled_set
  }
}
rm(list=ls()[c(1, 4:6,8,9,11:19,21,23:25)])

########### 3) End of Scale trees based on Cont. Trait ########

##### Discrete Trait Simulation ##########
# output objects setup for results
sim_results <- setNames(vector("list", length(scenario_names)), 
                        scenario_names)  # traits (now includes per-sf discrete)
# sf_trees will hold trees that have been scaled
sf_trees    <- setNames(vector("list", length(scenario_names)), 
                        scenario_names)  # scaled trees (sf1..sf10)
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
  # i will loop through all 5 sizes of trees
  for (i in seq_along(sizes)) {
    n <- sizes[i]
    tr_list <- bd_trees[[as.character(n)]]
    trait_list <- vector("list", length(tr_list))  # per original tree
    message(sprintf("[%s] Preparing %d trees (n=%d)...", scen, length(tr_list), n))
    # ti goes through all 100 trees of a given size
    # ---------- 100 trees ----------
    for (ti in seq_along(tr_list)) {
      tr <- tr_list[[ti]]
      # (C) Discrete tips **on each SCALED tree** (keeps pipeline consistent)
      disc_by_sf <- setNames(vector("list", length(bd_trees_scaled[[1]][[1]])), 
                             names(bd_trees_scaled[[1]][[1]]))
      for (sf_idx in seq_along(disc_by_sf)) {
        sf_name <- names(disc_by_sf)[sf_idx]
        tr_scaled <- bd_trees_scaled[[i]][[ti]][[sf_idx]]
        
        # root per scenario
        if (scen == "uni") {
          root_state_num <- 1
        } else {
          root_state_num <- sample(c(1, 2), size = 1)
        }
        # ensure non-degenerate discrete tips per scaled tree
        max_attempts <- 1000
        attempts <- 0
        repeat {
          attempts <- attempts + 1
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
            ok <- min(tab) >=10
          }
          if (ok || attempts >= max_attempts) break
        }
        disc_by_sf[[sf_idx]] <- setNames(as.character(disc_vec_num), tr_scaled$tip.label)
      }
      # (D) Store per-tree record
      if (scen == "uni") {
        # Keep backward-compat field `disc_trait` as sf1's vector
        trait_list[[ti]] <- list(
          tree         = tr,
          cont_trait   = cont.traits[[i]][[ti]],
          disc_by_sf   = disc_by_sf            # per-sf discrete vectors ("sf1".. "sf10")
        )
      }
      if(scen == "bi"){
        trait_list[[ti]] <- list(
          tree         = tr,
          cont_trait   = cont.traits[[i]][[ti]],
          disc_by_sf   = disc_by_sf            
        )
      }
    }
    per_scen_results[[as.character(n)]] <- trait_list
  }
  sim_results[[scen]] <- per_scen_results
}
##### End Discrete Trait Simulation ##########
# make name clearer
sim_data <- sim_results
save(sim_data,
     file = "../data/sim_data.RData")
message("Saved simulated data to /data/sim_data.RData")


