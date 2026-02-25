# ------------------------------------------------------------
# Run AncCond on selected scaling factors using sim_data,RData
# ------------------------------------------------------------

# Packages
library(ape)
library(phytools)
library(foreach)
library(doParallel)

num_cores <- detectCores() - 10
# 2. Setup the cluster (This works on Mac, Linux, and Windows)
# This creates a set of copies of R running in the background
my_cluster <- makeCluster(num_cores)
# 3. Register the cluster so 'foreach' knows to use it
registerDoParallel(my_cluster)

print(paste("Cluster registered with", getDoParWorkers(), "cores"))


# ---- Path to your saved data ----
# ---- Use the working directory for ALL I/O ----

# ---- AncCond settings ----
iter    <- 1000      # mc reps inside AncCond
verbose <- TRUE


# We’ll store results per scenario, per size, per tree, per sf
result.list <- replicate(
  2,  # level 1
  replicate(
    5,  # level 2
    replicate(
      1000,  # level 3
      as.list(rep(NA, 10)),  # level 4: list of 10 elements
      simplify = FALSE
    ),
    simplify = FALSE
  ),
  simplify = FALSE
)
load("../data/sim_data_UL.RData")  # loads 'sim_data'
names(result.list) <- names(sim_data)
names(sim_data[[1]])
for(i in 1:2){
  names(result.list[[i]]) <- names(sim_data[[1]])
  for(j in 1:5){
    names(result.list[[i]][[j]]) <- 1:1000
    for(k in 1:1000){
      names(result.list[[i]][[j]][[k]]) <- names(sim_data[[i]][[j]][[k]]$disc_by_sf)
    }
  }
}
scenario_names <- names(sim_data)

# ============================================================
# Running AncCond
# ============================================================
for(scen in seq_along(sim_data)){
  for (si in seq_along(sim_data[[scen]])) {
    size_name <- names(sim_data[[scen]])[si]
    if (verbose) message(sprintf("working on %s %s tip trees",  
                                 scenario_names[scen], size_name))
    parallel_results <- foreach(reps = seq_along(sim_data[[scen]][[si]])) %dopar% {
                                  library(ape)
                                  library(phytools)
                                  source("AncCond.R", keep.source = TRUE)
                                  load("../data/sim_data_UL.RData")  # loads 'sim_data'
                                  temp.res <- list()
                                    cur_dat <- sim_data[[scen]][[si]][[reps]]           # traits per original tree
                                    for(sf.strength in seq_along(sim_data[[scen]][[si]][[reps]]$disc_by_sf)){
                                      df <- data.frame(names(cur_dat$cont_trait),
                                                       cur_dat$cont_trait,
                                                       cur_dat$disc_by_sf[[sf.strength]])
                                      if(scen == 1){
                                        res <- AncCond(tree = cur_dat$tree,
                                                       data = df,
                                                       mat = c(0,0,1,0),
                                                       mc = iter)
                                      }
                                      if(scen == 2){
                                        res <- AncCond(tree = cur_dat$tree,
                                                       data = df,
                                                       mat = c(0,2,1,0),
                                                       mc = iter)
                                      }
                                      temp.res[[sf.strength]] <- res
                                    }
                                    temp.res
                                }
    result.list[[scen]][[si]] <- parallel_results
  }
}




  sim.results <- result.list
  save(result.list,
       file = ("../results/sim.results_UL.RData"))
  
  
  
  
