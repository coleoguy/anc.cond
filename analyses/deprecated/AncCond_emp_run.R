# ============================================================
# FULL PIPELINE: Prune -> Simulate -> AncCond -> Save
# ============================================================

# 1. Load Libraries
library(geiger)
library(ape)
library(phytools)
library(foreach)
library(doParallel)

# 2. Setup Parallel Cluster
num_cores <- detectCores() - 2
my_cluster <- makeCluster(num_cores)
registerDoParallel(my_cluster)

# 3. Load & Rescale Tree (Root-to-Tip = 1.0)
tr_raw <- read.tree("whales.tre")
tr_scaled <- tr_raw
tr_scaled$edge.length <- tr_scaled$edge.length / max(nodeHeights(tr_raw))

# 4. Load & Clean Continuous Trait Data
cont_raw <- read.csv("whale_sizes.csv") 
cont_raw[,1] <- gsub(" ", "_", cont_raw[,1]) # Standardize names

# 5. Prune Tree to Match CSV Data
# We do this BEFORE simulation so the frequency target is accurate
common_species <- intersect(tr_scaled$tip.label, cont_raw[,1])
tr_final <- keep.tip(tr_scaled, common_species)

# Align continuous data to the pruned tree
matching_idx <- match(tr_final$tip.label, cont_raw[,1])
cont_final   <- cont_raw[matching_idx, 2] 
names(cont_final) <- tr_final$tip.label

# ------------------------------------------------------------
# 6. Simulate 100 Discrete Datasets (5-95% Frequency)
# ------------------------------------------------------------
forward_rate <- 0.2   # Same forward rate as AncCond
nsim <- 100
simulated_datasets <- vector("list", nsim)

q_matrix <- matrix(c(-forward_rate, forward_rate, 0, 0), nrow = 2, byrow = TRUE)
dimnames(q_matrix) <- list(c("1","2"), c("1","2"))

message(sprintf("Simulating on %d species...", length(tr_final$tip.label)))

for (ti in 1:nsim) {
  repeat {
    sim_out <- geiger::sim.char(tr_final, par = q_matrix, model = "discrete", root = 1, nsim = 1)
    tips <- as.character(as.integer(sim_out[,,1]))
    names(tips) <- tr_final$tip.label
    
    prop_state2 <- sum(tips == "2") / length(tips)
    if (prop_state2 >= 0.05 && prop_state2 <= 0.95) {
      simulated_datasets[[ti]] <- tips
      break
    }
  }
}

# ------------------------------------------------------------
# 7. Run AncCond in Parallel
# ------------------------------------------------------------
# Define absolute path for workers to find the script in 'analyses'
anc_script_path <- normalizePath("../analyses/AncCond.R")
iter <- 1000

message("Running AncCond simulations...")

result.list <- foreach(i = seq_along(simulated_datasets), 
                       .packages = c("ape", "phytools"),
                       .export = c("tr_final", "cont_final", "iter", "anc_script_path")) %dopar% {
                         
                         source(anc_script_path, local = TRUE)
                         
                         # Align current simulation to tree
                         disc_final <- simulated_datasets[[i]]
                         
                         df <- data.frame(
                           species = tr_final$tip.label,
                           cont    = cont_final,
                           disc    = as.numeric(disc_final)
                         )
                         
                         # Run AncCond
                         res <- AncCond(tree = tr_final, data = df, mat = c(0, 0, 1, 0), mc = iter)
                         return(res)
                       }

# 8. Save and Cleanup
# Using the new unique filename to avoid overwriting original sim results
if(!dir.exists("../results")) dir.create("../results")
save(result.list, file = "../results/empirical_sim_results0.05.RData")

stopCluster(my_cluster)
message("Pipeline Complete! Results saved to ../results/empirical_sim_results.RData")