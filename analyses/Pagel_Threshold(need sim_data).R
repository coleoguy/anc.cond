library(phytools)
library(coda)
library(dplyr)
library(doParallel)
library(foreach)

# 1. Setup Parallel Backend
num_cores <- parallel::detectCores() - 1
cl <- makeCluster(num_cores)
registerDoParallel(cl)

message(sprintf("Running in parallel using %d cores...", num_cores))

# ... [Keep your helper functions: rescale_tree_on_the_fly, etc. here] ...

# 2. Parallel Loop
# We combine the directions and sf_vals into a single 'task list' to feed the cores
results_list <- foreach(dir = c("uni", "bi"), .combine = 'rbind') %:%
  foreach(r = 1:100, .combine = 'rbind', .packages = c("phytools", "coda", "dplyr")) %dopar% {
    
    # Locate data for this specific replicate
    rep_obj   <- sim_data[[dir]][["200"]][[r]]
    base_tree <- get_base_tree(rep_obj)
    cont_all  <- get_cont_trait(rep_obj)
    disc_list <- tryCatch(rep_obj$disc_by_sf, error = function(e) NULL)
    
    # We'll store the sf=1 and sf=5 results for this specific replicate
    rep_results <- data.frame()
    
    for (sf in c(1, 5)) {
      # 1. Scale tree
      current_tree <- rescale_tree_on_the_fly(base_tree, cont_all, sf)
      
      # 2. Align data
      cont_tip <- align_to_tips(cont_all, current_tree)
      d_trait  <- if(!is.null(disc_list)) disc_list[[sf]] else rep_obj$disc_trait
      disc_tip <- align_to_tips(d_trait, current_tree)
      
      # 3. Run Tests
      pagel_p <- run_pagel(current_tree, disc_tip, cont_tip)
      thr <- run_threshold(current_tree, disc_tip, cont_tip, ngen=50000, sample=1000)
      
      # 4. Bind row
      tmp <- data.frame(
        dir = dir, taxa = 200, rep = r, sf = sf,
        pagel_p = pagel_p, thresh_sig = thr$sig,
        r_hpd_lo = thr$r_lo, r_hpd_hi = thr$r_hi,
        stringsAsFactors = FALSE
      )
      rep_results <- rbind(rep_results, tmp)
    }
    return(rep_results)
  }

# 3. Shutdown Cluster
stopCluster(cl)

# 4. Process Results
res_df <- results_list %>%
  mutate(pagel_sig = !is.na(pagel_p) & pagel_p <= 0.05)

print(res_df %>% group_by(dir, sf) %>% 
        summarise(pagel_rate = mean(pagel_sig, na.rm=T), 
                  thresh_rate = mean(thresh_sig, na.rm=T), .groups="drop"))