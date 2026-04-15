library(phytools); library(coda); library(dplyr); library(doParallel); library(foreach)

# 1. SETUP
load("../data/sim_data_UL.RData") 
if(!exists("sim_data") && exists("sim_data_UL")) sim_data <- sim_data_UL

num_cores <- parallel::detectCores() - 1
cl <- makeCluster(num_cores)
registerDoParallel(cl)

message(sprintf("Launching 1000-rep Whole-Tree Analysis on %d cores...", num_cores))

# 2. MAIN LOOP
results_list <- foreach(dir = c("uni", "bi"), .combine = 'rbind') %:%
  foreach(r = 1:1000, 
          .combine = 'rbind', 
          .packages = c("phytools", "coda", "dplyr"), 
          .export = "sim_data") %dopar% {
            
            # --- A. Data Retrieval & Cleaning ---
            rep_obj   <- sim_data[[dir]][["200"]][[r]]
            base_tree <- rep_obj$tree
            n_taxa    <- length(base_tree$tip.label)
            
            # Force tips to a named vector
            cont_tip <- as.vector(rep_obj$cont_trait)
            names(cont_tip) <- base_tree$tip.label
            
            rep_results <- data.frame()
            
            # --- B. Scaling and Testing Loop ---
            for (sf in c(1, 5)) {
              current_tree <- base_tree
              
              if (sf > 1) {
                # 1. RECONSTRUCT ANCESTRAL STATES (since cont_trait_all is missing)
                asr <- try(phytools::anc.ML(base_tree, cont_tip, model = "BM"), silent = TRUE)
                
                if(!inherits(asr, "try-error")) {
                  # Combine Tips (1:200) and Reconstructed Nodes (201:399)
                  all_traits <- c(cont_tip, asr$ace)
                  
                  # 2. CALCULATE BRANCH MEANS FOR WHOLE TREE
                  branch_means <- numeric(nrow(base_tree$edge))
                  for(j in 1:nrow(base_tree$edge)){
                    p_node <- base_tree$edge[j,1]
                    c_node <- base_tree$edge[j,2]
                    branch_means[j] <- mean(c(all_traits[p_node], all_traits[c_node]))
                  }
                  
                  # 3. APPLY SCALING
                  qs <- quantile(branch_means, probs = c(0.25, 0.75))
                  current_tree$edge.length[branch_means < qs[1]] <- base_tree$edge.length[branch_means < qs[1]] / sf
                  current_tree$edge.length[branch_means > qs[2]] <- base_tree$edge.length[branch_means > qs[2]] * sf
                }
              }
              
              # --- C. ALIGN & ANALYZE ---
              cont_bin <- as.factor(ifelse(cont_tip > median(cont_tip), "high", "low"))
              names(cont_bin) <- names(cont_tip)
              
              sf_key   <- paste0("sf", sf)
              d_trait  <- if(!is.null(rep_obj$disc_by_sf)) rep_obj$disc_by_sf[[sf_key]] else rep_obj$disc_trait
              disc_tip <- as.factor(as.vector(d_trait))
              names(disc_tip) <- base_tree$tip.label
              
              # Test 1: Pagel
              p_val <- tryCatch({ fitPagel(current_tree, disc_tip, cont_bin, method='fitDiscrete')$P }, error = function(e) NA)
              
              # Test 2: Threshold (MCMC)
              X <- cbind(as.numeric(disc_tip) - 1, cont_tip)
              rownames(X) <- names(cont_tip)
              
              thr <- tryCatch({
                res <- threshBayes(current_tree, X, ngen=50000, control=list(sample=1000, quiet=TRUE))
                post_r <- res$par[(10000/1000 + 1):nrow(res$par), "r"] # 20% burn-in
                hpd <- HPDinterval(as.mcmc(post_r))
                list(sig = as.numeric(hpd[1] > 0 | hpd[2] < 0), lo = hpd[1], hi = hpd[2])
              }, error = function(e) list(sig=NA, lo=NA, hi=NA))
              
              rep_results <- rbind(rep_results, data.frame(
                scenario = dir, rep = r, sf = sf, pagel_p = p_val, 
                thresh_sig = thr$sig, r_hpd_lo = thr$lo, r_hpd_hi = thr$hi, stringsAsFactors = FALSE))
            }
            return(rep_results)
          }

summary_stats <- results_list %>%
  group_by(scenario, sf) %>%
  summarise(
    n_reps = n(),
    pagel_sig_rate = mean(pagel_p < 0.05, na.rm = TRUE),
    thresh_sig_rate = mean(thresh_sig, na.rm = TRUE),
    avg_r = mean((r_hpd_lo + r_hpd_hi)/2, na.rm = TRUE)
  )

print(summary_stats)

# 3. CLEANUP
stopCluster(cl)
save(results_list, file = "Final_Results_WholeTree.RData")
message("Run Complete. Results saved to Final_Results_WholeTree.RData")