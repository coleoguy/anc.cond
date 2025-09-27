#!/usr/bin/env Rscript
# Bidirectional scaling simulation described in the manuscript.

#!/usr/bin/env Rscript
# Bidirectional scaling simulation — parallel + cross-platform + fastBM

suppressPackageStartupMessages({
  library(phytools)
  library(diversitree)
  library(geiger)
  library(future)
  library(future.apply)
})

# ---------- Project root ----------
locate_project_root <- function() {
  file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(file_arg) == 1) {
    return(normalizePath(file.path(dirname(sub("^--file=", "", file_arg)), "..")))
  }
  if (file.exists(file.path("R", "anc_cond.R"))) return(normalizePath("."))
  normalizePath("..")
}
project_root <- locate_project_root()
source(file.path(project_root, "R", "anc_cond.R"))

set.seed(123)

# ---------- Cross-platform parallel plan ----------
# Windows uses multisession (processes); macOS/Linux can use multicore (forks)
n_workers <- max(1, future::availableCores() - 1)  # keep one core free
if (.Platform$OS.type == "windows") {
  plan(multisession, workers = n_workers)
} else {
  plan(multicore, workers = n_workers)
}
on.exit(plan(sequential), add = TRUE)   # restore when the script exits

# Some conservative safety options
options(
  future.globals.maxSize = 2 * 1024^3,  # 2 GB per worker
  future.wait.timeout   = 30*60,        # 30 min per future
  future.rng.onMisuse   = "ignore"
)

# ---------- Helpers ----------
branch_means <- function(tree, cont_trait) {
  cont_asr <- phytools::anc.ML(tree, cont_trait, model = "BM")
  means <- numeric(nrow(tree$edge))
  for (i in seq_len(nrow(tree$edge))) {
    parent <- tree$edge[i, 1]
    child  <- tree$edge[i, 2]
    parent_val <- if (parent <= length(tree$tip.label)) cont_trait[parent] else cont_asr$ace[as.character(parent)]
    child_val  <- if (child  <= length(tree$tip.label)) cont_trait[child]  else cont_asr$ace[as.character(child)]
    means[i] <- (parent_val + child_val) / 2
  }
  list(means = means, asr = cont_asr)
}

scale_tree_by_branch_means <- function(tree, branch_info, scale_factor) {
  scaled <- tree
  qs <- stats::quantile(branch_info$means, probs = c(0.25, 0.75))
  below <- branch_info$means < qs[[1]]
  above <- branch_info$means > qs[[2]]
  scaled$edge.length[below] <- scaled$edge.length[below] / scale_factor
  scaled$edge.length[above] <- scaled$edge.length[above] * scale_factor
  scaled
}

simulate_dataset <- function(n_taxa, cont_sigma) {
  base_tree <- sim.bdtree(b = 3, d = 1, stop = c("taxa"), n = n_taxa, extinct = FALSE)
  base_tree$edge.length <- base_tree$edge.length / max(ape::branching.times(base_tree))
  # Faster continuous trait simulation
  cont_trait <- phytools::fastBM(base_tree, sig2 = cont_sigma)
  names(cont_trait) <- base_tree$tip.label
  bm <- branch_means(base_tree, cont_trait)
  list(tree = base_tree, cont_trait = cont_trait, branch_info = bm)
}

# ---------- Main runner (parallel over trees for each scale) ----------
run_scaling_bidirectional <- function(n_trees = 10,
                                      n_taxa = 200,
                                      scaling_factors = c(1) #, 2, 4, 8),
                                      rate = 0.6,
                                      cont_sigma = 0.2,
                                      nsim = 10,
                                      iter = 200,
                                      verbose = FALSE,
                                      checkpoint_dir = file.path(project_root, "results", "checkpoints")) {
  dir.create(checkpoint_dir, showWarnings = FALSE, recursive = TRUE)
  
  # Two-state bidirectional Q (rows sum to zero)
  rate_matrix <- matrix(c(-rate, rate,
                          rate, -rate),
                        nrow = 2, byrow = TRUE)
  
  results <- vector("list", length(scaling_factors))
  names(results) <- paste0("scale_", scaling_factors)
  
  for (sf_idx in seq_along(scaling_factors)) {
    current_sf <- scaling_factors[sf_idx]
    if (verbose) message(sprintf("Scale %s using %d workers", current_sf, n_workers))
    
    # Parallelize the per-tree work
    tree_results <- future_lapply(seq_len(n_trees), function(t) {
      # Simulate base tree + continuous trait
      dataset <- simulate_dataset(n_taxa, cont_sigma)
      # Scale edges by branch means
      scaled_tree <- scale_tree_by_branch_means(dataset$tree, dataset$branch_info, current_sf)
      
      # Simulate discrete trait (retry if essentially all tips in one state)
      good_sim <- FALSE
      attempts <- 0L
      max_attempts <- 1000L
      disc_trait <- NULL
      
      while (!good_sim && attempts < max_attempts) {
        attempts <- attempts + 1L
        sim_out <- geiger::sim.char(
          scaled_tree,
          par   = rate_matrix,
          model = "discrete",
          root  = sample(1:2, 1),
          nsim  = 1
        )
        # sim.char can return matrix (tips x 1) or vector—coerce to simple vector
        if (is.matrix(sim_out)) disc_trait <- sim_out[, 1, drop = TRUE] else disc_trait <- sim_out
        disc_trait <- as.integer(disc_trait)
        
        freq <- mean(disc_trait == min(disc_trait), na.rm = TRUE)
        good_sim <- is.finite(freq) && freq > 0.05 && freq < 0.95
      }
      
      if (!good_sim) {
        # Could not obtain a balanced mapping; keep run moving
        return(c(p12 = NA_real_, p21 = NA_real_))
      }
      
      df <- data.frame(
        taxon = scaled_tree$tip.label,
        cont  = dataset$cont_trait,
        disc  = disc_trait
      )
      
      anc <- AncCond(
        tree = scaled_tree,
        data = df,
        # For bidirectional we typically *don’t* drop a state:
        # drop.state = NULL,
        nsim = nsim,
        iter = iter,
        mat = c(0, 2, 1, 0),  # allow 1<->2
        message = FALSE
      )
      
      c(p12 = anc$pvals[["12"]], p21 = anc$pvals[["21"]])
    }, future.seed = TRUE)  # reproducible parallel RNG
    
    sf_mat <- do.call(rbind, tree_results)
    results[[sf_idx]] <- sf_mat
    
    # Optional per-scale checkpoint
    saveRDS(sf_mat, file = file.path(checkpoint_dir, sprintf("tmp_scale_%s.rds", current_sf)))
    
    if (verbose) {
      message(sprintf("Finished scale %s: %d/%d trees", current_sf, nrow(sf_mat), n_trees))
    }
  }
  
  results
}

# ---------- Script entry ----------
if (sys.nframe() == 0) {
  dir.create(file.path(project_root, "results"), showWarnings = FALSE, recursive = TRUE)
  
  scaling_results <- run_scaling_bidirectional(
    n_trees = 10,
    n_taxa  = 200,
    scaling_factors = c(1), # 2, 4, 8),  # <- run all four scales
    rate = 0.6,
    cont_sigma = 0.2,
    nsim = 10,
    iter = 200,
    verbose = TRUE
  )
  
  # Save compact object
  rds_path <- file.path(project_root, "results", "scaling_bidirectional_results.rds")
  saveRDS(scaling_results, file = rds_path)
  
  # Save full snapshot for interactive work/resume
  rdata_path <- file.path(project_root, "results", "scaling_bidirectional_results.RData")
  sess <- sessionInfo()
  save(scaling_results, sess, file = rdata_path)
  
  message("Saved results to:\n  ", rds_path, "\n  ", rdata_path)
}
