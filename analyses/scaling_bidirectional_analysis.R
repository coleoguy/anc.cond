#!/usr/bin/env Rscript
# Scaling bidirectional analysis — minimal edits to ensure saving results

suppressPackageStartupMessages({
  library(phytools)
  library(diversitree)
  library(geiger)
  library(future)
  library(future.apply)
  library(progressr)
})

# ---------- Project root ----------
locate_project_root <- function() {
  file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(file_arg) == 1) return(normalizePath(file.path(dirname(sub("^--file=", "", file_arg)), "..")))
  if (file.exists(file.path("R", "anc_cond.R"))) return(normalizePath("."))
  normalizePath("..")
}
project_root <- locate_project_root()
source(file.path(project_root, "R", "anc_cond.R"))

set.seed(123)

# Run even when sourcing (set FALSE to restore “Rscript-only” behavior)
force_run <- TRUE

# ---------- Parallel plan ----------
n_workers <- max(1, future::availableCores() - 1)
if (.Platform$OS.type == "windows") {
  plan(multisession, workers = n_workers)
} else {
  plan(multicore, workers = n_workers)
}
on.exit(plan(sequential), add = TRUE)

options(
  future.globals.maxSize = 2 * 1024^3,
  future.wait.timeout   = 30*60,
  future.rng.onMisuse   = "ignore"
)
handlers(global = TRUE)

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
  quantiles <- stats::quantile(branch_info$means, probs = c(0.25, 0.75))
  below <- branch_info$means < quantiles[[1]]
  above <- branch_info$means > quantiles[[2]]
  scaled$edge.length[below] <- scaled$edge.length[below] / scale_factor
  scaled$edge.length[above] <- scaled$edge.length[above] * scale_factor
  scaled
}

simulate_dataset <- function(n_taxa, cont_sigma) {
  base_tree <- sim.bdtree(b = 3, d = 1, stop = c("taxa"), n = n_taxa, extinct = FALSE)
  base_tree$edge.length <- base_tree$edge.length / max(ape::branching.times(base_tree))
  # swap to fastBM for speed; keep sigma naming the same
  cont_trait <- phytools::fastBM(base_tree, sig2 = cont_sigma)
  names(cont_trait) <- base_tree$tip.label
  bm <- branch_means(base_tree, cont_trait)
  list(tree = base_tree, cont_trait = cont_trait, branch_info = bm)
}

# ---------- Main ----------
run_scaling_bidirectional <- function(n_trees = 10,
                                      n_taxa = 200,
                                      scaling_factors = c(1, 2, 4, 8),
                                      rate = 0.6,
                                      cont_sigma = 0.2,
                                      nsim = 10,
                                      iter = 200,
                                      verbose = TRUE,
                                      checkpoint_dir = file.path(project_root, "results", "scale_bi_checkpoints")) {
  
  dir.create(checkpoint_dir, showWarnings = FALSE, recursive = TRUE)
  
  rate_matrix <- matrix(c(-rate, rate, rate, -rate), nrow = 2, byrow = TRUE)
  results <- vector("list", length(scaling_factors))
  names(results) <- paste0("scale_", scaling_factors)
  
  with_progress({
    p <- progressor(steps = length(scaling_factors))
    
    for (sf in seq_along(scaling_factors)) {
      current_sf <- scaling_factors[sf]
      if (verbose) message(sprintf("Scale %s using %d workers", current_sf, n_workers))
      
      sf_results <- future_lapply(seq_len(n_trees), function(t) {
        ds <- simulate_dataset(n_taxa, cont_sigma)
        scaled_tree <- scale_tree_by_branch_means(ds$tree, ds$branch_info, current_sf)
        
        sim_out <- geiger::sim.char(
          scaled_tree,
          par = rate_matrix,
          model = "discrete",
          root = sample(1:2, 1),
          nsim = 1
        )
        disc_trait <- if (is.matrix(sim_out)) sim_out[, 1, drop = TRUE] else sim_out
        disc_trait <- as.integer(disc_trait)
        
        df <- data.frame(
          taxon = ds$tree$tip.label,
          cont  = ds$cont_trait,
          disc  = disc_trait
        )
        
        anc <- AncCond(
          tree = scaled_tree,
          data = df,
          nsim = nsim,
          iter = iter,
          mat  = c(0, 2, 1, 0),
          message = FALSE
        )
        
        c(p12 = anc$pvals[["12"]], p21 = anc$pvals[["21"]])
      }, future.seed = TRUE)
      
      sf_mat <- do.call(rbind, sf_results)
      results[[sf]] <- sf_mat
      
      # checkpoint per scale
      saveRDS(sf_mat, file = file.path(checkpoint_dir, sprintf("tmp_scale_bi_%s.rds", current_sf)))
      if (verbose) message(sprintf("Finished scale %s: %d/%d trees", current_sf, nrow(sf_mat), n_trees))
      
      p()
    }
  })
  
  results
}

# ---------- Entry point (ensures saving) ----------
if (force_run || sys.nframe() == 0) {
  dir.create(file.path(project_root, "results"), showWarnings = FALSE, recursive = TRUE)
  
  scaling_results <- run_scaling_bidirectional(
    n_trees = 10,
    n_taxa = 200,
    scaling_factors = c(1, 2, 4, 8),
    rate = 0.6,
    cont_sigma = 0.2,
    nsim = 10,
    iter = 200,
    verbose = TRUE
  )
  
  rds_path   <- file.path(project_root, "results", "scaling_bidirectional_results.rds")
  rdata_path <- file.path(project_root, "results", "scaling_bidirectional_results.RData")
  
  saveRDS(scaling_results, file = rds_path)
  sess <- sessionInfo()
  save(scaling_results, sess, file = rdata_path)
  
  message("Saved results to:\n  ", rds_path, "\n  ", rdata_path)
}
