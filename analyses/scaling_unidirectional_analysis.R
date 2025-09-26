#!/usr/bin/env Rscript
# Unidirectional scaling simulation focusing on gain of the derived state.
# Parallel + cross-platform + fastBM version

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
  if (length(file_arg) == 1) {
    return(normalizePath(file.path(dirname(sub("^--file=", "", file_arg)), "..")))
  }
  if (file.exists(file.path("R", "anc_cond.R"))) {
    return(normalizePath("."))
  }
  normalizePath("..")
}
project_root <- locate_project_root()
source(file.path(project_root, "R", "anc_cond.R"))

set.seed(321)

# ---------- Cross-platform parallel plan ----------
# Windows: multisession (process-based). macOS/Linux: multicore (fork).
n_workers <- max(1, future::availableCores() - 1)  # keep 1 core free
if (.Platform$OS.type == "windows") {
  plan(multisession, workers = n_workers)
} else {
  plan(multicore, workers = n_workers)
}
on.exit(plan(sequential), add = TRUE)   # restore when we exit

# Reasonable safety options
options(
  future.globals.maxSize = 2 * 1024^3,   # 2 GB per worker
  future.wait.timeout   = 30*60,         # 30 min wait per future
  future.rng.onMisuse   = "ignore"
)
handlers(global = TRUE)  # progressr

# ---------- Helpers ----------
simulate_unidirectional_tree <- function(n_taxa, cont_sigma) {
  tree <- sim.bdtree(b = 3, d = 1, stop = c("taxa"), n = n_taxa, extinct = FALSE)
  tree$edge.length <- tree$edge.length / max(ape::branching.times(tree))
  # fastBM is much faster than geiger::sim.char(..., model="BM")
  cont_trait <- phytools::fastBM(tree, sig2 = cont_sigma)
  # Names already set by fastBM; keep for safety:
  names(cont_trait) <- tree$tip.label
  list(tree = tree, cont_trait = cont_trait)
}

compute_branch_means <- function(tree, cont_trait) {
  cont_asr <- phytools::anc.ML(tree, cont_trait, model = "BM")
  means <- numeric(nrow(tree$edge))
  for (i in seq_len(nrow(tree$edge))) {
    parent <- tree$edge[i, 1]
    child  <- tree$edge[i, 2]
    parent_val <- if (parent <= length(tree$tip.label)) cont_trait[parent] else cont_asr$ace[as.character(parent)]
    child_val  <- if (child  <= length(tree$tip.label)) cont_trait[child]  else cont_asr$ace[as.character(child)]
    means[i] <- (parent_val + child_val) / 2
  }
  means
}

scale_edges <- function(tree, branch_means, scale_factor) {
  scaled <- tree
  qs <- stats::quantile(branch_means, probs = c(0.25, 0.75))
  scaled$edge.length[branch_means < qs[[1]]] <- scaled$edge.length[branch_means < qs[[1]]] / scale_factor
  scaled$edge.length[branch_means > qs[[2]]] <- scaled$edge.length[branch_means > qs[[2]]] * scale_factor
  scaled
}

# ---------- Main runner (parallel over trees) ----------
run_scaling_unidirectional <- function(n_trees = 10,
                                       n_taxa = 100,
                                       scaling_factors = c(1),   # e.g., c(1, 2, 4, 8)
                                       forward_rate = 0.2,
                                       reverse_rate = 1e-04,
                                       cont_sigma   = 0.2,
                                       nsim = 10,
                                       iter = 200,
                                       verbose = FALSE,
                                       checkpoint_dir = file.path(project_root, "results", "checkpoints")) {
  
  dir.create(checkpoint_dir, showWarnings = FALSE, recursive = TRUE)
  
  # Q with rows summing to zero
  q_matrix <- matrix(c(-forward_rate, forward_rate,
                       reverse_rate, -reverse_rate),
                     nrow = 2, byrow = TRUE)
  
  results <- vector("list", length(scaling_factors))
  names(results) <- paste0("scale_", scaling_factors)
  
  with_progress({
    p <- progressor(along = scaling_factors)
    
    for (sf_idx in seq_along(scaling_factors)) {
      sf <- scaling_factors[sf_idx]
      if (verbose) message(sprintf("Scale %s using %d workers", sf, n_workers))
      
      # Parallelize across trees
      tree_results <- future_lapply(seq_len(n_trees), function(t) {
        # ---- per-tree body (was the inner for-loop) ----
        dataset <- simulate_unidirectional_tree(n_taxa, cont_sigma)
        
        bm <- compute_branch_means(dataset$tree, dataset$cont_trait)
        scaled_tree <- scale_edges(dataset$tree, bm, sf)
        
        # Draw discrete trait until not ~fixed at a boundary
        good_sim <- FALSE
        attempts <- 0L
        max_attempts <- 1000L
        disc_trait <- NULL
        
        while (!good_sim && attempts < max_attempts) {
          attempts <- attempts + 1L
          # geiger::sim.char returns a matrix when nsim=1; coerce to vector of states
          sim_out <- geiger::sim.char(
            scaled_tree,
            par = q_matrix,
            model = "discrete",
            root = 1,
            nsim = 1
          )
          
          # Convert to a simple integer vector of states at the tips
          # sim_out can be matrix (tips x 1) or vector; handle both:
          if (is.matrix(sim_out)) disc_trait <- sim_out[, 1, drop = TRUE] else disc_trait <- sim_out
          disc_trait <- as.integer(disc_trait)
          
          freq <- mean(disc_trait == min(disc_trait), na.rm = TRUE)
          good_sim <- is.finite(freq) && freq > 0.05 && freq < 0.95
        }
        
        if (!good_sim) {
          # If we really can't get a good draw, return NAs (keeps run moving)
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
          drop.state = 2,
          nsim = nsim,
          iter = iter,
          mat = c(0, 0, 1, 0),
          message = FALSE
        )
        
        c(p12 = anc$pvals[["12"]], p21 = anc$pvals[["21"]])
      }, future.seed = TRUE)  # reproducible parallel RNG
      
      sf_mat <- do.call(rbind, tree_results)
      results[[sf_idx]] <- sf_mat
      
      # checkpoint each scaling factor
      saveRDS(sf_mat, file = file.path(checkpoint_dir, sprintf("tmp_scale_%s.rds", sf)))
      
      if (verbose) {
        msg <- sprintf("Finished scale %s: %d/%d trees", sf, nrow(sf_mat), n_trees)
        message(msg)
      }
      p()
    }
  })
  
  results
}

# ---------- Script entry point ----------
if (sys.nframe() == 0) {
  dir.create(file.path(project_root, "results"), showWarnings = FALSE, recursive = TRUE)
  
  scaling_results <- run_scaling_unidirectional(
    n_trees = 10,
    n_taxa  = 100,
    scaling_factors = c(1),  # change as needed
    forward_rate = 0.2,
    reverse_rate = 1e-04,
    cont_sigma   = 0.2,
    nsim = 10,
    iter = 200,
    verbose = TRUE
  )
  
  # 1) Compact, single-object file (recommended for programmatic use)
  rds_path <- file.path(project_root, "results", "scaling_unidirectional_results.rds")
  saveRDS(scaling_results, file = rds_path)
  
  # 2) Full RData snapshot (handy for interactive work/resuming)
  rdata_path <- file.path(project_root, "results", "scaling_unidirectional_results.RData")
  # Save exactly the key objects you’ll want later; add others if useful
  sess <- sessionInfo()
  save(scaling_results, sess, file = rdata_path)  # version=2 optional
  
  message("Saved results to:\n  ", rds_path, "\n  ", rdata_path)
}
View(as.data.frame(scaling_results))
