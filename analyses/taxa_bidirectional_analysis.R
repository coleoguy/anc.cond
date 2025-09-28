#!/usr/bin/env Rscript
# Power analysis varying the number of taxa under a bidirectional model.
# Parallel + cross-platform + fastBM + progress bar + checkpoints

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

set.seed(456)

# ---------- Cross-platform parallel plan ----------
n_workers <- max(1, future::availableCores() - 1)  # keep 1 core free
if (.Platform$OS.type == "windows") {
  plan(multisession, workers = n_workers)
} else {
  plan(multicore, workers = n_workers)
}
on.exit(plan(sequential), add = TRUE)

options(
  future.globals.maxSize = 2 * 1024^3,  # 2 GB per worker
  future.wait.timeout   = 30*60,
  future.rng.onMisuse   = "ignore"
)
handlers(global = TRUE)  # enable progress bars

# ---------- Helpers ----------
simulate_continuous <- function(tree, sigma) {
  cont_trait <- phytools::fastBM(tree, sig2 = sigma)
  names(cont_trait) <- tree$tip.label
  cont_trait
}

simulate_bidirectional_dataset <- function(n_taxa, sigma, rate) {
  tree <- sim.bdtree(b = 3, d = 1, stop = c("taxa"), n = n_taxa, extinct = FALSE)
  tree$edge.length <- tree$edge.length / max(ape::branching.times(tree))
  cont_trait <- simulate_continuous(tree, sigma)
  
  # simulate discrete; reject near-fixed outcomes
  rate_matrix <- matrix(c(-rate, rate, rate, -rate), nrow = 2, byrow = TRUE)
  good_sim <- FALSE; attempts <- 0L; max_attempts <- 1000L
  disc_trait <- NULL
  while (!good_sim && attempts < max_attempts) {
    attempts <- attempts + 1L
    sim_out <- geiger::sim.char(
      tree,
      par = rate_matrix,
      model = "discrete",
      root = sample(1:2, 1),
      nsim = 1
    )
    if (is.matrix(sim_out)) disc_trait <- sim_out[, 1, drop = TRUE] else disc_trait <- sim_out
    disc_trait <- as.integer(disc_trait)
    freq <- mean(disc_trait == min(disc_trait), na.rm = TRUE)
    good_sim <- is.finite(freq) && freq > 0.05 && freq < 0.95
  }
  if (!good_sim) disc_trait <- rep(NA_integer_, length(tree$tip.label))
  
  list(tree = tree, cont_trait = cont_trait, disc_trait = disc_trait)
}

# ---------- Main runner (parallel over trees WITH progress bar) ----------
run_taxa_bidirectional <- function(n_trees = 20,
                                   taxa_grid = c(20),   # e.g., seq(20, 200, length.out = 5)
                                   rate = 0.6,
                                   sigma = 0.2,
                                   nsim = 10,
                                   iter = 200,
                                   verbose = TRUE,
                                   checkpoint_dir = file.path(project_root, "results", "taxa_checkpoints")) {
  dir.create(checkpoint_dir, showWarnings = FALSE, recursive = TRUE)
  
  results <- vector("list", length(taxa_grid))
  names(results) <- paste0("taxa_", taxa_grid)
  
  for (n_taxa in taxa_grid) {
    if (verbose) message(sprintf("Taxa %s using %d workers", n_taxa, n_workers))
    
    with_progress({
      p <- progressor(steps = n_trees)  # one tick per tree at this taxa level
      
      taxa_results <- future_lapply(
        seq_len(n_trees),
        function(t) {
          p(sprintf("taxa %s: tree %d/%d", n_taxa, t, n_trees))
          
          ds <- simulate_bidirectional_dataset(n_taxa, sigma, rate)
          
          if (anyNA(ds$disc_trait)) {
            return(c(p12 = NA_real_, p21 = NA_real_))
          }
          
          df <- data.frame(
            taxon = ds$tree$tip.label,
            cont  = ds$cont_trait,
            disc  = ds$disc_trait
          )
          
          anc <- AncCond(
            tree = ds$tree,
            data = df,
            nsim = nsim,
            iter = iter,
            mat = c(0, 2, 1, 0),
            message = FALSE
          )
          
          c(p12 = anc$pvals[["12"]], p21 = anc$pvals[["21"]])
        },
        future.seed = TRUE
      )
      
      mat <- do.call(rbind, taxa_results)
      results[[paste0("taxa_", n_taxa)]] <- mat
      
      # checkpoint after each taxa level
      saveRDS(mat, file = file.path(checkpoint_dir, sprintf("tmp_taxa_%s.rds", n_taxa)))
      
      if (verbose) message(sprintf("Finished taxa %s: %d/%d trees", n_taxa, nrow(mat), n_trees))
    })
  }
  
  results
}

# ---------- Script entry ----------
if (sys.nframe() == 0) {
  dir.create(file.path(project_root, "results"), showWarnings = FALSE, recursive = TRUE)
  
  taxa_results <- run_taxa_bidirectional(
    n_trees = 20,
    taxa_grid = c(20),  # or: seq(20, 200, length.out = 5)
    rate = 0.6,
    sigma = 0.2,
    nsim = 10,
    iter = 200,
    verbose = TRUE
  )
  
  rds_path   <- file.path(project_root, "results", "taxa_bidirectional_results.rds")
  rdata_path <- file.path(project_root, "results", "taxa_bidirectional_results.RData")
  
  saveRDS(taxa_results, file = rds_path)
  sess <- sessionInfo()
  save(taxa_results, sess, file = rdata_path)
  
  message("Saved results to:\n  ", rds_path, "\n  ", rdata_path)
}
