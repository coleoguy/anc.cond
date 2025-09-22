#!/usr/bin/env Rscript
# Bidirectional scaling simulation described in the manuscript.

suppressPackageStartupMessages({
  library(phytools)
  library(diversitree)
  library(geiger)
})

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

set.seed(123)

branch_means <- function(tree, cont_trait) {
  cont_asr <- phytools::anc.ML(tree, cont_trait, model = "BM")
  means <- numeric(nrow(tree$edge))
  for (i in seq_len(nrow(tree$edge))) {
    parent <- tree$edge[i, 1]
    child <- tree$edge[i, 2]
    parent_val <- if (parent <= length(tree$tip.label)) cont_trait[parent] else cont_asr$ace[as.character(parent)]
    child_val <- if (child <= length(tree$tip.label)) cont_trait[child] else cont_asr$ace[as.character(child)]
    means[i] <- mean(c(parent_val, child_val))
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

simulate_dataset <- function(n_taxa, rate_matrix) {
  base_tree <- phytools::trees(
    pars = c(3, 1), type = "bd", n = 1, max.taxa = n_taxa, include.extinct = FALSE
  )[[1]]
  base_tree$edge.length <- base_tree$edge.length / max(phytools::branching.times(base_tree))
  cont_trait <- as.numeric(geiger::sim.char(base_tree, par = 0.2, model = "BM")[, 1])
  names(cont_trait) <- base_tree$tip.label
  bm <- branch_means(base_tree, cont_trait)
  list(tree = base_tree, cont_trait = cont_trait, branch_info = bm)
}

run_scaling_bidirectional <- function(n_trees = 10,
                                      n_taxa = 200,
                                      scaling_factors = c(1, 2, 4, 8),
                                      rate = 0.6,
                                      nsim = 10,
                                      iter = 200,
                                      verbose = FALSE) {
  rate_matrix <- matrix(c(-rate, rate, rate, -rate), nrow = 2, byrow = TRUE)
  results <- vector("list", length(scaling_factors))
  names(results) <- paste0("scale_", scaling_factors)
  for (sf in seq_along(scaling_factors)) {
    current_sf <- scaling_factors[sf]
    sf_results <- vector("list", n_trees)
    for (t in seq_len(n_trees)) {
      if (verbose) {
        message("Tree ", t, " of ", n_trees, ", scale factor ", current_sf)
      }
      dataset <- simulate_dataset(n_taxa, rate_matrix)
      scaled_tree <- scale_tree_by_branch_means(dataset$tree, dataset$branch_info, current_sf)
      disc_trait <- geiger::sim.char(
        scaled_tree,
        par = rate_matrix,
        model = "discrete",
        root = sample(1:2, 1)
      )[, 1]
      data_frame <- data.frame(
        taxon = dataset$tree$tip.label,
        cont = dataset$cont_trait,
        disc = disc_trait
      )
      anc_result <- AncCond(
        tree = scaled_tree,
        data = data_frame,
        nsim = nsim,
        iter = iter,
        mat = c(0, 2, 1, 0),
        message = FALSE
      )
      sf_results[[t]] <- c(p12 = anc_result$pvals[["12"]], p21 = anc_result$pvals[["21"]])
    }
    sf_results <- do.call(rbind, sf_results)
    results[[sf]] <- sf_results
  }
  results
}

if (sys.nframe() == 0) {
  dir.create(file.path(project_root, "results"), showWarnings = FALSE)
  scaling_results <- run_scaling_bidirectional()
  saveRDS(scaling_results, file = file.path(project_root, "results", "scaling_bidirectional_results.rds"))
}
