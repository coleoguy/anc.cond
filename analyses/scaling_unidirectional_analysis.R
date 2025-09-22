#!/usr/bin/env Rscript
# Unidirectional scaling simulation focusing on gain of the derived state.

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

set.seed(321)

simulate_unidirectional_tree <- function(n_taxa, cont_sigma) {
  tree <- phytools::trees(
    pars = c(3, 1), type = "bd", n = 1, max.taxa = n_taxa, include.extinct = FALSE
  )[[1]]
  tree$edge.length <- tree$edge.length / max(phytools::branching.times(tree))
  cont_trait <- as.numeric(geiger::sim.char(tree, par = cont_sigma, model = "BM")[, 1])
  names(cont_trait) <- tree$tip.label
  list(tree = tree, cont_trait = cont_trait)
}

compute_branch_means <- function(tree, cont_trait) {
  cont_asr <- phytools::anc.ML(tree, cont_trait, model = "BM")
  means <- numeric(nrow(tree$edge))
  for (i in seq_len(nrow(tree$edge))) {
    parent <- tree$edge[i, 1]
    child <- tree$edge[i, 2]
    parent_val <- if (parent <= length(tree$tip.label)) cont_trait[parent] else cont_asr$ace[as.character(parent)]
    child_val <- if (child <= length(tree$tip.label)) cont_trait[child] else cont_asr$ace[as.character(child)]
    means[i] <- mean(c(parent_val, child_val))
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

run_scaling_unidirectional <- function(n_trees = 10,
                                       n_taxa = 100,
                                       scaling_factors = c(1, 2, 4, 8),
                                       forward_rate = 0.2,
                                       reverse_rate = 1e-04,
                                       cont_sigma = 0.2,
                                       nsim = 10,
                                       iter = 200,
                                       verbose = FALSE) {
  q_matrix <- matrix(c(-forward_rate, forward_rate, reverse_rate, -reverse_rate), nrow = 2)
  results <- vector("list", length(scaling_factors))
  names(results) <- paste0("scale_", scaling_factors)
  for (sf in seq_along(scaling_factors)) {
    sf_results <- vector("list", n_trees)
    for (t in seq_len(n_trees)) {
      if (verbose) {
        message("Tree ", t, " of ", n_trees, ", scale factor ", scaling_factors[sf])
      }
      dataset <- simulate_unidirectional_tree(n_taxa, cont_sigma)
      branch_means <- compute_branch_means(dataset$tree, dataset$cont_trait)
      scaled_tree <- scale_edges(dataset$tree, branch_means, scaling_factors[sf])
      good_sim <- FALSE
      while (!good_sim) {
        disc_trait <- geiger::sim.char(
          scaled_tree,
          par = q_matrix,
          model = "discrete",
          root = 1
        )[, 1]
        freq <- mean(disc_trait == min(disc_trait))
        good_sim <- freq > 0.05 && freq < 0.95
      }
      data_frame <- data.frame(
        taxon = scaled_tree$tip.label,
        cont = dataset$cont_trait,
        disc = disc_trait
      )
      anc_result <- AncCond(
        tree = scaled_tree,
        data = data_frame,
        drop.state = 2,
        nsim = nsim,
        iter = iter,
        mat = c(0, 0, 1, 0),
        message = FALSE
      )
      sf_results[[t]] <- c(p12 = anc_result$pvals[["12"]], p21 = anc_result$pvals[["21"]])
    }
    results[[sf]] <- do.call(rbind, sf_results)
  }
  results
}

if (sys.nframe() == 0) {
  dir.create(file.path(project_root, "results"), showWarnings = FALSE)
  scaling_results <- run_scaling_unidirectional()
  saveRDS(scaling_results, file = file.path(project_root, "results", "scaling_unidirectional_results.rds"))
}
