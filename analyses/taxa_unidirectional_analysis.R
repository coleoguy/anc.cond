#!/usr/bin/env Rscript
# Power analysis varying the number of taxa under a nearly unidirectional model.

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

set.seed(654)

simulate_unidirectional_dataset <- function(n_taxa, sigma, forward_rate, reverse_rate) {
  tree <- sim.bdtree(b = 3, d = 1, stop = c("taxa"), n = n_taxa, extinct = FALSE)
  tree$edge.length <- tree$edge.length / max(ape::branching.times(tree))
  cont_trait <- as.numeric(geiger::sim.char(tree, par = sigma, model = "BM", nsim = 1))
  names(cont_trait) <- tree$tip.label
  q_matrix <- matrix(c(-forward_rate, forward_rate, reverse_rate, -reverse_rate), nrow = 2, byrow = TRUE)
  good_sim <- FALSE
  while (!good_sim) {
    disc_trait <- geiger::sim.char(
      tree,
      par = q_matrix,
      model = "discrete",
      root = 1,
      nsim = 1
    )
    freq <- mean(disc_trait == min(disc_trait))
    good_sim <- freq > 0.05 && freq < 0.95
  }
  list(tree = tree, cont_trait = cont_trait, disc_trait = disc_trait)
}

run_taxa_unidirectional <- function(n_trees = 20,
                                    #taxa_grid = seq(20, 200, length.out = 5),
                                    taxa_grid = c(20),
                                    forward_rate = 0.2,
                                    reverse_rate = 1e-04,
                                    sigma = 0.2,
                                    nsim = 10,
                                    iter = 200,
                                    verbose = FALSE) {
  results <- vector("list", length(taxa_grid))
  names(results) <- paste0("taxa_", taxa_grid)
  for (i in seq_along(taxa_grid)) {
    n_taxa <- taxa_grid[i]
    taxa_results <- vector("list", n_trees)
    for (t in seq_len(n_trees)) {
      if (verbose) {
        message("Tree ", t, " of ", n_trees, ", taxa ", n_taxa)
      }
      dataset <- simulate_unidirectional_dataset(n_taxa, sigma, forward_rate, reverse_rate)
      data_frame <- data.frame(
        taxon = dataset$tree$tip.label,
        cont = dataset$cont_trait,
        disc = dataset$disc_trait
      )
      anc_result <- AncCond(
        tree = dataset$tree,
        data = data_frame,
        drop.state = 2,
        nsim = nsim,
        iter = iter,
        mat = c(0, 0, 1, 0),
        message = FALSE
      )
      taxa_results[[t]] <- c(p12 = anc_result$pvals[["12"]], p21 = anc_result$pvals[["21"]])
    }
    results[[i]] <- do.call(rbind, taxa_results)
  }
  results
}

if (sys.nframe() == 0) {
  dir.create(file.path(project_root, "results"), showWarnings = FALSE)
  taxa_results <- run_taxa_unidirectional()
  saveRDS(taxa_results, file = file.path(project_root, "results", "taxa_unidirectional_results.rds"))
}
