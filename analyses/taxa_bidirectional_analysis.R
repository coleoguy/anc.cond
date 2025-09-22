#!/usr/bin/env Rscript
# Power analysis varying the number of taxa under a bidirectional model.

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

set.seed(456)

simulate_continuous <- function(tree, sigma) {
  cont_trait <- as.numeric(geiger::sim.char(tree, par = sigma, model = "BM")[, 1])
  names(cont_trait) <- tree$tip.label
  cont_trait
}

simulate_bidirectional_dataset <- function(n_taxa, sigma, rate) {
  tree <- phytools::trees(
    pars = c(3, 1), type = "bd", n = 1, max.taxa = n_taxa, include.extinct = FALSE
  )[[1]]
  tree$edge.length <- tree$edge.length / max(phytools::branching.times(tree))
  cont_trait <- simulate_continuous(tree, sigma)
  disc_trait <- geiger::sim.char(
    tree,
    par = matrix(c(-rate, rate, rate, -rate), 2),
    model = "discrete",
    root = sample(1:2, 1)
  )[, 1]
  list(tree = tree, cont_trait = cont_trait, disc_trait = disc_trait)
}

run_taxa_bidirectional <- function(n_trees = 20,
                                   taxa_grid = seq(20, 200, length.out = 5),
                                   rate = 0.6,
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
      dataset <- simulate_bidirectional_dataset(n_taxa, sigma, rate)
      data_frame <- data.frame(
        taxon = dataset$tree$tip.label,
        cont = dataset$cont_trait,
        disc = dataset$disc_trait
      )
      anc_result <- AncCond(
        tree = dataset$tree,
        data = data_frame,
        nsim = nsim,
        iter = iter,
        mat = c(0, 2, 1, 0),
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
  taxa_results <- run_taxa_bidirectional()
  saveRDS(taxa_results, file = file.path(project_root, "results", "taxa_bidirectional_results.rds"))
}
