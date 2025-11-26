#!/usr/bin/env Rscript
# Unidirectional scaling simulation (parallel version using parLapply on Windows)

library(phytools)
library(geiger)
library(parallel)
library(evobiR)

set.seed(321)
source("~/GitHub/evobir/R/AncCond.R", keep.source = TRUE)

## ---- parameters (were formerly function args) ----
n_trees         <- 10
n_taxa          <- 200
scaling_factors <- c(1)
forward_rate    <- 0.2
reverse_rate    <- 1e-04
cont_sigma      <- 0.2
iter            <- 200
verbose         <- TRUE

## ---- rate matrix ----
q_matrix <- matrix(c(-forward_rate, forward_rate,
                     reverse_rate, -reverse_rate),
                   nrow = 2, byrow = TRUE)

## ---- start cluster (Windows: PSOCK) ----
n_workers <- max(1L, parallel::detectCores() - 1L)
cl <- parallel::makeCluster(n_workers)
on.exit(parallel::stopCluster(cl), add = TRUE)

parallel::clusterEvalQ(cl, {
  library(phytools); library(geiger)
  source("~/GitHub/evobir/R/AncCond.R", local = TRUE)
  NULL
})

# export objects workers need
parallel::clusterExport(
  cl,
  varlist = c("q_matrix", "n_taxa", "cont_sigma", "iter"),
  envir = environment()
)

parallel::clusterSetRNGStream(cl, 321) # reproducible RNG

## ---- container for results ----
results <- vector("list", length(scaling_factors))
names(results) <- paste0("scale_", scaling_factors)

## ---- loop over scaling factors; parLapply across trees ----
for (sf_idx in seq_along(scaling_factors)) {
  sf <- scaling_factors[sf_idx]
  if (verbose) message(sprintf("Scale %s using %d workers", sf, n_workers))
  
  # >>> Option 1: export the current sf to workers each iteration <<<
  parallel::clusterExport(cl, varlist = "sf", envir = environment())
  
  uni_tree_results <- parallel::parLapply(cl, seq_len(n_trees), function(t) {
    ## simulate tree + continuous trait
    tree <- sim.bdtree(b = 3, d = 1, stop = c("taxa"), n = n_taxa, extinct = FALSE)
    tree$edge.length <- tree$edge.length / max(ape::branching.times(tree))
    cont_trait <- phytools::fastBM(tree, sig2 = cont_sigma)
    names(cont_trait) <- tree$tip.label
    
    ## compute branch means (ASR on continuous trait)
    cont_asr <- phytools::anc.ML(tree, cont_trait, model = "BM")
    bmeans <- numeric(nrow(tree$edge))
    for (i in seq_len(nrow(tree$edge))) {
      parent <- tree$edge[i, 1]
      child  <- tree$edge[i, 2]
      parent_val <- if (parent <= length(tree$tip.label)) cont_trait[parent] else cont_asr$ace[as.character(parent)]
      child_val  <- if (child  <= length(tree$tip.label)) cont_trait[child]  else cont_asr$ace[as.character(child)]
      bmeans[i] <- (parent_val + child_val) / 2
    }
    
    ## scale edges by branch-mean quantiles
    scaled_tree <- tree
    qs <- stats::quantile(bmeans, probs = c(0.25, 0.75))
    scaled_tree$edge.length[bmeans < qs[[1]]] <- scaled_tree$edge.length[bmeans < qs[[1]]] / sf
    scaled_tree$edge.length[bmeans > qs[[2]]] <- scaled_tree$edge.length[bmeans > qs[[2]]] * sf
    
    ## simulate discrete until not near-fixed
    good_sim <- FALSE
    attempts <- 0L
    max_attempts <- 1000L
    disc_trait <- NULL
    
    while (!good_sim && attempts < max_attempts) {
      attempts <- attempts + 1L
      sim_out <- geiger::sim.char(
        scaled_tree,
        par = q_matrix,
        model = "discrete",
        root = 1,
        nsim = 1
      )
      disc_trait <- sim_out[, , 1]
      disc_trait <- as.integer(disc_trait)
      freq <- mean(disc_trait == min(disc_trait), na.rm = TRUE)
      good_sim <- is.finite(freq) && freq > 0.05 && freq < 0.95
    }
    
    if (!good_sim) return(c(p12 = NA_real_, p21 = NA_real_))
    
    ## assemble data frame for AncCond
    df <- data.frame(
      taxon = scaled_tree$tip.label,
      cont  = cont_trait,
      disc  = disc_trait
    )
    
    ## AncCond (unidirectional)
    anc <- AncCond(
      tree = scaled_tree,
      data = df,
      mc   = as.integer(iter),
      drop.state = 2,
      mat  = c(0, 0, 1, 0),
      message = FALSE
    )
    
    ## return p-values in a consistent shape
    if (is.list(anc) && "pval1->2" %in% names(anc) && "pval2->1" %in% names(anc)) {
      c(p12 = anc[["pval1->2"]], p21 = anc[["pval2->1"]])
    } else if (is.list(anc) && "pval" %in% names(anc)) {
      c(p12 = anc[["pval"]], p21 = NA_real_)
    } else {
      c(p12 = NA_real_, p21 = NA_real_)
    }
  })
  
  results[[sf_idx]] <- do.call(rbind, uni_tree_results)
  colnames(results[[sf_idx]]) <- c("p12", "p21")
}

## ---- save ----
save(results, file='../../results/scaling_unidirectional_results.RData')

