## ===== AncCond runner (parallel; minimal helpers) =====
library(ape)
library(phytools)
library(evobiR)
library(parallel)

## 1) Load your patched AncCond from the repo and freeze it
source("~/GitHub/evobir/R/AncCond.R", keep.source = TRUE)
AncCond_fixed <- AncCond  # workers will call THIS

## --- user params ---
base_dir        <- "C:/Users/mcconnell.m.meghann/Documents/GitHub/anc.cond/results"
simdata_file    <- file.path(base_dir, "sim_tree_data.RData")
scaling_factors <- 1:10              # sf1..sf10
iter            <- 500               # MC reps inside AncCond
nperm           <- 100               # label permutations per tree
max_perm_retry  <- 50                # retries to avoid single-state permutes
verbose         <- TRUE

## --- load data ---
load(simdata_file)
stopifnot(exists("sim_tree_data"))
stopifnot(all(c("bd_trees","sim_results","sf_trees") %in% names(sim_tree_data)))
sim_results  <- sim_tree_data$sim_results
sf_trees     <- sim_tree_data$sf_trees
all_sf_names <- paste0("sf", scaling_factors)

## ====== parallel setup ======
n_workers <- max(1L, detectCores() - 1L)
cl <- makeCluster(n_workers, type = "PSOCK")
on.exit(stopCluster(cl), add = TRUE)

clusterEvalQ(cl, {
  suppressPackageStartupMessages({ library(evobiR); library(ape); library(phytools) })
  if (requireNamespace("RhpcBLASctl", quietly = TRUE)) {
    RhpcBLASctl::blas_set_num_threads(1L)
    RhpcBLASctl::omp_set_num_threads(1L)
  } else {
    Sys.setenv(OMP_NUM_THREADS = "1", MKL_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1")
  }
  NULL
})
clusterSetRNGStream(cl, 321)

## ====== SINGLE worker function (ti FIRST; no other helpers) ======
worker_fun <- function(ti, scen, size_name,
                       sim_results, sf_trees, all_sf_names,
                       iter, nperm, max_perm_retry) {
  bi <- identical(scen, "bi")
  rec <- sim_results[[scen]][[size_name]][[ti]]
  tr_scaled_set <- sf_trees[[scen]][[size_name]][[ti]]
  want <- intersect(names(tr_scaled_set), all_sf_names)
  if (!length(want)) return(NULL)
  
  out <- setNames(vector("list", length(want)), want)
  
  for (sf_name in want) {
    tr <- tr_scaled_set[[sf_name]]
    if (!inherits(tr, "phylo")) { out[[sf_name]] <- list(error="tree_not_phylo"); next }
    
    x   <- as.numeric(rec$cont_trait[tr$tip.label])
    y12 <- as.integer(rec$disc_by_sf[[sf_name]][tr$tip.label])   # 1/2 integers
    
    df3 <- data.frame(taxon = tr$tip.label, cont = x, disc = y12, stringsAsFactors = FALSE)
    
    ## Observed AncCond
    res_obs <- tryCatch(
      AncCond_fixed(
        tree       = tr,
        data       = df3,
        mc         = as.integer(iter),
        drop.state = if (bi) NULL else 2,
        mat        = if (bi) c(0,2,1,0) else c(0,0,1,0),
        message    = FALSE
      ),
      error = function(e) e
    )
    
    ## Extract observed p-values (handles bi/uni)
    p12 <- p21 <- NA_real_
    if (!inherits(res_obs, "error") && is.list(res_obs)) {
      if ("pval1->2" %in% names(res_obs)) p12 <- as.numeric(res_obs[["pval1->2"]])
      if ("pval2->1" %in% names(res_obs)) p21 <- as.numeric(res_obs[["pval2->1"]])
      if (is.na(p12) && "pval" %in% names(res_obs)) p12 <- as.numeric(res_obs[["pval"]])
    }
    
    ## Nulls via label permutations (skip if nperm == 0)
    p12_null <- p21_null <- NULL
    if (nperm > 0L) {
      p12_null <- rep(NA_real_, nperm)
      p21_null <- rep(NA_real_, nperm)
      for (pi in seq_len(nperm)) {
        tries <- 0L; ok <- FALSE; y_perm <- y12
        while (!ok && tries < max_perm_retry) {
          tries  <- tries + 1L
          y_perm <- sample(y12, replace = FALSE)
          ok <- (length(unique(y_perm)) == 2L)
        }
        if (!ok) next
        dfp3 <- data.frame(taxon = tr$tip.label, cont = x, disc = as.integer(y_perm), stringsAsFactors = FALSE)
        res_null <- tryCatch(
          AncCond_fixed(
            tree       = tr,
            data       = dfp3,
            mc         = as.integer(iter),
            drop.state = if (bi) NULL else 2,
            mat        = if (bi) c(0,2,1,0) else c(0,0,1,0),
            message    = FALSE
          ),
          error = function(e) e
        )
        if (!inherits(res_null, "error") && is.list(res_null)) {
          p12_null[pi] <- if ("pval1->2" %in% names(res_null)) as.numeric(res_null[["pval1->2"]])
          else if ("pval" %in% names(res_null)) as.numeric(res_null[["pval"]]) else NA_real_
          p21_null[pi] <- if ("pval2->1" %in% names(res_null)) as.numeric(res_null[["pval2->1"]]) else NA_real_
        }
      }
    }
    
    v12 <- if (length(p12_null)) p12_null[is.finite(p12_null)] else numeric(0)
    v21 <- if (length(p21_null)) p21_null[is.finite(p21_null)] else numeric(0)
    emp12 <- if (is.finite(p12) && length(v12)) (sum(v12 <= p12) + 1) / (length(v12) + 1) else NA_real_
    emp21 <- if (is.finite(p21) && length(v21)) (sum(v21 <= p21) + 1) / (length(v21) + 1) else NA_real_
    
    out[[sf_name]] <- list(
      sf_name   = sf_name,
      anccond   = res_obs,
      obs_p12   = p12,
      obs_p21   = p21,
      emp_p12   = emp12,
      emp_p21   = emp21,
      null_p12  = p12_null,
      null_p21  = p21_null
    )
  }
  
  out
}

## Export EVERYTHING the worker needs (including worker_fun!)
clusterExport(
  cl,
  varlist = c("AncCond_fixed","sim_results","sf_trees","all_sf_names",
              "iter","nperm","max_perm_retry","worker_fun"),
  envir = environment()
)

## ====== run UNI in parallel (per tree) ======
stopifnot("uni" %in% names(sim_results))
unidirectional_results <- setNames(vector("list", length(sim_results$uni)),
                                   names(sim_results$uni))

for (size_name in names(sim_results$uni)) {
  trees_n <- length(sim_results$uni[[size_name]])
  if (verbose) message(sprintf("[UNI] %s tips: %d trees (workers=%d)", size_name, trees_n, n_workers))
  per_size <- parLapplyLB(
    cl, seq_len(trees_n), worker_fun,
    scen = "uni", size_name = size_name,
    sim_results = sim_results, sf_trees = sf_trees, all_sf_names = all_sf_names,
    iter = iter, nperm = nperm, max_perm_retry = max_perm_retry
  )
  names(per_size) <- as.character(seq_len(trees_n))
  unidirectional_results[[size_name]] <- per_size
}
save(unidirectional_results, file = file.path(base_dir, "unidirectional_results_with_nulls.RData"))
if (verbose) message("Saved unidirectional_results_with_nulls.RData")

## ====== run BI in parallel (per tree) ======
stopifnot("bi" %in% names(sim_results))
bidirectional_results <- setNames(vector("list", length(sim_results$bi)),
                                  names(sim_results$bi))

for (size_name in names(sim_results$bi)) {
  trees_n <- length(sim_results$bi[[size_name]])
  if (verbose) message(sprintf("[BI ] %s tips: %d trees (workers=%d)", size_name, trees_n, n_workers))
  per_size <- parLapplyLB(
    cl, seq_len(trees_n), worker_fun,
    scen = "bi", size_name = size_name,
    sim_results = sim_results, sf_trees = sf_trees, all_sf_names = all_sf_names,
    iter = iter, nperm = nperm, max_perm_retry = max_perm_retry
  )
  names(per_size) <- as.character(seq_len(trees_n))
  bidirectional_results[[size_name]] <- per_size
}
save(bidirectional_results, file = file.path(base_dir, "bidirectional_results_with_nulls.RData"))
message("Done. Results saved to: ", base_dir)