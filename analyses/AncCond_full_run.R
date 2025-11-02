# ------------------------------------------------------------
# Run AncCond on selected scaling factors using sim_tree_data
# ------------------------------------------------------------

# Packages
library(evobiR)
library(ape)
library(phytools)

# ---- Path to your saved data ----
# ---- Use the working directory for ALL I/O ----
base_dir <- getwd()
if (!file.exists(file.path(base_dir, "sim_tree_data.RData"))) {
  stop("sim_tree_data.RData not found in working directory: ", base_dir)
}
load(file.path(base_dir, "sim_tree_data.RData"))  # loads 'sim_tree_data'

if (!exists("sim_tree_data")) stop("sim_tree_data not found after load().")
if (!all(c("bd_trees","sim_results","sf_trees") %in% names(sim_tree_data))) {
  stop("sim_tree_data must contain bd_trees, sim_results, and sf_trees.")
}

bd_trees    <- sim_tree_data$bd_trees
sim_results <- sim_tree_data$sim_results
sf_trees    <- sim_tree_data$sf_trees

# ---- Choose exactly which scaling factors to run ----
# e.g., c(1, 3, 5, 10) or just 1:10
scaling_factors <- c(1, 2, 8)

# ---- AncCond settings ----
iter    <- 1000      # mc reps inside AncCond
verbose <- TRUE

# If you use a custom AncCond, you can source it (optional):
# source("~/GitHub/evobir/R/AncCond.R", keep.source = TRUE)

# We’ll store results per scenario, per size, per tree, per sf
scenario_names <- c("uni", "bi")

# Containers for outputs
unidirectional_results <- setNames(vector("list", length(sim_results[["uni"]])),
                                   names(sim_results[["uni"]]))
bidirectional_results  <- setNames(vector("list", length(sim_results[["bi"]])),
                                   names(sim_results[["bi"]]))

# ============================================================
# Helper to ensure requested sf exist in data (names are "sf1".."sf10")
# ============================================================
all_sf_names <- paste0("sf", scaling_factors)

# ============================================================
# Unidirectional
# ============================================================
if (!("uni" %in% names(sim_results))) stop("No 'uni' scenario found in sim_results.")

for (si in seq_along(unidirectional_results)) {
  size_name <- names(unidirectional_results)[si]
  
  # Per-size lists from sim_tree_data
  uni_runs <- sim_results[["uni"]][[size_name]]           # traits per original tree
  sf_set   <- sf_trees   [["uni"]][[size_name]]           # scaled trees per original tree
  
  if (verbose) message(sprintf("[UNI] %s tips: %d trees", size_name, length(uni_runs)))
  
  per_size_out <- vector("list", length(uni_runs))
  
  for (ti in seq_along(uni_runs)) {
    rec <- uni_runs[[ti]]
    tr_scaled_set <- sf_set[[ti]]        # list: sf1..sf10 trees
    
    # Prepare per-tree container indexed by the sf names we want
    available_sfs <- intersect(names(tr_scaled_set), all_sf_names)
    if (length(available_sfs) == 0L) {
      stop(sprintf("[UNI] No requested scaling factors found for n=%s, tree %d", size_name, ti))
    }
    per_tree_out <- setNames(vector("list", length(available_sfs)), available_sfs)
    
    # Pull continuous once (we’ll align to each scaled tree's tip labels)
    cont_trait <- rec$cont_trait
    
    # Discrete traits per sf (from your new pipeline)
    disc_by_sf <- rec$disc_by_sf
    if (is.null(disc_by_sf)) stop("[UNI] disc_by_sf missing; run the updated sim_tree_data builder first.")
    
    for (sfi in seq_along(available_sfs)) {
      sf_name <- available_sfs[sfi]
      tr      <- tr_scaled_set[[sf_name]]          # scaled phylo
      disc    <- disc_by_sf[[sf_name]]             # "1"/"2" per tip for this scaled tree
      
      # Align vectors to tree tip order
      x <- cont_trait[tr$tip.label]                # numeric
      y <- as.integer(disc[tr$tip.label])          # 1/2 integers to match your AncCond call
      # (drop.state = 2 expects state "2" to be droppable)
      
      # Construct data frame in the format your example uses
      df <- data.frame(
        taxon = tr$tip.label,
        cont  = x,
        disc  = y
      )
      
      # AncCond
      ac_res <- tryCatch(
        evobiR::AncCond(
          tree = tr,
          data = df,
          mc   = as.integer(iter),
          drop.state = 2,
          mat  = c(0, 0, 1, 0),
          message = FALSE
        ),
        error = function(e) e
      )
      
      per_tree_out[[sfi]] <- list(
        sf_name   = sf_name,
        anccond   = ac_res,
        n_tips    = length(tr$tip.label),
        cont_trait= x,
        disc_trait= y,
        Q         = rec$Q,              # store the uni matrix used in sim
        root_info = rec$root_by_sf[[sf_name]],
        iter      = iter
      )
    }
    
    per_size_out[[ti]] <- per_tree_out
  }
  
  unidirectional_results[[si]] <- per_size_out
}

save(unidirectional_results,
     file = file.path(base_dir, "unidirectional_results.RData"))
if (verbose) message("Saved unidirectional_results.RData")

# ============================================================
# Bidirectional
# ============================================================
if (!("bi" %in% names(sim_results))) stop("No 'bi' scenario found in sim_results.")

for (si in seq_along(bidirectional_results)) {
  size_name <- names(bidirectional_results)[si]
  
  bi_runs <- sim_results[["bi"]][[size_name]]
  sf_set  <- sf_trees   [["bi"]][[size_name]]
  
  if (verbose) message(sprintf("[BI ] %s tips: %d trees", size_name, length(bi_runs)))
  
  per_size_out <- vector("list", length(bi_runs))
  
  for (ti in seq_along(bi_runs)) {
    rec <- bi_runs[[ti]]
    tr_scaled_set <- sf_set[[ti]]
    
    available_sfs <- intersect(names(tr_scaled_set), all_sf_names)
    if (length(available_sfs) == 0L) {
      stop(sprintf("[BI ] No requested scaling factors found for n=%s, tree %d", size_name, ti))
    }
    per_tree_out <- setNames(vector("list", length(available_sfs)), available_sfs)
    
    cont_trait <- rec$cont_trait
    disc_by_sf <- rec$disc_by_sf
    if (is.null(disc_by_sf)) stop("[BI ] disc_by_sf missing; run the updated sim_tree_data builder first.")
    
    for (sfi in seq_along(available_sfs)) {
      sf_name <- available_sfs[sfi]
      tr      <- tr_scaled_set[[sf_name]]
      disc    <- disc_by_sf[[sf_name]]
      
      x <- cont_trait[tr$tip.label]
      y <- as.integer(disc[tr$tip.label])   # 1/2 integers
      
      df <- data.frame(
        taxon = tr$tip.label,
        cont  = x,
        disc  = y
      )
      
      ac_res <- tryCatch(
        evobiR::AncCond(
          tree = tr,
          data = df,
          mc   = as.integer(iter),
          drop.state = 2,
          mat  = c(0, 0, 1, 0),
          message = FALSE
        ),
        error = function(e) e
      )
      
      per_tree_out[[sfi]] <- list(
        sf_name     = sf_name,
        anccond     = ac_res,
        n_tips      = length(tr$tip.label),
        cont_trait  = x,
        disc_trait  = y,
        rate_matrix = rec$rate_matrix,         # store the bi matrix used in sim
        root_info   = rec$root_by_sf[[sf_name]],
        iter        = iter
      )
    }
    
    per_size_out[[ti]] <- per_tree_out
  }
  
  bidirectional_results[[si]] <- per_size_out
}

save(bidirectional_results,
     file = file.path(base_dir, "bidirectional_results.RData"))
if (verbose) message("Saved bidirectional_results.RData")

# ------------------------------------------------------------
# Optional quick peek:
str(unidirectional_results[["25"]][[1]][["sf1"]]$anccond)
str(bidirectional_results [["100"]][[10]][["sf3"]]$anccond)