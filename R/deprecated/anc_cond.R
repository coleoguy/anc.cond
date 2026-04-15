#' AncCond core functions
#'
#' This module contains the implementation of the ancestral condition test
#' and helper utilities used across the simulation and empirical analyses.
#' The code consolidates the functionality that was previously scattered
#' across several scripts into a single, documented location.

#' Run the ancestral condition test.
#'
#' @param tree A phylo object.
#' @param data A data.frame with three columns: tip labels, continuous trait,
#'   and discrete trait (coded as 1/2).
#' @param drop.state Optional discrete state (1 or 2) to drop from the
#'   continuous trait reconstruction.
#' @param mat Vector representation of the 2x2 transition matrix supplied to
#'   phytools::make.simmap.
#' @param pi Root state prior passed to phytools::make.simmap.
#' @param n.tails Number of tails (1 or 2) for the p-value calculation.
#' @param nsim Number of stochastic maps to draw.
#' @param iter Number of null datasets generated per stochastic map.
#' @param message Logical; print progress messages and summaries when TRUE.
#'   Ignored when ncores > 1.
#' @param ncores Number of cores to use for parallel processing (default 1).
#'   Uses parallel::parLapply (PSOCK clusters) which works on Windows, Mac,
#'   and Linux.
#'
#' @return An object of class "AncCond" containing observed and null
#'   distributions, p-values, and the mean number of transitions.
AncCond <- function(tree,
                    data,
                    drop.state = NULL,
                    mat = c(0, 2, 1, 0),
                    pi = "estimated",
                    n.tails = 1,
                    nsim = 100,
                    iter = 100,
                    message = FALSE,
                    ncores = 1L) {
  checkPackages()
  InputTesting(tree, data, drop.state, mat, pi, n.tails, nsim, iter)

  unpacked <- UnpackData(data, drop.state)
  dt.vec <- unpacked[[1]]
  ct.vec <- unpacked[[2]]

  if (isTRUE(message)) {
    cat("Estimating ancestral states for the continuous trait\n")
  }
  anc.states.cont.trait <- phytools::anc.ML(tree, ct.vec, model = "BM")

  if (isTRUE(message)) {
    cat("Simulating stochastic mappings:\n")
  }
  anc.state.dt <- phytools::make.simmap(
    tree,
    dt.vec,
    model = matrix(mat, 2),
    nsim = nsim,
    pi = pi,
    Q = "mcmc",
    message = message
  )

  # Worker function for one stochastic map
  process_one_map <- function(j) {
    current.map <- anc.state.dt[[j]]
    obs <- exctractAncestral(
      current.map = current.map,
      anc.states.cont.trait = anc.states.cont.trait,
      count = TRUE
    )
    ntrans <- obs$ntrans
    obs$ntrans <- NULL

    null <- CreateNull(
      tree = tree,
      iter = iter,
      current.map = current.map,
      anc.states.cont.trait = anc.states.cont.trait,
      message = FALSE,
      j = j,
      nsim = nsim
    )
    list(obs = obs, null = null, ntrans = ntrans)
  }

  if (ncores > 1L) {
    cl <- parallel::makeCluster(ncores)
    on.exit(parallel::stopCluster(cl), add = TRUE)
    parallel::clusterExport(cl, varlist = c(
      "anc.state.dt", "anc.states.cont.trait", "tree", "iter", "nsim",
      "exctractAncestral", "CreateNull", ".count_transitions_fast"
    ), envir = environment())
    map_results <- parallel::parLapply(cl, seq_len(nsim), process_one_map)
  } else {
    map_results <- lapply(seq_len(nsim), function(j) {
      if (isTRUE(message)) cat("\014Analyzing map:", j, "of", nsim, "\n")
      process_one_map(j)
    })
  }

  observed.anc.cond <- lapply(map_results, `[[`, "obs")
  null.anc.cond <- lapply(map_results, `[[`, "null")
  meantrans <- Reduce(`+`, lapply(map_results, `[[`, "ntrans")) / nsim

  obs.dist <- ProcessObserved(observed.anc.cond)
  null.dist <- ProcessNull(null.anc.cond, iter)
  results <- list(observed = obs.dist, null = null.dist)
  pvals <- CalcPVal(results, n.tails)

  results$pvals <- pvals
  results$`mean n trans` <- meantrans
  class(results) <- "AncCond"

  if (isTRUE(message)) summary(results)
  results
}

# Fast transition counter -- avoids calling summary.simmap which is very slow.
# Counts total number of state changes across all edges by checking how many
# segments each edge has (transitions = segments - 1).
.count_transitions_fast <- function(simmap_obj) {
  sum(vapply(simmap_obj$maps, function(m) length(m) - 1L, integer(1)))
}

# Helper to ensure required packages are installed.
checkPackages <- function() {
  required <- c("phytools", "diversitree")
  missing_pkgs <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing_pkgs) > 0) {
    stop(
      "The following packages are required but not installed: ",
      paste(missing_pkgs, collapse = ", "),
      call. = FALSE
    )
  }
}

InputTesting <- function(tree,
                         data,
                         drop.state,
                         mat,
                         pi,
                         n.tails,
                         nsim,
                         iter) {
  if (!inherits(tree, "phylo")) {
    stop("tree must be class phylo")
  }
  if (!is.data.frame(data) || ncol(data) != 3) {
    stop("data should be a data.frame with three columns (tip labels, continuous trait, discrete trait)")
  }
  if (!is.null(drop.state) && !drop.state %in% c(1, 2)) {
    stop("drop.state must be NULL, 1, or 2")
  }
  valid_mats <- rbind(c(0, 0, 1, 0), c(0, 1, 1, 0), c(0, 2, 1, 0))
  if (!any(apply(valid_mats, 1, identical, y = mat))) {
    stop("mat must be one of c(0,0,1,0), c(0,1,1,0), or c(0,2,1,0)")
  }
  if (!pi[1] %in% c("equal", "estimated")) {
    if (!is.numeric(pi)) stop("pi must be 'equal', 'estimated', or a numeric vector of length 2")
    if (length(pi) != 2 || abs(sum(pi) - 1) > .Machine$double.eps^0.5) {
      stop("numeric pi must have length 2 and sum to 1")
    }
  }
  if (!n.tails %in% c(1, 2)) {
    stop("n.tails should be 1 or 2")
  }
  if (!is.numeric(nsim) || nsim < 1) {
    stop("nsim should be a positive numeric value")
  }
  if (!is.numeric(iter) || iter < 100) {
    stop("iter should be numeric and at least 100 to ensure stable null distributions")
  }
}

UnpackData <- function(data, drop.state) {
  dt.vec <- data[[3]]
  names(dt.vec) <- data[[1]]

  if (!is.null(drop.state)) {
    ct.data <- data[data[[3]] != drop.state, , drop = FALSE]
  } else {
    ct.data <- data
  }
  ct.vec <- as.numeric(ct.data[[2]])
  names(ct.vec) <- ct.data[[1]]

  if (anyNA(ct.vec) || anyNA(dt.vec)) {
    stop(
      "Missing trait data detected. Remove taxa with incomplete data from the tree before running AncCond."
    )
  }
  list(dt.vec, ct.vec)
}

CreateNull <- function(tree,
                       iter,
                       current.map,
                       anc.states.cont.trait,
                       message,
                       j,
                       nsim) {
  current.Q <- current.map$Q
  if (any(current.Q == 0)) {
    zero_idx <- which(current.Q == 0, arr.ind = TRUE)
    for (k in seq_len(nrow(zero_idx))) {
      r <- zero_idx[k, 1]
      cc <- zero_idx[k, 2]
      current.Q[r, cc] <- if (r == cc) -1e-25 else 1e-25
    }
    if (isTRUE(message)) {
      cat("\nSetting zero transition rates to +/- 1e-25 for simulation stability.\n")
    }
  }
  root.state <- c(`1` = 0, `2` = 0)
  root.state[names(root.state) == names(current.map$maps[[1]])[1]] <- 1

  # Cache empirical transition count outside the loop
  current_trans <- .count_transitions_fast(current.map)
  if (current_trans > 5) {
    lo <- 0.8 * current_trans
    hi <- 1.2 * current_trans
  } else {
    lo <- current_trans - 1
    hi <- current_trans + 1
  }

  nulldist <- vector("list", iter)
  batch_size <- 50L
  filled <- 0L
  total_attempts <- 0L
  max_attempts <- iter * 10000L

  while (filled < iter && total_attempts < max_attempts) {
    n_need <- min(batch_size, (iter - filled) * 20L)
    sims <- phytools::sim.history(
      tree = tree,
      Q = current.Q,
      nsim = n_need,
      message = FALSE,
      anc = root.state
    )
    if (!inherits(sims, "multiSimmap")) sims <- list(sims)

    for (s in seq_along(sims)) {
      total_attempts <- total_attempts + 1L
      sim_trans <- .count_transitions_fast(sims[[s]])
      if (sim_trans >= lo && sim_trans <= hi) {
        filled <- filled + 1L
        nulldist[[filled]] <- exctractAncestral(
          current.map = sims[[s]],
          anc.states.cont.trait = anc.states.cont.trait,
          count = FALSE
        )
        if (isTRUE(message)) {
          cat(
            "\014Analyzing map:", j, "of", nsim, "\n",
            "Number of transitions:\n",
            " Empirical map:\n", current_trans,
            "\n Null simulation:\n", sim_trans,
            "\n"
          )
        }
        if (filled >= iter) break
      }
    }
  }

  # Fill any remaining slots with NA if we timed out
  if (filled < iter) {
    if (isTRUE(message)) {
      warning("Unable to simulate a null with similar behavior to the observed map.")
    }
    for (n in (filled + 1L):iter) {
      nulldist[[n]] <- list(`12` = NA_real_, `21` = NA_real_)
    }
  }
  nulldist
}

exctractAncestral <- function(current.map,
                              anc.states.cont.trait,
                              count = FALSE) {
  me <- current.map$mapped.edge
  col1 <- if ("1" %in% colnames(me)) me[, "1", drop = TRUE] else rep(0, nrow(me))
  col2 <- if ("2" %in% colnames(me)) me[, "2", drop = TRUE] else rep(0, nrow(me))
  ss_nodes <- (col1 > 0) & (col2 > 0)
  
  # define wanted_nodes before using it
  rn <- rownames(me); if (is.null(rn)) rn <- as.character(seq_len(nrow(me)))
  wanted_nodes <- rn[ss_nodes]
  
  trans.maps <- current.map$maps[ss_nodes]
  wanted_nodes <- gsub(",.*", "", wanted_nodes)

  # Vectorized: get first state name for each transition edge
  first_states <- vapply(trans.maps, function(m) names(m)[1L], character(1))

  producing.nodes12 <- unique(wanted_nodes[first_states == "1"])
  producing.nodes21 <- unique(wanted_nodes[first_states == "2"])

  ace <- anc.states.cont.trait$ace
  ace_names <- names(ace)
  res <- list(
    `12` = ace[ace_names %in% producing.nodes12],
    `21` = ace[ace_names %in% producing.nodes21]
  )
  if (isTRUE(count)) {
    res$ntrans <- c(`12` = length(producing.nodes12), `21` = length(producing.nodes21))
  }
  res
}


ProcessObserved <- function(observed.anc.cond) {
  vals12 <- vapply(observed.anc.cond, function(x) mean(x$`12`, na.rm = TRUE), numeric(1))
  vals21 <- vapply(observed.anc.cond, function(x) mean(x$`21`, na.rm = TRUE), numeric(1))
  c(`12` = mean(vals12, na.rm = TRUE), `21` = mean(vals21, na.rm = TRUE))
}

ProcessNull <- function(null.anc.cond, iter) {
  n_maps <- length(null.anc.cond)
  # Vectorized: build matrices using vapply instead of nested for-loops
  mat12 <- vapply(null.anc.cond, function(map_nulls) {
    vapply(map_nulls, function(x) mean(x$`12`, na.rm = TRUE), numeric(1))
  }, numeric(iter))  # each column = one map
  mat21 <- vapply(null.anc.cond, function(map_nulls) {
    vapply(map_nulls, function(x) mean(x$`21`, na.rm = TRUE), numeric(1))
  }, numeric(iter))
  list(`12` = rowMeans(mat12, na.rm = TRUE), `21` = rowMeans(mat21, na.rm = TRUE))
}

CalcPVal <- function(results, n.tails) {
  calc_tail <- function(null_vals, observed) {
    null_vals <- null_vals[!is.na(null_vals)]
    if (length(null_vals) == 0 || is.na(observed)) return(NA_real_)
    bigger <- mean(null_vals >= observed)
    smaller <- mean(null_vals < observed)
    p <- min(bigger, smaller)
    if (n.tails == 2) {
      p <- min(1, 2 * p)
    }
    p
  }
  c(
    `12` = calc_tail(results$null$`12`, results$observed["12"]),
    `21` = calc_tail(results$null$`21`, results$observed["21"])
  )
}

summary.AncCond <- function(results) {
  cat("\nMean value for the continuous trait at 1 -> 2 transitions:",
      round(results$observed["12"], 4), "\n")
  cat("Mean value for the continuous trait at 2 -> 1 transitions:",
      round(results$observed["21"], 4), "\n\n")
  cat("Mean number of 1 -> 2 transitions:",
      round(results$`mean n trans`["12"], 4), "\n")
  cat("Mean number of 2 -> 1 transitions:",
      round(results$`mean n trans`["21"], 4), "\n\n")
  cat("Mean of null dist 1 -> 2:",
      round(mean(results$null$`12`, na.rm = TRUE), 4), "\n")
  cat("Mean of null dist 2 -> 1:",
      round(mean(results$null$`21`, na.rm = TRUE), 4), "\n\n")
  cat("SD of null dist 1 -> 2:",
      round(sd(results$null$`12`, na.rm = TRUE), 4), "\n")
  cat("SD of null dist 2 -> 1:",
      round(sd(results$null$`21`, na.rm = TRUE), 4), "\n\n")
  cat("pvalue 1 -> 2:", round(results$pvals["12"], 4), "\n")
  cat("pvalue 2 -> 1:", round(results$pvals["21"], 4), "\n\n")
  if (is.na(results$pvals["12"])) {
    cat("NA values indicate that no 1 -> 2 transitions occurred.\n\n")
  }
  if (is.na(results$pvals["21"])) {
    cat("NA values indicate that no 2 -> 1 transitions occurred.\n\n")
  }
}

plot.AncCond <- function(results) {
  op <- par(no.readonly = TRUE)
  on.exit(par(op))
  par(mfrow = c(1, 2))
  if (!is.na(results$pvals["12"])) {
    plot(density(results$null$`12`, na.rm = TRUE),
         main = "1 -> 2",
         xlim = range(c(results$null$`12`, results$observed["12"]), na.rm = TRUE),
         xlab = "Ancestral condition",
         ylab = "Frequency")
    abline(v = results$observed["12"], col = "red", lwd = 2)
    legend("topright", legend = c("Observed", "Null"), col = c("red", "black"), lwd = 2)
  }
  if (!is.na(results$pvals["21"])) {
    plot(density(results$null$`21`, na.rm = TRUE),
         main = "2 -> 1",
         xlim = range(c(results$null$`21`, results$observed["21"]), na.rm = TRUE),
         xlab = "Ancestral condition",
         ylab = "Frequency")
    abline(v = results$observed["21"], col = "red", lwd = 2)
    legend("topright", legend = c("Observed", "Null"), col = c("red", "black"), lwd = 2)
  }
}
