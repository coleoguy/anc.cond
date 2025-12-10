#############################
## Shared setup + helpers  ##
#############################
load("~/GitHub/anc.cond/results/sim.results.RData")

alpha     <- 0.025
taxa_vals <- as.integer(names(result.list$uni))     # 25,50,75,100,200
n_rep     <- length(result.list$uni[[1]])           # 100 trees
n_sf      <- length(result.list$uni[[1]][[1]])      # 10 scaling factors

# Unidirectional: single p-value
get_p_uni <- function(obj) {
  p <- as.numeric(obj$pval)[1]
  if (!is.finite(p)) return(NA_real_)
  p
}

# Bidirectional: return *all* p-values (handles $pval and any "pval*" elements)
get_p_bi_vec <- function(obj) {
  if (is.null(obj)) return(numeric(0))
  
  vals <- numeric(0)
  
  if (!is.null(obj$pval))
    vals <- c(vals, as.numeric(obj$pval))
  
  nm <- names(obj)
  if (!is.null(nm)) {
    idx <- grep("pval", nm)
    if (length(idx))
      vals <- c(vals, as.numeric(unlist(obj[idx])))
  }
  
  vals <- vals[is.finite(vals)]
  vals
}

##############################################
## Block 1: Percent significant (power/FPR) ##
##############################################

out <- data.frame(
  taxa_size = integer(),
  sf        = integer(),
  dir       = character(),
  percent   = numeric(),
  stringsAsFactors = FALSE
)

sig_uni <- function(obj) {
  p <- get_p_uni(obj)
  if (is.na(p)) return(NA)
  p <= alpha
}

sig_bi <- function(obj) {
  vals <- get_p_bi_vec(obj)
  if (!length(vals)) return(NA)
  any(vals <= alpha)
}

for (j in seq_along(taxa_vals)) {
  for (k in seq_len(n_sf)) {
    
    uni_sig <- vapply(
      seq_len(n_rep),
      function(r) sig_uni(result.list$uni[[j]][[r]][[k]]),
      logical(1)
    )
    
    bi_sig <- vapply(
      seq_len(n_rep),
      function(r) sig_bi(result.list$bi[[j]][[r]][[k]]),
      logical(1)
    )
    
    out <- rbind(
      out,
      data.frame(
        taxa_size = taxa_vals[j], sf = k,
        dir = "unidirectional",
        percent = mean(uni_sig, na.rm = TRUE)
      ),
      data.frame(
        taxa_size = taxa_vals[j], sf = k,
        dir = "bidirectional",
        percent = mean(bi_sig,  na.rm = TRUE)
      )
    )
  }
}

## Plot: percent significant
par(mfrow = c(2, 3), mar = c(5, 4, 3, 1))

for (i in seq_along(taxa_vals)) {
  tmp <- subset(out, taxa_size == taxa_vals[i])
  if (nrow(tmp) == 0) next
  
  x_sf <- sort(unique(tmp$sf))
  
  y_bi  <- tmp$percent[tmp$dir == "bidirectional" ][
    order(tmp$sf[tmp$dir == "bidirectional" ])
  ]
  y_uni <- tmp$percent[tmp$dir == "unidirectional"][
    order(tmp$sf[tmp$dir == "unidirectional"])
  ]
  
  ymax <- max(tmp$percent, 0.8, na.rm = TRUE)
  if (!is.finite(ymax)) ymax <- 1
  
  plot(x_sf, y_bi,
       type = "b", pch = 16, col = "blue",
       ylim = c(0, ymax),
       xlab = "scaling factor",
       ylab = "percent significant",
       main = paste(taxa_vals[i], "taxa"))
  
  points(x_sf, y_uni, pch = 16, col = "red")
  lines(x_sf, y_uni, col = "red")
  
  points(1, y_bi [x_sf == 1], pch = 0, cex = 1.3)  # bi FPR
  points(1, y_uni[x_sf == 1], pch = 4, cex = 1.3)  # uni FPR
  
  abline(h = alpha, lty = 2)
  
  legend("topleft",
         legend = c("bidirectional power", "unidirectional power",
                    "bidirectional false positive", "unidirectional false positive"),
         pch    = c(16, 16, 0, 4),
         col    = c("blue", "red", "black", "black"),
         bty    = "n", cex = 0.8)
}