#### ----- DIMENSIONS -----
#How many of each we have in each list
n_taxa <- length(unidirectional_results)            # 5
n_rep  <- length(unidirectional_results[[1]])       # 100
n_sf   <- length(unidirectional_results[[1]][[1]])  # 10 scaling factors

#### ----- STORAGE DATA FRAME (long format) -----
#Creates an empty data frame to later store the results
out <- data.frame(
  taxa_size = integer(),   # taxa scenario 1-5
  sf        = integer(),   # scaling factor index 1..10
  dir       = character(), # "uni" or "bi"
  percent   = numeric(),   # proportion of reps significant
  stringsAsFactors = FALSE
)

#### ----- FUNCTION TO EXTRACT MIN P-VALUE -----
get_p <- function(obj) {
  min(unlist(obj$anccond[4:5])) / 100   # two-tailed p
}

#### ----- LOOP: taxa × scaling factor × direction -----
for (i in 1:n_taxa) {
  for (k in 1:n_sf) {
    
    ## --- UNIDIRECTIONAL: vector of 100 p-values ---
    p_uni <- sapply(1:n_rep, function(j)
      get_p(unidirectional_results[[i]][[j]][[k]])
    )
    
    ## --- BIDIRECTIONAL: vector of 100 p-values ---
    p_bi <- sapply(1:n_rep, function(j)
      get_p(bidirectional_results[[i]][[j]][[k]])
    )
    
    ## --- percent/proportion significant (p <= 0.025) ---
    # Out of 100 replicates, what fraction had p ≤ 0.025?
    uni_pct <- mean(p_uni <= 0.025, na.rm = TRUE)
    bi_pct  <- mean(p_bi  <= 0.025, na.rm = TRUE)
    
    ## --- append to table ---
    out <- rbind(out,
                 data.frame(taxa_size = i, sf = k,
                            dir = "unidirectional", percent = uni_pct),
                 data.frame(taxa_size = i, sf = k,
                            dir = "bidirectional",  percent = bi_pct))
  }
}

#PLOT
#2x3 panel layout for 5 taxa sizes
par(mfrow = c(2, 3), mar = c(5, 4, 3, 1))

taxa_labels <- c("25 taxa", "50 taxa", "75 taxa", "100 taxa", "200 taxa")

for (i in 1:n_taxa) {
  tmp <- subset(out, taxa_size == i)
  
  # convenience vectors
  x_sf  <- sort(unique(tmp$sf))
  y_bi  <- tmp$percent[tmp$dir == "bidirectional"][order(tmp$sf[tmp$dir == "bidirectional"])]
  y_uni <- tmp$percent[tmp$dir == "unidirectional"][order(tmp$sf[tmp$dir == "unidirectional"])]
  
  # base plot with bidirectional line
  plot(x_sf, y_bi,
       type = "b", pch = 16, col = "blue",
       ylim = c(0, max(tmp$percent, 0.8)),
       xlab = "scaling factor",
       ylab = "percent significant",
       main = taxa_labels[i])
  
  # add unidirectional line
  points(x_sf, y_uni, pch = 16, col = "red")
  lines(x_sf, y_uni, col = "red")
  
  # mark sf = 1 as "false positive" points
  fp_bi  <- y_bi[x_sf == 1]
  fp_uni <- y_uni[x_sf == 1]
  points(1, fp_bi,  pch = 0, cex = 1.3)  # square
  points(1, fp_uni, pch = 4, cex = 1.3)  # X
  
  # dashed horizontal line at alpha
  abline(h = 0.025, lty = 2)
  
  legend("topleft",
         legend = c("bidirectional power", "unidirectional power",
                    "bidirectional false positive", "unidirectional false positive"),
         pch    = c(16, 16, 0, 4),
         col    = c("blue", "red", "black", "black"),
         bty    = "n", cex = 0.8)
}

#PLOT
# optional: 2x3 panel layout for 5 taxa sizes
par(mfrow = c(2, 3), mar = c(5, 4, 3, 1))

taxa_labels <- c("25 taxa", "50 taxa", "75 taxa", "100 taxa", "200 taxa")

for (i in 1:n_taxa) {
  tmp <- subset(out, taxa_size == i)
  
  # convenience vectors
  x_sf  <- sort(unique(tmp$sf))
  y_bi  <- tmp$percent[tmp$dir == "bidirectional"][order(tmp$sf[tmp$dir == "bidirectional"])]
  y_uni <- tmp$percent[tmp$dir == "unidirectional"][order(tmp$sf[tmp$dir == "unidirectional"])]
  
  # base plot with bidirectional line
  plot(x_sf, y_bi,
       type = "b", pch = 16, col = "blue",
       ylim = c(0, max(tmp$percent, 0.8)),
       xlab = "scaling factor",
       ylab = "percent significant",
       main = taxa_labels[i])
  
  # add unidirectional line
  points(x_sf, y_uni, pch = 16, col = "red")
  lines(x_sf, y_uni, col = "red")
  
  # mark sf = 1 as "false positive" points
  fp_bi  <- y_bi[x_sf == 1]
  fp_uni <- y_uni[x_sf == 1]
  points(1, fp_bi,  pch = 0, cex = 1.3)  # square
  points(1, fp_uni, pch = 4, cex = 1.3)  # X
  
  # dashed horizontal line at alpha
  abline(h = 0.025, lty = 2)
  
  legend("topleft",
         legend = c("bidirectional power", "unidirectional power",
                    "bidirectional false positive", "unidirectional false positive"),
         pch    = c(16, 16, 0, 4),
         col    = c("blue", "red", "black", "black"),
         bty    = "n", cex = 0.8)
}
