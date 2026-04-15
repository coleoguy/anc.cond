df <- read.csv("01_sim_power_results.csv")

anc_results <- subset(summary.df, method %in% c("anccond.01", "anccond.10"), 
                      select = c(n.tips, s, method, rate))

# 1. Isolate the AncCond methods
anc_subset <- summary.df[summary.df$method %in% c("anccond.01", "anccond.10"), ]

# 2. Reshape to 'wide' format so both rates are on one line
anc_wide <- reshape(anc_subset[, c("n.tips", "q.rate", "s", "method", "rate")], 
                    idvar = c("n.tips", "q.rate", "s"), 
                    timevar = "method", 
                    direction = "wide")

# 3. Calculate the difference
anc_wide$rate_diff <- anc_wide$rate.anccond.01 - anc_wide$rate.anccond.10

# 4. View the conditions with the largest discrepancies
anc_wide[order(abs(anc_wide$rate_diff), decreasing = TRUE), ]

# --- SAVE TO PDF ---
pdf("AncCond_Symmetry_Check.pdf", width = 7, height = 7) # Open PDF device

plot(anc_wide$rate.anccond.01, anc_wide$rate.anccond.10, 
     xlab = "Rate (0 to 1)", ylab = "Rate (1 to 0)",
     main = "Symmetry Check", pch = 19,
     col = rgb(0, 0, 0, 0.5)) # Added transparency to see overlapping points
abline(0, 1, col = "red", lty = 2) # The "Perfect Symmetry" line

dev.off() # Close the PDF device and save the file
# -------------------