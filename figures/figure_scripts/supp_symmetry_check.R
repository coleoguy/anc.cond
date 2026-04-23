df <- read.csv("01_sim_power_results.csv")

# 1. Isolate and Reshape (Original Logic)
anc_subset <- df[df$method %in% c("anccond.01", "anccond.10"), ]
anc_wide <- reshape(anc_subset[, c("n.tips", "q.rate", "s", "method", "rate")], 
                    idvar = c("n.tips", "q.rate", "s"), 
                    timevar = "method", 
                    direction = "wide")

# 2. Calculate CCC
x <- anc_wide$rate.anccond.01
y <- anc_wide$rate.anccond.10
ccc_val <- (2 * cov(x, y)) / (var(x) + var(y) + (mean(x) - mean(y))^2)

# 3. View discrepancies
anc_wide$rate_diff <- x - y
anc_wide[order(abs(anc_wide$rate_diff), decreasing = TRUE), ]

# --- SAVE TO PDF ---
pdf("AncCond_Symmetry_Check.pdf", width = 7, height = 7)

plot(x, y, 
     xlab = "Rate (0 to 1)", ylab = "Rate (1 to 0)",
     main = "Symmetry Check (Agreement Test)", 
     pch = 19, col = rgb(0, 0, 0, 0.5),
     xlim = c(0, 1), ylim = c(0, 1))

abline(0, 1, col = "red", lty = 2, lwd = 2)

# Add CCC label to the top-left of the plot
text(0.05, 0.95, labels = paste("CCC =", round(ccc_val, 4)), 
     adj = 0, font = 2, cex = 1.1)

dev.off()