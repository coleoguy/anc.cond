library(phytools)
library(plotfunctions)

# --- 0. INDIVIDUAL THICKNESS SETTINGS ---
thick1 <- 2.0  
thick2 <- 3.5  
thick3 <- 2.0  

# --- 1. Data Simulation ---
set.seed(80)
n_tips <- 20
tree <- pbtree(n = n_tips, scale = 1)
trait_all <- fastBM(tree, a = 1, sig2 = 0.5, internal = TRUE)
trait_tips <- trait_all[1:n_tips]

tree_labeled <- tree
tree_labeled$tip.label <- as.character(round(trait_tips, 2))
names(trait_tips) <- tree_labeled$tip.label

# --- 2. Scaling Logic for Panel 3 ---
tree.scaled <- tree
s <- 3 
branch.avg <- apply(tree$edge, 1, function(nodes) mean(trait_all[nodes]))
direction <- sample(c(1L, -1L), 1)
q25 <- quantile(branch.avg, 0.25); q75 <- quantile(branch.avg, 0.75)

if (s != 1) {
  if (direction == 1L) {
    tree.scaled$edge.length[branch.avg >= q75] <- tree.scaled$edge.length[branch.avg >= q75] * s
    tree.scaled$edge.length[branch.avg <= q25] <- tree.scaled$edge.length[branch.avg <= q25] / s
  } else {
    tree.scaled$edge.length[branch.avg <= q25] <- tree.scaled$edge.length[branch.avg <= q25] * s
    tree.scaled$edge.length[branch.avg >= q75] <- tree.scaled$edge.length[branch.avg >= q75] / s
  }
}

# --- 3. Plotting ---
common_mar <- c(7, 4, 3, 4) 
par(mfrow = c(1, 3)) 

# --- PANEL 1 ---
par(mar = common_mar)
plot(tree, show.tip.label = FALSE, type = "phylogram", main = "", edge.width = thick1) 
axisPhylo(backward = FALSE)

# --- PANEL 2 (Heatmap) ---
lims <- range(trait_tips)
smp <- contMap(tree_labeled, trait_tips, plot = FALSE, lims = lims)
n_cols <- length(smp$cols)
my_cols <- rainbow(n_cols, end = 4/6)

# SWITCH THE COLORS ON THE TREE HERE:
smp$cols[1:n_cols] <- rev(my_cols) 

plot(smp, 
     mar = common_mar, 
     legend = FALSE, 
     main = "", 
     ftype = 'reg', 
     fsize = 1.0, 
     lwd = thick2, 
     outline = TRUE, 
     axes = FALSE)

mtext("Cont. trait value", side = 1, line = 3.2, cex = 0.8)

# Legend remains switched to match the new tree colors
gradientLegend(depth = 0.015, 
               valRange = lims, 
               pos = 0.12, 
               side = 1, 
               color = rev(my_cols), 
               labels = rev(round(lims, 2)), 
               length = 0.5,
               n.ticks = 2,
               border = "black")

# --- PANEL 3 ---
par(mar = common_mar)
plot(tree.scaled, show.tip.label = FALSE, main = "", edge.width = thick3) 
axisPhylo(backward = FALSE)