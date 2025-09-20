library(R.utils)
library(phytools)
library(diversitree)
library(geiger)
library(doSNOW)
library(foreach)
cl<-makeCluster(3, type="SOCK")
on.exit(stopCluster(cl))
opts <- list(preschedule = FALSE)
registerDoSNOW(cl)

n.trees <- 100
n.taxa <- 200
message <- T
source('AncCond.R', local = TRUE)
percent.vec <- array(dim = c(100,10))
for(i in 1:100){
  
  trees <- trees(pars = c(3,1),
                 type = "bd",
                 n = 1,
                 max.taxa = n.taxa,
                 include.extinct = F)[[1]]
  trees$edge.length <- trees$edge.length / max(branching.times(trees))
  
  # we then simulate the continious character
  cont.trait <- sim.char(trees, 0.2, model = 'BM')
  names(cont.trait) <- trees$tip.label # this line somehow makes anc.ML work????
  
  # identifying which branch had a mean cont trait value in the upper and lower quartiles
  # we do this by 1st doing an ASR for the continious trait
  cont.trait.AC <- anc.ML(trees, cont.trait, model = "BM")
  # this will hold all of the branch means in the same order they are given in trees
  branch.means <- c()
  # branch names is essentially paste(rootward node, tipward node)
  branch.names <- c()
  # then for each branch we go through and calculate the name and mean
  for(j in 1:nrow(trees$edge)){
    # we first find the cont trait value at the rootward node
    node.o.int <- trees$edge[j,1]
    # we have to look in two different places for cont trait values, either in the cont.trait vector 
    # (if the node is a tip) or in the ASR if it is an interior node
    if(node.o.int <= n.taxa){
      one <- cont.trait[node.o.int]
    }else{
      one <- cont.trait.AC$ace[names(cont.trait.AC$ace) == as.character(node.o.int)]
    }
    # we do the same for the tipward node
    node.o.int <- trees$edge[j,2]
    if(node.o.int <= n.taxa){
      two <- cont.trait[node.o.int]
    }else{
      two <- cont.trait.AC$ace[names(cont.trait.AC$ace) == as.character(node.o.int)]
    }
    # to find the mean we avg the rootward and the tipward cont trait values
    branch.means <- c(branch.means, mean(one, two))
    # we create branch names by pasting the rootwward and tipward node labels together
    branch.names <- c(branch.names, paste(as.character(trees$edge[j,1]),as.character(trees$edge[j,2])))
  }
  # we name the branch names for nice bookkeeping
  names(branch.means) <- branch.names
  rm(branch.names)
  # finding upper and lower quartiles
  upper <- summary(branch.means)[[5]]
  lower <- summary(branch.means)[[2]]
  
  for(s in 1:10){
    total.sims <- 0
    total.bad <- 0
    scale.factor <- s
    # we leave the original trees un altered 
    alt.tree <- trees 
    
    # we then manipulate the branch lengths of those branches whose cont trait means are in the upper or lower quartiles
    for(j in 1:length(branch.means)){
      if(branch.means[j] < lower){alt.tree$edge.length[j] <- alt.tree$edge.length[j] / scale.factor}
      if(branch.means[j] > upper){alt.tree$edge.length[j] <- alt.tree$edge.length[j] * scale.factor}
    }
    # next we simulated a discrete trait on this altered tree
    # while loop is set up to make sure sufficient transitions occur on the tree
    good.sim <- F
    bad.count <- 0
    rate <- .1
    while(good.sim == F){
      disc.trait <- sim.char(phy = alt.tree,
                             par = matrix(c(-rate, 0, rate, 0), 2),
                             model = 'discrete',
                             root = 1)
      if((0.05 * n.taxa) < sum(disc.trait == min(disc.trait)) && 
         sum(disc.trait == min(disc.trait)) < (.95 * n.taxa)){
        good.sim <- T
      }else{bad.count <- bad.count + 1}
      
    }
    percent.vec[i, s] <- bad.count
  }
  
}
for(i in 1:10){
  total.sims[i] <- sum(percent.vec[,i] + 1)
  total.bad[i] <- sum(percent.vec[,i])
}
final.prop <- total.bad / total.sims
# save(final.prop, file = 'PercentDiscardedScaling.RData')

n.trees <- 100
scale.factor <- 5
n.taxa <- seq(20, 200, length.out = 10)
source('AncCond.R', local = TRUE)
message <- T
percent.vec <- array(dim = c(100,10))

for(s in 1:10){
  for(i in 1:100){
    total.sims <- 0
    total.bad <- 0
    trees <- trees(pars = c(3,1),
                   type = "bd",
                   n = 1,
                   max.taxa = n.taxa[s],
                   include.extinct = F)[[1]]
    trees$edge.length <- trees$edge.length / max(branching.times(trees))
    
    
    
    
    # we then simulate the continious character
    cont.trait <- sim.char(trees, 0.2, model = 'BM')
    names(cont.trait) <- trees$tip.label # this line somehow makes anc.ML work????
    
    # identifying which branch had a mean cont trait value in the upper and lower quartiles
    # we do this by 1st doing an ASR for the continious trait
    cont.trait.AC <- anc.ML(trees, cont.trait, model = "BM")
    # this will hold all of the branch means in the same order they are given in trees
    branch.means <- c()
    # branch names is essentially paste(rootward node, tipward node)
    branch.names <- c()
    # then for each branch we go through and calculate the name and mean
    for(j in 1:nrow(trees$edge)){
      # we first find the cont trait value at the rootward node
      node.o.int <- trees$edge[j,1]
      # we have to look in two different places for cont trait values, either in the cont.trait vector 
      # (if the node is a tip) or in the ASR if it is an interior node
      if(node.o.int <= n.taxa[s]){
        one <- cont.trait[node.o.int]
      }else{
        one <- cont.trait.AC$ace[names(cont.trait.AC$ace) == as.character(node.o.int)]
      }
      # we do the same for the tipward node
      node.o.int <- trees$edge[j,2]
      if(node.o.int <= n.taxa[s]){
        two <- cont.trait[node.o.int]
      }else{
        two <- cont.trait.AC$ace[names(cont.trait.AC$ace) == as.character(node.o.int)]
      }
      # to find the mean we avg the rootward and the tipward cont trait values
      branch.means <- c(branch.means, mean(one, two))
      # we create branch names by pasting the rootwward and tipward node labels together
      branch.names <- c(branch.names, paste(as.character(trees$edge[j,1]),as.character(trees$edge[j,2])))
    }
    # we name the branch names for nice bookkeeping
    names(branch.means) <- branch.names
    rm(branch.names)
    # finding upper and lower quartiles
    upper <- summary(branch.means)[[5]]
    lower <- summary(branch.means)[[2]]
    
    alt.tree <- trees
    # we then manipulate the branch lengths of those branches whose cont trait means are in the upper or lower quartiles
    for(j in 1:length(branch.means)){
      if(branch.means[j] < lower){alt.tree$edge.length[j] <- alt.tree$edge.length[j] / scale.factor}
      if(branch.means[j] > upper){alt.tree$edge.length[j] <- alt.tree$edge.length[j] * scale.factor}
    }
    # next we simulated a discrete trait on this altered tree
    # while loop is set up to make sure sufficient transitions occur on the tree
    good.sim <- F
    # count <- 0
    bad.count <- 0
    rate <- .1
    # withTimeout({
    while(good.sim == F){
      disc.trait <- sim.char(phy = alt.tree, 
                             par = matrix(c(-rate, 0, rate, 0), 2), 
                             model = 'discrete', 
                             root = 1)
      if(5 < sum(disc.trait == min(disc.trait)) && 
         sum(disc.trait == min(disc.trait)) < (n.taxa[s] - 5)){
        good.sim <- T
      }else{bad.count <- bad.count + 1}
      
    }
    percent.vec[i, s] <- bad.count
  }
}
for(i in 1:10){
  total.sims[i] <- sum(percent.vec[,i] + 1)
  total.bad[i] <- sum(percent.vec[,i])
}
final.prop <- total.bad / total.sims
# save(final.prop, file = 'PercentDiscardedTaxa.RData')