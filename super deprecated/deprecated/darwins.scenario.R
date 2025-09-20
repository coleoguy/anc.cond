# make data
library(diversitree)
tree <- trees(pars= c(3, 1), max.taxa = 200, type="bd", include.extinct = F)[[1]]
tree$edge.length <- tree$edge.length/max(branching.times(tree))
library(geiger)
cont.trait <- sim.char(tree, 0.2)[,,1]

# reconstruct values
library(phytools)
anc.states <- anc.ML(tree, cont.trait)
node.vals <- anc.states$ace

# order the reconstructions
node.vals <- node.vals[order(node.vals)]
searching <- T
i <- 0
while(searching){
  i <- i + 1
  curr.count <- sum(getDescendants(tree, node=names(node.vals)[i])<=200)
  if(curr.count>=15){
    searching <- F
  }
}

# this is the key node 
names(node.vals)[i]

# get the tips that belong to it
node.plus.tips <- getDescendants(tree, node=names(node.vals)[i])
tips <- node.plus.tips[node.plus.tips<=200]
disc.char <- rep("A", 200)
names(disc.char) <- tree$tip.label

disc.char[names(disc.char) %in% tree$tip.label[tips]] <- "B"

fitDiscrete(tree, disc.char, model = "ARD")
