
##### Bidirectional #####
library(phytools)
library(geiger)
library(diversitree)
ntaxa <- 200
scale.factor <- 1
rate <- .4

tree <- trees(pars = c(3,1),
              type = "bd",
              n = 1,
              max.taxa = ntaxa,
              include.extinct = F)[[1]]
tree$edge.length <- tree$edge.length / max(branching.times(tree))

# we then simulate the continious character
cont.trait <- sim.char(tree, 0.2, model = 'BM')[,,1]

# identifying which branch had a mean cont trait value in the upper and lower quartiles
# we do this by 1st doing an ASR for the continious trait
cont.trait.AC <- anc.ML(tree, cont.trait, model = "BM")
# this will hold all of the branch means in the same order they are given in tree
branch.means <- c()
# branch names is essentially paste(rootward node, tipward node)
branch.names <- c()
# then for each branch we go through and calculate the name and mean
for(j in 1:nrow(tree$edge)){
  # we first find the cont trait value at the rootward node
  node.o.int <- tree$edge[j,1]
  # we have to look in two different places for cont trait values, either in the cont.trait vector
  # (if the node is a tip) or in the ASR if it is an interior node
  if(node.o.int <= ntaxa){
    one <- cont.trait[node.o.int]
  }else{
    one <- cont.trait.AC$ace[names(cont.trait.AC$ace) == as.character(node.o.int)]
  }
  # we do the same for the tipward node
  node.o.int <- tree$edge[j,2]
  if(node.o.int <= ntaxa){
    two <- cont.trait[node.o.int]
  }else{
    two <- cont.trait.AC$ace[names(cont.trait.AC$ace) == as.character(node.o.int)]
  }
  # to find the mean we avg the rootward and the tipward cont trait values
  branch.means <- c(branch.means, mean(one, two))
  # we create branch names by pasting the rootwward and tipward node labels together
  branch.names <- c(branch.names, paste(as.character(tree$edge[j,1]),as.character(tree$edge[j,2])))
}
# we name the branch names for nice bookkeeping
names(branch.means) <- branch.names
# finding upper and lower quartiles
upper <- summary(branch.means)[[5]]
lower <- summary(branch.means)[[2]]

# next we perform the following analysis on this tree for each of the scaling factors


# we leave the original tree un altered
alt.tree <- tree

# we then manipulate the branch lengths of those branches whose cont trait means are in the upper or lower quartiles
for(j in 1:length(branch.means)){
  if(branch.means[j] < lower){alt.tree$edge.length[j] <- alt.tree$edge.length[j] / scale.factor}
  if(branch.means[j] > upper){alt.tree$edge.length[j] <- alt.tree$edge.length[j] * scale.factor}
}
# next we simulated a discrete trait on this altered tree
# while loop is set up to make sure sufficient transitions occur on the tree
good.sim <- F
while(good.sim == F){
  disc.trait <- sim.char(phy = alt.tree,
                         par = matrix(c(-rate, rate, rate, -rate), 2),
                         model = 'discrete',
                         root = sample(c(1,2),1))
  if((0.05 * ntaxa) < sum(disc.trait == min(disc.trait)) &&
     sum(disc.trait == min(disc.trait)) < (.95 * ntaxa)){
    good.sim <- T
  }
}
# we now apply the AncCond test to our simulated data and record its result
data <- data.frame(alt.tree$tip.label, cont.trait, disc.trait)

nsim <- 5
iter <- 100
drop.state <- NULL
mat <- c(0,2,1,0)
pi <- 'estimated'
n.tails <- 2
message <- T
make.plot <- T


rm(list=ls()[-c(24,6,8,14,21,16,18,10,15,14,13)])
##### Unidirectional #####
library(phytools)
library(geiger)
library(diversitree)
ntaxa <- 200
scale.factor <- 1
rate <- .2
tree <- trees(pars = c(3,1),
              type = "bd",
              n = 1,
              max.taxa = ntaxa,
              include.extinct = F)[[1]]
tree$edge.length <- tree$edge.length / max(branching.times(tree))

# we then simulate the continious character
cont.trait <- sim.char(tree, 0.2, model = 'BM')[,,1]

# identifying which branch had a mean cont trait value in the upper and lower quartiles
# we do this by 1st doing an ASR for the continious trait
cont.trait.AC <- anc.ML(tree, cont.trait, model = "BM")
# this will hold all of the branch means in the same order they are given in tree
branch.means <- c()
# branch names is essentially paste(rootward node, tipward node)
branch.names <- c()
# then for each branch we go through and calculate the name and mean
for(j in 1:nrow(tree$edge)){
  # we first find the cont trait value at the rootward node
  node.o.int <- tree$edge[j,1]
  # we have to look in two different places for cont trait values, either in the cont.trait vector
  # (if the node is a tip) or in the ASR if it is an interior node
  if(node.o.int <= ntaxa){
    one <- cont.trait[node.o.int]
  }else{
    one <- cont.trait.AC$ace[names(cont.trait.AC$ace) == as.character(node.o.int)]
  }
  # we do the same for the tipward node
  node.o.int <- tree$edge[j,2]
  if(node.o.int <= ntaxa){
    two <- cont.trait[node.o.int]
  }else{
    two <- cont.trait.AC$ace[names(cont.trait.AC$ace) == as.character(node.o.int)]
  }
  # to find the mean we avg the rootward and the tipward cont trait values
  branch.means <- c(branch.means, mean(one, two))
  # we create branch names by pasting the rootwward and tipward node labels together
  branch.names <- c(branch.names, paste(as.character(tree$edge[j,1]),as.character(tree$edge[j,2])))
}
# we name the branch names for nice bookkeeping
names(branch.means) <- branch.names
# finding upper and lower quartiles
upper <- summary(branch.means)[[5]]
lower <- summary(branch.means)[[2]]

# next we perform the following analysis on this tree for each of the scaling factors


# we leave the original tree un altered
alt.tree <- tree

# we then manipulate the branch lengths of those branches whose cont trait means are in the upper or lower quartiles
for(j in 1:length(branch.means)){
  if(branch.means[j] < lower){alt.tree$edge.length[j] <- alt.tree$edge.length[j] / scale.factor}
  if(branch.means[j] > upper){alt.tree$edge.length[j] <- alt.tree$edge.length[j] * scale.factor}
}
# next we simulated a discrete trait on this altered tree
# while loop is set up to make sure sufficient transitions occur on the tree
good.sim <- F
while(good.sim == F){
  disc.trait <- sim.char(phy = alt.tree,
                         par = matrix(c(-rate, 0, rate, 0), 2),
                         model = 'discrete',
                         root = 1)
  if((0.05 * ntaxa) < sum(disc.trait == min(disc.trait)) &&
     sum(disc.trait == min(disc.trait)) < (.95 * ntaxa)){
    good.sim <- T
  }
}
# we now apply the AncCond test to our simulated data and record its result
data <- data.frame(alt.tree$tip.label, cont.trait, disc.trait)

drop.state = NULL
mat = c(0,0,1,0)
pi = c(1,0)
nsim <- 100
iter <- 100
n.tails <- 1
message <- T
make.plot <- T


rm(list=ls()[-c(24,6,8,14,21,16,18,10,15,14,13)])

##### running #####
source('AncCond2.R')
start <- Sys.time()
profvis({
testing <-  AncCond(tree,
                    data,
                    drop.state = 2,
                    mat,
                    pi,
                    n.tails,
                    nsim,
                    iter,
                    message)
})
end <- Sys.time()
end - start
# 3.155765 mins
# 13.79768 mins
plot(testing)

