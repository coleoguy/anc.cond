##### Create testing data ######
library(phytools)
library(geiger)
library(diversitree)
tree <- trees(pars = c(3,1),
              type = "bd",
              n = 1,
              max.taxa = 100,
              include.extinct = F)[[1]]
tree$edge.length <- tree$edge.length / max(branching.times(tree))

# we then simulate the continious character
cont.trait <- sim.char(tree, 0.2, model = 'BM')[,,1]

# identifying which branch had a mean cont trait value in the upper and lower quartiles
# we do this by 1st doing an ASR for the continious trait
cont.trait.AC <- anc.ML(tree, cont.trait, model = "BM")
# this will hold all of the branch means in the same order they are given in trees
branch.means <- c()
# branch names is essentially paste(rootward node, tipward node)
branch.names <- c()
# then for each branch we go through and calculate the name and mean
for(j in 1:nrow(tree$edge)){
  # we first find the cont trait value at the rootward node
  node.o.int <- tree$edge[j,1]
  # we have to look in two different places for cont trait values, either in the cont.trait vector
  # (if the node is a tip) or in the ASR if it is an interior node
  if(node.o.int <= 100){
    one <- cont.trait[node.o.int]
  }else{
    one <- cont.trait.AC$ace[names(cont.trait.AC$ace) == as.character(node.o.int)]
  }
  # we do the same for the tipward node
  node.o.int <- tree$edge[j,2]
  if(node.o.int <= 100){
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

scale.factor <- 5
# we leave the original trees un altered
alt.tree <- tree

# we then manipulate the branch lengths of those branches whose cont trait means are in the upper or lower quartiles
for(j in 1:length(branch.means)){
  if(branch.means[j] < lower){alt.tree$edge.length[j] <- alt.tree$edge.length[j] / scale.factor}
  if(branch.means[j] > upper){alt.tree$edge.length[j] <- alt.tree$edge.length[j] * scale.factor}
}
# next we simulated a discrete trait on this altered tree
# while loop is set up to make sure sufficient transitions occur on the tree
good.sim <- F
rate <- .3
while(good.sim == F){
  disc.trait <- sim.char(phy = alt.tree,
                         par = matrix(c(-rate, rate, rate, -rate), 2),
                         model = 'discrete',
                         root = sample(c(1,2),1))
  if((0.05 * 100) < sum(disc.trait == min(disc.trait)) &&
     sum(disc.trait == min(disc.trait)) < (.95 * 100)){
    good.sim <- T
  }
}
# we now apply the AncCond test to our simulated data and record its result
data <- data.frame(alt.tree$tip.label, cont.trait, disc.trait)
##### create named vector for disc trait for all taxa #####
dt.vec <- data[, 3]
names(dt.vec) <- data[, 1]
ct.data <- data
ct.vec <- as.numeric(ct.data[, 2])
names(ct.vec) <- ct.data[, 1]
## ASR for the continuous trait
anc.states.cont.trait <- cont.trait.AC
## ASR for discrete trait
## using stochastic mappings to nail down specific transition points
anc.state.dt <- make.simmap(tree, dt.vec,
                            model = matrix(c(0,2,1,0), 2),
                            nsim = 1,
                            pi = "estimated",
                            message = F)
iter <- 100
mat <- c(0,2,1,0)


rm(list=ls()[-c(2,3,22,14,17)]) 
######







CreateNull <- function(tree,                     # a tree type phylo
                       iter,                     # number of simulations for null
                       anc.state.dt,             # for Q-matrix
                       mat,                      # model specification matrix
                       anc.states.cont.trait){   # ancestral state reconstruction for continuous
  if(sum(mat) > 1){
    nulldist <- matrix(NA,iter,2)
    colnames(nulldist) <- c('12','21')
  }else{
    nulldist <- vector(length = iter)
  }
  for(n in 1:iter){
    # while loop is set up to make sure sufficient transitions occur on the tree
    good.sim <- F
    while(good.sim == F){
      if(sum(mat) > 1){
        null.disc.trait <- sim.char(phy = tree,
                                    par = anc.state.dt$Q,
                                    model = 'discrete',
                                    root = sample(c(1,2),1))
      }else{
        null.disc.trait <- sim.char(phy = tree,
                                    par = matrix(c(-rate, 0, rate, 0), 2),
                                    model = 'discrete',
                                    root = sample(c(1,2),1))
      }
      if(5 < sum(null.disc.trait == min(null.disc.trait)) &&
         sum(null.disc.trait == min(null.disc.trait)) < (length(tree$tip.label) - 5)){
        good.sim <- T
      }
    }
    nullnames <- rownames(null.disc.trait)
    null.disc.trait <- as.factor(null.disc.trait)
    names(null.disc.trait) <- nullnames
    anc.state.dt <- make.simmap(tree, as.factor(null.disc.trait),
                                model = matrix(mat, 2),
                                nsim = 1,
                                pi = pi,
                                message = F)
    
    ## Parse simmap to get producing nodes
    # the mapped edge object has time spent in a state in
    # two columns so only branches with a change have an entry
    # in both columns
    ss_nodes <- anc.state.dt$mapped.edge[, 1] > 0 &
      anc.state.dt$mapped.edge[, 2] > 0
    
    # this returns the node pairs describing a branch with origins
    wanted_branches <- ss_nodes[ss_nodes == T]
    wanted_nodes <- names(wanted_branches)
    if(sum(mat) > 1){
      # for the general model we partition the producing nodes for 1->2 and 1<-2 transitions
      producing.nodes12 <- c()
      producing.nodes21 <-c()
      trans.maps <- anc.state.dt$maps[ss_nodes == T]
      # now we take the rootward node of each branch and get rid of duplicates
      wanted_nodes <- gsub(",.*", "", wanted_nodes)
      ##### Just realized we can do this with describe.simmap :(
      ##### But i dont want to change it, it would require match function
      for(i in 1:length(wanted_nodes)){
        if(names(trans.maps[[i]])[1] == '1'){
          producing.nodes12 <- c(producing.nodes12, wanted_nodes[i])
        }else if(names(trans.maps[[i]])[1] == '2'){
          producing.nodes21 <- c(producing.nodes21, wanted_nodes[i])
        }
      }
      producing.nodes12 <- unique(producing.nodes12)
      producing.nodes21 <- unique(producing.nodes21)
      ## get the mean ancestral value for the cont trait
      ## at nodes producing the derived state marginalizing across tree
      anc.states <- anc.states.cont.trait
      orig.val12 <- mean(anc.states$ace[names(anc.states$ace) %in%
                                          producing.nodes12])
      orig.val21 <- mean(anc.states$ace[names(anc.states$ace) %in%
                                          producing.nodes21])
      nulldist[n,] <- c(orig.val12,orig.val21)
    }else{
      # now we take the rootward node of each branch and get rid of duplicates
      wanted_nodes <- gsub(",.*", "", wanted_nodes)
      producing.nodes <- unique(wanted_nodes)
      ## get the mean ancestral value for the cont trait
      ## at nodes producing the derived state marginalizing across tree
      anc.states <- anc.states.cont.trait
      orig.val <- mean(anc.states$ace[names(anc.states$ace) %in%
                                        producing.nodes])
      nulldist[n] <- orig.val12
    }
    cat('\014')
    cat(n)
  }
  return(nulldist)
}
