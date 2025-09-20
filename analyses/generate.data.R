

GenerateData <- function(n.trees = 100,
             n.taxa = 100,
             sfac = c(1,2),
             drate = c(0.0001, .2),
             crate = .2){
  trees <- traits <- list()
  for(i in 1:n.trees){
    print(paste("working on tree", i))
    tree <- trees(pars = c(3,1),
                  type = "bd",
                  n = 1,
                  max.taxa = n.taxa,
                  include.extinct = F)[[1]]
    tree$edge.length <- tree$edge.length / max(branching.times(tree))
    trees[[i]] <- tree
    # we then simulate the continious character
    cont.trait <- sim.char(tree, crate, model = 'BM')[,,1]

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
      if(node.o.int <= n.taxa){
        one <- cont.trait[node.o.int]
      }else{
        one <- cont.trait.AC$ace[names(cont.trait.AC$ace) == as.character(node.o.int)]
      }
      # we do the same for the tipward node
      node.o.int <- tree$edge[j,2]
      if(node.o.int <= n.taxa){
        two <- cont.trait[node.o.int]
      }else{
        two <- cont.trait.AC$ace[names(cont.trait.AC$ace) == as.character(node.o.int)]
      }
      # to find the mean we avg the rootward and the tipward cont trait values
      branch.means <- c(branch.means, mean(one, two))
      # we create branch names by pasting the rootwward and tipward node labels together
      branch.names <- c(branch.names, paste(as.character(tree$edge[j,1]),
                                            as.character(tree$edge[j,2])))
    }
    # we name the branch names for nice bookkeeping
    names(branch.means) <- branch.names
    rm(branch.names)
    # finding upper and lower quartiles
    upper <- summary(branch.means)[[5]]
    lower <- summary(branch.means)[[2]]
      # we leave the original trees un altered
      alt.tree <- tree

      # we then manipulate the branch lengths of those branches whose cont trait means are in the upper or lower quartiles
      for(j in 1:length(branch.means)){
        if(branch.means[j] < lower){
          alt.tree$edge.length[j] <- alt.tree$edge.length[j] / sfac
        }
        if(branch.means[j] > upper){
          alt.tree$edge.length[j] <- alt.tree$edge.length[j] * sfac
        }
      }
      # next we simulated a discrete trait on this altered tree
      # while loop is set up to make sure sufficient transitions occur on the tree
      good.sim <- F
      while(good.sim == F){
        disc.trait <- sim.char(phy = alt.tree,
                               par = matrix(c(-drate[2], drate[1],
                                              drate[2], -drate[1]), 2),
                               model = 'discrete',
                               root = 1)
        if((0.05 * n.taxa) < sum(disc.trait == min(disc.trait)) &&
           sum(disc.trait == min(disc.trait)) < (.95 * n.taxa)){
          good.sim <- T
        }
      }
      traits[[i]] <- data.frame(alt.tree$tip.label, cont.trait, disc.trait)
      colnames(traits[[i]]) <- c("sp","cont","disc")
  }
  results <- list(trees, traits)
  return(results)
}
