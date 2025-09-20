
# install.packages("phytools")
# install.packages("diversitree")
# install.packages("geiger")
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
scale.factor <- 5
n.taxa <- seq(20, 200, length.out = 10)
message <- T
source('AncCond.R', local = TRUE)
rate <- .2

## this will hold the p.val for each of 100 tests for the 10 tree sizes
p.val.array <- array(dim = c(n.trees, 10))

# this time we vary the size of the tree

p.val.array <-foreach(s = 1:length(n.taxa), .options.multicore=opts, .combine = 'cbind', 
                      .packages=c("phytools","diversitree","geiger")) %dopar%{
                        
                        p.val.vec <- c()
                        # for each tree size we repeat the same analysis for the same number of trees
                        for(t in 1:n.trees){
                          # We begin with a single tree and test it at every scaling factor then move to the next tree
                          # first the tree
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
                          while(good.sim == F){
                            disc.trait <- sim.char(phy = alt.tree, 
                                                   par = matrix(c(-rate, 0, rate, 0), 2), 
                                                   model = 'discrete', 
                                                   root = 1)
                            if(5 < sum(disc.trait == min(disc.trait)) && 
                               sum(disc.trait == min(disc.trait)) < (n.taxa[s] - 5)){
                              good.sim <- T
                            }
                          }
                          if(message == T){cat('\n')}
                          # we now apply the AncCond test to our simulated data and record its result
                          dat <- data.frame(alt.tree$tip.label, cont.trait, disc.trait)
                          rslt <- AncCond(tree = trees, 
                                          data = dat, 
                                          drop.state = 2, 
                                          mat = c(0,0,1,0), 
                                          pi = c(1,0), 
                                          message = T)
                          
                          if(message == T){cat('\n')}
                          if(message == T){
                            cat('\n')
                            cat(' t = ', t)
                          }
                          p.val.vec[t] <- rslt$pvals[1]
                        }
                        if(message == T){
                          cat('\n')
                          cat(' s = ', s)
                        }
                        p.val.vec
                        # end <- sys.time
                      }
taxa.uni.power.results <- p.val.array
save(taxa.uni.power.results, file = 'UnidirectionalTaxaPowerResults.RData')
