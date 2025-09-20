
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

message <- T
source('AncCond2.R', local = TRUE)
source("generate.data.R")
n.trees = 10
n.taxa = 100
scaling.factors = c(1,2)
drate = c(0.0001, .2)
crate = .2
dat <- GenerateData(n.trees = n.trees,
                    n.taxa = n.taxa,
                    scaling.factors = scaling.factors,
                    drate = drate,
                    crate = crate)


# we do the following for each of 200 trees
# this will hold the p.val for each of 200 tests for the 10 scaling factors
# dont need this for mc
p.val.array <- array(dim = c(n.trees, 5))



p.val.array <-foreach(t = 1:n.trees, .options.multicore=opts, .combine = 'rbind',
                      .packages=c("phytools","diversitree","geiger")) %dopar%{

                        p.val.vec <- c()
                          rslt <- AncCond(tree = dat[[1]][[t]],
                                          data = dat[[2]][[t]],
                                          drop.state = NULL,
                                          mat = c(0,2,1,0),
                                          pi = "estimated",
                                          message = T)
                          p.val.vec[t] <- rslt$pval
                          if(message == T){cat(' s = ', s)}
                        }
                        if(message == T){cat('\n\n',' t = ', t)
                        }
                        p.val.vec
                      }

scaling.uni.results <- p.val.array
save(scaling.uni.results, file = 'UnidirectionalScalingAnalysisResults.RData')
