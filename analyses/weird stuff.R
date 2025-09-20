mytree <- rcoal(100)
Q <- matrix(c(-1,0,1,0),2,2)
root <- c(1,0)
names(root) <- c("1","2")
sim.history(tree=mytree, Q=Q, 
            nsim=1, 
            anc = root)
