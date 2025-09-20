
# just checks all given arguments to make sure that they make sense
# and gives good warning messages if not
InputTesting <- function(tree,
                         data,
                         drop.state,
                         mat,
                         pi,
                         n.tails,
                         nsim,
                         iter){
  ##### testing inputs #####

  if(class(tree) != 'phylo') {stop('tree must be class phylo')}
  if(!is.data.frame(data) & ncol(data) == 3){stop('data should be a dataframe with 3 columns\n(tip labels, cont data, discrete data)')}
  if(!is.null(drop.state)) if(!drop.state %in% c(1,2)){stop('drop.state must be NULL, or numeric 1 or 2')}
  if(!sum(mat == c(0,0,1,0)) == 4 & !sum(mat == c(0,1,1,0)) == 4 & !sum(mat == c(0,2,1,0)) == 4){
    stop('mat must be a vector of the form c(0,0,1,0), c(0,1,1,0), or c(0,2,1,0)')
  }
  if((!pi %in% c('equal', 'estimated'))[1]){
    if(!is.numeric(pi)) stop('pi must be equal, estimated or a vector of length 2\nwith probabilities for the state of the discrete character at the root')
    if(length(pi) != 2 | sum(pi) != 1) stop('pi must be equal, estimated or a vector of length 2\nwith probabilities for the state of the discrete character at the root')
  }
  if(n.tails != 1 & n.tails != 2){stop('n.tails should be numeric 1 or 2')}
  if(!is.numeric(nsim)){stop('nsim should be numeric')}
  if(!is.numeric(iter)){stop('iter should be numeric')}
  if(iter < 100){stop('iter is the number of null data points. It should be
                      greater or equal to 100 to ensure an accurate pvalue')}
}

quiet <- function(x) {
  sink(tempfile())
  on.exit(sink())
  invisible(force(x))
}

CreateNull <- function(tree,                     # a tree type phylo
                       iter,                     # number of simulations for null
                       current.map,             # for Q-matrix
                       anc.states.cont.trait,   # ancestral state reconstruction for continuous
                       dt.vec,
                       message,
                       j,
                       nsim){
  current.Q <- current.map$Q
  if(sum(current.Q == 0)>0){
    current.Q[current.Q == 0] <- c(10^(-25),-10^(-25))
    cat('\n')
    print("Your estimated transition matrix has a rate of zero for some parameters these are being set to 10e-25 for simulation purposes")
  }
  root.state <- c(0,0)
  names(root.state) <- 1:2
  root.state[names(root.state) == names(current.map$maps[[1]])[1]] <- 1
  nulldist <- vector(length=iter, mode="list")
  for(n in 1:iter){
    # while loop is set up to make sure sufficient transitions occur on the tree
    good.sim <- F
    sim.count <- 0
    TO <- F
    while(good.sim == F && TO == F){
      sim.count <- sim.count + 1
      sim.anc.state.dt <- sim.history(tree=tree, Q=current.Q,
                                      nsim=1, message = F,
                                      anc = root.state)
      if(summary(current.map)$N > 5){
        if(summary(sim.anc.state.dt)$N >= .8 * summary(current.map)$N &&
           summary(sim.anc.state.dt)$N <= 1.2 * summary(current.map)$N){
          good.sim <- T
        }
      }else if(summary(current.map)$N <= 5){
        if(summary(sim.anc.state.dt)$N >= summary(current.map)$N - 1 &&
           summary(sim.anc.state.dt)$N <= summary(current.map)$N + 1){
          good.sim <- T
        }
      }
      if(good.sim && message){
        cat('\014')
        cat('Analyzing map: ',j,' of ', nsim,'\n')
          cat('Number of transitions:\n')
          cat(' Emperical Map:\n')
          cat(summary(current.map)$N)
          cat('\n Null Simulation:\n')
          cat(summary(sim.anc.state.dt)$N)
      }
      if(sim.count > 10000){
        if(message){
          cat('Unable to simulate a null with similar behavior to the observed.\n')
          # return('SIG')
          warning('Unable to simulate a null with similar values to the observed.\n')
        }
        TO <- T
      }
    }
    if(good.sim == T){
    nulldist[[n]] <-  exctractAncestral(current.map = sim.anc.state.dt,
                                        anc.states.cont.trait = anc.states.cont.trait)
    }else if(TO == T){
      nulldist[[n]] <- list()
    }
  }
  return(nulldist)
}


# this takes a stochastic map and continuous trait
# and returns observed.anc.cond list of ancestral states at transitions
exctractAncestral <- function(current.map,
                              anc.states.cont.trait,
                              count = F){
  #### Parse simmap to get producing nodes ####
  # the mapped edge object has time spent in a state in
  # two columns so only branches with a change have an entry
  # in both columns
  #######
  # gets branches with transitions
  ss_nodes <- current.map$mapped.edge[, 1] > 0 &
    current.map$mapped.edge[, 2] > 0

  # this returns the node pairs describing a branch with transitions
  wanted_branches <- ss_nodes[ss_nodes == T]
  wanted_nodes <- names(wanted_branches)

  # for the general model we partition the producing nodes for 1->2 and 1<-2 transitions
  producing.nodes12 <- c()
  producing.nodes21 <-c()
  trans.maps <- current.map$maps[ss_nodes == T]
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
  ntrans <- c(length(producing.nodes12), length(producing.nodes21))
  names(ntrans) <- c('12','21')
  producing.nodes12 <- unique(producing.nodes12)
  producing.nodes21 <- unique(producing.nodes21)


  ##### get estimated ancestral conditions ######
  if(count == T){
    observed.anc.cond <- list('12' = anc.states.cont.trait$ace[names(anc.states.cont.trait$ace) %in%
                                                                 producing.nodes12],
                              '21' = anc.states.cont.trait$ace[names(anc.states.cont.trait$ace) %in%
                                                                 producing.nodes21],
                              'ntrans' = ntrans)
  }else{
    observed.anc.cond <- list('12' = anc.states.cont.trait$ace[names(anc.states.cont.trait$ace) %in%
                                                                 producing.nodes12],
                              '21' = anc.states.cont.trait$ace[names(anc.states.cont.trait$ace) %in%
                                                                 producing.nodes21])}
  return(observed.anc.cond)
}

# takes data and any state to be dropped and generates the
# dt.vec, ct.vec that are used by other functions
UnpackData <- function(data, drop.state){
  ##### create named vector for disc trait for all taxa #####
  dt.vec <- data[, 3]
  names(dt.vec) <- data[, 1]


  ##### create named vector for cont trait taxa not in derived state #####
  if(!is.null(drop.state)){
    ct.data <- data[(data[, 3] != drop.state),]
    ct.vec <- as.numeric(ct.data[, 2])
    names(ct.vec) <- ct.data[, 1]
  }else{
    ct.data <- data
    ct.vec <- as.numeric(ct.data[, 2])
    names(ct.vec) <- ct.data[, 1]
  }
  if(sum(is.na(ct.vec)) > 0 | sum(is.na(dt.vec)) > 0){
    stop('There exists missing trait data for some species in the phylogeny.\n
         Please remove such taxa from the tree.')
  }
  return(list(dt.vec, ct.vec))
}

ProcessObserved <- function(observed.anc.cond){
  vals12 <- vector(length = length(observed.anc.cond))
  vals21 <- vector(length = length(observed.anc.cond))
  for(i in 1:length(observed.anc.cond)){
    vals12[i] <- mean(observed.anc.cond[[i]]$'12', na.rm = T)
    vals21[i] <- mean(observed.anc.cond[[i]]$'21', na.rm = T)
  }
  res <- c(mean(vals12, na.rm = T), mean(vals21, na.rm = T))
  names(res) <- c("12", "21")
  return(res)
}
ProcessObservedOneMean <- function(observed.anc.cond){
  vals12 <- vector(length = length(observed.anc.cond))
  vals21 <- vector(length = length(observed.anc.cond))
  for(i in 1:length(observed.anc.cond)){
    vals12[i] <- mean(observed.anc.cond[[i]]$'12', na.rm = T)
    vals21[i] <- mean(observed.anc.cond[[i]]$'21', na.rm = T)
  }

  res <- list('12' = vals12,
              '21' = vals21)
  return(res)
}
ProcessObservedNoMean <- function(observed.anc.cond){
  vals12 <- vector(length = length(observed.anc.cond))
  vals21 <- vector(length = length(observed.anc.cond))
  for(i in 1:length(observed.anc.cond)){
    vals12 <- c(vals12,observed.anc.cond[[i]]$'12')
    vals21 <- c(vals21,observed.anc.cond[[i]]$'21')
  }

  res <- list('12' = vals12,
              '21' = vals21)
  return(res)
}

ProcessNull <- function(null.anc.cond, iter){
  vals12 <- vector(length = length(null.anc.cond))
  vals21 <- vector(length = length(null.anc.cond))
  for(j in 1:iter){
    cur.sim12 <- cur.sim21 <- c()
    for(i in 1:length(null.anc.cond)){
      cur.sim12[i] <- mean(null.anc.cond[[i]][[j]]$'12', na.rm = T)
      cur.sim21[i] <- mean(null.anc.cond[[i]][[j]]$'21', na.rm = T)
    }
    vals12[j] <- mean(cur.sim12, na.rm = T)
    vals21[j] <- mean(cur.sim21, na.rm = T)
  }
  return(list('12' = vals12, '21' = vals21))
}
ProcessNullOneMean <- function(null.anc.cond, iter){
  vals12 <- vector(length = length(null.anc.cond))
  vals21 <- vector(length = length(null.anc.cond))
  for(j in 1:iter){
    cur.sim12 <- cur.sim21 <- c()
    for(i in 1:length(null.anc.cond)){
      cur.sim12[i] <- mean(null.anc.cond[[i]][[j]]$'12', na.rm = T)
      cur.sim21[i] <- mean(null.anc.cond[[i]][[j]]$'21', na.rm = T)
    }
    vals12 <- c(vals12, cur.sim12)
    vals21 <- c(vals21, cur.sim21)
  }
  return(list('12' = vals12, '21' = vals21))
}
ProcessNullNoMean <- function(null.anc.cond, iter){
  vals12 <- vector(length = length(null.anc.cond))
  vals21 <- vector(length = length(null.anc.cond))
  for(j in 1:iter){
    cur.sim12 <- cur.sim21 <- c()
    for(i in 1:length(null.anc.cond)){
      cur.sim12 <- c(cur.sim12, null.anc.cond[[i]][[j]]$'12')
      cur.sim21 <- c(cur.sim21, null.anc.cond[[i]][[j]]$'21')
    }
    vals12 <- c(vals12, cur.sim12)
    vals21 <- c(vals21, cur.sim21)
  }
  return(list('12' = vals12, '21' = vals21))
}

CalcPVal <- function(results,n.tails){
  bigger12 <- sum(results$null$`12` >= results$observed[1], na.rm = T) / length(results$null$`12`)
  smaller12 <- sum(results$null$`12` < results$observed[1], na.rm = T) / length(results$null$`12`)
  if(sum(is.na(results$null$`12`)) < length(results$null$`12`)){
    if (bigger12 <= smaller12){pval12 <- bigger12}
    if (smaller12 < bigger12){pval12 <- smaller12}
    if (n.tails == 2){pval12 <- 2 * pval12}
  }else{pval12 <- NA}
  bigger21 <- sum(results$null$`21` >= results$observed[2], na.rm = T) / length(results$null$`21`)
  smaller21 <- sum(results$null$`21` < results$observed[2], na.rm = T) / length(results$null$`21`)
  if(sum(is.na(results$null$`21`)) < length(results$null$`21`)){
    if (bigger21 <= smaller21){pval21 <- bigger21}
    if (smaller21 < bigger21){pval21 <- smaller21}
    if (n.tails == 2){pval21 <- 2 * pval21}
  }else{pval21 <- NA}
  res <- c(pval12, pval21)
  names(res) <- c('12','21')
  return(res)
}

CountTrans <- function(current.map){
  ntrans <- sum(current.map$mapped.edge[, 1] > 0 &
                  current.map$mapped.edge[, 2] > 0)
  cat(ntrans)
}




summary.AncCond <- function(results){
  ## print results to terminal
  cat('\n')
  cat(paste(
    "Mean value for the continuous trait at 1 - > 2 transitions:",
    round(results$observed[1], digits = 4),
    "\n"
  ))
  cat(paste(
    "Mean value for the continuous trait at 2 - > 1 transitions:",
    round(results$observed[2], digits = 4),
    "\n\n"
  ))
  cat(paste('Mean number of 1 -> 2 transitions:', round(results$`mean n trans`[1], digits = 4), '\n'))
  cat(paste('Mean number of 2 -> 1 transitions:', round(results$`mean n trans`[2], digits = 4), '\n\n'))
  # cat(paste("Number of producing nodes 1->2:", round(mean(number.of.trans12),
  #                                                    digits = 4), "\n"))
  # cat(paste("Number of producing nodes 2->1:", round(mean(number.of.trans21),
  #                                                    digits = 4), "\n"))
  cat(paste("Mean of null dist 1->2:", round(mean(results$null$`12`, na.rm = T),
                                             digits = 4), "\n"))
  cat(paste("Mean of null dist 2->1:", round(mean(results$null$`21`, na.rm = T),
                                             digits = 4), "\n\n"))
  cat(paste("SD of null dist 1->2:", round(sd(results$null$`12`, na.rm = T), digits = 4), "\n"))
  cat(paste("SD of null dist 2->1:", round(sd(results$null$`21`, na.rm = T), digits = 4), "\n\n"))

  cat(paste("pvalue 1->2:", round(results$pvals[1], digits = 4), "\n"))
  cat(paste("pvalue 2->1:", round(results$pvals[2], digits = 4), "\n\n\n"))

  if(is.na(results$pvals[1])){
    cat('NA and NaN values are produced when no transitions of a type have occured. \n\n')
  }
  if(is.na(results$pvals[2])){
    cat('NA and NaN values are produced when no transitions of a type have occured. \n\n')
  }
}
plot.AncCond <- function(results){
  if(!is.na(results$pvals[1])){
    plot(density(results$null$`12`, na.rm = T),
         main = '1 -> 2',
         xlim= c(min(c(results$null$`12`,
                       results$observed[1])),
                 max(c(results$null$`12`,
                       results$observed[1]))),
         xlab = 'Ancestral Condition',
         ylab = 'Frequency')
    abline(v=results$observed[1], col = 'red')
    legend(x = 'topright', legend = c('Observed', 'Null'), col  = c('red', 'black'), lwd = 2)
  }
  if(!is.na(results$pvals[2])){
    plot(density(results$null$`21`, na.rm = T),
         main = '2 -> 1',
         xlim= c(min(c(results$null$`21`,
                       results$observed[2]), na.rm = T),
                 max(c(results$null$`21`,
                       results$observed[2]), na.rm = T)),
         xlab = 'Ancestral Condition',
         ylab = 'Frequency')
    abline(v=results$observed[2], col = 'red')
    legend(x = 'topright', legend = c('Observed', 'Null'), col  = c('red', 'black'), lwd = 2)
  }
}

