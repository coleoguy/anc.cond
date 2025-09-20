

plot(density(null.dist$`12`))
lines(density(obs.dist$`12`, adjust = 1), col = 'red')

plot(density(null.dist$`21`))
lines(density(obs.dist$`21`, adjust = 1), col = 'red')

# mean approach
obs.dist <- ProcessObserved2(observed.anc.cond)
null.dist <- ProcessNull2(null.anc.cond)

plot(density(null.dist$`12`), main = '12 Mean Approach', xlim=range(ct.vec))
abline(v=obs.dist[1], col = 'red')

plot(density(null.dist$`21`[!is.na(null.dist$`21`)]), main = '21 Mean Approach')
lines(density(obs.dist$`21`[!is.na(obs.dist$`21`)], adjust = 1), col = 'red')
}


##### DEPRECIEATED ######
# ## get the mean ancestral value for the cont trait
# ## at nodes producing the derived state marginalizing across tree
# orig.val12 <- mean(anc.states.cont.trait$ace[names(anc.states.cont.trait$ace) %in%
#                                                producing.nodes12])
# orig.val21 <- mean(anc.states.cont.trait$ace[names(anc.states.cont.trait$ace) %in%
#                                                producing.nodes21])
# ## Produce the null distribution of nodes in ancestral cond
# null.orig.val12 <- vector(length = mc)
# null.orig.val21 <- vector(length = mc)
# number.of.trans12 <- length(producing.nodes12)
# number.of.trans21 <- length(producing.nodes21)
# anc.dt <- anc.state.dt
# anc.ct <- anc.states.cont.trait
# node.states <- describe.simmap(anc.dt)$states
# anc.cond.nodes12 <- anc.ct$ace[names(anc.ct$ace) %in%
#                                  names(node.states)[node.states != '2']]
# anc.cond.nodes21 <- anc.ct$ace[names(anc.ct$ace) %in%
#                                  names(node.states)[node.states != '1']]
# ##### creating null distribution #
# ##### this is where the create null function will go
# for (j in 1:mc){
#   # set.seed(j)
#   null.orig.val12[j] <- mean(sample(anc.cond.nodes12,
#                                     length(producing.nodes12)))
#   null.orig.val21[j] <- mean(sample(anc.cond.nodes21,
#                                     length(producing.nodes21)))
# }
# ## how many more extreme
###### DEPRECIATED ######

###### MOVED TO END #####
# bigger12 <- (sum(null.orig.val12 >= orig.val12) / mc)
# smaller12 <- (sum(null.orig.val12 <= orig.val12) / mc)
# if(!is.null(producing.nodes12)){
#   if (bigger12 <= smaller12){pval12 <- bigger12}
#   if (smaller12 < bigger12){pval12 <- smaller12}
#   if (n.tails == 2){pval12 <- 2 * pval12}
# }else{
#   pval12 <- NA
# }
#
# bigger21 <- (sum(null.orig.val21 >= orig.val21) / mc)
# smaller21 <- (sum(null.orig.val21 <= orig.val21) / mc)
# if(!is.null(producing.nodes21)){
#   if (bigger21 <= smaller21){pval21 <- bigger21}
#   if (smaller21 < bigger21){pval21 <- smaller21}
#   if (n.tails == 2){pval21 <- 2 * pval21}
# }else{
#   pval21 <- NA
# }
#
# ## print results to terminal
# if (message == T){
#   cat(paste(
#     "Mean value for the continuous trait at origin oftrait 2:",
#     round(orig.val12, digits = 4),
#     "\n"
#   ))
#   cat(paste(
#     "Mean value for the continuous trait at origin of trait 1:",
#     round(orig.val21, digits = 4),
#     "\n"
#   ))
#   cat(paste("Number of producing nodes 1->2:", round(mean(number.of.trans12),
#                                                      digits = 4), "\n"))
#   cat(paste("Number of producing nodes 2->1:", round(mean(number.of.trans21),
#                                                      digits = 4), "\n"))
#   cat(paste("Mean of null dist 1->2:", round(mean(null.orig.val12),
#                                              digits = 4), "\n"))
#   cat(paste("Mean of null dist 2->1:", round(mean(null.orig.val21),
#                                              digits = 4), "\n"))
#   cat(paste("SD of null dist 1->2:", round(sd(null.orig.val12), digits = 4), "\n"))
#   cat(paste("SD of null dist 2->1:", round(sd(null.orig.val21), digits = 4), "\n"))
#
#   cat(paste("pvalue 1->2:", round(pval12, digits = 4), "\n"))
#   cat(paste("pvalue 2->1:", round(pval21, digits = 4), "\n\n"))
#   if(is.null(producing.nodes12)){cat('No 1 -> 2 transitions occured NA and NaN values produced.')}
#   if(is.null(producing.nodes21)){cat('No 2 -> 1 transitions occured NA and NaN values produced.')}
# }
#
# ## return results to user
# results <- list()
# results[[1]] <- orig.val12
# results[[2]] <- number.of.trans12
# results[[3]] <- null.orig.val12
# results[[4]] <- pval12
# results[[5]] <- orig.val21
# results[[6]] <- number.of.trans21
# results[[7]] <- null.orig.val21
# results[[8]] <- pval21
# names(results) <- c("OriginatingVals1->2", "NTrans1->2",
#                     "NullDist1->2", "pval1->2","OriginatingVals2->1", "NTrans2->1",
#                     "NullDist2->1", "pval2->1")
##### MOVED TO END #####
}else{
  # now we take the rootward node of each branch and get rid of duplicates
  wanted_nodes <- gsub(",.*", "", wanted_nodes)
  producing.nodes <- unique(wanted_nodes)
  ## get the mean ancestral value for the cont trait
  ## at nodes producing the derived state marginalizing across tree
  observed.anc.cond[[j]] <-  list(anc.states.cont.trait$ace[names(anc.states.cont.trait$ace) %in%
                                                              producing.nodes])
  null.anc.cond[[j]] <- CreateNull(tree = tree,
                                   iter = 10,
                                   anc.state.dt = current.map,
                                   mat = mat,
                                   anc.states.cont.trait = anc.states.cont.trait)
  ##### DEPRECATED #####
  # anc.states <- anc.states.cont.trait
  # orig.val <- mean(anc.states$ace[names(anc.states$ace) %in%
  #                                   producing.nodes])
  #
  # ## Produce the null distribution of nodes in ancestral cond
  # null.orig.val <- vector(length = mc)
  # number.of.trans <- length(producing.nodes)
  # anc.dt <- anc.state.dt
  # anc.ct <- anc.states.cont.trait
  # node.states <- describe.simmap(anc.dt)$states
  # anc.cond.nodes <- anc.ct$ace[names(anc.ct$ace) %in%
  #                                names(node.states)[node.states != '2']]
  #
  # for (j in 1:mc){
  #   # set.seed(j)
  #   null.orig.val[j] <- mean(sample(anc.cond.nodes,
  #                                   length(producing.nodes)))
  # }
  ##### DEPRECIATED ######

  ##### MOVED TO END #####
  # ## how many more extreme
  #
  # bigger <- (sum(null.orig.val >= orig.val) / mc)
  # smaller <- (sum(null.orig.val <= orig.val) / mc)
  # if (bigger <= smaller){pval <- bigger}
  # if (smaller < bigger){pval <- smaller}
  # if (n.tails == 2){
  #   pval <- 2 * pval
  # }
  # ## print results to terminal
  # if(message == T){cat(paste(
  #   "Mean value for the continuous trait at origin of derived trait:",
  #   round(orig.val, digits = 4),
  #   "\n"
  # ))
  #   cat(paste("Number of producing nodes:", round(mean(number.of.trans),
  #                                                 digits = 4), "\n"))
  #   cat(paste("Mean of null dist:", round(mean(null.orig.val),
  #                                         digits = 4), "\n"))
  #   cat(paste("SD of null dist:", round(sd(null.orig.val), digits = 4), "\n"))
  #
  #   cat(paste("pvalue:", round(pval, digits = 4), "\n\n"))
  # }
  #
  # ## return results to user
  # results <- list()
  # results[[1]] <- orig.val
  # results[[2]] <- number.of.trans
  # results[[3]] <- null.orig.val
  # results[[4]] <- pval
  # names(results) <- c("OriginatingVals", "NTrans",
  #                     "NullDist", "pval")
  ##### MOVED TO END ######
}
}
###### MOVED TO END #####
# bigger12 <- (sum(null.orig.val12 >= orig.val12) / mc)
# smaller12 <- (sum(null.orig.val12 <= orig.val12) / mc)
# if(!is.null(producing.nodes12)){
#   if (bigger12 <= smaller12){pval12 <- bigger12}
#   if (smaller12 < bigger12){pval12 <- smaller12}
#   if (n.tails == 2){pval12 <- 2 * pval12}
# }else{
#   pval12 <- NA
# }
#
# bigger21 <- (sum(null.orig.val21 >= orig.val21) / mc)
# smaller21 <- (sum(null.orig.val21 <= orig.val21) / mc)
# if(!is.null(producing.nodes21)){
#   if (bigger21 <= smaller21){pval21 <- bigger21}
#   if (smaller21 < bigger21){pval21 <- smaller21}
#   if (n.tails == 2){pval21 <- 2 * pval21}
# }else{
#   pval21 <- NA
# }
#
# ## print results to terminal
# if (message == T){
#   cat(paste(
#     "Mean value for the continuous trait at origin oftrait 2:",
#     round(orig.val12, digits = 4),
#     "\n"
#   ))
#   cat(paste(
#     "Mean value for the continuous trait at origin of trait 1:",
#     round(orig.val21, digits = 4),
#     "\n"
#   ))
#   cat(paste("Number of producing nodes 1->2:", round(mean(number.of.trans12),
#                                                      digits = 4), "\n"))
#   cat(paste("Number of producing nodes 2->1:", round(mean(number.of.trans21),
#                                                      digits = 4), "\n"))
#   cat(paste("Mean of null dist 1->2:", round(mean(null.orig.val12),
#                                              digits = 4), "\n"))
#   cat(paste("Mean of null dist 2->1:", round(mean(null.orig.val21),
#                                              digits = 4), "\n"))
#   cat(paste("SD of null dist 1->2:", round(sd(null.orig.val12), digits = 4), "\n"))
#   cat(paste("SD of null dist 2->1:", round(sd(null.orig.val21), digits = 4), "\n"))
#
#   cat(paste("pvalue 1->2:", round(pval12, digits = 4), "\n"))
#   cat(paste("pvalue 2->1:", round(pval21, digits = 4), "\n\n"))
#   if(is.null(producing.nodes12)){cat('No 1 -> 2 transitions occured NA and NaN values produced.')}
#   if(is.null(producing.nodes21)){cat('No 2 -> 1 transitions occured NA and NaN values produced.')}
# }
#
# ## return results to user
# results <- list()
# results[[1]] <- orig.val12
# results[[2]] <- number.of.trans12
# results[[3]] <- null.orig.val12
# results[[4]] <- pval12
# results[[5]] <- orig.val21
# results[[6]] <- number.of.trans21
# results[[7]] <- null.orig.val21
# results[[8]] <- pval21
# names(results) <- c("OriginatingVals1->2", "NTrans1->2",
#                     "NullDist1->2", "pval1->2","OriginatingVals2->1", "NTrans2->1",
#                     "NullDist2->1", "pval2->1")
##### MOVED TO END #####
##### MOVED TO END #####
# ## how many more extreme
#
# bigger <- (sum(null.orig.val >= orig.val) / mc)
# smaller <- (sum(null.orig.val <= orig.val) / mc)
# if (bigger <= smaller){pval <- bigger}
# if (smaller < bigger){pval <- smaller}
# if (n.tails == 2){
#   pval <- 2 * pval
# }
# ## print results to terminal
# if(message == T){cat(paste(
#   "Mean value for the continuous trait at origin of derived trait:",
#   round(orig.val, digits = 4),
#   "\n"
# ))
#   cat(paste("Number of producing nodes:", round(mean(number.of.trans),
#                                                 digits = 4), "\n"))
#   cat(paste("Mean of null dist:", round(mean(null.orig.val),
#                                         digits = 4), "\n"))
#   cat(paste("SD of null dist:", round(sd(null.orig.val), digits = 4), "\n"))
#
#   cat(paste("pvalue:", round(pval, digits = 4), "\n\n"))
# }
#
# ## return results to user
# results <- list()
# results[[1]] <- orig.val
# results[[2]] <- number.of.trans
# results[[3]] <- null.orig.val
# results[[4]] <- pval
# names(results) <- c("OriginatingVals", "NTrans",
#                     "NullDist", "pval")
##### MOVED TO END ######
return(results)
}



























