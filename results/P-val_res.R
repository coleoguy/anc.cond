

# i itterates through number taxa
# j itterates through 100 replicates
# k is scaling factor onyl one value for now
res.fp <- as.data.frame(matrix(,100,5))
colnames(res.fp) <- names(unidirectional_results)
for(i in 1:5){ # loop through taxa sizes and columns of res
  for(j in 1:100){
    for(k in 1:length(unidirectional_results[[i]][[j]])){
      res.fp[j,i] <- min(unlist(unidirectional_results[[i]][[j]][[k]]$anccond[4:5]))/100
    }
  }
}
sum(res.fp[,1]<=.025)/100

