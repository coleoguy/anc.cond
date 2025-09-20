
par(mfrow = c(1,1), mar = c(5,4,4,2) + .1)
load('../../results/UnidirectionalScalingAnalysisResults.RData')


probs <- vector()
for(i in 1:10){
  probs[i] <- sum(scaling.uni.results[1:100, i] <= .05) / 100
}
x <- rep(1:10)
load('../../results/BidirectionalScalingAnalysisResults.RData')


y <- vector()
for(j in 1:10){
  for(i in 1:100){
    y[2 * (100 * (j - 1) + i - 1) + 1] <- as.numeric(gsub(",.*", "",scaling.bi.results[i,j]))
    y[2 * (100 * (j - 1) + i - 1) + 2] <- as.numeric(substr(scaling.bi.results[i,j], (nchar(gsub(",.*", "",scaling.bi.results[i,j])) + 2), 
                                                            nchar(scaling.bi.results[i,j])))
  }
}
probs2 <- vector()
for(i in 1:10){
  probs2[i] <- round(sum(y[(200 * (i - 1) + 1):(200 * i)] <= .025, na.rm = T) / 
                       sum(!is.na(y[(200 * (i - 1) + 1):(200 * i)])),digits = 2)
}

results <- as.data.frame(matrix(,20,3))
colnames(results) <- c("value", "scaling.factor", "type")
results[, 1] <- c(probs2,probs)
results[, 2] <- c(1:10,1:10)
x <- factor(c("bidirectional false positive",
       rep("bidirectional power",9),
       "unidirectional false positive",
       rep("unidirectional power", 9)))
x <- factor(x, levels(x)[c(2,4,1,3)])
results[, 3] <- x
ggplot(results, aes(y=value, x=scaling.factor, group=type, colour=type, shape=type)) + 
  geom_point(stat="identity", size=4, alpha=.8) + 
  geom_line() +
  theme_bw() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        text=element_text(family="sans", face="plain", color="#000000", size=15, hjust=0.5, vjust=0.5)) + 
  scale_size(range=c(1, 4)) + 
  xlab("scaling factor") + 
  ylab("percent significant") +
  scale_colour_manual(labels = c("bidirectional power", 
                                 "unidirectional power",
                                 "bidirectional false positive",
                                 "unidirectional false positive"),
                      values=c("#377eb8","#e41a1c","#0d0101","#0d0101")) +
  scale_shape_manual(labels = c("bidirectional power", 
                                "unidirectional power",
                                "bidirectional false positive",
                                "unidirectional false positive"),
                     values=c(16,16,0,4))
