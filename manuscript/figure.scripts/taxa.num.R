
##### Fig 3 #####

load('../../results/UnidirectionalTaxaFPResults.RData')
load('../../results/UnidirectionalTaxaPowerResults.RData')

x <- seq(from=20, to=200, by=20)
y <- vector()
for(i in 1:10){
  y <- c(y, taxa.uni.power.results[1:100, i])
}
probs <- vector()
for(i in 1:10){
  probs[i] <- sum(taxa.uni.power.results[1:100, i] <= .05)
}
probsfp <- vector()
for(i in 1:10){
  probsfp[i] <- sum(taxa.uni.fp.results[1:100, i] <= .05)
}

load('../../results/BidirectionalTaxaPowerResults.RData')
load('../../results/BidirectionalTaxaFPResults.RData')

y <- vector()
for(i in 1:10){
  for(j in 1:10){
    for(i in 1:100){
      y[2 * (100 * (j - 1) + i - 1) + 1] <- as.numeric(gsub(",.*", "",taxa.bi.power.results[i,j]))
      y[2 * (100 * (j - 1) + i - 1) + 2] <- as.numeric(substr(taxa.bi.power.results[i,j], (nchar(gsub(",.*", "",taxa.bi.power.results[i,j])) + 2), 
                                                              nchar(taxa.bi.power.results[i,j])))
    }
  }
}
biprobs <- vector()
for(i in 1:10){
  biprobs[i] <- 
      round(100 * sum(y[(200 * (i - 1) + 1):(200 * i)] <= .025, na.rm = T) / sum(!is.na(y[(200 * (i - 1) + 1):(200 * i)])),
            digits = 0)
}

y <- vector()
for(i in 1:10){
  for(j in 1:10){
    for(i in 1:100){
      y[2 * (100 * (j - 1) + i - 1) + 1] <- as.numeric(gsub(",.*", "",taxa.bi.fp.results[i,j]))
      y[2 * (100 * (j - 1) + i - 1) + 2] <- as.numeric(substr(taxa.bi.fp.results[i,j], (nchar(gsub(",.*", "",taxa.bi.fp.results[i,j])) + 2), 
                                                              nchar(taxa.bi.fp.results[i,j])))
    }
  }
}
biprobsfp <- vector()
for(i in 1:10){
  biprobsfp[i] <- round(100 * sum(y[(200 * (i - 1) + 1):(200 * i)] <= .025, na.rm = T) / 
                          sum(!is.na(y[(200 * (i - 1) + 1):(200 * i)])),
            digits = 0)
  
}



up <- probs/100
bp <- biprobs/100
ufp <- probsfp/100
bfp <- biprobsfp/100



results <- as.data.frame(matrix(,40,3))
colnames(results) <- c("value", "taxa.number", "type")
results[, 1] <- c(up,bp,ufp,bfp)
results[, 2] <- rep(c(20,40,60,80,100,120,140,160,180,200), 4)
x <- factor(c(rep("unidirectional power", 10),
              rep("bidirectional power",10),
              rep("unidirectional false positive",10),
              rep("bidirectional false positive", 10)))
x <- factor(x, levels(x)[c(2,4,1,3)])
results[, 3] <- x

ggplot(results, aes(y=value, x=taxa.number, group=type, colour=type, shape=type)) + 
  geom_hline(yintercept = .05, linetype="dashed") +
geom_point(stat="identity", size=4, alpha=.8) + 
  geom_line() +
  theme_bw() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        text=element_text(family="sans", face="plain", color="#000000", size=15, 
                          hjust=0.5, vjust=0.5), legend.text=element_text(size=rel(0.6))) + 
  scale_size(range=c(1, 4)) + 
  xlab("number of taxa") + 
  ylab("percent significant") +
  scale_colour_manual(labels = c("bidirectional power", 
                                 "unidirectional power",
                                 "bidirectional false positive",
                                 "unidirectional false positive"),
                      values=c("#377eb8","#e41a1c","#377eb8","#e41a1c")) +
  scale_shape_manual(labels = c("bidirectional power", 
                                "unidirectional power",
                                "bidirectional false positive",
                                "unidirectional false positive"),
                     values=c(16,16,1,1)) 
# export at 6.5x4