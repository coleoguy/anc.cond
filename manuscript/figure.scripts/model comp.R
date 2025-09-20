results <- as.data.frame(matrix(,2,3))
row.names(results) <- c("power","false.positive")
colnames(results) <- c("pagel","thresh","ancond")
load("../../results/PagelThreshPower.RData")
results[1,1:3] <- c(colSums(pval.array), 50)

load("../../results/PagelThreshFP.RData")
results[2,1:3] <- c(colSums(pval.array), 6)

res <- as.data.frame(cbind(unlist(results),
      rep(c("power","false positive"), times=3),
      rep(c("Pagel's","Threshold","Ancestral condition"), each=2)))
colnames(res) <- c("pos","type","test")
res$pos <- as.numeric(as.character(res$pos))
library(ggraptR)
#ggraptR(res)

ggplot(res, aes(y=pos, x=as.factor(test))) + 
  geom_bar(fill=rep(c("#e41a1c", "#4f8cc0", "#60b45c"),times=2), stat="identity",alpha=.8) + 
  facet_grid(. ~ type) + 
  theme_bw() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        text=element_text(family="sans", face="plain", color="#000000", size=15, hjust=0.5, vjust=0.5)) +
  xlab("") + 
  ylab("rate")
  
# export 7x4
