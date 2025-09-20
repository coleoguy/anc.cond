fig_label <- function(text, region="figure", pos="topleft", cex=NULL, ...) {
  
  region <- match.arg(region, c("figure", "plot", "device"))
  pos <- match.arg(pos, c("topleft", "top", "topright", 
                          "left", "center", "right", 
                          "bottomleft", "bottom", "bottomright"))
  
  if(region %in% c("figure", "device")) {
    ds <- dev.size("in")
    # xy coordinates of device corners in user coordinates
    x <- grconvertX(c(0, ds[1]), from="in", to="user")
    y <- grconvertY(c(0, ds[2]), from="in", to="user")
    
    # fragment of the device we use to plot
    if(region == "figure") {
      # account for the fragment of the device that 
      # the figure is using
      fig <- par("fig")
      dx <- (x[2] - x[1])
      dy <- (y[2] - y[1])
      x <- x[1] + dx * fig[1:2]
      y <- y[1] + dy * fig[3:4]
    } 
  }
  
  # much simpler if in plotting region
  if(region == "plot") {
    u <- par("usr")
    x <- u[1:2]
    y <- u[3:4]
  }
  
  sw <- strwidth(text, cex=cex) * 60/100
  sh <- strheight(text, cex=cex) * 60/100
  
  x1 <- switch(pos,
               topleft     =x[1] + sw, 
               left        =x[1] + sw,
               bottomleft  =x[1] + sw,
               top         =(x[1] + x[2])/2,
               center      =(x[1] + x[2])/2,
               bottom      =(x[1] + x[2])/2,
               topright    =x[2] - sw,
               right       =x[2] - sw,
               bottomright =x[2] - sw)
  
  y1 <- switch(pos,
               topleft     =y[2] - sh,
               top         =y[2] - sh,
               topright    =y[2] - sh,
               left        =(y[1] + y[2])/2,
               center      =(y[1] + y[2])/2,
               right       =(y[1] + y[2])/2,
               bottomleft  =y[1] + sh,
               bottom      =y[1] + sh,
               bottomright =y[1] + sh)
  
  old.par <- par(xpd=NA)
  on.exit(par(old.par))
  
  text(x1, y1, text, cex=cex, ...)
  return(invisible(c(x,y)))
}
rep.row<-function(x,n){
  matrix(rep(x,each=n),nrow=n)
}
library(plotfunctions)
library(R.utils)
library(phytools)
library(diversitree)
library(geiger)
##### Fig 1 #####
par(mfrow = c(2,2), mar = c(4,4,6,6) + .1)
load('../../data/ data/Fig1ExampleTree.RData')
load('../../data/Fig1ExampleContTrait.RData')

smp <- contMap(Fig1.example.tree,Fig1.example.cont.trait, ftype = 'off', legend = F, lims = c(.24,2), plot = F)
n<-length(smp$cols)
smp$cols[1:n]<-rainbow(n, end = 4/6)
plot(smp, legend = F,ftype = 'off')
gradientLegend(depth = .03, valRange = c(.24,2), side = 1, pos = .17, color = rainbow(n, end = 4/6))
legend(x = 'bottomleft', legend = '', title = '           Cont. trait value', bg="transparent", bty = 'n')
fig_label('A',cex = 2.5)


load('../../data/Fig1ExampleDiscreteSimmap.RData')
plotSimmap(Fig1.example.disc.simmap, lwd = 4, ftype = 'off')
legend(x = 'bottomleft', legend = c('Ancestral','Derived'), col = c('black', 'red'), pch = 15, bty = 'n')
fig_label('B',cex = 2.5)


pies <- array(dim = c(Fig1.example.disc.simmap$Nnode, 3))
pies[1:4,] <- rep.row(c(1,0,0),4)
pies[5,] <- t(c(0,0,1))
pies[6:7,] <- rep.row(c(0,1,0),2)
pies[8:11,] <- rep.row(c(1,0,0),4)
pies[12,] <- t(c(0,0,1))
pies[13:20,] <- rep.row(c(1,0,0),8)
pies[21,] <- t(c(0,0,1))
pies[22:23,] <- rep.row(c(0,1,0),2)
pies[24:29,] <- rep.row(c(1,0,0),6)
plotSimmap(Fig1.example.disc.simmap, lwd = 4, ftype = 'off')
nodelabels(pie = pies, piecol = c('blue','green', 'red'),cex = .8)
legend(x = 'bottomleft', legend = c('Ancestral','Producing (Ancestral)','Derived'),
       col = c('blue', 'red','green'), pch = 16, bg="transparent", bty = 'n')
fig_label('C',cex = 2.5)

ss_nodes <- Fig1.example.disc.simmap$mapped.edge[, 1] > 0 &
  Fig1.example.disc.simmap$mapped.edge[, 2] > 0
wanted_branches <- ss_nodes[ss_nodes == T]
wanted_nodes <- names(wanted_branches)
wanted_nodes <- gsub(",.*", "", wanted_nodes)
producing.nodes <- unique(wanted_nodes)
anc.states <- anc.ML(Fig1.example.tree, Fig1.example.cont.trait, model = "BM")
orig.val <- mean(anc.states$ace[names(anc.states$ace) %in% producing.nodes])
null.orig.val <- vector(length = 1000)
number.of.trans <- length(producing.nodes)
anc.dt <- Fig1.example.disc.simmap
anc.ct <- anc.states
node.states <- describe.simmap(anc.dt)$states
anc.cond.nodes <- anc.ct$ace[names(anc.ct$ace) %in% names(node.states)[node.states != '2']]
for (j in 1:1000){
  # set.seed(j)
  null.orig.val[j] <- mean(sample(anc.cond.nodes,
                                  length(producing.nodes)))
}
par(mar = c(5,5,1,1) + .1)
plot(density(null.orig.val, bw = .025), ylab = 'Frequency', xlab = 'Mean coninuous trait', main = '')
abline(v = orig.val, col = 'red')
legend(x = 1.01, y=3.4, legend = c('Empirical','Null'), col = c('red', 'black'), pch = 15, bty = 'n')
fig_label('D',cex = 2.5)
