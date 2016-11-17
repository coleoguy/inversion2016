
###############################################################
## This file will run iterations to look at the effect of
## genetic architecture
##############################################################
setwd("~/Desktop/Dropbox/projects/fragileY/scripts")
source("functions/simFragileY.R")
#             XfA, Xfa, XfAi, Xfai, XmA,  Xma,  XmAi, Xmai,  YA,    Ya,   YAi,  Yai
genotypes <- c(.5, .5,   0,     0,  .25,  .25,    0,   0,   .25,   .25,    0,    0) 
sz <- 100
iter <- 6000
s.vec.m <- seq(from = 0, to = .3, length.out = sz)
s.vec.f <- seq(from = 0, to = -.3, length.out = sz)
h.vec <- seq(from = 0, to = 1, length.out = sz)
u <- .01
r<-.1

# These are the two beneficial inversions under male benefit female cost
#           XmAi,   Yai
inv.ind.m <- c(7,     12)
#           Xmai,   YAi
inv.ind.f <- c(8,     11)
init.freq <- .001


results <- list()  
counter <- 1

  # setup result container
  results.mw <- as.data.frame(matrix(,sz,sz))
  colnames(results.mw) <- h.vec
  row.names(results.mw) <- s.vec.m
  results.fa <- results.fw <- results.ma <- results.mw
  for(ix in 1:sz){ #across recombination rates
    cat(ix)
    for(iy in 1:sz){ #across selection coefficients
      
      ## MALE ANTAGONISM
      
      # let system equilibrate
      equi <- simFragileY(genotypes=genotypes, h = h.vec[ix], u = u, 
                          s = s.vec.m[iy], r = r, report = "FATE", 
                          criterion = "STABLE", reporting=1)
      equi.use <- equi
      # introduce mutant chromosome
      equi.use[inv.ind.m[1]] <- init.freq
      #cat("\nitterating")
      results.mw[iy,ix] <- simFragileY(genotypes=equi.use, h = h.vec[ix], u = u, 
                                      s = s.vec.m[iy], r = r, 
                                      report = "FATE", criterion = "GEN", iter = iter,
                                      reporting=1)[inv.ind.m[1]]
      # introduce mutant chromosome
      equi[inv.ind.m[2]] <- init.freq
      results.ma[iy,ix] <- simFragileY(genotypes=equi, h = h.vec[ix], u = u, 
                                      s = s.vec.m[iy], r = r, 
                                      report = "FATE", criterion = "GEN", iter = iter,
                                      reporting=1)[inv.ind.m[2]]
      
      ## FEMALE ANTAGONISM
      
      # let system equilibrate
      equi <- simFragileY(genotypes=genotypes, h = h.vec[ix], u = u, 
                          s = s.vec.f[iy], r = r, report = "FATE", 
                          criterion = "STABLE", reporting=1)
      equi.use <- equi
      # introduce mutant chromosome
      equi.use[inv.ind.f[1]] <- init.freq
      #cat("\nitterating")
      results.fa[iy,ix] <- simFragileY(genotypes=equi.use, h = h.vec[ix], u = u, 
                                       s = s.vec.f[iy], r = r, 
                                       report = "FATE", criterion = "GEN", iter = iter,
                                       reporting=1)[inv.ind.f[1]]
      # introduce mutant chromosome
      equi[inv.ind.f[2]] <- init.freq
      results.fw[iy,ix] <- simFragileY(genotypes=equi, h = h.vec[ix], u = u, 
                                       s = s.vec.f[iy], r = r, 
                                       report = "FATE", criterion = "GEN", iter = iter,
                                       reporting=1)[inv.ind.f[2]]
    }
  }
  library(viridis)
  library(lattice)
  results.mw2 <- results.mw
  for(i in 1:sz) results.mw2[i, results.mw2[i, ]<.001] <- 0
  levelplot(t(results.mw2), col.regions=viridis(100),zlim=c(0,.5))
 
  
  layout(mat=matrix(c(1,2,3,4,5,5),2,3), widths=c(5,5,1))
  
  image(t(results.mw),zlim=c(0,.5), col=viridis(100),main="XA", axes=F, xlab="dominance factor", ylab="selection coefficient")
  axis(side=1,at=seq(0,1,length.out=6), labels=seq(0,1,length.out=6))
  axis(side=2,at=seq(0,1,length.out=6), labels=seq(0,.3,length.out=6))
  
  image(t(results.ma),zlim=c(0,.5), col=viridis(100),main="Ya", axes=F, xlab="dominance factor", ylab="selection coefficient")
  axis(side=1,at=seq(0,1,length.out=6), labels=seq(0,1,length.out=6))
  axis(side=2,at=seq(0,1,length.out=6), labels=seq(0,.3,length.out=6))
  
  image(t(results.fa),zlim=c(0,.5), col=viridis(100),main="Xa", axes=F, xlab="dominance factor", ylab="selection coefficient")
  axis(side=1,at=seq(0,1,length.out=6), labels=seq(0,1,length.out=6))
  axis(side=2,at=seq(0,1,length.out=6), labels=seq(0,.3,length.out=6))
  
  image(t(results.fw),zlim=c(0,.5), col=viridis(100),main="YA", axes=F, xlab="dominance factor", ylab="selection coefficient")
  axis(side=1,at=seq(0,1,length.out=6), labels=seq(0,1,length.out=6))
  axis(side=2,at=seq(0,1,length.out=6), labels=seq(0,.3,length.out=6))
  
  vals <- round(seq(from=0, to=1, length.out=4), digits=3)
  legend_image <- as.raster(matrix(viridis(100)[100:1], ncol=1))
  par(mar=c(0,0,0,0))
  plot(c(0,2),c(0,1),type = 'n', axes = F,xlab = '', ylab = '', main = "")
  text(x=.8, y = seq(from=.75, to=.95, length.out=4), labels = vals, pos=4)
  rasterImage(legend_image, .5, .75, .75,.95)
  #xleft, ybottom, xright, ytop
  
  
  results[[1]] <- results.mw
  results[[2]] <- results.ma
  results[[3]] <- results.fw
  results[[4]] <- results.fa
  names(results) <- c("m.wildtype u = .01", "m.antagonistic u = .01",
                      "f.wildtype u = .01", "f.antagonistic u = .01")


  