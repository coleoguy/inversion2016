source("simFragileY.R")
#             XfA, Xfa, XfAi, Xfai, XmA,  Xma,  XmAi, Xmai,  YA,    Ya,   YAi,  Yai
genotypes <- c(.5, .5,   0,     0,  .25,  .25,    0,   0,   .25,   .25,    0,    0) 
sz <- 100
u.vec <- seq(from = 0, to = .08, length.out = sz)
r <- .1
inv.ind <- c(7,    8,    11,  12)
#           XmAi, Xmai, YAi,  Yai

# iter.f <- 3000
# iter.e <- 10000

h <- .5
results <- list()  
for(j in 1:4){
  cat("\nadditive", j)
  s.vec <- seq(from = 0, to = .25, length.out = sz)
  results[[j]] <- as.data.frame(matrix(,sz,sz))
  colnames(results[[j]]) <- u.vec
  row.names(results[[j]]) <- s.vec
  if(j==2 | j==3) s.vec <- -1*s.vec
  for(ix in 1:sz){ #across aneuploidy rates
    if(ix %% 5 == 0) cat(", ix")
    for(iy in 1:sz){ #across selection coefficients
      # let system equilibrate
      #cat("\nequilibrate")
      equi <- simFragileY(genotypes=genotypes, h = h, u = u.vec[ix], s = s.vec[iy],
                          r = r, report = "FATE", criterion = "STABLE", reporting=1)
      # introduce rare mutation type equals 3,4,7,8,11,12
      # XfA, Xfa, XfAi, Xfai, XmA,  Xma,  XmAi, Xmai,  YA,    Ya,   YAi,  Yai
      equi[inv.ind[j]] <- .005
      #cat("\nitterating")
      results[[j]][iy,ix] <- simFragileY(genotypes=equi, h = h, u = u.vec[ix], s = s.vec[iy],
                                         r = r, report = "FATE", criterion = "STABLE", reporting=1)[inv.ind[j]]
    }
  }
}
names(results) <- c("XAi", "Xai", "YAi", "Yai")
results.add <- results



cat("\ndominance")
h <- 1
results <- list()  
for(j in 1:4){
  cat(j)
  s.vec <- seq(from = 0, to = .25, length.out = sz)
  results[[j]] <- as.data.frame(matrix(,sz,sz))
  colnames(results[[j]]) <- u.vec
  row.names(results[[j]]) <- s.vec
  if(j==2 | j==3) s.vec <- -1*s.vec
  for(ix in 1:sz){ #across aneuploidy rates
    for(iy in 1:sz){ #across selection coefficients
      # let system equilibrate
      #cat("\nequilibrate")
      equi <- simFragileY(genotypes=genotypes, h = h, u = u.vec[ix], s = s.vec[iy],
                          r = r, report = "FATE", criterion = "STABLE", reporting=1)
      # introduce rare mutation type equals 3,4,7,8,11,12
      # XfA, Xfa, XfAi, Xfai, XmA,  Xma,  XmAi, Xmai,  YA,    Ya,   YAi,  Yai
      equi[inv.ind[j]] <- .005
      #cat("\nitterating")
      results[[j]][iy,ix] <- simFragileY(genotypes=equi, h = h, u = u.vec[ix], s = s.vec[iy],
                                         r = r, report = "FATE", criterion = "STABLE", reporting=1)[inv.ind[j]]
    }
  }
}
names(results) <- c("XAi", "Xai", "YAi", "Yai")
results.dom <- results




cat("\nrecessive")
h <- 0
results <- list()  
for(j in 1:4){
  cat(j)
  s.vec <- seq(from = 0, to = .25, length.out = sz)
  results[[j]] <- as.data.frame(matrix(,sz,sz))
  colnames(results[[j]]) <- u.vec
  row.names(results[[j]]) <- s.vec
  if(j==2 | j==3) s.vec <- -1*s.vec
  for(ix in 1:sz){ #across aneuploidy rates
    for(iy in 1:sz){ #across selection coefficients
      # let system equilibrate
      #cat("\nequilibrate")
      equi <- simFragileY(genotypes=genotypes, h = h, u = u.vec[ix], s = s.vec[iy],
                          r = r, report = "FATE", criterion = "STABLE", reporting=1)
      # introduce rare mutation type equals 3,4,7,8,11,12
      # XfA, Xfa, XfAi, Xfai, XmA,  Xma,  XmAi, Xmai,  YA,    Ya,   YAi,  Yai
      equi[inv.ind[j]] <- .005
      #cat("\nitterating")
      results[[j]][iy,ix] <- simFragileY(genotypes=equi, h = h, u = u.vec[ix], s = s.vec[iy],
                                         r = r, report = "FATE", criterion = "STABLE", reporting=1)[inv.ind[j]]
    }
  }
}
names(results) <- c("XAi", "Xai", "YAi", "Yai")
results.rec <- results


image(t(results[[4]])) 

# save as par-space.RData

rm(has, result,cols,i,x,x.leg,xvals,y,y.leg,y2,ConvertMatrix, equi,genotypes,h,inv.ind, ix,iy,j,r,results,s.vec,u.vec,simFragileY)
