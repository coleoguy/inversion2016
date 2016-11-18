source("simFragileY.R")
#             XfA, Xfa, XfAi, Xfai, XmA,  Xma,  XmAi, Xmai,  YA,    Ya,   YAi,  Yai
genotypes <- c(.5, .5,   0,     0,  .25,  .25,    0,   0,   .25,   .25,    0,    0) 

# sz sets the size of the vectors of parameters that we will test
# so 100 means we will test 100 values from x to y of each variable
sz <- 100

# this will be our vector of aneuploidy rates associated with chromosome inversions
u.vec <- seq(from = 0, to = .08, length.out = sz)

# this will be the recombination distance between the SDR and the SA locus
r <- .1

# these are the indices for the chromosomes that we will want to track below
#           XmAi, Xmai, YAi,  Yai
inv.ind <- c(7,    8,    11,  12)

# this is the domminance factor that we will be setting 
# recessive=0, additive=.5, dominant=1
h <- .5

# results is just a list to hold results in
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
      cat("\nitterating")
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

# results.rec, results.add, and results.dom are plotted in figure 2 of the manuscript