source("simFragileY.R")
#             XfA, Xfa, XfAi, Xfai, XmA,  Xma,  XmAi, Xmai,  YA,    Ya,   YAi,  Yai
genotypes <- c(.5, .5,   0,     0,  .25,  .25,    0,   0,   .25,   .25,    0,    0) 

# the number of points to test
sz <- 200

# the number of itterations to run some of the setting we are evaluating can
# take as much as 3000 generations to fix/equilibrate
iter <- 3000

# vector of selection coefficients in males
s.vec.m <- seq(from = 0, to = .5, length.out = sz)

# vector of selection coefficients in females
s.vec.f <- seq(from = 0, to = -.5, length.out = sz)

# vector of recombination distances to evaluate
r.vec <- seq(from = 0, to = .5, length.out = sz)
u <- .01

# These are the two beneficial inversions under male benefit female cost
#           XmAi,   Yai
inv.ind.m <- c(7,     12)
#           Xmai,   YAi
inv.ind.f <- c(8,     11)

# the starting frequency to use
init.freq <- .001

# we will be limiting ourselves to additive genetic architecture
h <- .5

results <- list()  
counter <- 1

  # setup result container
  results.mw <- as.data.frame(matrix(, sz, sz))
  colnames(results.mw) <- r.vec
  row.names(results.mw) <- s.vec.m
  results.fa <- results.fw <- results.ma <- results.mw
  for(ix in 1:sz){ #across recombination rates
    cat(ix)
    for(iy in 1:sz){ #across selection coefficients
      
      ## MALE ANTAGONISM
      
      # let system equilibrate
      equi <- simFragileY(genotypes=genotypes, h = h, u = u, 
                          s = s.vec.m[iy], r = r.vec[ix], report = "FATE", 
                          criterion = "STABLE", reporting=1)
      equi.use <- equi
      # introduce mutant chromosome
      equi.use[inv.ind.m[1]] <- init.freq
      #cat("\nitterating")
      results.mw[iy,ix] <- simFragileY(genotypes=equi.use, h = h, u = u, 
                                      s = s.vec.m[iy], r = r.vec[ix], 
                                      report = "FATE", criterion = "GEN", iter = iter,
                                      reporting=1)[inv.ind.m[1]]
      # introduce mutant chromosome
      equi[inv.ind.m[2]] <- init.freq
      results.ma[iy,ix] <- simFragileY(genotypes=equi, h = h, u = u, 
                                      s = s.vec.m[iy], r = r.vec[ix], 
                                      report = "FATE", criterion = "GEN", iter = iter,
                                      reporting=1)[inv.ind.m[2]]
      
      ## FEMALE ANTAGONISM
      
      # let system equilibrate
      equi <- simFragileY(genotypes=genotypes, h = h, u = u, 
                          s = s.vec.f[iy], r = r.vec[ix], report = "FATE", 
                          criterion = "STABLE", reporting=1)
      equi.use <- equi
      # introduce mutant chromosome
      equi.use[inv.ind.f[1]] <- init.freq
      #cat("\nitterating")
      results.fa[iy,ix] <- simFragileY(genotypes=equi.use, h = h, u = u, 
                                       s = s.vec.f[iy], r = r.vec[ix], 
                                       report = "FATE", criterion = "GEN", iter = iter,
                                       reporting=1)[inv.ind.f[1]]
      # introduce mutant chromosome
      equi[inv.ind.f[2]] <- init.freq
      results.fw[iy,ix] <- simFragileY(genotypes=equi, h = h, u = u, 
                                       s = s.vec.f[iy], r = r.vec[ix], 
                                       report = "FATE", criterion = "GEN", iter = iter,
                                       reporting=1)[inv.ind.f[2]]
    }
  }
  results[[1]] <- results.mw
  results[[2]] <- results.ma
  results[[3]] <- results.fw
  results[[4]] <- results.fa
  names(results) <- c("m.wildtype u = .01", "m.antagonistic u = .01",
                      "f.wildtype u = .01", "f.antagonistic u = .01")

result <- as.data.frame(matrix(,4,200))
row.names(result) <- names(results)
colnames(result) <- as.numeric(colnames(results[[1]]))
result[1, ] <- ConvertMatrix(results[[1]],tol=init.freq/3)
result[2, ] <- ConvertMatrix(results[[2]],tol=init.freq)
result[3, ] <- ConvertMatrix(results[[3]],tol=init.freq)
result[4, ] <- ConvertMatrix(results[[4]],tol=init.freq/3)

# if no tested values fixed then ConvertMatrix returns NA
# NA since we are only plotting to .25 set these to .3 so that
# we dont indicate that we achieved fixation at any tested parameter
# space.

result[1,is.na(result[1,])] <- .3
result[2,is.na(result[2,])] <- .3
result[3,is.na(result[3,])] <- .3
result[4,is.na(result[4,])] <- .3


# just for plotting ease lets pull these out
mal.se <- as.numeric(result[1, ])
mal.sa <- as.numeric(result[2, ])
fem.se <- as.numeric(result[3, ])
fem.sa <- as.numeric(result[4, ])

x <- colnames(result)
colvec <- c("#92c5de","#0571b0","#ca0020","#f4a582")

# the variables mal.se, mal.sa, fem.se, and fem.sa are plotted in figure 4 of the manuscript