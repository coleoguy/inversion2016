source("functions/simFragileY.R")
#             XfA, Xfa, XfAi, Xfai, XmA,  Xma,  XmAi, Xmai,  YA,    Ya,   YAi,  Yai
genotypes <- c(.5, .5,   0,     0,  .25,  .25,    0,   0,   .25,   .25,    0,    0) 

# the number of points to test
sz <- 100

# 3000 is usually safe but lets overkill this to just make sure we really are stable
iter <- 6000

# male and female selection vectors
s.vec.m <- seq(from = 0, to = .3, length.out = sz)
s.vec.f <- seq(from = 0, to = -.3, length.out = sz)

# dominance factor vector
h.vec <- seq(from = 0, to = 1, length.out = sz)

# aneuploidy rate associated with inversion
u <- .01

# recombination distance
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
# The variables results.fa, results.fw, results.ma, and results. mw are plotted
# in figure 3 of manuscript