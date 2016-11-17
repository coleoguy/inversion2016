simFragileY <- function(genotypes, h, s, r, u, iter=100, report="FULL", criterion="GEN", d=6, reporting=1){
  ## genotypes
  # vector of length 12 that contains starting frequencies
  
  ## h
  # dominance term
  
  ## s
  # selection coefficient
  
  ## r
  # probability of recombination event between sex determining locus and the 
  # sexually antagonistic locus
  
  ## u
  # aneuploidy rate in males that carry inversion
  
  ## iter
  # number of itterations to run for
  
  ## report
  # "FULL" gives matrix of states at every generation
  # "FATE" gives final frequencies
  
  ## criterion
  # "GEN" runs for set number of generations
  # "STABLE" runs till allele frequincies are unchanging
  
  ## d
  # digits to look at to determine whether allele freq is changing
  
  ## report
  # 1  minimal reporting
  # 2 moderte reporting
  # 3 troubleshooting
  
  ## Set up internal variables
  # egg frequencies
  XfA <-  genotypes[1]
  Xfa <-  genotypes[2]
  XfAi <- genotypes[3]
  Xfai <- genotypes[4]
  # Sperm frequencies
  XmA <-  genotypes[5]
  Xma <-  genotypes[6]
  XmAi <- genotypes[7]
  Xmai <- genotypes[8]
  YA <-   genotypes[9]
  Ya <-   genotypes[10]
  YAi <-  genotypes[11]
  Yai <-  genotypes[12]
  # selection coefficients
  w1 <- 1
  w2 <- 1 - h * s
  w3 <- 1 - s
  w4 <- 1 + h * s
  w5 <- 1 + s
  
  ## setup data for output at end
  if(report=="FULL"){
    results <- as.data.frame(matrix(genotypes, 12, 1))
  }
  if(report == "FATE"){
    results <- vector()
  }
  
  continue <- TRUE
  counter <- 1
  if(reporting > 1) cat("\nCalculating change in generation: ")
  while(continue == TRUE){
    if(reporting > 1) cat(counter, " ")
    ## Eggs
    FD <- XfA*XmA*w1 + XfA*Xma*w2 + Xfa*XmA*w2 + Xfa*Xma*w3 + XfAi*XmA*w1 + 
      XfA*XmAi*w1 + XfAi*XmAi*w1 + XfAi*Xma*w2 + XfA*Xmai*w2 + Xfai*XmA*w2 + 
      Xfa*XmAi*w2 + Xfai*Xma*w3 + Xfa*Xmai*w3 + Xfai*Xmai*w3
    
    XfA2 <- (XfA*XmA*w1 + .5*XfA*Xma*w2 + .5*Xfa*XmA*w2 + .5*XfAi*XmA*w1 +
               .5*XfA*XmAi*w1 + .5*XfA*Xmai*w2 + .5*Xfai*XmA*w2) / FD
    
    Xfa2 <- (.5*XfA*Xma*w2 + .5*Xfa*XmA*w2 + Xfa*Xma*w3 + .5*XfAi*Xma*w2 +
               .5*Xfa*XmAi*w2 + .5*Xfai*Xma*w3 + .5*Xfa*Xmai*w3) / FD
    
    XfAi2 <- (.5*XfAi*XmA*w1 + .5*XfA*XmAi*w1 + XfAi*XmAi*w1 + .5*XfAi*Xma*w2 +
                .5*Xfa*XmAi*w2) / FD
    
    Xfai2 <- (.5*XfA*Xmai*w2 + .5*Xfai*XmA*w2 + .5*Xfai*Xma*w3 +
                .5*Xfa*Xmai*w3 + Xfai*Xmai*w3) / FD
    
    ## Sperm
    MD <- XfA*YA*w1 + XfA*Ya*w4 + Xfa*YA*w4 + Xfa*Ya*w5 + XfAi*YA*(w1-u) + 
      XfA*YAi*(w1-u) + XfAi*Ya*(w4-u) + XfA*Yai*(w4-u) + Xfai*YA*(w4-u) +
      Xfa*YAi*(w4-u) + Xfai*Ya*(w5-u) + Xfa*Yai*(w5-u)
    
    XmA2 <- (.5*XfA*YA*w1 + .5*(1-r)*XfA*Ya*w4 + .5*r*Xfa*YA*w4 + 
               .5*XfA*YAi*(w1-u) + .5*XfA*Yai*(w4-u)) / MD
    
    Xma2 <- (.5*r*XfA*Ya*w4 + .5*(1-r)*Xfa*YA*w4 + .5*Xfa*Ya*w5 + 
               .5*Xfa*YAi*(w4-u) + .5*Xfa*Yai*(w5-u)) / MD
    
    XmAi2 <- (.5*XfAi*YA*(w1-u) + .5*XfAi*Ya*(w4-u)) / MD
    
    Xmai2 <- (.5*Xfai*YA*(w4-u) + .5*Xfai*Ya*(w5-u)) / MD
    
    YA2 <- (.5*XfA *YA*w1 + .5*r*XfA*Ya*w4 + .5*(1-r)*Xfa*YA*w4 + 
              .5*XfAi*YA*(w1-u) + .5*Xfai*YA*(w4-u)) / MD
    
    Ya2 <- (.5*(1-r)*XfA*Ya*w4 + .5*r*Xfa*YA*w4 + .5*Xfa*Ya*w5 + 
              .5*XfAi*Ya*(w4-u) + .5*Xfai*Ya*(w5-u)) / MD
    
    YAi2 <- (.5*XfA*YAi*(w1-u) + .5*Xfa*YAi*(w4-u)) / MD
    
    Yai2 <- (.5*XfA*Yai*(w4-u) + .5*Xfa*Yai*(w5-u)) / MD
    
     ## record results if every generation is requested
    if(report=="FULL"){
      results[1:12, counter] <- c(XfA, Xfa, XfAi, Xfai, XmA, Xma,
                                  XmAi, Xmai, YA, Ya, YAi, Yai)
    }
# iterated for generations
    counter <- counter + 1
    # check to see if stopping conditions have been met
    if(criterion == "GEN" & counter > iter) continue <- FALSE
    if(criterion == "STABLE"){
      if(round(XfA, digits=d) == round(XfA2, digits = d) & 
         round(Xfa, digits=d) == round(Xfa2, digits = d) & 
         round(XfAi, digits=d) == round(XfAi2, digits = d) & 
         round(Xfai, digits=d) == round(Xfai2, digits = d) &
         round(XmA, digits=d) == round(XmA2, digits = d) & 
         round(Xma, digits=d) == round(Xma2, digits = d) & 
         round(XmAi, digits=d) == round(XmAi2, digits = d) & 
         round(Xmai, digits=d) == round(Xmai2, digits = d) &
         round(YA, digits=d) == round(YA2, digits = d) & 
         round(Ya, digits=d) == round(Ya2, digits = d) & 
         round(YAi, digits=d) == round(YAi2, digits = d) & 
         round(Yai, digits=d) == round(Yai2, digits = d)){
        continue <- FALSE
      }
    }    
    
    ## pass data to next generation
    ## Eggs
    XfA <- XfA2
    Xfa <- Xfa2
    XfAi <- XfAi2
    Xfai <- Xfai2
    ## Sperm
    XmA <- XmA2 
    Xma <- Xma2
    XmAi <- XmAi2 
    Xmai <- Xmai2 
    YA <- YA2 
    Ya <- Ya2 
    YAi <- YAi2 
    Yai <- Yai2
    
    
    
  }
  ## setup data for output at end
  if(report=="FATE"){
    results <- c(XfA2, Xfa2, XfAi, Xfai2, XmA2, Xma2,
                 XmAi2, Xmai2, YA2, Ya2, YAi2, Yai2)
    names(results) <- c("XfA", "Xfa", "XfAi", "Xfai", "XmA", "Xma",
                        "XmAi", "Xmai", "YA", "Ya", "YAi", "Yai")
    
  }
  if(report=="FULL"){
    row.names(results) <- c("XfA", "Xfa", "XfAi", "Xfai", "XmA", "Xma",
                                "XmAi", "Xmai", "YA", "Ya", "YAi", "Yai")
  }
  
  return(results)
}



# This function lets us process the output from above

ConvertMatrix <- function(in.matrix, tol){
  Slide <- function (FUN, data, window, step){
    total <- length(data)
    spots <- seq(from = 1, to = (total - (window - 1)), by = step)
    result <- vector(length = length(spots))
    for (i in 1:length(spots)) {
      result[i] <- match.fun(FUN)(data[spots[i]:(spots[i] + window - 1)])
    }
    return(result)
  }
  result <- vector(mode="numeric", length=ncol(in.matrix))
  from <- as.numeric(row.names(in.matrix))[1]
  to <- as.numeric(row.names(in.matrix))[nrow(in.matrix)]
  
  x <- seq(from, to, length.out=nrow(in.matrix))
  vals <-Slide(mean, x, 2, 1)
  for(i in 1:length(result)){
    # this gives the ones that increased
    increased <- in.matrix[,i]>tol
    
    # all selection coefficients caused increase so min = 0
    if(sum(increased) == nrow(in.matrix)) result[i] <- 0
    
    # no selection coefficient was sufficient to allow invasion
    if(sum(increased) == 0) result[i] <- NA
    
    # there is a mix
    if(sum(increased) != 0 & sum(increased) != nrow(in.matrix)){
      # gets the value where they flip
      result[i] <- vals[max(which(in.matrix[,i] <= tol))]
    }
  }
  result
}
