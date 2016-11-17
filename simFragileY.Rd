\name{simFragileY}
\alias{simFragileY}
\title{Itteration of sex chromosome model}
\description{This function allows for the definition and itteration of a set of sex chromosome models based on the fragile Y hypothesis.}
\usage{simFragileY(genotypes, h, s, r, u, criterion="GEN", iter=100, report="FULL", d=6, reporting=1}
\arguments{
  \item{genotypes}{a vector of length 12 with starting chromosome frequencies. Frequencies should be in the order: XfA, Xfa, XfAi, Xfai, XmA, Xma, XmAi, Xmai, YA, Ya, YAi, Yai}
  \item{h}{dominance factor 0 to 1}
  \item{s}{selection coefficient associated with the A locus male fitness is 1, 1+hs, 1+s, and female fitness is 1, 1-hs, 1-s for genotypes AA, Aa, and aa respectively.}
  \item{r}{recombination distance between SDR and A locus 0 to .5}
  \item{u}{aneuploidy rate in males that carry an inversion}
  \item{criterion}{indicates if equations should be itterated for a specified number of generations "GEN" or should be itterated tills stable "STABLE" see d below.}
  \item{iter}{number of generations to run if criterion is set to "GEN"}
  \item{report}{if "FULL" then a matrix is returned with each chromosome frequency in every generation if "FATE" then vector of length 12 is returned with frequencies in the final generation.}
  \item{d}{threshold for detecting change in chromosome frequency. Frequencies are rounded to the nearest 10-d}
  \item{reporting}{level of messages to return default is 1 can be set to 2 or 3 for more detailed status and troublshooting}
  }
\details{provides either a matrix or vector of chromosome frequencies in gametes.
}
\author{
Heath Blackmon}
\examples{
}
