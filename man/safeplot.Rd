\name{safeplot}
\alias{safeplot}
\title{ SAFE plot }
\description{
  A SAFE plot for a given category displays the empirical distribution
  function for the ranked local statistics of a given category.
}
\synopsis{
safeplot(safe = NULL, cat.name = NULL, c.vec = NULL, local.stats = NULL, 
         p.val = NULL, one.sided = NA, limits = c(-Inf,Inf), extreme = NA,
         italic = FALSE, x.label = "Ranked local statistic")
}
\usage{
safeplot(safe)
safeplot(safe , cat.name)
safeplot(c.vec=, local.stats= , p.val=, one.sided=, limits=,
         extreme=, italic =, x.label=)
}
\arguments{
  \item{safe}{ Object of class \code{SAFE}. }
  \item{cat.name}{ Name of the category to be plotted. If omitted, the
    most significant category is plotted. }
  \item{c.vec}{ Optional logical vector specifying membership to a gene category. }
  \item{local.stats}{ Optional numeric vector of local statistics. Gene names should 
    be provided as \code{names(local.stats)}.}
  \item{p.val}{ Optional numeric value of the category's empirical p-value}
  \item{one.sided}{ Optional logical value indicating if local statistics are
    one-sided. }
  \item{limits}{ Limits of the shaded region in the plot on the unranked scale. }
  \item{extreme}{ Optional logical value whether only genes in the shaded region
    should be labeled. }
  \item{italic}{ Optional logical value whether gene names should be italic. }
  \item{x.label}{ Character string for the x-axis label. }
}
\details{
  SAFE-plots are suggested as appropriate for visualizing the differential expression
  in a given category relative to the complementary set of genes. The empirical
  cumulative distribution is plotted for the ranked local statistics in the category.
  Tick marks are drawn along the top of the graph to indicate each gene's positions, and
  labeled when sufficient space permits. In this manner, genes with the
  most extreme local statistics can be identified as contributing to a
  categories significance.
}
\references{ W. T. Barry, A. B. Nobel and F.A. Wright, 2005, \emph{Significance Analysis
    of functional categories in gene expression studies: a structured permutation approach},
    \emph{Bioinformatics} \bold{21}(9) 1943--1949. 

  See also the vignette included with this package. }
\author{ William T. Barry: \email{bill.barry@duke.edu} }

\seealso{\code{\link{safe}}.}
\examples{
## Simulate a dataset with 1000 genes and 20 arrays in a 2-sample design.
## The top 100 genes will be differentially expressed at varying levels

g.alt <- 100
g.null <- 900
n <- 20

data<-matrix(rnorm(n*(g.alt+g.null)),g.alt+g.null,n)
data[1:g.alt,1:(n/2)] <- data[1:g.alt,1:(n/2)] + 
                         seq(2,2/g.alt,length=g.alt)
dimnames(data) <- list(c(paste("Alt",1:g.alt),
                         paste("Null",1:g.null)),
                       paste("Array",1:n))

## A treatment vector 
trt <- rep(c("Trt","Ctr"),each=n/2)

## 2 alt. categories and 18 null categories of size 50

C.matrix <- kronecker(diag(20),rep(1,50))
dimnames(C.matrix) <- list(dimnames(data)[[1]],
    c(paste("TrueCat",1:2),paste("NullCat",1:18)))
dim(C.matrix)

results <- safe(data,trt,C.matrix,Pi.mat = 100)
results

## SAFE-plot made for the first category
if (interactive()) { 
safeplot(results,"TrueCat 1")
}
}
\keyword{ hplot }

