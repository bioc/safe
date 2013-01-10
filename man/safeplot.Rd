\name{safeplot}
\alias{safeplot}
\title{ SAFE plot }
\description{
  A SAFE plot for a given category displays the empirical distribution function for the ranked (or unranked) local statistics of a given category.
}
\usage{
  safeplot(safe = NULL, cat.name = "", limits = NULL, 
           c.vec = NULL, local.stats = NULL, gene.names = NULL,
           rank = TRUE, x.limits = NULL, c.thresh = 0, 
           colors = NULL, x.ticks = NULL, t.cex = 1, 
           p.val = NULL, cat.desc = NULL, title = "", ...)
}
\arguments{
  \item{safe}{ Object of class \code{SAFE} }
  \item{cat.name}{ Name of the category to be plotted. If omitted, the most significant category is plotted. }
  \item{limits}{ Limits of the shaded region in the plot on the unranked scale }
  \item{c.vec}{ Logical vector specifying membership to a gene category }
  \item{local.stats}{ Numeric vector of local statistics }
  \item{gene.names}{ Optional character vector to replace
    \code{names(local.stats)} in labels }
  \item{rank}{ Logical to plotted raned (TRUE) or unranked (FALSE) local
  statistics on the x-axis }
  \item{x.limits}{ Optional limits of the x-axis. By default will be
    range(local.stats) }
  \item{c.thresh}{ Optional threshold for plotting tickmarks for a weighted
    (``soft'') gene category }
  \item{colors}{ Optional vector specificy colors for gene labels }
  \item{x.ticks}{ Optional location of x-axis tick marks on the ranked scale }
  \item{t.cex}{ Text size for gene labels }
  \item{p.val}{ Optional numeric value of the category's empirical p-value }
  \item{cat.desc}{ Optional character string as a sub-title beneath the
    category name }
  \item{title}{ Optional title to the plot }
  \item{\dots}{ Allows arguments from version 2.0 to be ignored }
}
\details{
  SAFE-plots display the differential expression in a given category relative to the complementary set of genes. The empirical cumulative distribution is plotted for local statistics in the category, on either a ranked or unranked scale. Tick marks are drawn along the top of the graph to indicate each gene's positions, and labeled when sufficient space permits. In this manner, genes with the most extreme local statistics can be identified as contributing to the category's significance.
  
  Typical usages are
  \preformatted{
  safeplot(safe)
  safeplot(safe, cat.name)
  safeplot(c.vec, local.stats, p.val, limits)
}
}
\references{ W. T. Barry, A. B. Nobel and F.A. Wright, 2005, \emph{Significance Analysis of functional categories in gene expression studies: a structured permutation approach}, \emph{Bioinformatics} \bold{21}(9) 1943--1949. 

  See also the vignette included with this package. }
\author{ William T. Barry: \email{bbarry@jimmy.harvard.edu} }

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

results <- safe(data,trt,C.mat = C.matrix,Pi.mat = 100)
results

## SAFE-plot made for the first category
if (interactive()) { 
safeplot(results,"TrueCat 1")
}
}
\keyword{ hplot }