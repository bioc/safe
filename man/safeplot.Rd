\name{safeplot}
\alias{safeplot}
\title{ SAFE plot }
\description{
  A SAFE plot for a given category displays the distribution of ranked
  local statistics for the genes in the category against their complement.
}
\synopsis{
safeplot(safe = NULL, cat.name = NULL, c.vec = NULL, local.stats = NULL, 
         p.val = NULL, one.sided = FALSE, limits = c(-Inf,Inf), extreme = NA,
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
  \item{cat.name}{ Name of the category to be plotted. }
  \item{c.vec}{ Optional numeric vector specifying the gene category. }
  \item{local.stats}{ Optional numeric vector of local statistics. Gene names should 
    be provided as \code{names(local.stats)}.}
  \item{p.val}{ Optional numeric value of the category's empirical p-value}
  \item{one.sided}{ Optional logical value indicating if local statistics are
    one-sided. }
  \item{limits}{ Limits of the shaded region in the plot on the unranked scale. }
  \item{extreme}{ Optional logical value whether only genes in the shaded region
    should be labeled. }
  \item{italic}{ Optional logical value whether gene names should be italic. }
  \item{x.label}{ Optional character string for the x-axis label. }
}
\details{
  SAFE-plots are suggeted as appropriate for visualizing the differential expression
  in a given category relative to the complementary set of genes. The empirical
  cumulative distribution is plotted for the ranked local statistics in the category.
  Tick marks are drawn along the top of the graph to indicate the genes' positions, and
  labeled. In this manner, genes with the most extreme local statistics can be
  identified in the plot.
}
\references{ W. T. Barry, A. B. Nobel and F.A. Wright, 2004, \emph{Significance Analysis
    of functional categories in gene expression studies: a structured permutation approach},
    \emph{Bioinformatics} In press. 

  See also the vignette included with this package. }
\author{ William T. Barry: \email{wbarry@bios.unc.edu} }

\seealso{{\code{\link{safe}}.}}
\examples{
## Consider a dataset with 1000 genes and 20 arrays in a 2-sample design.
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

## A treatment vector is also made
trt <- rep(c("Trt","Ctr"),each=n/2)
trt

## 2 alternative catagories and  18 null categories
## will be made of 50 null genes. 

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

