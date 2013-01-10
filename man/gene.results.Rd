\name{gene.results}
\alias{gene.results}
\title{ Gene-specific results from SAFE }
\description{
  Provides gene-specific local statistics and resampling-based p-values for every feature in the category of interest. Features are ordered by the degree and direction of differential expression.
}
\usage{
gene.results(object, cat.name = NULL, error = "none", 
             print.it = TRUE, gene.names = NULL)
}
\arguments{
  \item{object}{ Object of class \code{SAFE}. }
  \item{cat.name}{ Name of the category to be displayed. If omitted, the most significant category is displayed. }
  \item{error}{ Specifies a non-resampling based method for adjusting the empirical p-values.  A Bonferroni, ("FWER.Bonf"), Holm's step-up ("FWER.Holm"), and Benjamini-Hochberg step down ("FDR.BH") adjustment can be selected.  By default ("none") no error rates are computed. }
  \item{print.it}{ Logical determining whether results are printed to screen or returned as a list of data.frames for up- and down-regulated genes. }
  \item{gene.names}{ Optional character vector of gene names to append to the SAFE output. }
}
\references{ W. T. Barry, A. B. Nobel and F.A. Wright, 2005, \emph{Significance Analysis of functional categories in gene expression studies: a structured permutation approach}, \emph{Bioinformatics}\bold{21}(9) 1943--1949. 

  See also the vignette included with this package. }
\author{ William T. Barry: \email{bbarry@jimmy.harvard.edu} }

\seealso{\code{\link{safe}}.}
\keyword{ htest }