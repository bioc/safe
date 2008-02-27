\name{gene.results}
\alias{gene.results}
\title{ Gene-specific results from SAFE }
\description{
  Prints gene-specific local statistics and resampling-based p-values
  for every probeset in the gene category of interest.  Probesets
  are ordered by the degree and direction of differential expression.
}
\usage{
gene.results(object = NULL, cat.name = NULL, error = "none", print.it = TRUE)
}
\arguments{
  \item{object}{ Object of class \code{SAFE}. }
  \item{cat.name}{ Name of the category to be plotted. If omitted, the
    most significant category is plotted. }
  \item{error}{ Specifies a non-resampling based method for adjusting
    the empirical p-values.  A Bonferroni, ("FWER.Bonf"), Holm's step-up
    ("FWER.Holm"), and Benjamini-Hochberg step down ("FDR.BH")
    adjustment can be selected. By default ("none") no error rates are
    computed. }
  \item{print.it}{ Logical determining whether results are printed to
    screen or returned as a list of two data.frames. }
}
\references{ W. T. Barry, A. B. Nobel and F.A. Wright, 2005, \emph{Significance Analysis
    of functional categories in gene expression studies: a structured permutation approach},
    \emph{Bioinformatics} {\bf 21}(9) 1943--1949. 

  See also the vignette included with this package. }
\author{ William T. Barry: \email{bill.barry@duke.edu} }

\seealso{{\code{\link{safe}}.}}
\keyword{ hplot }
