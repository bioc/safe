\name{safe.toptable}
\alias{safe.toptable}
\title{ Category-specific results from SAFE }
\description{
  Prints annotated results from \code{SAFE} as a data.frame for categories with the strongest association to response/phenotype. This includes (a) category name; (b) category size; (c) global statistic; (d) nominal empirical p-values; (e) adjusted p-values; and (f) descriptions if available.
}
\usage{
safe.toptable(safe, number = 10, 
              pretty = TRUE, description = TRUE)
}
\arguments{
  \item{safe}{ Object of class \code{SAFE}. }
  \item{number}{ Number of categories to be printed. If omitted, the first 10 categories are reported. }
  \item{pretty}{ Logical determining whether p-values smaller than 10^-4 are output in scientific notation, rather than as decimals. By default \code{pretty = TRUE}. }
  \item{description}{ Logical determining whether category descriptions are appended to printed output. By default \code{description = TRUE}. }
}
\references{ W. T. Barry, A. B. Nobel and F.A. Wright, 2005, \emph{Significance Analysis of functional categories in gene expression studies: a structured permutation approach}, \emph{Bioinformatics}\bold{21}(9) 1943--1949. 

  See also the vignette included with this package. }
\author{ William T. Barry: \email{bbarry@jimmy.harvard.edu} }

\seealso{\code{\link{safe}}.}
\keyword{ htest }
