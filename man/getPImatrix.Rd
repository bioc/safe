\name{getPImatrix}
\alias{getPImatrix}
\title{ Generation of a Pi matrix }
\description{
  A function to generate random permutations of the data. Permutations are stored as
  a matrix with each row containing the permuted order as the reordered integers from
  1 to n.
}
\synopsis{
getPImatrix(n = NULL, y.vec = NULL, block.vec = NULL, K, method="random")
}
\usage{
getPImatrix(n, K= ,method=)
getPImatrix(y.vec= , K=)
getPImatrix(block.vec=, K=)
}
\arguments{
  \item{n}{ Number of samples or length of the response vector }
  \item{y.vec}{ Response vector for a SAFE analysis }
  \item{block.vec}{ A character or numeric vector containing the blocking assignments
    (see Details:).}
  \item{K}{ Number of permutations to generate }
  \item{method}{ For \code{method} = "random" each permutation is a random reordering of
    1 through n (with or without blocks). Future versions will allow \code{method} = "exaust"
    to give all unique permuations for a moderate number of arrays. }
}
\details{
  \code{getPImatrix} is provided as a separate function from \code{safe}, so that a
  matrix of permutation assignments can be created and stored to reproduce identical
  results. In some experimental designs it may be appropriate to generate blocked 
  permutaions of the samples, that hold another effect constant. By specifying
  \code{block.vec}, only the column numbers sharing a unique value in the vector will
  be permuted. 
}
\value{
  A matrix is returned with each row containing a permutation of the integers from 1
  to n. The first row is restricted to the ordered integers to represent the 
  observed data.
}
\references{ W. T. Barry, A. B. Nobel and F.A. Wright, 2004, \emph{Significance Analysis
  of functional categories in gene expression studies: a structured permutation approach},
  \emph{Bioinformatics} In press. 

  See also the vignette included with this package. }
\author{ William T. Barry: \email{wbarry@bios.unc.edu} }

\seealso{{\code{\link{safe}}.}}
\examples{
## Random permutations for 8 samples
getPImatrix( n = 8, K = 5)

## Random permutations within blocks 
block <- rep(0:1,each=4)
block

getPImatrix(block.vec = block, K = 5)
}
\keyword{ htest }
