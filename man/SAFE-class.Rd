\name{SAFE-class}
\docType{class}
\alias{SAFE-class}
\alias{show,SAFE-method}
\alias{[,SAFE-method}
\title{Class SAFE}
\description{
  The class SAFE is the output from the function \code{\link{safe}}. It is also the input to the plotting function \code{\link{safeplot}}.
}
\section{Slots}{
  \describe{
    \item{\code{local}:}{Object of class \code{"character"} indicating the local statistic used.}
    \item{\code{local.stat}:}{Object of class \code{"numeric"} containing the (unsorted) observed local statistics for genes.}
    \item{\code{local.pval}:}{Object of class \code{"numeric"} containing the (unsorted) empirical p-values for genes}
    \item{\code{global}:}{Object of class \code{"character"} indicating the global statistic used.}
    \item{\code{global.stat}:}{Object of class \code{"numeric"} containing the (unsorted) observed global statistics for categories.}
    \item{\code{global.pval}:}{Object of class \code{"numeric"} containing the (unsorted) empirical p-values for categories.}
    \item{\code{error}:}{Object of class \code{"character"} indicating the method used to estimate error rates across multiple comparisons.}
    \item{\code{global.error}:}{Object of class \code{"numeric"} containing the (unsorted) error rates (adjusted p-values) for categories. If not computed, it will be set to \code{NA}.}
    \item{\code{C.mat}:}{Object of class \code{"SparseM"} containing the category assignments.}
    \item{\code{alpha}:}{Object of class \code{"numeric"} containing the alpha level used in returning significant results.}
    \item{\code{method}:}{Object of class \code{"character"} indicating the resampling method used in \code{safe}.}
  }
}
\section{Methods}{ 
  \describe{
    \item{show}{(safe.results): Summarizes the test results of significant categories.}
    \item{[}{(safe.results): Returns a SAFE object for categories indicated by integer or character strings.}
    \item{safeplot}{(safe.results): \code{\link{safeplot}} produces CDFs of the association of expression to phenotype in a category relative to its complement.}
  }
}

\author{William T Barry: \email{bbarry@jimmy.harvard.edu }}

\seealso{\code{\link{safe}}, \code{\link{safeplot}}.}

\keyword{methods}
