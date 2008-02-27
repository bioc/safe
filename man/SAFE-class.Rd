\name{SAFE-class}
\docType{class}
\alias{SAFE-class}
\alias{show,SAFE-method}
\alias{[,SAFE-method}
\title{Class SAFE}
\description{
  The class SAFE is the output from the function \code{\link{safe}}. It is also the input   
  to the plotting function \code{\link{safeplot}}.
}
\section{Slots}{
  \describe{
    \item{\code{local}:}{Object of class \code{"character"} describing the local
      statistic used.}
    \item{\code{local.stat}:}{Object of class \code{"numeric"} containing
       the (unsorted) observed local statistics for genes.}
    \item{\code{local.pval}:}{Object of class \code{"numeric"} containing
       the (unsorted) empirical p-values for genes}
    \item{\code{global}:}{Object of class \code{"character"} describing the local
      statistic used.}
    \item{\code{global.stat}:}{Object of class \code{"numeric"} containing
       the (unsorted) observed global statistics for categories.}
    \item{\code{global.pval}:}{Object of class \code{"numeric"} containing
       the (unsorted) empirical p-values for categories.}
    \item{\code{error}:}{Object of class \code{"character"} describing the method
      used to estimate error rates across multiple comparisons.}
    \item{\code{global.error}:}{Object of class \code{"numeric"} containing
       the (unsorted) error rates associated with the p-values for categories. If not computed, 
       it will be set to \code{NA}.}
    \item{\code{C.mat}:}{Object of class \code{"matrix"} containing the category
      assignments. Each row corresponds to a gene, and each column a category.}
    \item{\code{alpha}:}{Object of class \code{"numeric"} containing the alpha level
      for significance of a category.}
    \item{\code{method}:}{Object of class \code{"character"} describing the resampling 
      method used in \code{safe}.}
  }
}
\section{Methods}{ 
  \describe{
    \item{show}{(gt.result): Summarizes the test results of significant categories.}
    \item{[}{(gt.result): Returns a SAFE object for categories indicated
      by integer of character strings.}
    \item{safeplot}{(gt.result): The \code{\link{safeplot}} produces a plot of the
      relative association of expression in a category of genes relative to their
      complement.}
  }
}

\author{William T Barry: \email{bill.barry@duke.edu}}

\seealso{\code{\link{safe}}, \code{\link{safeplot}}.}

\keyword{methods}

