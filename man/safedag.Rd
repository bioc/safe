\name{safedag}
\alias{safedag}
\title{ SAFE results displayed in Gene Ontology }
\description{
  SAFE results are displayed on the directed acyclic graph for one of
  the ontologies under investigation.  Category-wide significance
  displayed by node color.
}
\synopsis{
safedag(object = NULL, ontology = NULL, top = NULL, file = NULL,
        color.cutoffs = c(0.1, 0.01, 0.001), filter = 0, max.GOnames = 200)
}
\arguments{
  \item{object}{ Object of class \code{SAFE} }
  \item{ontology}{ Gene Ontology of interest. Character strings of
     "GO.CC", "GO.BP", and "GO.MF" accepted.}
  \item{top}{ Optional character string giving the node name from which to draw
     a subgraph of the tree}
  \item{file}{Optional filename for a post-script of the graph }
  \item{color.cutoffs}{ Numeric vector of length 3 for the cutoffs for
    coloring significant nodes. Nodes with unadjusted p-values less than
    \code{color.cutoff[3]} are drawn in blue; less than
    \code{color.cutoff[2]} are drawn in green; less than
    \code{color.cutoff[1]} are drawn in red.}
  \item{filter}{ Optional integer (1,2,3) to only include branches that
    contain at least one node as significant as the respective color.cutoff. }
  \item{max.GOnames}{ Maximum size of DAG to include category names as labels.}
}
\details{
  DAG-plots are suggested as a means for visualizing the extent of differential expression
  in Gene Ontology categories. The relatedness of significant categories
  suggests whether similar or disparate biological findings are
  identified.
}
\references{ W. T. Barry, A. B. Nobel and F.A. Wright, 2005, \emph{Significance Analysis
    of functional categories in gene expression studies: a structured permutation approach},
    \emph{Bioinformatics} {\bf 21}(9) 1943--1949. 

  See also the vignette included with this package. }
\author{ William T. Barry: \email{bill.barry@duke.edu} }

\seealso{{\code{\link{safe}}.}}
\keyword{ hplot }
