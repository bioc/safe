\name{getCmatrix}
\alias{getCmatrix}
\title{ Generation of a C matrix }
\description{
  This function will convert a text file, vector or list of gene annotation into
  a C matrix. Size constraints, and present/absent calls can be set to filter 
  categories and genes accordingly. 
}
\synopsis{
  getCmatrix<-function(list = NULL,vector = NULL, file = NULL, delimiter = ",", 
            present.genes = NULL, min.size = 0, max.size = Inf, ...)
}
\usage{
  getCmatrix(list, present.genes=, min.size=, max.size=)
  getCmatrix(vector=, delimiter=)
  getCmatrix(file=, delimiter=, \dots)
}
\arguments{
  \item{list}{ A list containing character vectors for each gene that specify the
    functional categories it belongs to by keywords. }
  \item{vector}{ A character vector of gene annotation with a specified delimiter
    between category keywords. }
  \item{file}{ A file containing the character vector of gene annotation. }
  \item{delimiter}{ Delimiter used between category keywords when provided as a vector
    or file. }
  \item{present.genes}{ An optional vector used to filter genes in the C matrix. Can be 
     provided as an unordered character vector of gene names that match \code{names(list)},
     or as an ordered vector of presence (\code{1}) and absence (\code{0}) calls.  }
  \item{min.size}{ Optional minimum category size to be considered. }
  \item{max.size}{ Optional maximum category size to be considered. }
  \item{...}{Any extra arguments will be forwarded to the read.table function when
    category assigments are given as a file.}
}
\value{
  A matrix is returned with the rows corresponding to the genes as labeled by
  \code{names(list)} and each column corresponding to a category.
}
\references{ W. T. Barry, A. B. Nobel and F.A. Wright, 2004, \emph{Significance Analysis
    of functional categories in gene expression studies: a structured permutation approach},
    \emph{Bioinformatics} In press. 

    See also the vignette included with this package. }
\author{ William T. Barry: \email{wbarry@bios.unc.edu} }

\seealso{{\code{\link{safe}}, \code{\link{safeplot}}, 
  \code{\link{getPImatrix}},.}}
\examples{
## A simple illustration
anno <- c("Keyword1","Keyword2;Keyword3","",
          "Keyword3;Keyword1","Keyword3")
names(anno) <- paste("Gene",1:5)

getCmatrix(vector = anno, delimiter = ";")
}
\keyword{ htest }
