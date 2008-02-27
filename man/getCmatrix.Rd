\name{getCmatrix}
\alias{getCmatrix}
\title{ Generation of a C matrix }
\description{
  This function will convert a list, vector or file of gene annotation into
  a C matrix. Size constraints, and present/absent calls can be set to filter 
  categories and genes accordingly. 
}
\synopsis{
  getCmatrix(keyword.list = NULL, gene.list = NULL, vector = NULL, file = NULL, 
             delimiter = ",", present.genes = NULL, GO.ont = NULL, min.size = 0, 
             max.size = Inf, as.matrix = FALSE, ...)
}
\usage{
  getCmatrix(keyword.list, present.genes=, GO.ont=)
  getCmatrix(gene.list, present.genes=, min.size=, max.size=)
  getCmatrix(vector=, delimiter=, as.matrix=)
  getCmatrix(file=, delimiter=, \dots)
}
\arguments{
  \item{keyword.list}{ A list containing character vectors for each keyword that specify the
    gene members. }
  \item{gene.list}{ A list containing character vectors for each gene that specify the
    functional categories it belongs to. }
  \item{vector}{ A character vector of gene annotation with a specified delimiter
    between category keywords. }
  \item{file}{ A file containing the character vector of gene annotation. }
  \item{delimiter}{ Delimiter used between category keywords when provided as a vector
    or file. }
  \item{present.genes}{ An optional vector used to filter genes in the C matrix. Can be 
     provided as an unordered character vector of gene names that match \code{names(list)},
     or as an ordered vector of presence (\code{1}) and absence (\code{0}) calls.  }
  \item{GO.ont}{ "CC","BP",or "MF" specify the ontology to limit categories to. }
  \item{min.size}{ Optional minimum category size to be considered. }
  \item{max.size}{ Optional maximum category size to be considered. }
  \item{as.matrix}{ Optional argument to specify a matrix is returned rather than a matrix.csr.} 
  \item{...}{Any extra arguments will be forwarded to the read.table function when
    category assignments are given as a file.}
}
\value{
  \item{C.mat.csr}{ If \code{as.matrix=F} a sparse matrix is returned with the rows corresponding 
    to the genes and columns are categories}
  \item{row.names}{ Character vector of gene names}
  \item{col.names}{ Character vector of category names}
}
\references{ W. T. Barry, A. B. Nobel and F.A. Wright, 2005, \emph{Significance Analysis
    of functional categories in gene expression studies: a structured permutation approach},
    \emph{Bioinformatics} {\bf 21}(9) 1943-9. 

    See also the vignette included with this package. }
\author{ William T. Barry: \email{bill.barry@duke.edu} }

\seealso{{\code{\link{safe}}, \code{\link{safeplot}}, 
  \code{\link{getPImatrix}},.}}
\examples{
## A simple illustration
anno <- c("Keyword1","Keyword2;Keyword3","",
          "Keyword3;Keyword1","Keyword3")
names(anno) <- paste("Gene",1:5)

getCmatrix(vector = anno, delimiter = ";",
           as.matrix=TRUE)
}
\keyword{ htest }
