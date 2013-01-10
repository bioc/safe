\name{getCmatrix}
\alias{getCmatrix}
\title{ Generation of a C matrix }
\description{
  This function will construct a matrix of indicator variables for category membership from keyword or gene-indexed lists. Size constraints, the option to prune identical categories, and a vector of present genes can be defined to filter categories and order genes. New to version 3.0.0, replicate gene features can be provided such that each gene has equal weight in a category.
}
\usage{
  getCmatrix(keyword.list = NULL, gene.list = NULL, 
             present.genes = NULL, GO.ont = NULL, min.size = 2, 
             max.size = Inf, by.gene = FALSE, gene.names =  NULL, 
             prune = FALSE, as.matrix = FALSE, ...)
}
\details{
  Typical usages are
  \preformatted{
  getCmatrix(keyword.list, present.genes)
  getCmatrix(gene.list, present.genes)
  }
}
\arguments{
  \item{keyword.list}{ A list containing character vectors for each keyword that specify the gene members. }
  \item{gene.list}{ A list containing character vectors for each gene that specify the annotated functional categories. }
  \item{present.genes}{ An optional vector used to filter genes in the C matrix. Can be provided as an unordered character vector of gene names that match \code{names(list)}, or as an ordered vector of presence (\code{1}) and absence (\code{0}) calls.  }
  \item{GO.ont}{ "CC", "BP", or "MF" specify the ontology to limit categories to. }
  \item{min.size}{ Optional minimum category size to be considered. }
  \item{max.size}{ Optional maximum category size to be considered. }
  \item{by.gene}{ Optional logical to build 'soft' categories at the gene level, instead of the feature level. }
  \item{gene.names}{ Optional character vector of gene names for 'soft' categories. }
  \item{prune}{ Optional logical to remove duplicate categories. }
  \item{as.matrix}{ Optional argument to specify a matrix is returned rather than a matrix.csr.} 
  \item{\dots}{Any extra arguments will be forwarded to the \code{read.table} function when category assignments are given as a file.}
}
\value{
  \item{C.mat.csr}{ If \code{as.matrix=F} a sparse matrix is returned with the rows corresponding to the genes and columns are categories}
  \item{row.names}{ Character vector of gene names}
  \item{col.names}{ Character vector of category names}
}
\references{ W. T. Barry, A. B. Nobel and F.A. Wright, 2005, \emph{Significance Analysis of functional categories in gene expression studies: a structured permutation approach}, \emph{Bioinformatics} \bold{21}(9) 1943-9. 

    See also the vignette included with this package. }
\author{ William T. Barry: \email{bbarry@jimmy.harvard.edu} }

\seealso{\code{\link{safe}}, \code{\link{safeplot}}, 
  \code{\link{getPImatrix}}.}
\examples{
require(hgu133a.db)
genes <- unlist(as.list(hgu133aSYMBOL))
RS.list <- list(Genes21 = c("ACTB","RPLP0","MYBL2","BIRC5","BAG1",
                             "GUSB","CD68","BCL2","MMP11","AURKA",
                             "GSTM1","ESR1","TFRC","PGR","CTSL2",
                             "GRB7","ERBB2","MKI67","GAPDH","CCNB1",
                             "SCUBE2"),
                Genes16 = c("MYBL2","BIRC5","BAG1","CD68","BCL2",
                            "MMP11","AURKA","GSTM1","ESR1","PGR","CTSL2",
                            "GRB7","ERBB2","MKI67","CCNB1","SCUBE2"))
RS.list <- lapply(RS.list,function(x) return(names(genes[which( match(genes, x, nomatch = 0) > 0)])))
C. <- getCmatrix(keyword.list = RS.list)
}
\keyword{ htest }