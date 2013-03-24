\name{safe}
\alias{safe}
\title{Significance Analysis of Function and Expression}
\description{
  Performs a significance analysis of function and expression (SAFE) for a gene expression experiment and a set of functional categories specified by the user. SAFE is a two-stage resampling-based method that can be applied to a 2-sample, paired, multi-class, simple linear and right-censored linear regression models. Other experimental designs can also be accommodated through user-defined functions.
}
\usage{
safe(X.mat, y.vec, C.mat = NULL, Z.mat = NULL,
     method = "permutation", platform = NULL, 
     annotate = NULL, min.size = 2, max.size = Inf, 
     by.gene = FALSE, local = "default", global = "default", 
     args.local = NULL, args.global = list(one.sided = FALSE), 
     Pi.mat = NULL, error = "FDR.BH", parallel=FALSE, alpha = NA, 
     epsilon = 10^(-10), print.it = TRUE, ...)
}
\arguments{
  \item{X.mat}{ A matrix or data.frame of expression data of size m by n where each row corresponds to a gene feature and each column to a sample. Data should be properly normalized and cannot contain missing values. }
  \item{y.vec}{ A numeric, integer or character vector of length n containing the response of interest. For examples of the acceptable forms \code{y.vec} can take, see the vignette. }
  \item{C.mat}{ A matrix containing the gene category assignments. Each column represents a category and should be named accordingly. For each column, values of 1 (\code{TRUE}) and 0 (\code{FALSE}) indicate whether the features in the corresponding rows of \code{X.mat} are contained in the category. This can also be a list containing a sparse matrix and names as created by \code{getCmatrix}. }
  \item{Z.mat}{ A data.frame of size n by p, with p covariates as numeric or factors. }
  \item{method}{ Type of hypothesis test can be specified as "permutation", "bootstrap.t", and "bootstrap.q". "express" calls the dependent package \code{safeExpress}. See vignette for details.}
  \item{platform}{ If \code{C.mat} is unspecified, a character string of a Bioconductor annotation package can be used to build gene categories. See vignette for details and examples.}
  \item{annotate}{ If \code{C.mat} is unspecified, a character string to specify the type of gene categories to build from annotation packages. "GO.MF", "GO.BP", "GO.CC", and "GO.ALL" (default) specify one or all Gene Ontologies. "KEGG" specifies pathways, and "PFAM" homologous families from the respective sources. }
  \item{min.size}{ Optional minimum category size in building \code{C.mat}. }
  \item{max.size}{ Optional maximum category size in building \code{C.mat}. }
  \item{by.gene}{ Logical argument (default = \code{FALSE}) specifying whether multiple features to a single gene should be down-weighted. }
  \item{local}{ Specifies the gene-specific statistic from the following options: "t.Student", "t.Welch", and "t.paired", for 2-sample designs, "f.ANOVA" for 1-way ANOVAs, and "t.LM" for simple linear regressions. "default" will choose between "t.Student", "f.ANOVA", and "t.LM" based on the form of \code{y.vec}. User-defined local statistics can also be used; details are provided in the vignette. }
  \item{global}{ Specifies the global statistic for a gene categories. By default, the Wilcoxon rank sum ("Wilcoxon") is used. Else, a Fisher's Exact test statistic ("Fisher"), a Pearson's chi-squared type statistic ("Pearson") or t-statistic for average difference ("AveDiff") is available. User-defined global statistics can also be implemented. }
  \item{args.local}{ An optional list to be passed to user-defined local statistics that require additional arguments. By default \code{args.local = NULL}. }
  \item{args.global}{ An optional list to be passed to global statistics that require additional arguments. For two-sided local statistics, \code{args.global} = list(one.sided=F) allows bi-directional differential expression to be considered. }
  \item{Pi.mat}{ Either an integer, or a matrix or data.frame containing the permutations. See \code{getPImatrix} for the acceptable form of a matrix or data.frame. If \code{Pi.mat} is an integer, B, then \code{safe} will generate B resamples of \code{X.mat}. }
  \item{error}{ Specifies the method for computing error rate estimates. By default, Benjamini-Hochberg step down ("FDR.BH") FDR estimates are computed. A Bonferroni ("FWER.Bonf") and Holm's step-up ("FWER.Holm") adjustment can also be specified. Under permutation, "FDR.YB" computes the Yekutieli-Benjamini FDR estimate, and "FWER.WY" computes the Westfall-Young FWER estimate. The user can also specify "none" if no error rates are desired. }
  \item{parallel}{ Logical argument (default = \code{FALSE}) specifying whether hypothesis test of \code{method} should be conducted with parallel processing. Only compatible with \code{error = "none", "FWER.Bonf",} or \code{FDR.BH}. See vignette for details.} 
  \item{alpha}{ The threshold for significant results to return. By default, alpha will be 0.05 for nominal p-values (\code{error} = "none" ), and 0.1 for adjusted p-values. }
  \item{epsilon}{ Numeric argument sets the minimum difference for ranking local and global statistics, correcting a numerical precision issue when computing empirical p-values in small data sets (n < 15). The default value is 10^(-10). }
  \item{print.it}{ Logical argument (default = \code{TRUE}) specifying whether to print progress updates to the log for permutation and bootstrap calculations. }
  \item{\dots}{ Allows arguments from version 2.0 to be ignored. }
}
\details{
  \code{safe} utilizes a general framework for testing differential expression across gene categories that allows it to be used in various experimental designs. Through structured resampling of the data, \code{safe} accounts for the unknown correlation among genes, and enables proper estimation of error rates when testing multiple categories. \code{safe} also provides statistics and empirical p-values for the gene-specific differential expression.
}
\value{
  The function returns an object of class \code{SAFE}. See help for \code{SAFE-class} for more details.
}
\references{ W. T. Barry, A. B. Nobel and F.A. Wright, 2005, \emph{Significance Analysis of functional categories in gene expression studies: a structured permutation approach}, \emph{Bioinformatics} \bold{21}(9) 1943--1949. 

    See also the vignette included with this package. }
\author{ William T. Barry: \email{bbarry@jimmy.harvard.edu} }

\seealso{\code{\link{safeplot}}, \code{\link{safe.toptable}}, \code{\link{gene.results}}, \code{\link{getCmatrix}}, \code{\link{getPImatrix}}.}
\examples{
## Simulate a dataset with 1000 genes and 20 arrays in a 2-sample design.
## The top 100 genes will be differentially expressed at varying levels

g.alt <- 100
g.null <- 900
n <- 20

data<-matrix(rnorm(n*(g.alt+g.null)),g.alt+g.null,n)
data[1:g.alt,1:(n/2)] <- data[1:g.alt,1:(n/2)] + 
                         seq(2,2/g.alt,length=g.alt)
dimnames(data) <- list(c(paste("Alt",1:g.alt),
                         paste("Null",1:g.null)),
                       paste("Array",1:n))

## A treatment vector 
trt <- rep(c("Trt","Ctr"),each=n/2)

## 2 alt. categories and 18 null categories of size 50

C.matrix <- kronecker(diag(20),rep(1,50))
dimnames(C.matrix) <- list(dimnames(data)[[1]],
    c(paste("TrueCat",1:2),paste("NullCat",1:18)))
dim(C.matrix)

results <- safe(data,trt,C.mat = C.matrix,Pi.mat = 100)
results

## SAFE-plot made for the first category
if (interactive()) { 
safeplot(results,"TrueCat 1")
}
}
\keyword{ htest }