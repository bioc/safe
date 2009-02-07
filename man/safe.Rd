\name{safe}
\alias{safe}
\title{Significance Analysis of Function and Expression}
\description{
  Performs a significance analysis of function and expression (SAFE) for a given
  gene expression experiment and a given set of functional categories. SAFE is a
  two-stage permutation-based method that can be applied to a 2-sample,
  multi-class, simple linear regression, and other linear models. Other
  experimental designs can also be accommodated through user-defined functions.
}
\usage{
safe(X.mat, y.vec, C.mat = NULL, platform = NULL, annotate = NULL, Pi.mat = NULL, 
     local = "default", global = "Wilcoxon", args.local = NULL, 
     args.global = list(one.sided = FALSE), error = "none", alpha = NA, 
     method = "permutation", min.size = 2, max.size = Inf, ...)
}
\arguments{
  \item{X.mat}{ A matrix or data.frame of expression data; each row corresponds to a gene
    and each column to a sample. Data can also be given as the Bioconductor class
    \code{\link[Biobase:ExpressionSet-class]{ExpressionSet}}.
    Data should be properly normalized and may not contain missing values.}
  \item{y.vec}{a numeric, integer or character vector of length \code{ncol(X.mat)}
    containing the response of interest. If \code{X.mat} is an
    \code{\link[Biobase:ExpressionSet-class]{ExpressionSet}}, \code{y.vec} can also be the name or
    column number of a covariate in the \code{\link[Biobase:phenoData-class]{phenoData}}
    slot. For examples of the acceptable forms \code{y.vec} can take, see the vignette. }
  \item{C.mat}{ A matrix or data.frame containing the gene category assignments. Each column
    represents a category and should be named accordingly. For each column, values of
    1 (\code{TRUE}) and 0 (\code{FALSE}) indicate whether the genes in the corresponding rows of
    \code{X.mat} are contained in the category. This can also be a list containing a sparse 
    matrix and dimnames as created by \code{getCmatrix}}
  \item{platform}{ If \code{C.mat} is unspecified, a character string of a Bioconductor annotation
    package can be used to build gene categories. See vignette for details and examples.}
  \item{annotate}{ If \code{C.mat} is unspecified, a character string to specify the type of gene
    categories to build from annotation packages. "GO.MF", "GO.BP", "GO.CC", and "GO.ALL" (default) 
    specify one or all Gene Ontologies. "KEGG" specifies pathways, and "PFAM" homologous families
    from the respective sources.}
 \item{Pi.mat}{ Either a matrix or data.frame containing the permutations, or an integer. See
    \code{getPImatrix} for the acceptable form of a matrix or data.frame. If \code{Pi.mat} is
    an integer, then \code{safe} will
    automatically generate as many random permutations of \code{X.mat}. }
  \item{local}{ Specifies the gene-specific statistic from the following options: "t.Student",
    "t.Welch" and "t.SAM" for 2-sample designs, "f.ANOVA" for 1-way ANOVAs, "t.LM" for
    simple linear regressions, and "z.COXPH" for a Cox
    proportional hazards survival model.  "default" will choose
    between "t.Student" and "f.ANOVA", based on the form of \code{y.vec}. User-defined local statistics
    can also be used; details are provided in the vignette. }
  \item{global}{ Specifies the global statistic for a gene categories. By default, the Wilcoxon rank sum
    ("Wilcoxon") is used. Else, a Fisher's Exact test statistic ("Fisher") based on the hypergeometric
    dist'n, a chi-squared type Pearson's test ("Pearson") or t-test of average difference ("AveDiff")
    is available. User-defined global statistics can also be implemented. }
  \item{args.local}{ An optional list to be passed to user-defined local statistics that require
    additional arguments. By default \code{args.local = NULL}.  }
  \item{args.global}{ An optional list to be passed to global statistics that require
    additional arguments. For two-sided local statistics, \code{args.global} = list(one.sided=F) allows
    bi-directional differential expression to be considered. }
  \item{error}{ Specifies the method for computing error rate estimates. "FDR.YB" computes the
    Yekutieli-Benjamini FDR estimate, "FWER.WY" computes the
    Westfall-Young FWER estimate. A Bonferroni, ("FWER.Bonf"), Holm's step-up ("FWER.Holm"),
    and Benjamini-Hochberg step down ("FDR.BH") adjustment can also be
    specified. By default ("none") no error rates are computed. }
  \item{alpha}{ Allows the user to define the criterion for significance. By default, alpha will be
    0.05 for nominal p-values (\code{error} = "none" ), and 0.1 otherwise. }
  \item{method}{ Type of hypothesis test can be specified as "permutation", "bootstrap.t", and 
                 "bootstrap.q". See vignette for details}
  \item{min.size}{ Optional minimum category size to be considered. }
  \item{max.size}{ Optional maximum category size to be considered. }
  \item{\dots}{  Allows arguments from version 1.0 to be ignored }
}
\details{
  \code{safe} utilizes a general framework for testing differential expression across gene categories
  that allows it to be used in various experimental designs. Through structured resampling of the data,
  \code{safe} accounts for the unknown correlation among
  genes, and enables proper estimation of error rates when testing multiple categories. 
  \code{safe} also provides statistics and empirical p-values for the gene-specific 
  differential expression.
}
\value{
  The function returns an object of class \code{SAFE}. See help for \code{SAFE-class} for more details.
}
\references{ W. T. Barry, A. B. Nobel and F.A. Wright, 2005, \emph{Significance Analysis
    of functional categories in gene expression studies: a structured permutation approach},
    \emph{Bioinformatics} {\bf 21}(9) 1943--1949. 

    See also the vignette included with this package. }
\author{ William T. Barry: \email{bill.barry@duke.edu} }

\seealso{\code{\link{safeplot}}, \code{\link{getCmatrix}}, 
  \code{\link{getPImatrix}}.}}
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

results <- safe(data,trt,C.matrix,Pi.mat = 100)
results

## SAFE-plot made for the first category
if (interactive()) { 
safeplot(results,"TrueCat 1")
}
}
\keyword{ htest }
