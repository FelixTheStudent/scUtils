% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/variance.R
\name{colVars_spm}
\alias{colVars_spm}
\alias{rowVars_spm}
\title{Variance computation for sparse matrices}
\usage{
colVars_spm(spm)

rowVars_spm(spm)
}
\arguments{
\item{spm}{A sparse matrix of class dgCMatrix from the Matrix package.}
}
\value{
Vector with variances.
}
\description{
Compute variance for each column / each row of
a dgCMatrix (from Matrix package).
}
\details{
The only supported format currently is dgCMatrix. While the Matrix
package has other formats, this one is used for scRNAseq raw count data.
Function code written by Simon Anders.
}
\examples{
library(Matrix)
 mat <- as(matrix(rpois(900,1), ncol=3), "dgCMatrix")
 colVars_spm(mat)
}
\seealso{
\code{vignette("Intro2Matrix", package="Matrix")}
\link[Matrix]{CsparseMatrix-class}
}
