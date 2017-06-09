\name{readSamplesSpRfn}}
\alias{readSamplesSpRfn}
\title{hapRFN: Read Sparse Matrix Samples}
\description{

  \code{readSamplesSpRfn}: C implementation of  \code{readSamplesSpRfn}.
}
\usage{

readSamplesSpRfn(X, samples = 0, lowerB = 0.0, upperB = 1000.0)

}
\arguments{
  \item{X}{the file name of the sparse matrix in sparse format.}
  \item{samples}{vector of samples which should be read; default
    = 0 (all samples)}
  \item{lowerB}{lower bound for filtering the inputs columns, the
    minimal column sum; default = 0.0.}
  \item{upperB}{upper bound for filtering the inputs columns, the
    maximal column sum; default = 1000.0.}
}
\details{


  The data matrix is
  directly scanned by the C-code and must be in sparse matrix format.

  Sparse matrix format:
  *first line: number of nonzero elements (nnz)
  *second line: number of rows (nrow)
  *third line: row pointers, a list of integers of length nrow+1, 
    where the ith element is the number of nonzero elements 
    until the i-1th row.
  *fourth line: column indices, a list of integers of lentgth nnz, 
    where the ith element is the (zero based) column index 
    of the ith nonzero element
  *fifth line: values, the list of nonzero elements 

  The code is implemented in C.

}
\value{
  \code{X} (data of the given samples)
  }
\seealso{
\code{\link{fabia}},
}
\author{Balazs Bencze}
\examples{

}


\references{
  S. Hochreiter et al.,
  \sQuote{FABIA: Factor Analysis for Bicluster Acquisition},
  Bioinformatics 26(12):1520-1527, 2010.
  http://bioinformatics.oxfordjournals.org/cgi/content/abstract/btq227
  

}
\keyword{}
\concept{}
