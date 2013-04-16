\name{scaler}
\alias{scaler}
\title{Row and Column Scaling of a Data Matrix}
\description{
The function provides some popular (and effective) standardization
methods for community ecologists.
}

\usage{
scaler(x, col = c("mean1", "max1", "min0", "ssq1", "range01",
	"zsc", "pa", "rank")[2], row = c("mean1", "max1", "min0",
	"ssq1", "range01", "zsc", "pa", "rank")[1])
}
\arguments{
  \item{x}{ Data matrix.}
  \item{col}{ Character vector of column standardizations.}
  \item{row}{ Character vector of row standardizations.}
}
\details{
  The function offers following data matrix standardizations:
    \itemize{
    \item \code{mean1 }: scale to mean of 1.
    \item \code{max1 }: scale to maximum of 1.
    \item \code{ssq1 }: scale to sums of sqaures equal 1.
    \item \code{range01 }: scale range to 0-1.
    \item \code{zsc }: standardize to z-scores (mean=0, sd=1).
    \item \code{pa }: scale to presence/absence scale (0/1).
    \item \code{rank }: scale to rank order (1=lowest).
}


  Standardizations are performed first on columns then on rows.
  "pa" applies to the whole matrix and can be specied using row or col.
}
\value{
  Returns the standardized matrix.
}
\author{Jari Oksanen -- modified Glenn De'ath (Dec 03)}
\note{Common transformations can be made with standard \R functions.}


\examples{
data(spider)
spid.data <- scaler(spider, col = "max", row="mean1")
}
\keyword{ multivariate}%-- one or more ...
\keyword{ manip }

